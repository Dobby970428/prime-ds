"""
PRIME-DS ICH M7 Hybrid Ensemble Matrix — Layers A (DB) + B (SMARTS) + C (QSAR) + D (read-across).

Automated consensus rescoring with conservative toxicity profile and
seamless fallback when deep-learning inference is unavailable.

2026-09-08 추가 — Layer D (`p69_read_across_arbitration`)
    Layer B(전문가 규칙)와 Layer C(통계 모델)의 판정이 **엇갈릴 때에 한하여**
    실측 풀에서 구조 유사체를 찾아 제3의 근거로 삼습니다. ICH M7(R2)이
    «두 방법론이 상충하면 전문가 검토»를 요구하는데, 종전 구현은 상충 자체를
    정보로 쓰지 않고 임계값만 보았습니다.

    유사체가 양성으로 몰리면 승급을 제안하고(승급 전용 원칙 유지),
    음성으로 몰리면 **자동으로 내리지 않고** 하향 후보와 근거만 기록한 뒤
    전문가 검토를 요구합니다.
"""

from __future__ import annotations

import sqlite3
import time
from dataclasses import dataclass, field
from typing import Any

from rdkit import Chem

from p64_gnn_mutagenicity_predictor import GnnMutagenicityResult, predict_mutagenicity_probability
from p66_ich_m7 import (
    _hit_from_db_query,
    _shared_alerts,
    route_m7_consensus_class_detailed,
    scan_alerts_for_smiles,
)
from p60_m7_empirical_db import EmpiricalM7Hit, get_default_db_connection, query_empirical_toxicology_db
from p69_read_across_arbitration import ReadAcrossArbitration, arbitrate as _arbitrate_read_across

# 2026-09-11 추가 — 적용범위(applicability domain) 보고.
# 판정 결과를 바꾸지 않는 부가 정보이므로, 모듈이 없거나 산출에 실패해도
# 기존 동작이 그대로 유지되어야 한다.
try:
    from p73_applicability_domain import applicability_domain as _applicability_domain

    _AD_AVAILABLE = True
except Exception:  # noqa: BLE001
    _applicability_domain = None  # type: ignore[assignment]
    _AD_AVAILABLE = False

# Layer B — enzymatic / metabolization toxicophore extensions
METABOLIZATION_FRAGMENT_SMARTS: tuple[tuple[str, str], ...] = (
    ("Aromatic N-oxide (CYP metabolite)", "[n+][O-]"),
    ("N-hydroxy intermediate", "[#7]-[#8]"),
    ("Epoxide bioactivation", "[OX2r3]"),
    ("Aldehyde oxidase substrate", "[CH1](=O)[#6]"),
)

AI_INFERENCE_THRESHOLD = 0.70
AI_ESCALATION_THRESHOLD = 0.85
AI_INFERRED_LABEL = "Class 3-Alert (AI Inferred Liability)"


@dataclass
class EnsembleLayerA:
    """Exact-match relational record (Tox21 / Hansen / ISSCAN matrix)."""

    is_in_db: bool
    ames_result: str | None
    carc_result: str | None
    empirical: EmpiricalM7Hit
    latency_ms: float = 0.0


@dataclass
class EnsembleLayerB:
    """Deterministic SMARTS + metabolization fragment scan."""

    alert_present: bool
    alerts: list[str]
    metabolization_fragments: list[str] = field(default_factory=list)

    @property
    def combined_alerts(self) -> list[str]:
        return list(dict.fromkeys(self.alerts + self.metabolization_fragments))


@dataclass
class EnsembleLayerC:
    """AttentiveFP / graph deep-learning P(Mutagenic)."""

    gnn: GnnMutagenicityResult

    @property
    def p_mutagenic(self) -> float | None:
        return self.gnn.p_mutagenic


@dataclass
class EnsembleMatrixResult:
    layer_a: EnsembleLayerA
    layer_b: EnsembleLayerB
    layer_c: EnsembleLayerC
    # 2026-09-08: Layer D — B/C 상충 시에만 동작하는 read-across 중재층
    layer_d: ReadAcrossArbitration
    rule_class: int
    final_class: int
    class_label: str
    ai_inferred_liability: bool
    description: str
    preliminary: bool
    ensemble_log: list[str] = field(default_factory=list)
    shared_alerts: list[str] = field(default_factory=list)
    # 2026-09-03 수정: 규칙 기반 Class 4와 통계 모델이 충돌해도 기록이 남지 않던 문제.
    expert_review_required: bool = False
    expert_review_reason: str = ""
    # 2026-09-11 추가: Layer C 적용범위. ApplicabilityDomainResult 또는 None.
    # 등급 부여에 관여하지 않고, 「이 예측을 얼마나 믿을 수 있는가」만 보고한다.
    ad: Any = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "Ensemble Layer A (DB Hit)": self.layer_a.is_in_db,
            "Ensemble Layer B (Alert Present)": self.layer_b.alert_present,
            "Ensemble Layer B (Metabolization Fragments)": ", ".join(
                self.layer_b.metabolization_fragments
            )
            or "None",
            "Ensemble Layer C P(Mutagenic)": self.layer_c.p_mutagenic,
            "Ensemble Layer C Model": self.layer_c.gnn.model_name,
            "Rule-Based Class": f"Class {self.rule_class}",
            "Final Assigned Class": self.class_label,
            "AI Inferred Liability": self.ai_inferred_liability,
            "Ensemble Log": " | ".join(self.ensemble_log),
            # 2026-09-03 수정: evaluate_ich_m7_class가 이 dict로 최종 페이로드를
            # 덮어쓰므로, 누락돼 있던 shared_alerts·description·preliminary·
            # Layer C fallback 사유를 함께 내보냅니다.
            "Shared Alerts With Parent API": ", ".join(self.shared_alerts)
            if self.shared_alerts
            else "None",
            "Ensemble Description": self.description,
            "Is Preliminary": self.preliminary,
            "Ensemble Layer C Fallback Reason": self.layer_c.gnn.fallback_reason or "",
            "Expert Review Required": self.expert_review_required,
            "Expert Review Reason": self.expert_review_reason,
            **self.layer_d.to_dict(),
            # 2026-09-11 추가: 적용범위. 산출 실패 시 키 자체가 나타나지 않는다.
            **(self.ad.to_dict() if self.ad is not None else {}),
        }


def _scan_metabolization_fragments(mol: Chem.Mol) -> list[str]:
    hits: list[str] = []
    for name, smarts in METABOLIZATION_FRAGMENT_SMARTS:
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat and mol.HasSubstructMatch(pat):
                hits.append(name)
        except Exception:
            continue
    return hits


def _run_layer_a(
    conn: sqlite3.Connection,
    canonical_smiles: str,
) -> EnsembleLayerA:
    t0 = time.perf_counter()
    is_in_db, ames, carc, meta = query_empirical_toxicology_db(conn, canonical_smiles)
    latency = (time.perf_counter() - t0) * 1000
    empirical = _hit_from_db_query(
        canonical_smiles, is_in_db, ames, carc, meta, latency_ms=latency
    )
    return EnsembleLayerA(
        is_in_db=is_in_db,
        ames_result=ames,
        carc_result=carc,
        empirical=empirical,
        latency_ms=latency,
    )


def _run_layer_b(
    mol: Chem.Mol,
    canonical_smiles: str,
    *,
    precomputed_alerts: list[str] | None,
) -> EnsembleLayerB:
    alert_present, alerts = scan_alerts_for_smiles(
        canonical_smiles, mol=mol, precomputed=precomputed_alerts
    )
    metab = _scan_metabolization_fragments(mol)
    combined_present = alert_present or bool(metab)
    return EnsembleLayerB(
        alert_present=combined_present,
        alerts=alerts,
        metabolization_fragments=metab,
    )


def _run_layer_c(canonical_smiles: str) -> EnsembleLayerC:
    try:
        gnn = predict_mutagenicity_probability(canonical_smiles)
    except Exception as exc:
        from p64_gnn_mutagenicity_predictor import GnnMutagenicityResult as _GnnResult

        gnn = _GnnResult(
            p_mutagenic=None,
            model_name="ECFP-2048 + RandomForest (Tox21 SR-p53 fallback)",
            available=False,
            fallback_used=True,
            fallback_reason=f"Layer C error — RDKit/DB standalone: {exc}"[:240],
        )
    return EnsembleLayerC(gnn=gnn)


def _dl_suggested_class(p: float | None) -> int | None:
    if p is None:
        return None
    if p >= AI_ESCALATION_THRESHOLD:
        return 3
    if p >= AI_INFERENCE_THRESHOLD:
        return 3
    if p >= 0.50:
        return 3
    return 5


def rescoring_consensus_matrix(
    *,
    layer_a: EnsembleLayerA,
    layer_b: EnsembleLayerB,
    layer_c: EnsembleLayerC,
    rule_class: int,
    rule_description: str,
    rule_preliminary: bool,
    shared_alerts: list[str],
) -> tuple[int, str, bool, bool, str, list[str]]:
    """
    Conservative ensemble rescoring — DL signal can escalate but not downgrade
    verified empirical Class 1/2 outcomes.

    하위 호환 래퍼 — 6-튜플을 그대로 반환합니다. 전문가 검토 플래그까지 필요하면
    rescoring_consensus_matrix_detailed()를 쓰세요.
    """
    (
        final_class,
        class_label,
        ai_inferred,
        preliminary,
        description,
        log,
        _review,
        _reason,
    ) = rescoring_consensus_matrix_detailed(
        layer_a=layer_a,
        layer_b=layer_b,
        layer_c=layer_c,
        rule_class=rule_class,
        rule_description=rule_description,
        rule_preliminary=rule_preliminary,
        shared_alerts=shared_alerts,
    )
    return final_class, class_label, ai_inferred, preliminary, description, log


def rescoring_consensus_matrix_detailed(
    *,
    layer_a: EnsembleLayerA,
    layer_b: EnsembleLayerB,
    layer_c: EnsembleLayerC,
    rule_class: int,
    rule_description: str,
    rule_preliminary: bool,
    shared_alerts: list[str],
    layer_d: ReadAcrossArbitration | None = None,
) -> tuple[int, str, bool, bool, str, list[str], bool, str]:
    """
    2026-09-03 수정: elif 사슬을 «각 분기가 제안한 class의 최솟값» 방식으로 재구성.

    이전 구조의 문제
      · `max(rule_class, dl_class)` 때문에 `final_class > rule_class`가 결코
        참이 되지 않아 승급 분기 전체가 죽은 코드였습니다(ICH M7은 숫자가
        작을수록 위험하므로 승급 = 최솟값).
      · 그 elif가 뒤 분기를 가려 P=0.60은 Class 3, P=0.99는 Class 5가 되는
        비단조 거동이 있었습니다.
    이제 모든 분기가 후보를 «제안»하고 최솟값을 취하므로 분기 순서와 무관하며,
    P가 커질수록 후보 집합이 커지기만 하므로 단조입니다.

    단, Class 4(모약물 실측 Ames 음성 근거)는 통계 모델만으로 승급시키지 않고
    preliminary + 전문가 검토 로그로 남깁니다(ICH M7의 실측·전문가 판단 우선).
    """
    log: list[str] = [
        f"Layer A: DB hit={layer_a.is_in_db} ames={layer_a.ames_result} carc={layer_a.carc_result}",
        f"Layer B: alert={layer_b.alert_present} alerts={len(layer_b.combined_alerts)}",
    ]
    # 2026-09-03 수정: 인자로 받고도 본문에서 전혀 쓰이지 않던 shared_alerts를
    # 판정 로그에 반영합니다(Class 4 근거의 출처).
    if shared_alerts:
        log.append(f"Shared alerts with parent API: {', '.join(shared_alerts)}")

    p = layer_c.p_mutagenic
    if p is not None:
        log.append(
            f"Layer C: P(Mutagenic)={p:.3f} model={layer_c.gnn.model_name}"
            + (" [fallback]" if layer_c.gnn.fallback_used else "")
        )
    else:
        log.append(f"Layer C: unavailable ({layer_c.gnn.fallback_reason or 'N/A'})")

    if layer_d is not None and layer_d.triggered:
        log.append(
            f"Layer D: disagreement={layer_d.disagreement} verdict={layer_d.verdict} "
            f"analogues={layer_d.n_analogues}"
            + (
                f" weighted_positive_fraction={layer_d.weighted_positive_fraction:.2f}"
                if layer_d.weighted_positive_fraction is not None
                else ""
            )
        )
        log.append(f"Layer D rationale: {layer_d.rationale}")

    expert_review_required = False
    expert_review_reason = ""

    if rule_class in (1, 2):
        log.append(f"Empirical Class {rule_class} locked — no DL downgrade.")
        return (
            rule_class,
            f"Class {rule_class}",
            False,
            rule_preliminary,
            rule_description,
            log,
            expert_review_required,
            expert_review_reason,
        )

    ai_inferred = False
    preliminary = rule_preliminary

    # (제안 class, 라벨, 설명, 로그, AI 추론 여부) — 순서는 동점 시 우선순위일 뿐
    proposals: list[tuple[int, str, str, str, bool]] = []

    if rule_class == 4:
        # Class 4는 실측(모약물 Ames 음성) 근거이므로 통계 모델만으로 승급하지 않음.
        if p is not None and p >= 0.50:
            preliminary = True
            expert_review_required = True
            expert_review_reason = (
                "expert review required — statistical model conflicts with "
                f"rule-based Class 4 (P(Mutagenic)={p:.2f}"
                + (f"; shared alerts: {', '.join(shared_alerts)}" if shared_alerts else "")
                + ")"
            )
            log.append(expert_review_reason)
    else:
        if (
            not layer_b.alert_present
            and not layer_a.is_in_db
            and p is not None
            and p > AI_INFERENCE_THRESHOLD
            and rule_class >= 5
        ):
            proposals.append(
                (
                    3,
                    AI_INFERRED_LABEL,
                    # 2026-09-03 수정: 실제 모델은 RandomForest이므로 "graph deep
                    # learning"이라는 사용자 노출 문구를 사실대로 고쳤습니다.
                    f"Hidden mutagenic liability inferred by the statistical (Q)SAR model "
                    f"(P(Mutagenic)={p:.2f} > {AI_INFERENCE_THRESHOLD:.2f}); "
                    "no SMARTS alert or DB hit — conservative Class 3 escalation.",
                    "AI override: Class 3-Alert (AI Inferred Liability).",
                    True,
                )
            )
        if p is not None and p >= AI_ESCALATION_THRESHOLD and rule_class >= 5:
            dl_class = _dl_suggested_class(p) or rule_class
            proposals.append(
                (
                    dl_class,
                    f"Class {dl_class}",
                    f"{rule_description} Statistical (Q)SAR consensus escalation "
                    f"(P(Mutagenic)={p:.2f} ≥ {AI_ESCALATION_THRESHOLD:.2f}).",
                    f"(Q)SAR escalation: rule Class {rule_class} → Class {dl_class}.",
                    False,
                )
            )
        if p is not None and p >= 0.50 and layer_b.alert_present and rule_class >= 5:
            proposals.append(
                (
                    3,
                    "Class 3",
                    f"{rule_description} Reinforced by statistical (Q)SAR signal "
                    f"P(Mutagenic)={p:.2f} with structural alerts.",
                    "(Q)SAR reinforcement: alert + P≥0.50 → Class 3.",
                    False,
                )
            )

    # ------------------------------------------------------------------
    # Layer D — read-across 중재 (B/C 상충 시에만 triggered)
    # ------------------------------------------------------------------
    if layer_d is not None and layer_d.triggered and rule_class not in (1, 2, 4):
        if layer_d.suggested_class is not None and layer_d.suggested_class < rule_class:
            proposals.append(
                (
                    layer_d.suggested_class,
                    f"Class {layer_d.suggested_class}",
                    f"{rule_description} Read-across arbitration: {layer_d.rationale}",
                    f"Read-across escalation: rule Class {rule_class} → "
                    f"Class {layer_d.suggested_class}.",
                    False,
                )
            )
        if layer_d.expert_review_required:
            expert_review_required = True
            if not expert_review_reason:
                expert_review_reason = (
                    "expert review required — expert rule and statistical layers "
                    f"disagree ({layer_d.disagreement}); {layer_d.rationale}"
                )
    elif layer_d is not None and layer_d.triggered and rule_class == 4:
        # Class 4 는 실측 근거이므로 유사체로 흔들지 않되, 상충 사실은 남깁니다.
        expert_review_required = True
        if not expert_review_reason:
            expert_review_reason = (
                "expert review required — rule-based Class 4 with a Layer B / Layer C "
                f"disagreement; read-across verdict: {layer_d.verdict}"
            )

    final_class = min([rule_class] + [c for c, _l, _d, _g, _ai in proposals])
    class_label = f"Class {rule_class}"
    description = rule_description

    if final_class < rule_class:
        # 최솟값에 도달한 첫 제안을 채택 — 목록 순서가 고정이므로 결과는 결정적입니다.
        for cls, label, desc, entry, ai_flag in proposals:
            if cls == final_class:
                class_label = label
                description = desc
                ai_inferred = ai_flag
                log.append(entry)
                break
        preliminary = True

    return (
        final_class,
        class_label,
        ai_inferred,
        preliminary,
        description,
        log,
        expert_review_required,
        expert_review_reason,
    )


def run_m7_ensemble_matrix(
    smiles: str,
    db_connection: sqlite3.Connection | None = None,
    *,
    drug_smiles: str | None = None,
    alerts: list[str] | None = None,
    drug_has_negative_mutagenicity_data: bool = False,
) -> EnsembleMatrixResult:
    """
    Three-pronged ensemble routing for a single canonical structure.
    Thread-safe DB + GNN layers; preserves RDKit-only fallback path.
    """
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        empty = EmpiricalM7Hit()
        gnn = GnnMutagenicityResult(
            p_mutagenic=None,
            model_name="none",
            available=False,
            fallback_used=True,
            fallback_reason="Invalid SMILES",
        )
        return EnsembleMatrixResult(
            layer_a=EnsembleLayerA(False, None, None, empty),
            layer_b=EnsembleLayerB(False, [], []),
            layer_c=EnsembleLayerC(gnn=gnn),
            layer_d=ReadAcrossArbitration(),
            rule_class=5,
            final_class=5,
            class_label="Class 5",
            ai_inferred_liability=False,
            description="Invalid SMILES",
            preliminary=True,
            ensemble_log=["Invalid input"],
        )

    canon = Chem.MolToSmiles(mol, canonical=True)
    conn = db_connection or get_default_db_connection()

    layer_a = _run_layer_a(conn, canon)
    layer_b = _run_layer_b(mol, canon, precomputed_alerts=alerts)
    layer_c = _run_layer_c(canon)

    drug_mol = None
    drug_empirical: EmpiricalM7Hit | None = None
    if drug_smiles:
        try:
            drug_mol = Chem.MolFromSmiles(drug_smiles.strip())
            if drug_mol is not None:
                d_can = Chem.MolToSmiles(drug_mol, canonical=True)
                d_in, d_ames, d_carc, d_meta = query_empirical_toxicology_db(conn, d_can)
                drug_empirical = _hit_from_db_query(d_can, d_in, d_ames, d_carc, d_meta)
        except Exception:
            drug_mol = None

    shared = _shared_alerts(mol, drug_mol)
    combined = layer_b.combined_alerts

    # 2026-09-03 수정: 전문가 검토 플래그까지 받기 위해 상세판 라우터를 사용합니다.
    decision = route_m7_consensus_class_detailed(
        alert_present=layer_b.alert_present,
        alerts=combined,
        is_in_db=layer_a.is_in_db,
        ames_result=layer_a.ames_result,
        carc_result=layer_a.carc_result,
        empirical=layer_a.empirical,
        drug_empirical=drug_empirical,
        shared_alerts=shared,
        drug_has_negative_mutagenicity_data=drug_has_negative_mutagenicity_data,
    )
    rule_class = decision.class_number
    rule_desc = decision.rationale
    rule_prelim = decision.preliminary

    # Layer D — Layer B 와 Layer C 가 엇갈릴 때에만 실측 유사체를 찾습니다.
    layer_d = _arbitrate_read_across(
        canon,
        alert_present=layer_b.alert_present,
        p_mutagenic=layer_c.p_mutagenic,
        layer_a_hit=layer_a.is_in_db,
        rule_class=rule_class,
    )

    (
        final_class,
        class_label,
        ai_inferred,
        preliminary,
        description,
        log,
        expert_review,
        expert_reason,
    ) = rescoring_consensus_matrix_detailed(
        layer_a=layer_a,
        layer_b=layer_b,
        layer_c=layer_c,
        rule_class=rule_class,
        rule_description=rule_desc,
        rule_preliminary=rule_prelim,
        shared_alerts=shared,
        layer_d=layer_d,
    )
    # 규칙 기반 판정 단계에서 이미 전문가 검토가 필요하다고 나온 경우도 합칩니다.
    if decision.expert_review_required:
        expert_review = True
        expert_reason = expert_reason or decision.expert_review_reason
        log.append(f"Rule-based expert review flag: {decision.expert_review_reason}")

    # 2026-09-11 추가 — Layer C 적용범위.
    # ICH M7(R2)은 통계 모델의 적용범위를 벗어나는 화합물을 전문가 검토로
    # 회부할 것을 요구한다. 등급 자체는 건드리지 않고 플래그만 세운다.
    ad_result = None
    if _AD_AVAILABLE:
        try:
            ad_result = _applicability_domain(canon)
            if ad_result.available:
                log.append(f"Layer C AD: {ad_result.summary}")
                if ad_result.in_domain is False:
                    expert_review = True
                    ad_reason = (
                        f"Outside Layer C applicability domain "
                        f"(d_max = {ad_result.d_max:.3f} < {ad_result.threshold:.2f}); "
                        f"ICH M7(R2) requires expert review."
                    )
                    expert_reason = (
                        f"{expert_reason} | {ad_reason}" if expert_reason else ad_reason
                    )
                    log.append("AD flag: out of domain → expert review required.")
            else:
                log.append(f"Layer C AD: not evaluated ({ad_result.reason})")
        except Exception as exc:  # noqa: BLE001
            ad_result = None
            log.append(f"Layer C AD: skipped ({type(exc).__name__})")

    return EnsembleMatrixResult(
        layer_a=layer_a,
        layer_b=layer_b,
        layer_c=layer_c,
        layer_d=layer_d,
        rule_class=rule_class,
        final_class=final_class,
        class_label=class_label,
        ai_inferred_liability=ai_inferred,
        description=description,
        preliminary=preliminary,
        ensemble_log=log,
        shared_alerts=shared,
        expert_review_required=expert_review,
        expert_review_reason=expert_reason,
        ad=ad_result,
    )
