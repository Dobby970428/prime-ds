"""
ICH M7(R2) Consensus Classifier — empirical DB + Ashby-Tennant alert routing.

Integrates Hansen / ISSCAN / Tox21 seed records (SQLite) with structural alerts
for deterministic Class 1–5 assignment per ICH M7 consensus framework.
"""

from __future__ import annotations

import sqlite3
import time
from dataclasses import dataclass, field
from typing import Any

from rdkit import Chem

from p60_m7_empirical_db import (
    EmpiricalM7Hit,
    get_default_db_connection,
    lookup_empirical_m7,
    query_empirical_toxicology_db,
)
from p62_toxicity_alerts import alert_exception_applies, scan_ashby_tennant_alerts

# 2026-09-03 수정: ICH M7 용어 단일 출처(p16)를 실제로 참조하도록 import.
from p16_ich_m7_terms import (
    COHORT_OF_CONCERN_AI_TEXT,
    COHORT_OF_CONCERN_AI_TEXT_KO,
    class_action,
    class_definition,
    class_label,
    is_cohort_of_concern_class,
)

# Core M7 alert patterns (canonical SMARTS) — Operation A
# 2026-09-03 수정: 두 패턴이 의도한 작용기를 잡지 못했습니다.
#   - "N-nitroso": 재귀 배제 !$(N[*]=O) 가 니트로사민의 아민 질소 자기 자신을
#     제외해 NDMA(CN(C)N=O)에 매치되지 않았습니다. 바로 아래 이름만 다른
#     중복 항목이 가려주고 있었을 뿐이므로 중복을 없애고 하나로 합쳤습니다.
#   - "Aromatic amine": [c][N;H2,H1;!$(NC=O)] 는 설폰아미드·히드라진 질소까지
#     방향족 아민으로 잡았습니다. p62_toxicity_alerts 의 교정 패턴과 맞췄습니다.
CORE_M7_ALERT_SMARTS: tuple[tuple[str, str], ...] = (
    ("Aromatic amine", "[NX3;H2,H1;!$(N[C,S,P]=[O,S]);!$(N-[N,O])]-c"),
    ("Nitro group", "[N+](=O)[O-]"),
    ("Alkyl halide", "[Cl,Br,I][CH2;!$(C=O)]"),
    ("Epoxide", "[OX2r3]"),
    ("N-nitroso", "[#7;!$([#7]=O);!$([#7+][O-])]-[#7X2]=[#8]"),
)

# 2026-09-03 수정: Class별 기본 한도 문구. 상수 dict 조회로는 cohort of
# concern(§7.5) 판정을 반영할 수 없어 아래 m7_ttc_threshold() 함수로 감쌌습니다.
_M7_TTC_BASE_TEXT: dict[int, str] = {
    1: (
        "Compound-specific Acceptable Intake (AI) — TTC not applicable "
        "(Class 1 verified mutagenic carcinogen)"
    ),
    2: "TTC 1.5 µg/day (lifetime exposure) or staged TTC per ICH M7(R2)",
    3: "TTC 1.5 µg/day — confirmatory bacterial mutagenicity (Ames) recommended",
    4: (
        "ICH Q3A/Q3B qualification limits — impurity non-mutagenicity justified "
        "by parent API or empirical negative Ames lineage"
    ),
    5: "ICH Q3A/Q3B general impurity limits (non-mutagenic / sufficient negative data)",
}


# 2026-09-03 수정: NDMA/NDEA 같은 cohort of concern에 TTC 1.5 µg/day가
# 그대로 인쇄되던 문제. ICH M7(R2) §7.5는 CoC를 TTC 적용 대상에서 제외하고
# 화합물별 AI를 요구합니다(NDMA 실제 한도 96 ng/day = TTC의 1/15).
def m7_ttc_threshold(class_number: Any, alerts: list[str] | None = None) -> str:
    """(class_number, alerts) → 권고 한도 문구.

    구조 경보에 cohort of concern(aflatoxin-like · N-nitroso · alkyl-azoxy)이
    포함되면 Class와 무관하게 TTC 대신 화합물별 AI 문구를 반환합니다.
    """
    if is_cohort_of_concern_class(alerts):
        return COHORT_OF_CONCERN_AI_TEXT
    try:
        n = int(class_number)
    except (TypeError, ValueError):
        return _M7_TTC_BASE_TEXT[3]
    return _M7_TTC_BASE_TEXT.get(n, _M7_TTC_BASE_TEXT[3])


# 하위 호환 별칭 — dict가 아니라 호출 가능한 함수입니다(.get() 조회는 더 이상 불가).
M7_TTC_THRESHOLDS = m7_ttc_threshold

# 2026-09-03 수정: p16(ICH M7 용어 단일 출처)의 복제본이었던 dict를
# p16 헬퍼로부터 생성하도록 변경. 공개 키(name/name_ko/definition/control)는 유지.
_M7_CONTROL_DETAIL_KO: dict[int, str] = {
    1: "TTC 미적용. 화합물별 Acceptable Intake(AI) 수준으로 관리",
    2: "TTC 1.5 μg/일 (평생 노출 기준) 또는 단계별 TTC 적용",
    3: "TTC 1.5 μg/일 (평생 노출) — 추가 in silico / Ames 검토 권고",
    4: "ICH Q3A/Q3B 한도 적용",
    5: "ICH Q3A/Q3B 일반 불순물 한도 적용",
}

ICH_M7_CLASS_INFO: dict[int, dict[str, str]] = {
    n: {
        "name": class_label(n),
        "name_ko": class_label(n, lang="ko"),
        "definition": class_definition(n),
        "control": f"{class_action(n, lang='ko')} — {_M7_CONTROL_DETAIL_KO[n]}",
    }
    for n in (1, 2, 3, 4, 5)
}


@dataclass
class ICHM7Result:
    class_number: int
    class_name_ko: str
    definition: str
    control_strategy: str
    rationale: str
    alerts: list[str] = field(default_factory=list)
    shared_alerts_with_drug: list[str] = field(default_factory=list)
    is_preliminary: bool = True
    alert_present: bool = False
    empirical_db_match: str = "None"
    ames_evidence: str = "Not evaluated"
    carcinogenicity_evidence: str = "Not evaluated"
    recommended_ttc: str = ""
    consensus_sources: list[str] = field(default_factory=list)
    empirical_compound_name: str = ""
    classifier_engine: str = "ICH M7 Hybrid Ensemble (DB + SMARTS + GNN)"
    gnn_mutagenic_probability: float | None = None
    gnn_model_name: str = ""
    gnn_fallback_log: str = ""
    ai_inferred_liability: bool = False
    assigned_class_label: str = ""
    rule_based_class: int | None = None
    ensemble_log: list[str] = field(default_factory=list)
    # 2026-09-03 수정: CoC 판정과 전문가 검토 요구를 결과에 남기기 위한 필드 추가.
    cohort_of_concern: bool = False
    expert_review_required: bool = False
    expert_review_reason: str = ""

    @property
    def badge_color(self) -> str:
        return {1: "#b71c1c", 2: "#e65100", 3: "#f57c00", 4: "#1565c0", 5: "#2e7d32"}.get(
            self.class_number, "#546e7a"
        )

    def to_outcome_dict(self) -> dict[str, Any]:
        """Report serialization — Assigned Class + empirical evidence + TTC."""
        return {
            "Assigned Class": f"Class {self.class_number}",
            # 2026-09-03 수정: 미등록 Class로 KeyError가 나던 직접 인덱싱을 .get()으로 교체.
            "ICH M7 Classification": ICH_M7_CLASS_INFO.get(
                self.class_number, ICH_M7_CLASS_INFO[5]
            )["name"],
            # 2026-09-15 수정: p60 과 함께 «(Hansen/ISSCAN)» 고정 라벨을 뗐습니다.
            #   실측 DB 5,324건 중 5,274건은 TDC 가 컴파일한 Ames 학습 분할이고
            #   1차 실험 기록은 50건뿐인데, 라벨이 전부를 Hansen/ISSCAN 으로
            #   보이게 하고 있었습니다. 출처 문자열을 그대로 노출합니다.
            "Empirical DB Source Match": self.empirical_db_match,
            "Ames Evidence": self.ames_evidence,
            "Carcinogenicity Evidence": self.carcinogenicity_evidence,
            # 2026-09-03 수정: dict 조회 → CoC를 반영하는 m7_ttc_threshold() 호출.
            "Recommended Action Threshold (TTC Guidelines)": self.recommended_ttc
            or m7_ttc_threshold(self.class_number, self.alerts),
            "Cohort of Concern (ICH M7 §7.5)": self.cohort_of_concern,
            "Expert Review Required": self.expert_review_required,
            "Expert Review Reason": self.expert_review_reason,
            "Structural Alerts": ", ".join(self.alerts) if self.alerts else "None detected",
            "Alert Present": self.alert_present,
            "Shared Alerts With Parent API": ", ".join(self.shared_alerts_with_drug)
            if self.shared_alerts_with_drug
            else "None",
            "Consensus Engine": self.classifier_engine,
            "Is Preliminary": self.is_preliminary,
            "P(Mutagenic) GNN": self.gnn_mutagenic_probability,
            "GNN Model": self.gnn_model_name,
            "AI Inferred Liability": self.ai_inferred_liability,
            "Assigned Class Label": self.assigned_class_label or f"Class {self.class_number}",
            "Rule-Based Class": f"Class {self.rule_based_class}"
            if self.rule_based_class is not None
            else f"Class {self.class_number}",
        }


def _scan_core_alerts(mol: Chem.Mol) -> list[str]:
    """Operation A — canonical SMARTS alert scan."""
    ashby, _ = scan_ashby_tennant_alerts(mol)
    if ashby:
        return ashby
    hits: list[str] = []
    for name, smarts in CORE_M7_ALERT_SMARTS:
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat and mol.HasSubstructMatch(pat):
                # 2026-09-29: 이 목록은 위 Ashby-Tennant 스캔이 빈 경우의
                # 폴백이며 같은 독성단을 짧은 이름으로 한 번 더 들고 있다.
                # 큐레이션된 예외를 여기에도 적용하지 않으면, 예외로 꺼둔
                # 경보가 폴백에서 되살아나 예외가 무력화된다.
                if alert_exception_applies(mol, name):
                    continue
                hits.append(name)
        except Exception:
            continue
    return hits


def _shared_alerts(impurity_mol: Chem.Mol, drug_mol: Chem.Mol | None) -> list[str]:
    if drug_mol is None:
        return []
    imp_alerts = _scan_core_alerts(impurity_mol)
    drug_alerts = _scan_core_alerts(drug_mol)
    return sorted(set(imp_alerts) & set(drug_alerts))


def _build_result(
    class_num: int,
    *,
    rationale: str,
    alerts: list[str],
    shared: list[str],
    is_preliminary: bool,
    empirical: EmpiricalM7Hit,
    alert_present: bool,
    ensemble: Any | None = None,
    expert_review_required: bool = False,
    expert_review_reason: str = "",
) -> ICHM7Result:
    info = ICH_M7_CLASS_INFO.get(class_num, ICH_M7_CLASS_INFO[3])
    gnn_p = None
    gnn_model = ""
    gnn_fb = ""
    ai_inf = False
    assigned_label = f"Class {class_num}"
    rule_cls: int | None = None
    ens_log: list[str] = []

    if ensemble is not None:
        gnn_p = ensemble.layer_c.p_mutagenic
        gnn_model = ensemble.layer_c.gnn.model_name
        gnn_fb = ensemble.layer_c.gnn.fallback_reason or ""
        ai_inf = ensemble.ai_inferred_liability
        # 2026-09-03 수정: p16 헬퍼 class_label과 이름이 충돌하던 지역변수를 개명.
        assigned_label = ensemble.class_label
        rule_cls = ensemble.rule_class
        ens_log = list(ensemble.ensemble_log)
        if getattr(ensemble, "expert_review_required", False):
            expert_review_required = True
            expert_review_reason = (
                expert_review_reason or getattr(ensemble, "expert_review_reason", "")
            )

    # 2026-09-03 수정: cohort of concern이면 TTC 문구 대신 화합물별 AI 문구를 넣고
    # 판정 근거(rationale)에도 그 사실을 남깁니다 (ICH M7(R2) §7.5).
    coc = is_cohort_of_concern_class(alerts)
    ttc_text = m7_ttc_threshold(class_num, alerts)
    control_text = info["control"]
    if coc:
        control_text = f"{COHORT_OF_CONCERN_AI_TEXT_KO} (Class {class_num} 판정)"
        rationale = (
            f"{rationale} Cohort of concern structure detected "
            f"({', '.join(alerts) or 'N/A'}) — the TTC is explicitly excluded by "
            "ICH M7(R2) §7.5; a compound-specific Acceptable Intake (AI) is required."
        )
        ens_log = list(ens_log) + [
            "Cohort of concern (ICH M7 §7.5): TTC not applicable — compound-specific AI required."
        ]

    return ICHM7Result(
        class_number=class_num,
        class_name_ko=info["name_ko"] + (" — AI Inferred" if ai_inf else ""),
        definition=info["definition"],
        control_strategy=control_text,
        rationale=rationale,
        alerts=alerts,
        shared_alerts_with_drug=shared,
        is_preliminary=is_preliminary,
        alert_present=alert_present,
        empirical_db_match=empirical.source_label,
        ames_evidence=empirical.ames_result if empirical.matched else "Not in empirical DB",
        carcinogenicity_evidence=empirical.carcinogenicity_result
        if empirical.matched
        else "Not in empirical DB",
        recommended_ttc=ttc_text,
        consensus_sources=list(empirical.sources) if empirical.sources else [],
        empirical_compound_name=empirical.compound_name,
        gnn_mutagenic_probability=gnn_p,
        gnn_model_name=gnn_model,
        gnn_fallback_log=gnn_fb,
        ai_inferred_liability=ai_inf,
        assigned_class_label=assigned_label,
        rule_based_class=rule_cls,
        ensemble_log=ens_log,
        cohort_of_concern=coc,
        expert_review_required=expert_review_required,
        expert_review_reason=expert_review_reason,
    )


def scan_alerts_for_smiles(
    smiles: str,
    *,
    mol: Chem.Mol | None = None,
    precomputed: list[str] | None = None,
) -> tuple[bool, list[str]]:
    """Operation A — returns (alert_present, matched_alerts)."""
    if precomputed is not None:
        matched = list(precomputed)
        return bool(matched), matched
    try:
        m = mol or Chem.MolFromSmiles((smiles or "").strip())
        if m is None:
            return False, []
        matched = _scan_core_alerts(m)
        return bool(matched), matched
    except Exception:
        return False, []


def _hit_from_db_query(
    canonical_smiles: str,
    is_in_db: bool,
    ames_result: str | None,
    carc_result: str | None,
    meta: dict[str, str],
    *,
    latency_ms: float = 0.0,
) -> EmpiricalM7Hit:
    if not is_in_db:
        return EmpiricalM7Hit(
            matched=False,
            canonical_smiles=canonical_smiles,
            match_type="none",
            latency_ms=latency_ms,
        )
    source = meta.get("source", "")
    sources = tuple(source.split("/")) if source else ()
    return EmpiricalM7Hit(
        matched=True,
        canonical_smiles=canonical_smiles,
        compound_name=meta.get("compound_name", ""),
        sources=sources,
        ames_result=ames_result or "unknown",
        carcinogenicity_result=carc_result or "unknown",
        match_type="exact",
        latency_ms=latency_ms,
    )


CARCINOGEN_AMES_NEGATIVE_RATIONALE = (
    "Carcinogenic but Ames-negative — non-mutagenic mechanism; "
    "Q3A/Q3B limits with expert justification."
)


@dataclass
class M7RoutingDecision:
    """2026-09-03 추가: 3-튜플로는 실을 수 없던 전문가 검토 플래그 전달용."""

    class_number: int
    rationale: str
    preliminary: bool
    expert_review_required: bool = False
    expert_review_reason: str = ""


def route_m7_consensus_class_detailed(
    *,
    alert_present: bool,
    alerts: list[str],
    is_in_db: bool,
    ames_result: str | None,
    carc_result: str | None,
    empirical: EmpiricalM7Hit,
    drug_empirical: EmpiricalM7Hit | None,
    shared_alerts: list[str],
    drug_has_negative_mutagenicity_data: bool = False,
) -> M7RoutingDecision:
    """
    Deterministic ICH M7 routing — mirrors evaluate_ich_m7_class pseudocode
    with PRIME-DS parent-API Class 4 extension.

    route_m7_consensus_class()의 상세판. 전문가 검토 필요 여부까지 반환합니다.
    """
    src = empirical.source_label

    # 2026-09-03 수정: Class 1의 정의는 *변이원성* 발암물질입니다. 발암성 양성만으로
    # Class 1을 주면 비유전독성 발암물질까지 Class 1/2 lock에 걸려 후속 보정이
    # 차단되므로, Ames 양성을 함께 요구하도록 조건을 강화했습니다.
    if is_in_db and carc_result == "positive" and ames_result == "positive":
        return M7RoutingDecision(
            1,
            "Known mutagenic carcinogen based on empirical ISSCAN/EPA data "
            "(Ames positive + carcinogenicity positive).",
            False,
        )

    # 2026-09-03 수정: 발암성 양성 + Ames 음성 → 비변이원성 기전. Class 1이 아니라
    # Class 5 + 전문가 검토 플래그로 분기합니다.
    if is_in_db and carc_result == "positive" and ames_result == "negative":
        return M7RoutingDecision(
            5,
            CARCINOGEN_AMES_NEGATIVE_RATIONALE,
            True,
            expert_review_required=True,
            expert_review_reason=(
                f"Empirical DB ({src}) reports positive carcinogenicity with negative "
                "bacterial mutagenicity — non-DNA-reactive mechanism outside the scope "
                "of the ICH M7 mutagenic-impurity classification; qualified expert "
                "assessment required."
            ),
        )

    # 2026-09-03 수정: 발암성 양성인데 Ames 자료가 없으면 변이원 기전을 확정할 수
    # 없으므로 Class 1이 아니라 보수적으로 Class 3(자료 없음 + Ames 시험 권고)입니다.
    if is_in_db and carc_result == "positive":
        return M7RoutingDecision(
            3,
            f"Carcinogenicity positive in empirical DB ({src}) but bacterial "
            f"mutagenicity data are unavailable (Ames {ames_result}); the mutagenic "
            "mechanism required for Class 1 is unconfirmed. Ames testing required.",
            True,
            expert_review_required=True,
            expert_review_reason=(
                "Positive carcinogenicity without bacterial mutagenicity data — "
                "expert review and confirmatory Ames testing required before "
                "assigning a final ICH M7 class."
            ),
        )

    if is_in_db and ames_result == "positive":
        return M7RoutingDecision(
            2,
            "Known mutagen with unknown carcinogenic potential (Hansen DB Positive).",
            False,
        )

    if alert_present:
        drug_ames_neg = bool(
            drug_has_negative_mutagenicity_data
            or (
                drug_empirical
                and drug_empirical.matched
                and drug_empirical.ames_result == "negative"
            )
        )
        if shared_alerts and drug_ames_neg:
            return M7RoutingDecision(
                4,
                f"Shared alert with parent API ({', '.join(shared_alerts)}); "
                "parent empirical Ames negative — Q3A/Q3B limits apply.",
                False,
            )
        if is_in_db and ames_result == "negative":
            return M7RoutingDecision(
                4,
                "Alerting structure proven to be non-mutagenic by empirical testing.",
                False,
            )
        if not is_in_db:
            return M7RoutingDecision(
                3,
                "Alerting structure with unknown mutagenic potential. Requires qualification.",
                True,
            )
        return M7RoutingDecision(
            3,
            f"Alerting structure ({', '.join(alerts)}); empirical DB ({src}) "
            f"inconclusive (Ames {ames_result}, Carc {carc_result}).",
            True,
        )

    return M7RoutingDecision(
        5,
        "No structural alerts identified, or verified non-mutagenic.",
        not (is_in_db and ames_result == "negative"),
    )


def route_m7_consensus_class(
    *,
    alert_present: bool,
    alerts: list[str],
    is_in_db: bool,
    ames_result: str | None,
    carc_result: str | None,
    empirical: EmpiricalM7Hit,
    drug_empirical: EmpiricalM7Hit | None,
    shared_alerts: list[str],
    drug_has_negative_mutagenicity_data: bool = False,
) -> tuple[int, str, bool]:
    """하위 호환 래퍼 — (class_number, rationale, is_preliminary) 3-튜플 반환."""
    d = route_m7_consensus_class_detailed(
        alert_present=alert_present,
        alerts=alerts,
        is_in_db=is_in_db,
        ames_result=ames_result,
        carc_result=carc_result,
        empirical=empirical,
        drug_empirical=drug_empirical,
        shared_alerts=shared_alerts,
        drug_has_negative_mutagenicity_data=drug_has_negative_mutagenicity_data,
    )
    return d.class_number, d.rationale, d.preliminary


def evaluate_ich_m7_class(
    smiles: str,
    db_connection: sqlite3.Connection | None = None,
    *,
    drug_smiles: str | None = None,
    alerts: list[str] | None = None,
) -> dict[str, Any]:
    """
    Hybrid ensemble API — Layers A (DB) + B (SMARTS) + C (AttentiveFP GNN)
    with automated consensus rescoring.
    """
    from p65_m7_ensemble_matrix import run_m7_ensemble_matrix

    t0 = time.perf_counter()
    ensemble = run_m7_ensemble_matrix(
        smiles,
        db_connection,
        drug_smiles=drug_smiles,
        alerts=alerts,
    )

    try:
        mol = Chem.MolFromSmiles((smiles or "").strip())
        canon_smiles = Chem.MolToSmiles(mol, canonical=True) if mol else ""
    except Exception:
        canon_smiles = ""

    result = _build_result(
        ensemble.final_class,
        rationale=ensemble.description,
        alerts=ensemble.layer_b.combined_alerts,
        shared=ensemble.shared_alerts,
        is_preliminary=ensemble.preliminary,
        empirical=ensemble.layer_a.empirical,
        alert_present=ensemble.layer_b.alert_present,
        ensemble=ensemble,
    )

    out = result.to_outcome_dict()
    out.update(
        {
            "class": ensemble.class_label,
            # 2026-09-03 수정: CoC 근거가 붙은 rationale을 description으로 노출.
            "description": result.rationale,
            "alerts": ensemble.layer_b.combined_alerts,
            "Canonical SMILES": canon_smiles,
            "alert_present": ensemble.layer_b.alert_present,
            "is_in_db": ensemble.layer_a.is_in_db,
            "ames_result": ensemble.layer_a.ames_result,
            "carcinogenicity_result": ensemble.layer_a.carc_result,
            "Empirical Lookup ms": f"{ensemble.layer_a.latency_ms:.2f}",
            "Ensemble Total ms": f"{(time.perf_counter() - t0) * 1000:.2f}",
        }
    )
    out.update(ensemble.to_dict())
    # 2026-09-03 수정: ensemble.to_dict()가 뒤에 오므로 CoC·전문가 검토 판정이
    # 덮이지 않도록 마지막에 다시 확정합니다.
    out["Cohort of Concern (ICH M7 §7.5)"] = result.cohort_of_concern
    out["Recommended Action Threshold (TTC Guidelines)"] = result.recommended_ttc
    return out


def classify_ich_m7(
    mol: Chem.Mol | None,
    alerts: list[str] | None = None,
    drug_mol: Chem.Mol | None = None,
    *,
    drug_has_negative_mutagenicity_data: bool = False,
) -> ICHM7Result:
    """
    Backward-compatible classifier — hybrid ensemble matrix (A+B+C).
    """
    if mol is None:
        info = ICH_M7_CLASS_INFO[5]
        return ICHM7Result(
            class_number=5,
            class_name_ko=info["name_ko"],
            definition=info["definition"],
            control_strategy=info["control"],
            rationale="유효하지 않은 구조 — M7 분류 불가",
            is_preliminary=True,
            recommended_ttc=m7_ttc_threshold(5),  # 2026-09-03 수정: dict 조회 → 함수 호출
        )

    try:
        smi = Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        info = ICH_M7_CLASS_INFO[5]
        return ICHM7Result(
            class_number=5,
            class_name_ko=info["name_ko"],
            definition=info["definition"],
            control_strategy=info["control"],
            rationale="SMILES generation failed",
            is_preliminary=True,
            recommended_ttc=m7_ttc_threshold(5),  # 2026-09-03 수정: dict 조회 → 함수 호출
        )

    drug_smi = None
    if drug_mol is not None:
        try:
            drug_smi = Chem.MolToSmiles(drug_mol, canonical=True)
        except Exception:
            drug_smi = None

    from p65_m7_ensemble_matrix import run_m7_ensemble_matrix

    ensemble = run_m7_ensemble_matrix(
        smi,
        drug_smiles=drug_smi,
        alerts=alerts,
        drug_has_negative_mutagenicity_data=drug_has_negative_mutagenicity_data,
    )

    return _build_result(
        ensemble.final_class,
        rationale=ensemble.description,
        alerts=ensemble.layer_b.combined_alerts,
        shared=ensemble.shared_alerts,
        is_preliminary=ensemble.preliminary,
        empirical=ensemble.layer_a.empirical,
        alert_present=ensemble.layer_b.alert_present,
        ensemble=ensemble,
    )


def ich_m7_html_badge(result: ICHM7Result) -> str:
    """Streamlit unsafe_allow_html용 M7 배지."""
    pre = " (추정)" if result.is_preliminary else ""
    emp = (
        f'<br><span style="color:#546e7a;font-size:0.76rem;">'
        f"Empirical: {result.empirical_db_match} · Ames {result.ames_evidence} · "
        f"Carc {result.carcinogenicity_evidence}</span>"
        if result.empirical_db_match != "None"
        else ""
    )
    gnn = ""
    if result.gnn_mutagenic_probability is not None:
        gnn = (
            f'<br><span style="color:#546e7a;font-size:0.76rem;">'
            f"GNN P(Mutagenic)={result.gnn_mutagenic_probability:.2f} · {result.gnn_model_name}</span>"
        )
    ai = ""
    if result.ai_inferred_liability:
        ai = (
            '<br><span style="color:#e65100;font-size:0.76rem;">'
            "⚠ Class 3-Alert (AI Inferred Liability)</span>"
        )
    # 2026-09-03 수정: CoC는 화면에서도 TTC가 아니라 AI가 필요함을 알려야 합니다.
    coc = ""
    if result.cohort_of_concern:
        coc = (
            '<br><span style="color:#b71c1c;font-size:0.78rem;">'
            "⚠ Cohort of concern — TTC 적용 불가, 화합물별 AI 필요 (ICH M7(R2) §7.5)</span>"
        )
    label = result.assigned_class_label or f"Class {result.class_number}"
    return (
        f'<div style="background:{result.badge_color}15;border-left:4px solid {result.badge_color};'
        f'padding:0.55rem 0.75rem;margin:0.35rem 0;border-radius:4px;font-size:0.84rem;">'
        f'<strong style="color:{result.badge_color};">ICH M7 {result.class_name_ko} ({label}){pre}</strong><br>'
        f'<span style="color:#37474f;">{result.definition}</span><br>'
        f'<span style="color:#546e7a;font-size:0.8rem;">관리: {result.control_strategy}</span><br>'
        f'<span style="color:#78909c;font-size:0.78rem;">근거: {result.rationale}</span>'
        f"{emp}{gnn}{ai}{coc}"
        f"</div>"
    )
