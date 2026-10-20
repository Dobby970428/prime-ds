"""
PRIME-DS Phase 7 — ICH M7 Class 1~5 자동 등급 판정 코어.

Structural Alert 스크리닝 + DeepChem/AttentiveFP Ames QSAR 앙상블 입력을 받아
공정 유연물질 후보에 ich_m7_class · regulatory_action_text 를 바인딩한다.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from rdkit import Chem

from p66_ich_m7 import ICH_M7_CLASS_INFO, _shared_alerts, scan_alerts_for_smiles
from p60_m7_empirical_db import EmpiricalM7Hit, lookup_empirical_m7

MUTAGENICITY_POSITIVE_THRESHOLD = 0.5

REGULATORY_ACTION_TEXT: dict[int, str] = {
    1: (
        "Class 1 — 알려진 변이원성 발암물질로 분류됩니다. "
        "TTC 미적용 · 화합물별 Acceptable Intake(AI) 수준 관리 및 공정 제거/대체가 필수입니다."
    ),
    2: (
        "Class 2 — 알려진 변이원성(발암성 미확인)으로 분류됩니다. "
        "TTC 1.5 µg/일(평생 노출) 기준 적용 · 공정 개선 및 Ames 확인 시험을 권고합니다."
    ),
    3: (
        "Class 3 — 구조적 알럿(Toxicophore) 감지 · 유전독성 리스크가 있습니다. "
        "TTC 1.5 µg/일 관리 · in silico/Ames 확인 및 공정 개선을 권고합니다."
    ),
    4: (
        "Class 4로 분류되어 TTC 기준 관리가 해제됩니다. "
        "구조적 알럿이 있으나 QSAR/ML Ames 앙상블 최종 음성 — ICH Q3A/Q3B 일반 불순물 한도 적용."
    ),
    5: (
        "Class 5 — 구조적 알럿 없음 · QSAR 음성. "
        "비유전독성 불순물로 ICH Q3A/Q3B 일반 한도 내 관리 가능합니다."
    ),
}

BADGE_COLORS: dict[int, str] = {
    1: "#b71c1c",
    2: "#e65100",
    3: "#f57c00",
    4: "#1565c0",
    5: "#2e7d32",
}


@dataclass
class ICHM7ClassificationResult:
    ich_m7_class: int
    regulatory_action_text: str
    alert_present: bool
    structural_alerts: list[str] = field(default_factory=list)
    ames_mutagenicity_prob: float | None = None
    ames_prediction_label: str = "Not evaluated"
    empirical_db_match: bool = False
    ames_evidence: str = "Not in empirical DB"
    carcinogenicity_evidence: str = "Not in empirical DB"
    rationale: str = ""

    @property
    def class_display(self) -> str:
        info = ICH_M7_CLASS_INFO.get(self.ich_m7_class, ICH_M7_CLASS_INFO[5])
        return f"Class {self.ich_m7_class} — {info['name_ko']}"

    @property
    def ich_m7_label(self) -> str:
        return f"ICH M7 Class {self.ich_m7_class}"

    @property
    def ich_m7_detail(self) -> str:
        if self.structural_alerts:
            return f"Alerts: {', '.join(self.structural_alerts)} · {self.ames_prediction_label}"
        return self.ames_prediction_label


class ICHM7Classifier:
    """
    ICH M7(R2) consensus Class 1~5 판정기.

    입력: Structural Alert 유무, empirical DB 매칭, QSAR/AI Ames 확률.
    """

    def __init__(self, *, mutagenicity_threshold: float = MUTAGENICITY_POSITIVE_THRESHOLD) -> None:
        self.mutagenicity_threshold = mutagenicity_threshold

    def classify(
        self,
        *,
        has_structural_alert: bool,
        structural_alerts: list[str] | None = None,
        empirical_hit: EmpiricalM7Hit | None = None,
        ames_mutagenicity_prob: float | None = None,
        parent_shared_alerts: list[str] | None = None,
        parent_ames_negative: bool = False,
    ) -> ICHM7ClassificationResult:
        """ICH M7(R2) Class 1~5 판정.

        2026-09-08 수정 — Class 4 판정 근거를 규정에 맞게 고쳤습니다.
        --------------------------------------------------------------
        종전 구현은 «구조 알럿이 있고 통계 모델이 음성이면 Class 4» 였습니다.
        그러나 ICH M7(R2) 에서 Class 4 는 «해당 알럿이 **원료의약품(모 화합물)과
        공유되고**, 그 모 화합물에 변이원성 음성 실측 자료가 있는 경우» 입니다.
        통계 모델의 음성 예측은 그 자체로 알럿을 해소하지 못하며, 전문가 규칙과
        통계 모델의 판정이 상충하는 상태이므로 Class 3 을 유지하고 전문가 검토
        대상으로 남기는 것이 규정에 부합합니다.

        따라서 이 판정에는 **생성 경로 정보(모 화합물 구조)** 가 필요하며,
        이는 후보 유연물질 생성부가 부모–자식 관계를 함께 전달할 때에만
        존재합니다.
        """
        alerts = list(structural_alerts or [])
        shared = list(parent_shared_alerts or [])
        hit = empirical_hit or EmpiricalM7Hit()
        qsar_label = self._ames_label(ames_mutagenicity_prob)

        if hit.matched and hit.carcinogenicity_result == "positive":
            return self._result(
                1,
                has_structural_alert=has_structural_alert or bool(alerts),
                alerts=alerts,
                hit=hit,
                ames_prob=ames_mutagenicity_prob,
                qsar_label=qsar_label,
                rationale="외부 DB 매칭 — 알려진 발암성·변이원성 물질",
            )

        if hit.matched and hit.ames_result == "positive":
            return self._result(
                2,
                has_structural_alert=has_structural_alert or bool(alerts),
                alerts=alerts,
                hit=hit,
                ames_prob=ames_mutagenicity_prob,
                qsar_label=qsar_label,
                rationale="외부 DB 매칭 — 알려진 변이원성(발암성 미확인)",
            )

        if has_structural_alert or alerts:
            # ---- Class 4 — 모 화합물과 알럿을 공유하고 모 화합물이 실측 음성 ----
            if shared and parent_ames_negative:
                return self._result(
                    4,
                    has_structural_alert=True,
                    alerts=alerts,
                    hit=hit,
                    ames_prob=ames_mutagenicity_prob,
                    qsar_label=qsar_label,
                    rationale=(
                        f"구조적 알럿({', '.join(shared)})이 원료의약품과 공유되며 "
                        "해당 원료의약품에 변이원성 음성 실측 자료가 존재 — "
                        "ICH M7 Class 4"
                    ),
                )
            # ---- Class 4 — 후보 자체가 실측 DB 에 있고 Ames 음성 ----
            if hit.matched and hit.ames_result == "negative":
                return self._result(
                    4,
                    has_structural_alert=True,
                    alerts=alerts,
                    hit=hit,
                    ames_prob=ames_mutagenicity_prob,
                    qsar_label=qsar_label,
                    rationale=(
                        f"구조적 알럿({', '.join(alerts) or 'detected'}) 존재하나 "
                        "해당 구조 자체에 변이원성 음성 실측 자료가 존재 — ICH M7 Class 4"
                    ),
                )
            # ---- 알럿 양성 · 통계 음성 = 두 방법론 상충 → Class 3 유지 + 전문가 검토 ----
            if self._is_qsar_negative(ames_mutagenicity_prob):
                return self._result(
                    3,
                    has_structural_alert=True,
                    alerts=alerts,
                    hit=hit,
                    ames_prob=ames_mutagenicity_prob,
                    qsar_label=qsar_label,
                    rationale=(
                        f"구조적 알럿({', '.join(alerts) or 'detected'}) 존재 · "
                        f"통계 (Q)SAR 음성(P={ames_mutagenicity_prob:.3f}) — "
                        "전문가 규칙과 통계 모델의 판정이 상충하므로 Class 3 유지, "
                        "ICH M7(R2)에 따라 전문가 검토 대상"
                    ),
                )
            return self._result(
                3,
                has_structural_alert=True,
                alerts=alerts,
                hit=hit,
                ames_prob=ames_mutagenicity_prob,
                qsar_label=qsar_label,
                rationale=(
                    f"구조적 알럿({', '.join(alerts) or 'detected'}) — "
                    "유전독성 리스크 평가 필요"
                ),
            )

        if self._is_qsar_positive(ames_mutagenicity_prob):
            return self._result(
                3,
                has_structural_alert=False,
                alerts=alerts,
                hit=hit,
                ames_prob=ames_mutagenicity_prob,
                qsar_label=qsar_label,
                rationale=f"알럿 없음 · QSAR Ames 양성 예측(P={ames_mutagenicity_prob:.3f})",
            )

        return self._result(
            5,
            has_structural_alert=False,
            alerts=alerts,
            hit=hit,
            ames_prob=ames_mutagenicity_prob,
            qsar_label=qsar_label,
            rationale="구조적 알럿 없음 · QSAR 음성 또는 평가 불가(보수적 Class 5)",
        )

    def classify_smiles(
        self,
        smiles: str,
        *,
        mol: Chem.Mol | None = None,
        run_qsar: bool = True,
        parent_smiles: str | None = None,
    ) -> ICHM7ClassificationResult:
        """Structural Alert + empirical DB + DeepChem Ames 파이프라인 일괄 실행."""
        smi = (smiles or "").strip()
        m = mol
        if m is None and smi:
            m = Chem.MolFromSmiles(smi)

        alert_present, alerts = scan_alerts_for_smiles(smi, mol=m)
        empirical = lookup_empirical_m7(smi) if smi else EmpiricalM7Hit()

        # ---- 생성 경로 정보(모 화합물)로부터 알럿 공유 여부와 모 화합물 실측을 확보 ----
        shared: list[str] = []
        parent_ames_negative = False
        parent = (parent_smiles or "").strip()
        if parent and m is not None:
            try:
                parent_mol = Chem.MolFromSmiles(parent)
                if parent_mol is not None:
                    shared = _shared_alerts(m, parent_mol)
                    parent_hit = lookup_empirical_m7(parent)
                    parent_ames_negative = bool(
                        parent_hit.matched and parent_hit.ames_result == "negative"
                    )
            except Exception:  # noqa: BLE001
                shared = []
                parent_ames_negative = False

        ames_prob: float | None = None
        if run_qsar and smi:
            try:
                from p64_gnn_mutagenicity_predictor import predict_mutagenicity_probability

                gnn = predict_mutagenicity_probability(smi)
                ames_prob = gnn.p_mutagenic
            except Exception:
                ames_prob = None

        return self.classify(
            has_structural_alert=alert_present,
            structural_alerts=alerts,
            empirical_hit=empirical,
            ames_mutagenicity_prob=ames_prob,
            parent_shared_alerts=shared,
            parent_ames_negative=parent_ames_negative,
        )

    def _is_qsar_negative(self, prob: float | None) -> bool:
        return prob is not None and prob < self.mutagenicity_threshold

    def _is_qsar_positive(self, prob: float | None) -> bool:
        return prob is not None and prob >= self.mutagenicity_threshold

    def _ames_label(self, prob: float | None) -> str:
        if prob is None:
            return "QSAR not evaluated"
        if prob >= self.mutagenicity_threshold:
            return f"QSAR Positive (P={prob:.3f})"
        return f"QSAR Negative (P={prob:.3f})"

    def _result(
        self,
        class_num: int,
        *,
        has_structural_alert: bool,
        alerts: list[str],
        hit: EmpiricalM7Hit,
        ames_prob: float | None,
        qsar_label: str,
        rationale: str,
    ) -> ICHM7ClassificationResult:
        return ICHM7ClassificationResult(
            ich_m7_class=class_num,
            regulatory_action_text=REGULATORY_ACTION_TEXT.get(class_num, ""),
            alert_present=has_structural_alert,
            structural_alerts=alerts,
            ames_mutagenicity_prob=ames_prob,
            ames_prediction_label=qsar_label,
            empirical_db_match=hit.matched,
            ames_evidence=hit.ames_result if hit.matched else "Not in empirical DB",
            carcinogenicity_evidence=(
                hit.carcinogenicity_result if hit.matched else "Not in empirical DB"
            ),
            rationale=rationale,
        )


_DEFAULT_CLASSIFIER = ICHM7Classifier()


def classify_process_impurity(
    smiles: str,
    *,
    mol: Chem.Mol | None = None,
    run_qsar: bool = True,
    parent_smiles: str | None = None,
) -> ICHM7ClassificationResult:
    """공정 유연물질 단일 SMILES — ICH M7 Class 1~5 판정.

    `parent_smiles` 는 이 유연물질이 유래한 원료의약품(모 화합물)의 구조입니다.
    Class 4 판정에 필수이며, 생성부가 부모–자식 관계를 전달할 때에만 존재합니다.
    """
    return _DEFAULT_CLASSIFIER.classify_smiles(
        smiles, mol=mol, run_qsar=run_qsar, parent_smiles=parent_smiles
    )


def bind_ich_m7_to_record(
    record: dict[str, Any],
    classification: ICHM7ClassificationResult,
) -> dict[str, Any]:
    """예측 레코드에 ICH M7 Phase 7 필드 바인딩."""
    record = dict(record)
    record["ich_m7_class"] = classification.ich_m7_class
    record["regulatory_action_text"] = classification.regulatory_action_text
    record["ich_m7"] = classification.ich_m7_label
    record["ich_m7_detail"] = classification.ich_m7_detail
    record["structural_alerts"] = classification.structural_alerts
    record["ames_mutagenicity_prob"] = classification.ames_mutagenicity_prob
    record["ames_prediction_label"] = classification.ames_prediction_label
    record["alert_present"] = classification.alert_present
    return record


# 레코드에서 모 화합물 구조를 찾을 때 순서대로 살펴보는 키
PARENT_SMILES_KEYS: tuple[str, ...] = (
    "parent_smiles",
    "api_smiles",
    "drug_smiles",
    "source_smiles",
    "precursor_smiles",
)


def resolve_parent_smiles(
    record: dict[str, Any], default_parent: str | None = None
) -> str | None:
    """레코드에 실린 생성 경로 정보에서 모 화합물 구조를 꺼냅니다.

    레코드가 자체 부모를 가지고 있으면(다세대 분해 트리의 중간 노드 등)
    그것을 우선하고, 없으면 배치 전체의 기본 모 화합물을 씁니다.
    """
    for key in PARENT_SMILES_KEYS:
        val = str(record.get(key) or "").strip()
        if val:
            return val
    return (default_parent or "").strip() or None


def enrich_impurity_records(
    records: list[dict[str, Any]],
    *,
    run_qsar: bool = True,
    parent_smiles: str | None = None,
) -> list[dict[str, Any]]:
    """공정 시뮬레이션 batch 결과에 ICH M7 Class 1~5 일괄 부여.

    2026-09-08 — `parent_smiles`(모 화합물 구조)를 함께 받습니다. 각 레코드가
    자체 부모 구조를 실어 보내면 그것을 우선합니다. 이 정보가 없으면 Class 4
    (모 화합물과 알럿을 공유하는 유연물질) 판정 분기가 열리지 않습니다.
    """
    enriched: list[dict[str, Any]] = []
    for rec in records:
        smi = str(rec.get("smiles") or "").strip()
        if not smi:
            enriched.append(rec)
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            enriched.append(rec)
            continue
        parent = resolve_parent_smiles(rec, parent_smiles)
        classification = classify_process_impurity(
            smi, mol=mol, run_qsar=run_qsar, parent_smiles=parent
        )
        out = bind_ich_m7_to_record(rec, classification)
        out["parent_smiles_used"] = parent or ""
        enriched.append(out)
    return enriched


def ich_m7_class_badge_html(class_num: int, *, detail: str = "") -> str:
    """Streamlit 테이블용 ICH M7 등급 뱃지 (Class 1~3 빨강/주황, 4 청색, 5 녹색)."""
    color = BADGE_COLORS.get(class_num, "#546e7a")
    risk = "High Risk" if class_num <= 3 else ("De-risked by QSAR" if class_num == 4 else "Safe")
    info = ICH_M7_CLASS_INFO.get(class_num, ICH_M7_CLASS_INFO[5])
    extra = f'<br><span style="color:#546e7a;font-size:0.75rem;">{detail}</span>' if detail else ""
    return (
        f'<span style="display:inline-block;background:{color}18;color:{color};'
        f"border:1px solid {color};border-radius:6px;padding:0.2rem 0.55rem;"
        f'font-weight:700;font-size:0.82rem;">'
        f"Class {class_num} · {risk}</span>{extra}"
    )
