"""ICH M7(R2) regulatory-style IN-SILICO TOXICITY ASSESSMENT REPORT (English)."""

# p80_formal_assessment_report - 본문 생성기

# ProductAssessment: 개별 화합물 평가 결과
# AssessmentContext: 보고서 전체 정보(물질명, 구조, 등급, 물성,,,)
# M7RegulatoryOptionSummary: ICH M7 판정 요약

# 보고서 목차(영문 규제 양식)
# DOCUMENT CONTROL / CONFIDENTIALITY
# 1. SUBSTANCE IDENTIFICATION
# 2. PROCESS ROUTE AND IMPURITY FORMATION
#    2.1 DRUG SUBSTANCE (FINAL PRODUCT)
# 3. COMPLEMENTARY HYBRID ASSESSMENT
#    3.1 Expert Rule-Based
#    3.4 ICH M7(R2) Regulatory Option Determination
# 4. PHYSICOCHEMICAL PROPERTIES AND EXPOSURE
# 5. FINAL ICH M7(R2) HAZARD CLASSIFICATION
# 6. RECOMMENDED ACTIONS AND ACCEPTABLE INTAKE
# 7. CONCLUSION AND LIMITATIONS
# APPENDIX A — ALL PREDICTED SPECIES
# ⚠️ 3.4절의 "Option 3/4/5" 표기가 ICH M7의 실제 정의와 다릅니다. 가이드라인의 Option 1~4는 평가 방법이 아니라 관리 전략이며 Option 5는 존재하지 않습니다. 수정이 필요한 부분입니다.

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime

from p66_ich_m7 import ICHM7Result

# 2026-09-03 수정: ICH M7 등급 라벨·정의·조치를 p16(단일 출처)에서 가져옵니다.
from p16_ich_m7_terms import class_action, class_definition, class_label
from p83_m7_executive_visuals import build_m7_executive_rows, format_m7_executive_text_block

PLATFORM_NAME = "PRIME-DS"
PLATFORM_FULL = (
    "Platform for Reaction-based Impurity and Mutagenicity Evaluation for Drug Substances"
)
DOCUMENT_VERSION = "1.0"
REPORT_TITLE_LINE1 = "Prediction of Reaction & Mutagenic Impurity"
REPORT_TITLE_LINE2 = "Assessment Report"
REPORT_TITLE = f"{REPORT_TITLE_LINE1} — {REPORT_TITLE_LINE2}"

# 2026-09-03 수정: 13040을 보고서 네 곳에 하드코딩해 두어 drug_index.db가 갱신되면
# 거짓 숫자가 인쇄됐습니다. p20에서 실제 카운트를 조회하고 실패 시에만 이 상수를
# fallback으로 쓰되, fallback이라는 사실을 값(라벨)에 표시합니다. (p20은 읽기 전용)
DEFAULT_INDEXED_DB_COUNT = 13040

# 2026-09-03 수정: p16의 복제본이던 등급 표를 p16 헬퍼로부터 생성합니다.
# name/definition은 ICH M7(R2) Table 1 원문(p16), 아래 dict는 보고서용 부연 설명.
_M7_DEFINITION_DETAIL: dict[int, str] = {
    1: (
        "Positive bacterial mutagenicity (Ames) and positive carcinogenicity data "
        "in experimental animals or humans."
    ),
    2: "Positive in the Ames test; carcinogenicity data are not available.",
    3: (
        "Structural alert(s) for mutagenicity identified; no Ames or carcinogenicity "
        "data available for the impurity."
    ),
    4: (
        "Alerting structure present in the impurity and drug substance (or metabolite); "
        "mutagenicity/carcinogenicity concern for the impurity is considered addressed "
        "by negative data for the drug substance."
    ),
    5: (
        "No structural alerts for mutagenicity, or adequate experimental data demonstrate "
        "absence of mutagenic potential."
    ),
}

_M7_CONTROL_DETAIL: dict[int, str] = {
    1: (
        "Compound-specific Acceptable Intake (AI) applies; Threshold of Toxicological "
        "Concern (TTC) is not applicable."
    ),
    2: (
        "Control at or below the TTC of 1.5 µg/day (lifetime exposure) or staged TTC "
        "as justified per ICH M7(R2)."
    ),
    3: (
        "Control at or below the TTC of 1.5 µg/day (lifetime exposure); "
        "confirmatory in silico and/or bacterial mutagenicity testing is recommended."
    ),
    4: "Manage as a non-mutagenic impurity per ICH Q3A/Q3B qualification limits.",
    5: "Manage per ICH Q3A/Q3B general impurity limits.",
}

ICH_M7_REGULATORY: dict[int, dict[str, str]] = {
    n: {
        "name": class_label(n),
        "definition": f"{class_definition(n)}. {_M7_DEFINITION_DETAIL[n]}",
        "control": f"{class_action(n)}. {_M7_CONTROL_DETAIL[n]}",
    }
    for n in (1, 2, 3, 4, 5)
}


@dataclass
class ProductAssessment:
    """Final product (API) — structure and hybrid toxicity evaluation."""

    smiles: str
    name: str
    mw: str
    alerts: list[str] = field(default_factory=list)
    m7: ICHM7Result | None = None
    pubchem_cid: str | None = None
    pubchem_status: str = ""
    logp: str = ""
    stability: str = ""
    stability_detail: str = ""
    expert_comment: str = ""
    db_match: str = ""


@dataclass
class AssessmentContext:
    report_id: str
    issue_date: str
    is_reactant: bool
    is_process_impurity: bool
    name: str
    smiles: str
    mw: str
    is_standalone_compound: bool = False
    reactant1: str = ""
    reactant2: str = ""
    reagent: str = ""
    product: str = ""
    product_assessment: ProductAssessment | None = None
    is_final_product: bool = False
    reference_structure: str = ""
    primary_impurity_name: str = ""
    primary_impurity_smiles: str = ""
    nitrite: bool = False
    acidic_condition: bool = False
    mechanism: str = ""
    predicted_impurity_smiles: str = ""
    expert_engine: str = "RDKit Substructure Search / Ashby-Tennant Expert Ruleset"
    alerts: list[str] = field(default_factory=list)
    expert_comment: str = ""
    # 2026-09-03 수정: DeepChem은 requirements에서 제거된 패키지이므로 언급 삭제.
    qsar_engine: str = "Scikit-Learn statistical (Q)SAR (Tox21/Ames) — integration pending"
    qsar_ames: str = "Not integrated; expert rule-based assessment applied per ICH M7(R2)"
    oecd_connection_status: str = "Not queried"
    oecd_hydrolysis: str = "N/A"
    oecd_biodegradability: str = "N/A"
    oecd_mutagenicity_alert: str = "N/A"
    pubchem_status: str = ""
    pubchem_cid: str | None = None
    db_match: str = ""
    db_ref: str = ""
    logp: str = ""
    logd_pka_bp: str = ""
    stability: str = ""
    stability_detail: str = ""
    exposure_comment: str = ""
    m7: ICHM7Result | None = None
    indexed_db_count: int = DEFAULT_INDEXED_DB_COUNT
    # 2026-09-03 수정: 위 건수가 실측인지 fallback 상수인지 보고서에 표시하기 위한 플래그.
    indexed_db_count_is_fallback: bool = True
    appendix_species: list[dict] = field(default_factory=list)


def resolve_indexed_db_count(
    fallback: int = DEFAULT_INDEXED_DB_COUNT,
) -> tuple[int, bool]:
    """(건수, is_fallback) — p20_db_index.count_indexed() 조회, 실패 시 상수 사용."""
    try:
        from p20_db_index import count_indexed

        n = int(count_indexed())
        if n > 0:
            return n, False
        return fallback, True
    except Exception:
        return fallback, True


def indexed_db_count_label(ctx: "AssessmentContext") -> str:
    """3.3절 제목에 넣을 건수 문구 — fallback이면 그 사실을 함께 적습니다."""
    n = ctx.indexed_db_count
    if getattr(ctx, "indexed_db_count_is_fallback", False):
        return f"{n:,} compounds — static reference count; live index query unavailable"
    return f"{n:,} compounds — live index count"


def new_report_id() -> str:
    return f"IN-SILICO-{datetime.now().strftime('%Y-%m%d')}"


def english_issue_date() -> str:
    return datetime.now().strftime("%d %B %Y")


def _chk(on: bool) -> str:
    return "[X]" if on else "[ ]"


def substance_category_label(ctx: AssessmentContext) -> str:
    if ctx.is_standalone_compound:
        return "Standalone compound (single-compound toxicity & property scan)"
    if ctx.is_final_product:
        return "Drug substance (final product): Yes | Reactant: No | Process impurity: No"
    return (
        f"Reactant: {'Yes' if ctx.is_reactant else 'No'} | "
        f"Process impurity: {'Yes' if ctx.is_process_impurity else 'No'}"
    )


def includes_process_route_section(ctx: AssessmentContext) -> bool:
    """화학공정 반응식·불순물 형성 경로 섹션 — 공정 모드에서만 포함."""
    if ctx.is_standalone_compound:
        return False
    if ctx.is_process_impurity or ctx.is_reactant or ctx.is_final_product:
        return True
    if ctx.product_assessment is not None:
        return True
    if any((ctx.reactant1, ctx.reactant2, ctx.reagent, ctx.product, ctx.predicted_impurity_smiles)):
        return True
    return False


def _translate_stability(label: str, detail: str) -> tuple[str, str]:
    """Map UI stability strings to regulatory English."""
    low = (label + " " + detail).lower()
    if "unstable" in low or "불안정" in label:
        en_label = "Unstable under typical aqueous process conditions"
        en_detail = detail
        if "가수분해" in detail or "hydrolysis" not in detail.lower():
            en_detail = (
                "Elevated risk of hydrolytic degradation under aqueous process conditions "
                "(ester or similarly labile functionality)."
            )
    elif "moderate" in low or "보통" in label:
        en_label = "Moderately stable"
        en_detail = (
            "Limited hydrolytic susceptibility under strongly acidic/basic or elevated "
            "temperature conditions cannot be excluded."
        )
    elif "stable" in low or "안정" in label:
        en_label = "Stable"
        en_detail = (
            "No major hydrolysis- or oxidation-prone functionalities identified; "
            "chemical behavior under typical process conditions is predicted to be stable."
        )
    else:
        en_label = label or "Not assessed"
        en_detail = detail or "Insufficient data for stability assessment."
    return en_label, en_detail


def _translate_mechanism(text: str) -> str:
    if "단일 화합물" in text:
        return "Standalone compound toxicity and physicochemical property assessment (no process route)."
    return text or "Not applicable — standalone assessment mode."


def m7_action_flags(
    class_number: int,
    *,
    cohort_of_concern: bool = False,
) -> tuple[bool, bool, str]:
    # 2026-09-03 수정: cohort of concern이면 Class와 무관하게 TTC를 적용할 수
    # 없습니다(ICH M7(R2) §7.5). 기존 호출부는 키워드 없이 그대로 동작합니다.
    if cohort_of_concern:
        return (
            True,
            False,
            (
                "Cohort of concern (aflatoxin-like, N-nitroso or alkyl-azoxy): the TTC is "
                "not applicable (ICH M7(R2) Section 7.5). Derive and apply a "
                "compound-specific Acceptable Intake (AI) for this impurity."
            ),
        )
    # 2026-09-03 수정: `class_number <= 3`이 Class 1까지 TTC로 보내는 바람에,
    # 같은 보고서 안에서 Class 1 정의("TTC is not applicable")와 정면으로
    # 모순되는 문장이 나란히 인쇄됐습니다. Class 1을 별도 분기로 분리합니다.
    if class_number == 1:
        return (
            True,
            False,
            (
                "Apply a compound-specific Acceptable Intake (AI) derived from "
                "carcinogenicity potency data for this impurity (ICH M7(R2) Section 7.2); "
                "the TTC is NOT applicable to Class 1 known mutagenic carcinogens."
            ),
        )
    if class_number <= 3:  # Class 2·3만 TTC 적용
        return (
            True,
            False,
            (
                "Apply the Threshold of Toxicological Concern (TTC): control the impurity "
                "to an acceptable intake of < 1.5 µg/day (lifetime exposure basis)."
            ),
        )
    if class_number == 4:
        return (
            False,
            True,
            (
                "Apply ICH Q3A/Q3B qualification limits for non-mutagenic impurities "
                "(Class 4 — impurity non-mutagenicity justified by drug substance data)."
            ),
        )
    return (
        False,
        True,
        "Apply ICH Q3A/Q3B general impurity limits (Class 5).",
    )


def _expert_rationale_en(m7: ICHM7Result | None, alerts: list[str]) -> str:
    if m7 is None:
        return "Classification could not be established."
    # 2026-09-03 수정: 미등록 Class로 KeyError가 나던 직접 인덱싱을 .get()으로 교체.
    info = ICH_M7_REGULATORY.get(m7.class_number, ICH_M7_REGULATORY[5])
    alert_part = f" Structural alerts: {', '.join(alerts)}." if alerts else ""
    emp_part = ""
    if getattr(m7, "empirical_db_match", "None") not in ("None", "", "Not evaluated"):
        emp_part = (
            f" Empirical DB ({m7.empirical_db_match}): Ames {getattr(m7, 'ames_evidence', 'N/A')}, "
            f"carcinogenicity {getattr(m7, 'carcinogenicity_evidence', 'N/A')}."
        )
    ttc = getattr(m7, "recommended_ttc", "") or info.get("control", "")
    prelim = "Preliminary in silico conclusion" if m7.is_preliminary else "ICH M7 consensus classification"
    return (
        f"{info['definition']}{alert_part}{emp_part} "
        f"Recommended control: {ttc}. "
        f"{prelim} — confirm with literature and experimental data."
    )


def hybrid_ensemble_report_rows(m7: ICHM7Result | None) -> list[tuple[str, str]]:
    """Section 3.2b — Hybrid Ensemble Matrix (Layers A/B/C) for PDF/DOCX/TXT."""
    if m7 is None:
        return [
            ("Ensemble engine", "ICH M7 Hybrid Ensemble (DB + SMARTS + GNN)"),
            ("Layer C P(Mutagenic)", "Not evaluated"),
        ]
    rows: list[tuple[str, str]] = [
        ("Ensemble engine", getattr(m7, "classifier_engine", "ICH M7 Hybrid Ensemble")),
        ("Layer A (Empirical DB)", m7.empirical_db_match),
        (
            "Layer B (SMARTS / metabolization)",
            "Positive" if m7.alert_present else "Negative",
        ),
    ]
    if m7.gnn_mutagenic_probability is not None:
        rows.append(
            (
                "Layer C P(Mutagenic)",
                f"{m7.gnn_mutagenic_probability:.3f} · {m7.gnn_model_name or 'GNN'}",
            )
        )
    else:
        rows.append(("Layer C P(Mutagenic)", "Unavailable (RDKit/DB standalone mode)"))
    if m7.gnn_fallback_log:
        rows.append(("Layer C fallback log", m7.gnn_fallback_log))
    label = m7.assigned_class_label or f"Class {m7.class_number}"
    if m7.rule_based_class is not None and m7.rule_based_class != m7.class_number:
        rows.append(
            (
                "Ensemble assigned class",
                f"{label} (rule Class {m7.rule_based_class} → final Class {m7.class_number})",
            )
        )
    else:
        rows.append(("Ensemble assigned class", label))
    if m7.ai_inferred_liability:
        rows.append(("AI inferred liability", "Class 3-Alert (AI Inferred Liability)"))
    # 2026-09-03 수정: CoC·전문가 검토 판정 근거를 보고서 표에도 남깁니다.
    if getattr(m7, "cohort_of_concern", False):
        rows.append(
            (
                "Cohort of concern",
                "Yes — TTC not applicable; compound-specific AI required "
                "(ICH M7(R2) Section 7.5)",
            )
        )
    if getattr(m7, "expert_review_required", False):
        rows.append(
            (
                "Expert review required",
                getattr(m7, "expert_review_reason", "") or "Yes",
            )
        )
    if m7.ensemble_log:
        rows.append(("Ensemble log", " | ".join(m7.ensemble_log[:6])))
    return rows


@dataclass
class M7RegulatoryOptionSummary:
    """ICH M7(R2) Option 3/4 적용 여부 및 Class 판정 근거."""
    primary_option: str
    applied_options: list[str]
    not_applied_options: list[str]
    supporting_data: list[str]
    evidence_chain: list[str]
    conclusion_paragraph: str


def _qsar_option_applied(ctx: AssessmentContext) -> bool:
    """Option 3(통계 (Q)SAR) 적용 여부.

    2026-09-03 수정: AMBIT 응답 문자열만 보던 탓에, Hybrid Ensemble Layer C가
    정상 동작해 P(Mutagenic)을 산출해도 False가 되어 "Option 3 미적용"과
    "Layer C 적용"이 한 보고서에 함께 실렸습니다. ICH M7 정의상 통계 (Q)SAR가
    곧 Option 3이므로, Layer C 확률이 있으면 적용으로 판정하고 AMBIT 판정은
    fallback으로 둡니다.
    """
    if getattr(ctx.m7, "gnn_mutagenic_probability", None) is not None:
        return True
    text = f"{ctx.qsar_engine} {ctx.qsar_ames} {ctx.oecd_connection_status}".lower()
    skip = ("not integrated", "pending", "not invoked", "not queried")
    return bool(ctx.qsar_ames.strip()) and not any(s in text for s in skip)


def _statistical_qsar_source(ctx: AssessmentContext) -> str:
    """Option 3의 실제 출처 문구 — Layer C 통계 모델 우선, 없으면 AMBIT."""
    if getattr(ctx.m7, "gnn_mutagenic_probability", None) is not None:
        model = getattr(ctx.m7, "gnn_model_name", "") or "statistical (Q)SAR model"
        return (
            f"Hybrid Ensemble Layer C statistical (Q)SAR — {model}; "
            f"P(Mutagenic)={ctx.m7.gnn_mutagenic_probability:.3f}"
        )
    return ctx.qsar_engine


def _oecd_context_kwargs(smiles: str) -> dict[str, str]:
    try:
        from p52_ambit_client import build_oecd_assessment

        return build_oecd_assessment(smiles)
    except Exception as exc:
        return {
            "qsar_engine": "OECD QSAR Toolbox (AMBIT2)",
            "qsar_ames": f"OECD QSAR query failed: {exc}",
            "oecd_connection_status": "Error",
            "oecd_hydrolysis": "N/A",
            "oecd_biodegradability": "N/A",
            "oecd_mutagenicity_alert": "N/A",
        }


def _qsar_context_fields(smiles: str) -> tuple[str, str]:
    oecd = _oecd_context_kwargs(smiles)
    return oecd["qsar_engine"], oecd["qsar_ames"]


def build_m7_regulatory_option_summary(ctx: AssessmentContext) -> M7RegulatoryOptionSummary:
    """
    ICH M7(R2) Option 3 vs Option 4 — 실제 적용 경로와 Class 판정 근거를 분리 기술.
    두 Option을 동시 적용한 것처럼 보이지 않도록 primary route를 명시합니다.
    """
    m7 = ctx.m7
    class_num = m7.class_number if m7 else 5
    # 2026-09-03 수정: KeyError 방지 — 미등록 Class는 Class 5 기본값.
    m7_info = ICH_M7_REGULATORY.get(class_num, ICH_M7_REGULATORY[5])
    alerts = list(ctx.alerts or [])
    alert_positive = bool(alerts)
    qsar_applied = _qsar_option_applied(ctx)

    applied: list[str] = []
    not_applied: list[str] = []
    supporting: list[str] = []
    evidence: list[str] = []

    applied.append(
        "Option 4 — Expert rule-based assessment (ICH M7(R2) structural alert evaluation)"
    )
    evidence.append(
        f"Expert rule engine: {ctx.expert_engine}."
    )
    evidence.append(
        f"Structural alert outcome: {'POSITIVE' if alert_positive else 'NEGATIVE'}"
        + (f" — {', '.join(alerts)}" if alerts else " — none detected.")
    )

    if qsar_applied:
        applied.append("Option 3 — (Q)SAR-based assessment (statistical Ames/mutagenicity models)")
        # 2026-09-03 수정: Option 3의 근거를 실제 통계 모델(Layer C)로 명시.
        evidence.append(f"(Q)SAR engine: {_statistical_qsar_source(ctx)}.")
        evidence.append(f"Ames / mutagenicity (Q)SAR output: {ctx.qsar_ames}.")
    else:
        not_applied.append(
            "Option 3 — (Q)SAR-based assessment (statistical bacterial mutagenicity models)"
        )
        evidence.append(
            "(Q)SAR Option 3 not applied: validated statistical Ames models are not integrated "
            f"in this release ({ctx.qsar_ames}). Classification does not rely on (Q)SAR output."
        )

    not_applied.extend(
        [
            "Option 1 — Literature search (manual; not performed in silico)",
            "Option 2 — Read-across assessment (not performed)",
            "Option 5 — Experimental bacterial mutagenicity (Ames) testing (not performed)",
        ]
    )

    db_low = (ctx.db_match or "").lower()
    if "validated" in db_low:
        supporting.append(f"In-house database (drug_index.db): {ctx.db_match}")
        evidence.append(
            f"Experimental reference in indexed database: {ctx.db_ref or 'matched entry'}."
        )
    else:
        supporting.append("In-house database: no exact experimental match (theoretical prediction)")
        evidence.append(
            "In-house database cross-validation: no validated experimental match for this SMILES."
        )

    if ctx.pubchem_cid or ctx.pubchem_status:
        pub_line = ctx.pubchem_status or f"PubChem CID {ctx.pubchem_cid}"
        supporting.append(f"Open-source structure mapping: {pub_line}")
        evidence.append(
            "PubChem open data used for structure identity and public BioAssay screening only; "
            "does not replace ICH M7 Option 3 (Q)SAR when (Q)SAR is not deployed."
        )

    if m7 and getattr(m7, "empirical_db_match", "None") not in ("None", "", "Not evaluated"):
        supporting.append(
            f"Empirical M7 index (Hansen/ISSCAN/Tox21): {m7.empirical_db_match}"
        )
        evidence.append(
            f"Empirical Ames evidence: {getattr(m7, 'ames_evidence', 'N/A')}; "
            f"carcinogenicity: {getattr(m7, 'carcinogenicity_evidence', 'N/A')}."
        )
        if getattr(m7, "recommended_ttc", ""):
            evidence.append(f"Recommended TTC threshold: {m7.recommended_ttc}.")

    if m7 and getattr(m7, "gnn_mutagenic_probability", None) is not None:
        # 2026-09-03 수정: 실제 사용 모델명을 그대로 적어 Option 3과의 관계를 명확히 함.
        applied.append(
            "Hybrid Ensemble Layer C — statistical (Q)SAR mutagenicity inference "
            f"({getattr(m7, 'gnn_model_name', '') or 'statistical model'}); "
            "this is the Option 3 model layer"
        )
        supporting.append(
            f"Layer C GNN ({getattr(m7, 'gnn_model_name', 'GNN')}): "
            f"P(Mutagenic)={m7.gnn_mutagenic_probability:.3f}"
        )
        evidence.append(
            f"Hybrid ensemble deep-learning signal: P(Mutagenic)={m7.gnn_mutagenic_probability:.3f}."
        )
        if getattr(m7, "ai_inferred_liability", False):
            evidence.append(
                "AI Inferred Liability override — Class 3-Alert (hidden mutagenic signal without SMARTS/DB hit)."
            )
        if getattr(m7, "rule_based_class", None) is not None and m7.rule_based_class != m7.class_number:
            evidence.append(
                f"Ensemble rescoring: rule Class {m7.rule_based_class} → final Class {m7.class_number}."
            )

    if alert_positive and class_num == 3:
        primary = "Option 4 — Expert rule-based assessment (PRIMARY)"
        class_driver = (
            "Class 3 assigned because a structural alert for mutagenicity was identified "
            f"({', '.join(alerts)}) and no bacterial mutagenicity or carcinogenicity data "
            "are available for the impurity (ICH M7 Class 3 definition)."
        )
    elif alert_positive and class_num == 4:
        primary = "Option 4 — Expert rule-based assessment (PRIMARY)"
        class_driver = (
            "Class 4 considered because alerting structure is shared with the drug substance "
            "and non-mutagenicity may be justified by drug-related data; expert review required."
        )
    elif not alert_positive and class_num >= 5:
        primary = "Option 4 — Expert rule-based assessment (PRIMARY)"
        class_driver = (
            "Class 5 assigned because no Ashby-Tennant / ICH M7 structural alerts were detected "
            "for mutagenicity under the expert rule set applied."
        )
    elif qsar_applied and not alert_positive:
        primary = "Option 3 — (Q)SAR-based assessment (PRIMARY)"
        class_driver = f"(Q)SAR prediction supports Class {class_num} with no conflicting expert alerts."
    else:
        primary = "Option 4 — Expert rule-based assessment (PRIMARY)"
        class_driver = f"Class {class_num} — {m7_info['definition']}"

    evidence.append(f"ICH M7 hazard class assignment logic: {class_driver}")

    applied_txt = "; ".join(applied) if applied else "None"
    excluded_txt = "; ".join(not_applied[:3]) + (
        f"; +{len(not_applied) - 3} other options not performed" if len(not_applied) > 3 else ""
    )

    conclusion = (
        f"PRIMARY REGULATORY ROUTE: {primary}. "
        f"This report applies {applied_txt}. "
        f"Option 3 and Option 4 were not used interchangeably: "
        f"{'both contributed' if qsar_applied and alert_positive else 'only the primary route above determined the hazard class'}. "
        f"Not applied in silico: {excluded_txt}. "
        "Confirmatory literature review, qualified expert sign-off, and experimental Ames testing "
        "remain required where ICH M7(R2) demands."
    )

    return M7RegulatoryOptionSummary(
        primary_option=primary,
        applied_options=applied,
        not_applied_options=not_applied,
        supporting_data=supporting,
        evidence_chain=evidence,
        conclusion_paragraph=conclusion,
    )


def generate_formal_assessment_report(ctx: AssessmentContext) -> str:
    m7 = ctx.m7
    class_num = m7.class_number if m7 else 5
    # 2026-09-03 수정: KeyError 방지 — 미등록 Class는 Class 5 기본값.
    m7_info = ICH_M7_REGULATORY.get(class_num, ICH_M7_REGULATORY[5])
    # 2026-09-03 수정: CoC 판정을 6절 Acceptable intake에도 반영.
    m7_control, q3_control, ai_limit = m7_action_flags(
        class_num, cohort_of_concern=bool(getattr(m7, "cohort_of_concern", False))
    )
    alert_positive = bool(ctx.alerts)
    alert_status = "POSITIVE" if alert_positive else "NEGATIVE"
    alerts_text = ", ".join(ctx.alerts) if ctx.alerts else "None detected"

    pubchem_line = ctx.pubchem_status or (
        f"PubChem CID {ctx.pubchem_cid} — structure mapped successfully"
        if ctx.pubchem_cid
        else "No PubChem lookup performed"
    )

    stab_label, stab_detail = _translate_stability(ctx.stability, ctx.stability_detail)
    mechanism = _translate_mechanism(ctx.mechanism)
    expert_comment = ctx.expert_comment or _expert_rationale_en(m7, ctx.alerts)
    option_summary = build_m7_regulatory_option_summary(ctx)

    lines = [
        "=" * 88,
        REPORT_TITLE_LINE1.center(88),
        REPORT_TITLE_LINE2.center(88),
        "ICH M7(R2) — Complementary (Q)SAR and Expert Rule-Based Evaluation".center(88),
        "=" * 88,
        "",
        "DOCUMENT CONTROL",
        "-" * 88,
        f"  Report Number      : {ctx.report_id}",
        f"  Document Version   : {DOCUMENT_VERSION}",
        f"  Issue Date         : {ctx.issue_date}",
        f"  Software Platform  : {PLATFORM_NAME} — {PLATFORM_FULL}",
        f"  Regulatory Basis   : ICH M7(R2); ICH Q3A(R2); ICH Q3B(R2)",
        f"  OECD QSAR status   : {ctx.oecd_connection_status}",
        f"  OECD hydrolysis    : {ctx.oecd_hydrolysis}",
        f"  OECD biodegradability : {ctx.oecd_biodegradability}",
        f"  OECD mutagenicity alert : {ctx.oecd_mutagenicity_alert}",
        "",
        "CONFIDENTIALITY STATEMENT",
        "  This document is intended for regulatory submission support and internal GMP use.",
        "  Computational predictions do not replace qualified expert review or confirmatory testing.",
        "",
        "=" * 88,
        "1. SUBSTANCE IDENTIFICATION",
        "=" * 88,
        f"  Substance category     {substance_category_label(ctx)}",
        f"  Substance name         {ctx.name}",
        f"  Chemical structure     {ctx.smiles}",
        f"  Molecular weight       {ctx.mw} g/mol (RDKit)",
        "",
    ]
    if includes_process_route_section(ctx):
        lines.extend(
            [
                "=" * 88,
                "2. PROCESS ROUTE AND IMPURITY FORMATION PATHWAY",
                "=" * 88,
                f"  Process inputs         Reactant 1: {ctx.reactant1 or 'N/A'}",
                f"                         Reactant 2: {ctx.reactant2 or 'N/A'}",
                f"                         Reagent: {ctx.reagent or 'N/A'}",
                f"  Drug substance (API)   {ctx.product or ctx.smiles or 'N/A'}",
            ]
        )
        if ctx.reference_structure:
            lines.append(
                f"  Reference structure    {ctx.reference_structure} "
                "(Class 4 comparison basis only — not the assessed substance in this report)"
            )
        if ctx.primary_impurity_smiles or ctx.primary_impurity_name:
            lines.append(
                f"  Representative impurity {ctx.primary_impurity_name or 'N/A'} — "
                f"{ctx.primary_impurity_smiles or 'N/A'} (see Appendix A)"
            )
        lines.extend(
            [
                f"  Process conditions     Nitrite source present {_chk(ctx.nitrite)}",
                f"                         Acidic conditions {_chk(ctx.acidic_condition)}",
                f"  Formation mechanism    {mechanism}",
                f"  Predicted impurity     {ctx.predicted_impurity_smiles or 'N/A'}",
                "",
            ]
        )
        if ctx.product_assessment and not ctx.is_final_product:
            pa = ctx.product_assessment
            pa_m7 = pa.m7.class_number if pa.m7 else 5
            # 2026-09-03 수정: KeyError 방지 — 미등록 Class는 Class 5 기본값.
            pa_info = ICH_M7_REGULATORY.get(pa_m7, ICH_M7_REGULATORY[5])
            pa_alert_status = "POSITIVE" if pa.alerts else "NEGATIVE"
            pa_alerts_text = ", ".join(pa.alerts) if pa.alerts else "None detected"
            pa_stab_label, pa_stab_detail = _translate_stability(pa.stability, pa.stability_detail)
            pa_pub = pa.pubchem_status or (
                f"PubChem CID {pa.pubchem_cid}" if pa.pubchem_cid else "No PubChem lookup performed"
            )
            lines.extend(
                [
                    "=" * 88,
                    "2.1 DRUG SUBSTANCE (FINAL PRODUCT) — STRUCTURE AND ASSESSMENT",
                    "=" * 88,
                    f"  Substance name         {pa.name}",
                    f"  Chemical structure     {pa.smiles}",
                    f"  Molecular weight       {pa.mw} g/mol (RDKit)",
                    f"  Structural alert scan  {pa_alert_status}",
                    f"  Alerts identified      {pa_alerts_text}",
                    # 2026-09-03 수정: name 자체가 "Class N — …" 형태라 등급이 두 번 찍혔습니다.
                    f"  ICH M7 hazard class    {pa_info['name']}",
                    f"  Assessment conclusion  {pa.expert_comment or _expert_rationale_en(pa.m7, pa.alerts)}",
                    f"  Lipophilicity (Log P)  {pa.logp}",
                    f"  Process stability      {pa_stab_label}",
                    f"  Stability rationale    {pa_stab_detail}",
                    f"  External database      {pa_pub}",
                    f"  In-house DB match      {pa.db_match or 'N/A'}",
                    "",
                ]
            )
    lines.extend(
        [
        "=" * 88,
        "3. COMPLEMENTARY HYBRID ASSESSMENT (ICH M7(R2))",
        "=" * 88,
        "",
        "3.1 Expert Rule-Based Assessment (Derek-style)",
        f"  Computational engine   {ctx.expert_engine}",
        f"  Structural alert scan  {alert_status}",
        f"  Alerts identified      {alerts_text}",
        f"  Assessment conclusion  {expert_comment}",
        "",
        "3.2 Statistical / (Q)SAR-Based Assessment (Sarah-style)",
        "3.2a OECD QSAR Toolbox (AMBIT2)",
        f"  OECD QSAR engine       {ctx.qsar_engine}",
        f"  AMBIT connection       {ctx.oecd_connection_status}",
        f"  OECD hydrolysis        {ctx.oecd_hydrolysis}",
        f"  OECD biodegradability  {ctx.oecd_biodegradability}",
        f"  OECD mutagenicity alert {ctx.oecd_mutagenicity_alert}",
        f"  Ames / (Q)SAR summary  {ctx.qsar_ames}",
        f"  External database      {pubchem_line}",
        "",
        "3.2b Hybrid Ensemble Matrix (Layers A / B / C)",
        ]
    )
    for label, val in hybrid_ensemble_report_rows(m7):
        lines.append(f"  {label:<24} {val}")
    lines.extend(
        [
        "",
        # 2026-09-03 수정: 실측/fallback 여부를 함께 표시.
        f"3.3 In-House Master Database Cross-Validation ({indexed_db_count_label(ctx)})",
        f"  Database match status  {ctx.db_match}",
        f"  Reference index        {ctx.db_ref or 'N/A'}",
        "",
        "3.4 ICH M7(R2) Regulatory Option Determination and Classification Rationale",
        f"  Primary assessment route   {option_summary.primary_option}",
        f"  Options applied in silico  {'; '.join(option_summary.applied_options)}",
        f"  Options NOT applied        {'; '.join(option_summary.not_applied_options[:4])}",
        "  Supporting open data       " + (
            "; ".join(option_summary.supporting_data) if option_summary.supporting_data else "N/A"
        ),
        "  Evidence chain for Class assignment:",
        ]
    )
    for i, ev in enumerate(option_summary.evidence_chain, 1):
        lines.append(f"    ({i}) {ev}")
    lines.extend(
        [
        "",
        "=" * 88,
        "4. PHYSICOCHEMICAL PROPERTIES AND EXPOSURE CONSIDERATIONS",
        "=" * 88,
        f"  Lipophilicity (Log P)  {ctx.logp}",
        f"  Log D / pKa / bp       {ctx.logd_pka_bp}",
        f"  OECD hydrolysis (AMBIT) {ctx.oecd_hydrolysis}",
        f"  OECD biodegradability  {ctx.oecd_biodegradability}",
        f"  Process stability      {stab_label}",
        f"  Stability rationale    {stab_detail}",
        f"  Exposure assessment    {ctx.exposure_comment}",
        "",
        "=" * 88,
        "5. FINAL ICH M7(R2) HAZARD CLASSIFICATION",
        "=" * 88,
        "",
    ]
    )
    class_basis = (
        option_summary.evidence_chain[-1]
        if option_summary.evidence_chain
        else "See Section 3.4"
    )
    lines.extend(
        format_m7_executive_text_block(
            class_num,
            build_m7_executive_rows(
                class_num,
                classification=m7_info["name"],
                definition=m7_info["definition"],
                primary_route=option_summary.primary_option,
                classification_basis=class_basis,
                recommended_control=m7_info["control"],
            ),
        )
    )
    lines.extend(
        [
        # 2026-09-03 수정: 하드코딩된 "Class N" 요약을 p16 class_label()로 교체(단일 출처).
        "  ICH M7(R2) Class definitions (summary):",
        *[f"    {class_label(_n)}" for _n in (1, 2, 3, 4, 5)],
        "",
        "=" * 88,
        "6. RECOMMENDED ACTIONS AND ACCEPTABLE INTAKE",
        "=" * 88,
        f"  Regulatory action      {_chk(m7_control)} Control per ICH M7 (mutagenic impurity limits)",
        f"                         {_chk(q3_control)} Manage per ICH Q3A/Q3B (general impurity limits)",
        f"  Acceptable intake      {ai_limit}",
        "",
        "=" * 88,
        "7. CONCLUSION AND LIMITATIONS",
        "=" * 88,
        f"  {option_summary.conclusion_paragraph}",
        "  Results remain preliminary until qualified toxicological review, literature search,",
        "  and where warranted, confirmatory bacterial mutagenicity (Ames) testing are completed.",
        "",
        f"  Report status          In silico assessment complete — expert sign-off required",
        f"  Generated by           {PLATFORM_NAME} v{DOCUMENT_VERSION}",
        "=" * 88,
        ]
    )
    return "\n".join(lines)


def build_appendix_all_species(rows: list[dict]) -> str:
    lines = [
        "",
        "=" * 88,
        "APPENDIX A — ALL PREDICTED PROCESS-RELATED SPECIES",
        "=" * 88,
        "",
    ]
    for i, row in enumerate(rows, 1):
        alerts = row.get("alerts") or []
        alert_text = ", ".join(alerts) if alerts else "None"
        lines.append(
            f"  [{i:03d}] {row.get('name', 'Unknown')} | ICH M7 Class {row.get('ich_m7_class', '—')} | "
            f"{row.get('classification', '—')}"
        )
        lines.append(f"       Log P: {row.get('logp', '—')} | Stability: {row.get('stability', '—')}")
        lines.append(f"       Alerts: {alert_text}")
        lines.append(f"       SMILES: {row.get('smiles', '—')}")
        lines.append("")
    return "\n".join(lines)


def build_single_compound_context(
    *,
    name: str,
    smiles: str,
    mw: str | float,
    properties: dict | None,
    advanced: dict | None,
    stability: tuple[str, str],
    m7: ICHM7Result,
    alerts: list[str],
    db_status: str,
    db_ref: str = "",
    pubchem_cid: str | int | None = None,
    # 2026-09-03 수정: 하드코딩 13040 → None이면 p20에서 실제 카운트를 조회합니다.
    indexed_count: int | None = None,
    reference_structure: str = "",
) -> AssessmentContext:
    # 2026-09-03 수정: 보고서에 찍히는 DB 건수를 실제 인덱스에서 조회(실패 시 상수).
    if indexed_count is None:
        resolved_count, count_is_fallback = resolve_indexed_db_count()
    else:
        resolved_count, count_is_fallback = int(indexed_count), False
    logp = str(properties.get("LogP", "N/A")) if properties else "N/A"
    adv = advanced or {}
    logd = adv.get("LogD (My Algo)", "N/A")
    pka = adv.get("pKa (My Algo)", "N/A")
    bp = adv.get("Boiling Point", "N/A")
    if logd in ("추후 업데이트 예정", "pending", ""):
        logd = "N/A"
    if pka in ("추후 업데이트 예정", "pending", ""):
        pka = "N/A"
    if bp in ("추후 업데이트 예정", "pending", ""):
        bp = "N/A"

    stab_en, stab_detail_en = _translate_stability(stability[0], stability[1])
    exposure = (
        f"ICH M7 Class {m7.class_number} assigned. {stab_detail_en} "
        "Residual impurity level and patient exposure must be evaluated against the recommended "
        "acceptable intake in the context of maximum daily dose and duration of treatment."
    )
    oecd_kwargs = _oecd_context_kwargs(smiles)
    return AssessmentContext(
        report_id=new_report_id(),
        issue_date=english_issue_date(),
        is_reactant=False,
        is_process_impurity=False,
        is_standalone_compound=True,
        name=name,
        smiles=smiles,
        mw=str(mw),
        mechanism="Standalone compound toxicity and physicochemical property assessment",
        expert_comment=_expert_rationale_en(m7, alerts),
        alerts=alerts,
        db_match=db_status,
        db_ref=db_ref or "N/A",
        pubchem_cid=str(pubchem_cid) if pubchem_cid else None,
        logp=f"{logp} (RDKit MolLogP)",
        logd_pka_bp=f"LogD {logd} | pKa {pka} | Boiling point {bp}",
        stability=stab_en,
        stability_detail=stab_detail_en,
        exposure_comment=exposure,
        m7=m7,
        indexed_db_count=resolved_count,
        indexed_db_count_is_fallback=count_is_fallback,
        qsar_engine=oecd_kwargs["qsar_engine"],
        qsar_ames=oecd_kwargs["qsar_ames"],
        oecd_connection_status=oecd_kwargs["oecd_connection_status"],
        oecd_hydrolysis=oecd_kwargs["oecd_hydrolysis"],
        oecd_biodegradability=oecd_kwargs["oecd_biodegradability"],
        oecd_mutagenicity_alert=oecd_kwargs["oecd_mutagenicity_alert"],
        reference_structure=reference_structure or "",
    )


def build_process_impurity_context(
    *,
    impurity_name: str,
    impurity_smiles: str,
    mw: str | float,
    reactant1: str,
    reactant2: str,
    reagent: str = "",
    product: str,
    nitrite: bool,
    mechanism: str,
    stability: tuple[str, str],
    m7: ICHM7Result,
    alerts: list[str],
    db_status: str,
    db_ref: str = "",
    pubchem_cid: str | int | None = None,
    properties: dict | None = None,
    advanced: dict | None = None,
    # 2026-09-03 수정: 하드코딩 13040 → None이면 p20에서 실제 카운트를 조회합니다.
    indexed_count: int | None = None,
    is_unreacted_sm: bool = False,
    product_assessment: ProductAssessment | None = None,
) -> AssessmentContext:
    # 2026-09-03 수정: 보고서에 찍히는 DB 건수를 실제 인덱스에서 조회(실패 시 상수).
    if indexed_count is None:
        resolved_count, count_is_fallback = resolve_indexed_db_count()
    else:
        resolved_count, count_is_fallback = int(indexed_count), False
    logp = str(properties.get("LogP", "N/A")) if properties else "N/A"
    adv = advanced or {}
    logd = adv.get("LogD (My Algo)", "N/A")
    pka = adv.get("pKa (My Algo)", "N/A")
    bp = adv.get("Boiling Point", "N/A")
    if logd in ("추후 업데이트 예정", "pending", ""):
        logd = "N/A"
    if pka in ("추후 업데이트 예정", "pending", ""):
        pka = "N/A"
    if bp in ("추후 업데이트 예정", "pending", ""):
        bp = "N/A"

    stab_en, stab_detail_en = _translate_stability(stability[0], stability[1])
    mech_en = mechanism or "RDKit rdChemReactions (SMIRKS) — side-reaction pathway prediction"
    if nitrite and "nitro" not in mech_en.lower():
        mech_en += "; nitrosamine formation pathway screened (nitrite source present)."

    oecd_kwargs = _oecd_context_kwargs(impurity_smiles)
    return AssessmentContext(
        report_id=new_report_id(),
        issue_date=english_issue_date(),
        is_reactant=is_unreacted_sm,
        is_process_impurity=not is_unreacted_sm,
        name=impurity_name,
        smiles=impurity_smiles,
        mw=str(mw),
        reactant1=reactant1,
        reactant2=reactant2,
        reagent=reagent or "",
        product=product,
        nitrite=nitrite,
        mechanism=mech_en,
        predicted_impurity_smiles=impurity_smiles,
        alerts=alerts,
        expert_comment=_expert_rationale_en(m7, alerts),
        db_match=db_status,
        db_ref=db_ref or "N/A",
        pubchem_cid=str(pubchem_cid) if pubchem_cid else None,
        logp=f"{logp} (RDKit MolLogP)",
        logd_pka_bp=f"LogD {logd} | pKa {pka} | Boiling point {bp}",
        stability=stab_en,
        stability_detail=stab_detail_en,
        exposure_comment=(
            f"Process-related impurity — ICH M7 Class {m7.class_number if m7 else 5}. {stab_detail_en} "
            "Carryover, purge, and fate in the synthetic route should be assessed."
        ),
        m7=m7,
        indexed_db_count=resolved_count,
        indexed_db_count_is_fallback=count_is_fallback,
        product_assessment=product_assessment,
        qsar_engine=oecd_kwargs["qsar_engine"],
        qsar_ames=oecd_kwargs["qsar_ames"],
        oecd_connection_status=oecd_kwargs["oecd_connection_status"],
        oecd_hydrolysis=oecd_kwargs["oecd_hydrolysis"],
        oecd_biodegradability=oecd_kwargs["oecd_biodegradability"],
        oecd_mutagenicity_alert=oecd_kwargs["oecd_mutagenicity_alert"],
    )


def build_process_product_context(
    *,
    product_name: str,
    product_smiles: str,
    mw: str | float,
    reactant1: str,
    reactant2: str,
    reagent: str = "",
    nitrite: bool,
    reference_structure: str = "",
    primary_impurity_name: str = "",
    primary_impurity_smiles: str = "",
    mechanism: str = "",
    stability: tuple[str, str],
    m7: ICHM7Result,
    alerts: list[str],
    db_status: str,
    db_ref: str = "",
    pubchem_cid: str | int | None = None,
    properties: dict | None = None,
    advanced: dict | None = None,
    # 2026-09-03 수정: 하드코딩 13040 → None이면 p20에서 실제 카운트를 조회합니다.
    indexed_count: int | None = None,
) -> AssessmentContext:
    """Process-mode regulatory report — assesses [4] final product (API), not reference structure."""
    # 2026-09-03 수정: 보고서에 찍히는 DB 건수를 실제 인덱스에서 조회(실패 시 상수).
    if indexed_count is None:
        resolved_count, count_is_fallback = resolve_indexed_db_count()
    else:
        resolved_count, count_is_fallback = int(indexed_count), False
    logp = str(properties.get("LogP", "N/A")) if properties else "N/A"
    adv = advanced or {}
    logd = adv.get("LogD (My Algo)", "N/A")
    pka = adv.get("pKa (My Algo)", "N/A")
    bp = adv.get("Boiling Point", "N/A")
    if logd in ("추후 업데이트 예정", "pending", ""):
        logd = "N/A"
    if pka in ("추후 업데이트 예정", "pending", ""):
        pka = "N/A"
    if bp in ("추후 업데이트 예정", "pending", ""):
        bp = "N/A"

    stab_en, stab_detail_en = _translate_stability(stability[0], stability[1])
    mech_en = mechanism or "RDKit rdChemReactions (SMIRKS) — process-related species prediction"
    if nitrite and "nitro" not in mech_en.lower():
        mech_en += "; nitrosamine formation pathway screened (nitrite source present)."

    oecd_kwargs = _oecd_context_kwargs(product_smiles)
    return AssessmentContext(
        report_id=new_report_id(),
        issue_date=english_issue_date(),
        is_reactant=False,
        is_process_impurity=False,
        is_final_product=True,
        name=product_name,
        smiles=product_smiles,
        mw=str(mw),
        reactant1=reactant1,
        reactant2=reactant2,
        reagent=reagent or "",
        product=product_smiles,
        reference_structure=reference_structure or "",
        primary_impurity_name=primary_impurity_name or "",
        primary_impurity_smiles=primary_impurity_smiles or "",
        nitrite=nitrite,
        mechanism=mech_en,
        predicted_impurity_smiles=primary_impurity_smiles or "",
        alerts=alerts,
        expert_comment=_expert_rationale_en(m7, alerts),
        db_match=db_status,
        db_ref=db_ref or "N/A",
        pubchem_cid=str(pubchem_cid) if pubchem_cid else None,
        logp=f"{logp} (RDKit MolLogP)",
        logd_pka_bp=f"LogD {logd} | pKa {pka} | Boiling point {bp}",
        stability=stab_en,
        stability_detail=stab_detail_en,
        exposure_comment=(
            f"Drug substance (final product) — ICH M7 Class {m7.class_number if m7 else 5}. "
            f"{stab_detail_en} Process-related impurity carryover is summarized in Section 2 "
            "and Appendix A."
        ),
        m7=m7,
        indexed_db_count=resolved_count,
        indexed_db_count_is_fallback=count_is_fallback,
        qsar_engine=oecd_kwargs["qsar_engine"],
        qsar_ames=oecd_kwargs["qsar_ames"],
        oecd_connection_status=oecd_kwargs["oecd_connection_status"],
        oecd_hydrolysis=oecd_kwargs["oecd_hydrolysis"],
        oecd_biodegradability=oecd_kwargs["oecd_biodegradability"],
        oecd_mutagenicity_alert=oecd_kwargs["oecd_mutagenicity_alert"],
    )
