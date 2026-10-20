"""ICH M7(R2) 공식 용어 — 단일 출처.

화면·보고서에서 ICH M7 용어를 쓸 때는 반드시 이 모듈의 상수를
가져다 쓰세요. 문자열을 직접 적으면 표현이 갈라집니다.

여기 담긴 영문은 전부 ICH M7(R2) 가이드라인 원문에서 가져온 것입니다.
근거와 한국어 대응은 `md파일/ICH_M7_공식용어집.md` 참조.

출처
    ICH M7(R2) — Assessment and Control of DNA Reactive (Mutagenic)
    Impurities in Pharmaceuticals to Limit Potential Carcinogenic Risk
    EMA/CHMP/ICH/83812/2013

2026-08-25 작성.
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "GUIDELINE_ID",
    "GUIDELINE_TITLE",
    "GUIDELINE_TITLE_KO",
    "CLASS_DEFINITIONS",
    "CLASS_ACTIONS",
    "CLASS_SHORT",
    "CLASS_SHORT_KO",
    "CONTROL_OPTIONS",
    "QSAR_TERMS",
    "COHORT_OF_CONCERN",
    "COHORT_OF_CONCERN_AI_TEXT",
    "TTC_DEFAULT_UG_DAY",
    "LTL_LIMITS_UG_DAY",
    "LTL_LIMIT_TABLE",
    "class_label",
    "class_definition",
    "class_action",
    "is_cohort_of_concern_class",
    "ltl_limit_ug_day",
    "ltl_limit_label",
]


# ---------------------------------------------------------------------------
# 가이드라인 식별
# ---------------------------------------------------------------------------

GUIDELINE_ID = "ICH M7(R2)"

GUIDELINE_TITLE = (
    "Assessment and Control of DNA Reactive (Mutagenic) Impurities in "
    "Pharmaceuticals to Limit Potential Carcinogenic Risk"
)

GUIDELINE_TITLE_KO = (
    "의약품 중 DNA 반응성(변이원성) 불순물의 평가 및 관리 — "
    "잠재적 발암 위해 제한"
)

GUIDELINE_REFERENCE = "EMA/CHMP/ICH/83812/2013"


# ---------------------------------------------------------------------------
# Table 1 — 불순물 분류
# 가이드라인 원문 그대로입니다. 수정하지 마세요.
# ---------------------------------------------------------------------------

CLASS_DEFINITIONS: dict[int, str] = {
    1: "Known mutagenic carcinogens",
    2: (
        "Known mutagens with unknown carcinogenic potential "
        "(bacterial mutagenicity positive, no rodent carcinogenicity data)"
    ),
    3: (
        "Alerting structure, unrelated to the structure of the drug substance; "
        "no mutagenicity data"
    ),
    4: (
        "Alerting structure, same alert in drug substance or compounds related "
        "to the drug substance (e.g., process intermediates) which have been "
        "tested and are non-mutagenic"
    ),
    5: (
        "No structural alerts, or alerting structure with sufficient data to "
        "demonstrate lack of mutagenicity or carcinogenicity"
    ),
}

CLASS_DEFINITIONS_KO: dict[int, str] = {
    1: "알려진 변이원성 발암물질",
    2: "알려진 변이원성 물질, 발암성 불명 (세균 변이원성 양성, 설치류 발암성 자료 없음)",
    3: "알러트 구조, 원료의약품 구조와 무관, 변이원성 데이터 없음",
    4: "알러트 구조, 원료의약품 또는 관련 화합물과 동일 알러트이며 시험 결과 비변이원성",
    5: "알러트 없음, 또는 비변이원성·비발암성을 입증할 충분한 데이터 보유",
}

CLASS_ACTIONS: dict[int, str] = {
    1: "Control at or below compound-specific acceptable limit",
    2: "Control at or below acceptable limits (appropriate TTC)",
    3: (
        "Control at or below acceptable limits (appropriate TTC) or conduct "
        "bacterial mutagenicity assay; if non-mutagenic = Class 5, "
        "if mutagenic = Class 2"
    ),
    4: "Treat as non-mutagenic impurity",
    5: "Treat as non-mutagenic impurity",
}

CLASS_ACTIONS_KO: dict[int, str] = {
    1: "화합물 특이적 허용 한도 이하로 관리",
    2: "허용 한도(해당 TTC) 이하로 관리",
    3: "허용 한도(해당 TTC) 이하로 관리하거나 세균 변이원성 시험 수행",
    4: "비변이원성 불순물로 취급",
    5: "비변이원성 불순물로 취급",
}

# 배지·카드처럼 좁은 자리에서 쓸 축약형
CLASS_SHORT: dict[int, str] = {
    1: "Known mutagenic carcinogen",
    2: "Known mutagen, carcinogenic potential unknown",
    3: "Alerting structure, no mutagenicity data",
    4: "Alerting structure shared with drug substance",
    5: "No structural alert / demonstrated non-mutagenic",
}

CLASS_SHORT_KO: dict[int, str] = {
    1: "알려진 변이원성 발암물질",
    2: "알려진 변이원성, 발암성 불명",
    3: "알러트 구조, 시험 데이터 없음",
    4: "원료와 동일한 알러트",
    5: "알러트 없음 · 비변이원성 확인",
}

CLASS_FOOTNOTE = (
    "Or other relevant positive mutagenicity data indicative of DNA-reactivity "
    "related induction of gene mutations (e.g., positive findings in in vivo "
    "gene mutation studies)"
)


# ---------------------------------------------------------------------------
# 8.1 Control of Process Related Impurities
#
# ⚠️ Options는 네 개뿐이고 전부 «관리 전략»입니다.
#    평가 방법이 아니며 Option 5는 존재하지 않습니다.
#    예측 플랫폼이 자동 배정할 수 있는 것이 아니므로, 보고서에서는
#    Option을 배정하지 말고 Class별 권고 조치만 적으세요.
# ---------------------------------------------------------------------------

CONTROL_OPTIONS: dict[int, str] = {
    1: (
        "Include a test for the impurity in the drug substance specification "
        "with an acceptance criterion at or below the acceptable limit using "
        "an appropriate analytical procedure."
    ),
    2: (
        "Include a test for the impurity in the specification for a raw "
        "material, starting material or intermediate, or as an in-process "
        "control, with an acceptance criterion at or below the acceptable "
        "limit using an appropriate analytical procedure."
    ),
    3: (
        "Include a test for the impurity in the specification for a raw "
        "material, starting material or intermediate, or as an in-process "
        "control, with an acceptance criterion above the acceptable limit of "
        "the impurity in the drug substance, coupled with demonstrated "
        "understanding of fate and purge and associated process controls."
    ),
    4: (
        "Understand process parameters and impact on residual impurity levels "
        "(including fate and purge knowledge) with sufficient confidence that "
        "the level of the impurity in the drug substance will be below the "
        "acceptable limit such that no analytical testing is recommended."
    ),
}

CONTROL_OPTIONS_SHORT: dict[int, str] = {
    1: "Test in drug substance specification",
    2: "Test in upstream specification or in-process control",
    3: "Upstream test above limit + fate and purge understanding",
    4: "Process understanding only — no analytical testing",
}

CONTROL_SECTION_TITLE = "Control of Process Related Impurities"


# ---------------------------------------------------------------------------
# 6. Hazard assessment — (Q)SAR
# ---------------------------------------------------------------------------

QSAR_TERMS: dict[str, str] = {
    "methodology": "(Q)SAR prediction methodology",
    "expert": "expert rule-based",
    "statistical": "statistical-based",
    "complementary": "complementary",
    "validation": "OECD validation principles",
    "expert_review": "expert knowledge",
    "assay": "bacterial mutagenicity assay",
    "assay_display": "Bacterial mutagenicity (Ames)",
    "computational": "computational toxicology assessment",
}

# 이 플랫폼의 4층 구조가 가이드라인 어디에 대응하는지
#   ⚠️ Layer A는 (Q)SAR가 아닙니다. 6절 첫 문장의
#      "database and literature searches"에 해당합니다.
LAYER_MAPPING: dict[str, dict[str, str]] = {
    "A": {
        "name": "Empirical database lookup",
        "name_ko": "실측 데이터베이스 조회",
        "guideline": "database and literature searches (Section 6)",
        "is_qsar": "false",
    },
    "B": {
        "name": "Expert rule-based alert scan",
        "name_ko": "전문가 규칙 기반 구조 경보",
        "guideline": "expert rule-based (Q)SAR methodology (Section 6)",
        "is_qsar": "true",
    },
    "C": {
        "name": "Statistical (Q)SAR model",
        "name_ko": "통계 기반 (Q)SAR 모델",
        "guideline": "statistical-based (Q)SAR methodology (Section 6)",
        "is_qsar": "true",
    },
    # 2026-10-06 추가: Layer D 가 빠져 있었다. B 와 C 가 상충할 때 발화하는
    #   read-across 중재 계층이며, ICH M7 6절의 "expert knowledge" 와
    #   7.3절의 전문가 검토 요건에 대응한다. (Q)SAR 방법론은 아니다.
    "D": {
        "name": "Read-across arbitration",
        "name_ko": "읽기교차 중재",
        "guideline": "expert knowledge / expert review (Sections 6 and 7.3)",
        "is_qsar": "false",
    },
}

METHODOLOGY_STATEMENT = (
    "Database search for carcinogenicity and bacterial mutagenicity data, "
    "followed by two complementary (Q)SAR methodologies "
    "(expert rule-based and statistical-based)."
)


# ---------------------------------------------------------------------------
# 7.5 Cohort of concern
# ---------------------------------------------------------------------------

COHORT_OF_CONCERN = ("aflatoxin-like", "N-nitroso", "alkyl-azoxy")

COHORT_STATEMENT = (
    'High potency mutagenic carcinogens referred to as the "cohort of '
    'concern", comprising aflatoxin-like-, N-nitroso-, and alkyl-azoxy '
    "compounds."
)

# 2026-09-03 수정: CoC는 §7.5에서 TTC 적용 대상에서 명시적으로 제외되므로
# TTC 문자열 대신 이 문장을 쓰도록 단일 출처 상수로 추가했습니다.
COHORT_OF_CONCERN_AI_TEXT = (
    "Cohort of concern — TTC not applicable; compound-specific Acceptable "
    "Intake (AI) required (ICH M7(R2) §7.5)"
)

COHORT_OF_CONCERN_AI_TEXT_KO = (
    "Cohort of concern — TTC 적용 불가. 화합물별 허용섭취량(AI) 산출 필요 "
    "(ICH M7(R2) §7.5)"
)


# ---------------------------------------------------------------------------
# 7. Risk characterization — 용량 · 한도
# ---------------------------------------------------------------------------

TTC_DEFAULT_UG_DAY = 1.5

TTC_STATEMENT = (
    "A TTC-based acceptable intake of a mutagenic impurity of 1.5 µg per "
    "person per day is considered to be associated with a negligible risk "
    "(theoretical excess cancer risk of <1 in 100,000 over a lifetime of "
    "exposure)."
)

# 2026-09-03 수정: DOSE_TERMS["ltl"]에 용어만 있고 Table 2의 실제 한도가
# 코드에 없어 LTL 산출이 불가능했습니다. 원문 표를 상수로 추가합니다.
#
# ICH M7(R2) Table 2 — Acceptable intakes for an individual impurity
#   (less-than-lifetime exposure, single mutagenic impurity)
#     ≤1개월          120 µg/day
#     >1–12개월        20 µg/day
#     >1–10년          10 µg/day
#     >10년 ~ 평생      1.5 µg/day
#
# (upper_bound_months, µg/day, 영문 구간명, 한국어 구간명)
# upper_bound_months = None 은 "평생"(상한 없음)을 뜻합니다.
LTL_LIMIT_TABLE: tuple[tuple[float | None, float, str, str], ...] = (
    (1.0, 120.0, "≤1 month", "1개월 이하"),
    (12.0, 20.0, ">1–12 months", "1개월 초과 12개월 이하"),
    (120.0, 10.0, ">1–10 years", "1년 초과 10년 이하"),
    (None, 1.5, ">10 years to lifetime", "10년 초과 ~ 평생"),
)

# 구간명 → 한도(µg/day) 조회용
LTL_LIMITS_UG_DAY: dict[str, float] = {row[2]: row[1] for row in LTL_LIMIT_TABLE}

LTL_STATEMENT = (
    "For less-than-lifetime (LTL) exposure, higher acceptable intakes than the "
    "lifetime TTC of 1.5 µg/day are justified based on the duration of "
    "treatment (ICH M7(R2), Section 7.3 and Table 2)."
)

DOSE_TERMS: dict[str, str] = {
    "ttc": "Threshold of Toxicological Concern (TTC)",
    "ai": "acceptable intake",
    "csai": "compound-specific acceptable intake (CSAI)",
    "limit": "acceptable limit",
    "pde": "Permissible Daily Exposure (PDE)",
    "noel": "No-Observed Effect Level (NOEL)",
    "ltl": "Less-Than-Lifetime (LTL)",
}


# ---------------------------------------------------------------------------
# 5. 불순물 종류
# ---------------------------------------------------------------------------

SUBSTANCE_TERMS: dict[str, str] = {
    "drug_substance": "drug substance",
    "drug_product": "drug product",
    "actual": "actual impurity",
    "potential": "potential impurity",
    "process_related": "process related impurity",
    "degradation": "degradation product",
    "starting_material": "starting material",
    "intermediate": "intermediate",
    "reagent": "reagent",
    "raw_material": "raw material",
    "excipient": "excipient",
}

ALERT_TERMS: dict[str, str] = {
    "alert": "structural alert",
    "alerting": "alerting structure",
    "none": "no structural alerts",
    "same": "same alert",
    "shared": "shared alert",
}

# NDSRI — ICH M7 본문에는 없고 EMA·FDA 니트로사민 가이던스 용어입니다.
NDSRI_FULL = "Nitrosamine Drug Substance-Related Impurity (NDSRI)"


# ---------------------------------------------------------------------------
# 헬퍼
# ---------------------------------------------------------------------------


def class_label(class_number: Any, *, lang: str = "en", short: bool = True) -> str:
    """`Class 3 — Alerting structure, no mutagenicity data` 형태 문자열."""
    try:
        n = int(class_number)
    except (TypeError, ValueError):
        return "Not classified" if lang == "en" else "미판정"
    if n not in CLASS_DEFINITIONS:
        return "Not classified" if lang == "en" else "미판정"

    if lang == "ko":
        table = CLASS_SHORT_KO if short else CLASS_DEFINITIONS_KO
    else:
        table = CLASS_SHORT if short else CLASS_DEFINITIONS
    return f"Class {n} — {table[n]}"


def class_definition(class_number: Any, *, lang: str = "en") -> str:
    """Table 1의 Definition 원문."""
    try:
        n = int(class_number)
    except (TypeError, ValueError):
        return ""
    table = CLASS_DEFINITIONS_KO if lang == "ko" else CLASS_DEFINITIONS
    return table.get(n, "")


def class_action(class_number: Any, *, lang: str = "en") -> str:
    """Table 1의 Proposed action for control 원문."""
    try:
        n = int(class_number)
    except (TypeError, ValueError):
        return ""
    table = CLASS_ACTIONS_KO if lang == "ko" else CLASS_ACTIONS
    return table.get(n, "")


# 2026-09-03 수정: 노출 기간(개월)으로 Table 2 한도를 얻는 함수가 없어
# LTL이 코드에서 쓰이지 못했습니다. 다른 모듈이 import해 쓰도록 추가합니다.
def ltl_limit_ug_day(duration_months: Any) -> float:
    """노출 기간(개월) → ICH M7(R2) Table 2 허용섭취량(µg/day).

    기간을 해석할 수 없거나 0 이하면 가장 보수적인 평생 한도(1.5 µg/day)를
    반환합니다. 경계값은 "이하" 포함입니다(예: 12개월 → 20 µg/day).
    """
    try:
        months = float(duration_months)
    except (TypeError, ValueError):
        return TTC_DEFAULT_UG_DAY
    if months <= 0:
        return TTC_DEFAULT_UG_DAY
    for upper, limit, _en, _ko in LTL_LIMIT_TABLE:
        if upper is None or months <= upper:
            return limit
    return TTC_DEFAULT_UG_DAY


def ltl_limit_label(duration_months: Any, *, lang: str = "en") -> str:
    """`20 µg/day (>1–12 months)` 형태의 표시용 문자열."""
    try:
        months = float(duration_months)
    except (TypeError, ValueError):
        months = -1.0
    for upper, limit, en, ko in LTL_LIMIT_TABLE:
        if months > 0 and (upper is None or months <= upper):
            name = ko if lang == "ko" else en
            return f"{limit:g} µg/day ({name})"
    name = LTL_LIMIT_TABLE[-1][3] if lang == "ko" else LTL_LIMIT_TABLE[-1][2]
    return f"{TTC_DEFAULT_UG_DAY:g} µg/day ({name})"


def is_cohort_of_concern_class(alerts: list[str] | None) -> bool:
    """구조 경보 목록에 cohort of concern 구조가 있는지."""
    text = " ".join(alerts or []).lower()
    return any(
        marker in text
        for marker in ("nitroso", "nitrosamine", "aflatoxin", "azoxy")
    )
