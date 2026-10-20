"""
가혹 조건 기반 분해산물(Degradant) 예측 — RDKit 고급 화학 엔진 + DB RAG 교차 검증.

Mock 데이터셋 없이 RDKit Substructure Matching·반응 SMARTS·ExactMolWt 연산과
PRIME-DS 1만 건+ 내부 검증 DB Tanimoto 매트릭스로 Confidence Score를 산출합니다.
"""

# 정리 + 등급 연결

# 분자량 큰 순 → 신뢰도 높은 순으로 정렬
# 각 산물에 ICH M7 등급 붙임 (p66 호출)
# 구조 경보 검사 (p62 호출)
# 화면·보고서용 표로 변환

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from rdkit import Chem

from p42_degradation_rag_engine import (
    _canonical_smiles,
    _generate_diverse_fallback_hits,
    compute_rag_confidence,
    deduplicate_hits,
    exact_mw_from_smiles,
    filter_valid_mechanism_hits,
    generate_mechanism_hits,
    parse_parent_structure,
)
from p44_degradation_shelf_life_filter import ShelfLifeChildFilter
from p66_ich_m7 import ICH_M7_CLASS_INFO, classify_ich_m7
from p62_toxicity_alerts import scan_ashby_tennant_alerts

# 2026-09-03 수정: ICH M7 등급 라벨·조치 문자열을 p16(단일 출처)에서 가져옵니다.
from p16_ich_m7_terms import class_action, class_label

# ── 가혹 조건 선택지 (UI selectbox와 1:1 매핑) ──────────────────────────────
STRESS_ENV_CHOICES: tuple[str, ...] = (
    "선택 안함",
    "산성 가혹 가수분해 (Acidic Hydrolysis)",
    "염기성 가혹 가수분해 (Alkaline Hydrolysis)",
    "수분 노출 (Neutral Hydrolysis)",
    "과산화 산화 조건 (Oxidative Stress)",
    "광분해 조건 (Photolytic Degradation)",
)

from p02_degradation_constants import ENGINE_CACHE_VERSION

@dataclass
class DegradantProduct:
    """예측 분해산물 1건 — 리포트·UI 공통 데이터 모델."""

    degradant_id: str
    pathway_name: str
    name: str
    smiles: str
    molecular_weight: float
    ich_m7_class: int
    ich_m7_label: str
    signal_tier: str  # danger | caution | safe
    signal_label_ko: str
    mechanism: str
    biotransformer_pathway: str
    relative_abundance: str
    confidence_score: float = 0.0
    confidence_label: str = ""
    db_match_detail: str = ""
    alerts: list[str] = field(default_factory=list)

    def to_report_row(self) -> dict[str, Any]:
        return {
            "ID": self.degradant_id,
            "분해 경로": self.pathway_name,
            "예측 분해산물 구조 (SMILES)": self.smiles,
            "실시간 연산 분자량 (g/mol)": self.molecular_weight,
            # 2026-09-03 수정: ich_m7_label 자체가 이미 "Class N — …" 형태여서
            # 등급이 두 번 찍혔습니다. p16 class_label()로 단일화.
            "ICH M7 독성 등급": class_label(self.ich_m7_class, lang="ko"),
            "예측 신뢰도 및 DB 매칭률 (%)": f"{self.confidence_score:.1f}% — {self.confidence_label}",
            "예측 분해산물명": self.name,
            "위험 신호": self.signal_label_ko,
            "예상 비율": self.relative_abundance,
            "반응 메커니즘": self.mechanism,
            "DB RAG 매칭": self.db_match_detail,
            "Ashby-Tennant 알러트": ", ".join(self.alerts) if self.alerts else "없음",
        }

    def to_report_row_en(self) -> dict[str, Any]:
        from p80_formal_assessment_report import ICH_M7_REGULATORY

        m7 = ICH_M7_REGULATORY.get(self.ich_m7_class, ICH_M7_REGULATORY[5])
        return {
            "ID": self.degradant_id,
            "Degradation pathway": translate_degradant_text(self.pathway_name),
            "Predicted degradant": translate_degradant_text(self.name),
            "SMILES": self.smiles,
            "Molecular weight (g/mol)": self.molecular_weight,
            "ICH M7 Class": self.ich_m7_class,
            "ICH M7 classification": m7["name"],
            "Confidence score (%)": self.confidence_score,
            "Confidence tier": self.confidence_label,
            "DB RAG match": self.db_match_detail,
            "Mutagenicity signal": ich_m7_signal_label_en(self.ich_m7_class),
            "Predicted abundance": abundance_en(self.relative_abundance),
            "Reaction mechanism": mechanism_en(self.mechanism),
            "BioTransformer route": self.biotransformer_pathway,
            "Ashby-Tennant alerts": ", ".join(self.alerts) if self.alerts else "None detected",
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "degradant_id": self.degradant_id,
            "pathway_name": self.pathway_name,
            "name": self.name,
            "smiles": self.smiles,
            "molecular_weight": self.molecular_weight,
            "ich_m7_class": self.ich_m7_class,
            "ich_m7_label": self.ich_m7_label,
            "signal_tier": self.signal_tier,
            "signal_label_ko": self.signal_label_ko,
            "mechanism": self.mechanism,
            "biotransformer_pathway": self.biotransformer_pathway,
            "relative_abundance": self.relative_abundance,
            "confidence_score": self.confidence_score,
            "confidence_label": self.confidence_label,
            "db_match_detail": self.db_match_detail,
            "alerts": list(self.alerts),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DegradantProduct:
        return cls(
            degradant_id=str(data["degradant_id"]),
            pathway_name=str(data["pathway_name"]),
            name=str(data["name"]),
            smiles=str(data["smiles"]),
            molecular_weight=float(data["molecular_weight"]),
            ich_m7_class=int(data["ich_m7_class"]),
            ich_m7_label=str(data["ich_m7_label"]),
            signal_tier=str(data["signal_tier"]),
            signal_label_ko=str(data["signal_label_ko"]),
            mechanism=str(data["mechanism"]),
            biotransformer_pathway=str(data["biotransformer_pathway"]),
            relative_abundance=str(data["relative_abundance"]),
            confidence_score=float(data.get("confidence_score", 0)),
            confidence_label=str(data.get("confidence_label", "")),
            db_match_detail=str(data.get("db_match_detail", "")),
            alerts=list(data.get("alerts") or []),
        )


def ich_m7_signal_tier(class_number: int) -> str:
    # 2026-09-03 수정: ICH M7 Table 1에서 Class 4의 조치는 Class 5와 동일한
    # "Treat as non-mutagenic impurity"(p16 CLASS_ACTIONS)입니다. 화면·영문
    # 보고서에만 "주의(잠재적 위험)"로 나가 규정과 어긋나 있어 safe로 옮깁니다.
    # (caution tier는 다른 호출부 호환을 위해 라벨 표에 남겨 둡니다.)
    if class_number <= 3:
        return "danger"
    return "safe"


def ich_m7_signal_label_ko(class_number: int) -> str:
    tier = ich_m7_signal_tier(class_number)
    # 2026-09-03 수정: Class 4 라벨을 p16 class_action("비변이원성 불순물로 취급")에 맞춤.
    if class_number == 4:
        return f"🟢 안전 ({class_action(4, lang='ko')} — 원료와 동일 알러트)"
    labels = {
        "danger": "🔴 위험 (유전독성 경고 구조 포함)",
        "caution": "🟡 주의 (잠재적 위험)",
        "safe": "🟢 안전 (독성 위험 낮음)",
    }
    return labels[tier]


def ich_m7_signal_label_en(class_number: int) -> str:
    tier = ich_m7_signal_tier(class_number)
    # 2026-09-03 수정: "MODERATE — Potential concern"은 p16의 Class 4 조치와
    # 모순이므로 가이드라인 조치 문구(Treat as non-mutagenic impurity)로 교체.
    if class_number == 4:
        return (
            f"LOW RISK — {class_action(4)} "
            "(alert shared with the drug substance; ICH Q3A/Q3B limits apply)"
        )
    labels = {
        "danger": "HIGH RISK — Mutagenicity alerting structure identified",
        "caution": "MODERATE — Potential concern (drug-related context)",
        "safe": "LOW RISK — No significant mutagenicity concern",
    }
    return labels[tier]


_DEGRADANT_TEXT_EN: dict[str, str] = {
    "에스테르 절단 — 카르복실산 잔기": "Ester cleavage — carboxylic acid fragment",
    "카르복실산 잔기": "Carboxylic acid fragment",
    "아민/아마인 잔기": "Amine fragment",
    "알코올/페놀릭 잔기": "Alcohol / phenolic fragment",
    "탈할로겐·탈수소 부산물 — 방향족 할로겐 제거체": "Dehalogenation byproduct — aryl halogen removed",
    "탈할로겐·탈수소 부산물": "Dehalogenation / elimination byproduct",
    "N-탈알킬화 2차 아민 대사체": "N-demethylated secondary amine metabolite",
    "N-옥사이드 대사체": "N-oxide metabolite",
    "에폭사이드 산화 분해체": "Epoxide oxidative degradant",
    "설폭시드/아울폭사이드형 산화체": "Sulfoxide / sulfone-type oxidant",
    "다성분 하이드록실화 잔기": "Poly-hydroxylated oxidative residue",
    "카르복실산 산화 분해체": "Carboxylic acid oxidation product",
    "아릴아민 산화 전구체": "Arylamine oxidative precursor",
    "🔥 실제 검출 이력과 구조 90% 이상 일치 (검증 완료)": "High confidence — literature/DB structural match verified (≥90%)",
    "⚠️ 화학 메커니즘 및 유사 작용기 유효성 확인": "Mechanistically plausible — similar functional group validated in DB",
    "💡 정보화학 모델 기반 가상 시뮬레이션 결과": "In-silico simulation — cheminformatics model prediction",
    "🔥 실제 검출 이력과 구조 90% 이상 일치 (신뢰도 높음)": "High confidence — literature/DB structural match (≥90%)",
    "⚠️ 이론적 메커니즘 일치 및 유사 작용기 DB 확인": "Mechanistically plausible — similar functional group in DB",
    "💡 이론적 예측 최적 모델 기반 도출": "Theoretical prediction — optimized in-silico model",
    "🔥 실제 검출 이력 일치 (매우 확실)": "High confidence — literature/DB detection match",
    "⚠️ 메커니즘적 유효": "Mechanistically plausible",
    "💡 이론적 가능성": "Theoretical possibility",
}


def translate_degradant_text(text: str) -> str:
    text = (text or "").strip()
    if not text:
        return "N/A"
    if text in _DEGRADANT_TEXT_EN:
        return _DEGRADANT_TEXT_EN[text]
    for match in re.finditer(r"\(([^)]+)\)", text):
        inner = match.group(1).strip()
        if re.search(r"[A-Za-z]", inner) and len(inner) > 2:
            return inner
    return text


def stress_environment_en(stress: str) -> str:
    text = (stress or "").strip()
    if not text or text == "선택 안함":
        return "Not selected"
    return translate_degradant_text(text)


def abundance_en(text: str) -> str:
    return translate_degradant_text(text)


def temp_stress_label_en(temp_c: float) -> str:
    if temp_c >= 60:
        return "Accelerated (high temperature)"
    if temp_c >= 40:
        return "Intermediate stress"
    return "Ambient / room temperature"


def mechanism_en(mechanism: str) -> str:
    s = mechanism or ""
    s = s.replace("·", " | ")
    for ko, en in _DEGRADANT_TEXT_EN.items():
        s = s.replace(ko, en)
    for stress in STRESS_ENV_CHOICES:
        if stress in s:
            s = s.replace(stress, stress_environment_en(stress))
    return s


def _temp_stress_label(temp_c: float) -> str:
    if temp_c >= 60:
        return "고온 가혹 (Accelerated)"
    if temp_c >= 40:
        return "중간 가혹 (Intermediate)"
    return "저온/실온 (Ambient)"


def sanitize_degradant_products(
    products: list[DegradantProduct],
    parent_canonical: str,
) -> list[DegradantProduct]:
    """구 Mock(DP-)·Parent 복사·탈할로겐 미적용·라디칼 중간체 결과를 최종 단계에서 제거."""
    from p42_degradation_rag_engine import (
        _canonical_smiles,
        _halogen_count,
        _has_aryl_halide_smiles,
    )
    from p44_degradation_shelf_life_filter import is_spin_stabilized_neutral_smiles

    parent = _canonical_smiles(parent_canonical)
    parent_h = _halogen_count(parent)
    parent_mw = exact_mw_from_smiles(parent) if parent else 0.0
    _transform_kw = ("탈할로", "탈수소", "탈알킬", "N-옥사이드", "에폭", "β-제거", "E2")
    out: list[DegradantProduct] = []
    for p in products:
        did = str(p.degradant_id or "")
        if did.startswith("DP-"):
            continue
        canon = _canonical_smiles(p.smiles)
        if not canon or canon == parent:
            continue
        if not is_spin_stabilized_neutral_smiles(canon):
            continue
        code = p.biotransformer_pathway or ""
        if "탈할로" in p.name or "DEHAL" in code:
            if parent_h > 0 and _halogen_count(canon) >= parent_h:
                continue
            if parent_h > 0 and _has_aryl_halide_smiles(canon):
                continue
        if parent_mw > 0 and abs(p.molecular_weight - parent_mw) < 0.05:
            if any(k in p.name for k in _transform_kw):
                continue
        out.append(p)
    return out


# 2026-09-03 추가: 결합 절단(cleavage) 기전 코드 — 생성물은 Parent보다 반드시 가볍습니다.
# 락탐 개환(LAC_HYD)은 절단이 아니라 물 첨가 개환이므로 제외합니다(Parent보다 무거움).
_CLEAVAGE_MECHANISM_CODES = frozenset({
    "EST_HYD", "AMD_HYD", "CAR_HYD",
    "SMIRKS_EST", "SMIRKS_AMD", "SMIRKS_CARB",
    "FB_EST", "FB_AMD",
})

# 개환 반응 표기 — 이름에 "잔기"가 들어가도 절단 산물이 아닙니다.
_RING_OPENING_KEYWORDS = ("개환", "Ring Opening", "ring opening")


def _is_cleavage_fragment(product: DegradantProduct) -> bool:
    """절단 기전으로 생성된 «모핵 잔기 / 파편»인지 판정."""
    name = product.name or ""
    if any(k in name for k in _RING_OPENING_KEYWORDS):
        return False
    # biotransformer_pathway = "PRIME-DS Smart Engine — {mechanism_code}"
    code = (product.biotransformer_pathway or "").rsplit("—", 1)[-1].strip()
    if code in _CLEAVAGE_MECHANISM_CODES:
        return True
    return ("잔기" in name) or ("파편" in name) or ("절단" in name)


def _finalize_degradant_products(
    products: list[DegradantProduct],
    parent_canonical: str,
    parent_exact_mw: float,
) -> list[DegradantProduct]:
    """
    최종 검증 — Parent 동일 SMILES 제거 + ExactMolWt 재연산.
    절단 기전으로 생성된 모핵 잔기·파편의 MW는 Parent보다 반드시 감소해야 합니다
    (개환·산화 부가처럼 질량이 늘어나는 기전은 이 검사에서 제외).
    """
    parent = _canonical_smiles(parent_canonical)
    out: list[DegradantProduct] = []
    from p44_degradation_shelf_life_filter import is_spin_stabilized_neutral_smiles

    for p in products:
        canon = _canonical_smiles(p.smiles)
        if not canon or canon == parent:
            continue
        if not is_spin_stabilized_neutral_smiles(canon):
            continue
        mw = exact_mw_from_smiles(canon)
        # 2026-09-03 수정: 310행에서 이미 `canon == parent`를 걸러내므로 여기 있던
        # `and canon == parent` 조건은 절대 참이 될 수 없는 죽은 코드였고, 그 결과
        # docstring이 선언한 "잔기 MW 감소" 검사가 실제로는 없었습니다. 이제
        # 절단 산물이 Parent보다 가볍지 않으면 실제로 버립니다.
        if parent_exact_mw > 0 and _is_cleavage_fragment(p) and mw >= parent_exact_mw - 0.02:
            continue
        out.append(
            DegradantProduct(
                degradant_id=p.degradant_id,
                pathway_name=p.pathway_name,
                name=p.name,
                smiles=canon,
                molecular_weight=mw,
                ich_m7_class=p.ich_m7_class,
                ich_m7_label=p.ich_m7_label,
                signal_tier=p.signal_tier,
                signal_label_ko=p.signal_label_ko,
                mechanism=p.mechanism,
                biotransformer_pathway=p.biotransformer_pathway,
                relative_abundance=p.relative_abundance,
                confidence_score=p.confidence_score,
                confidence_label=p.confidence_label,
                db_match_detail=p.db_match_detail,
                alerts=list(p.alerts),
            )
        )
    return out


def _hit_to_product(
    hit,
    *,
    stress_environment: str,
    temperature_c: float,
    parent_analysis,
) -> DegradantProduct:
    """MechanismHit → DegradantProduct (ICH M7 + RAG confidence + ExactMolWt)."""
    smiles = _canonical_smiles(hit.smiles) or hit.smiles
    mol = Chem.MolFromSmiles(smiles)
    alerts, _ = scan_ashby_tennant_alerts(mol) if mol else ([], [])
    if mol:
        m7 = classify_ich_m7(mol, alerts)
        ich_class = m7.class_number
        ich_label = m7.class_name_ko
    else:
        ich_class = 5
        ich_label = ICH_M7_CLASS_INFO[5]["name_ko"]

    conf_score, conf_label, db_detail = compute_rag_confidence(
        smiles,
        mechanism_base_score=hit.mechanism_base_score,
        stress_environment=stress_environment,
        temperature_c=temperature_c,
        mechanism_code=hit.mechanism_code,
    )

    mechanism = (
        f"{stress_environment_en(stress_environment)} | {_temp_stress_label(temperature_c)} "
        f"({temperature_c:.0f}°C) | Parent MW {parent_analysis.molecular_weight} g/mol "
        f"(Exact {parent_analysis.exact_molecular_weight}) | {hit.mechanism_detail}"
    )

    return DegradantProduct(
        degradant_id=hit.degradant_id,
        pathway_name=hit.pathway_name,
        name=hit.name,
        smiles=smiles,
        molecular_weight=exact_mw_from_smiles(smiles),
        ich_m7_class=ich_class,
        ich_m7_label=ich_label,
        signal_tier=ich_m7_signal_tier(ich_class),
        signal_label_ko=ich_m7_signal_label_ko(ich_class),
        mechanism=mechanism,
        biotransformer_pathway=f"PRIME-DS Smart Engine — {hit.mechanism_code}",
        relative_abundance=hit.relative_abundance,
        confidence_score=conf_score,
        confidence_label=conf_label,
        db_match_detail=db_detail,
        alerts=alerts,
    )


@dataclass
class DegradationPredictionResult:
    """predict_advanced_degradants 반환값 — 예측 목록 + 검증 메타데이터."""

    products: list[DegradantProduct] = field(default_factory=list)
    error_message: str | None = None
    parent_smiles: str = ""
    parent_canonical_smiles: str = ""
    parent_exact_mw: float = 0.0
    parent_formula: str = ""
    stress_condition: str = ""
    temperature_c: float = 25.0

    @property
    def ok(self) -> bool:
        return self.error_message is None and bool(self.products)


def validate_parent_smiles(input_smiles: str) -> tuple[Chem.Mol | None, str | None]:
    """
    Chem.MolFromSmiles + SanitizeMol로 Parent SMILES 검증.
    None 반환 시 UI에 한국어 오류 메시지를 표시합니다.
    """
    raw = (input_smiles or "").strip()
    if not raw:
        return None, (
            "Parent SMILES가 비어 있습니다.\n"
            "Ketcher 캔버스에 구조를 그리거나 SMILES를 입력한 뒤 **독성 스캔 시작**을 눌러 주세요."
        )

    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        return None, (
            f"**RDKit 구조 해석 실패** — 유효하지 않은 SMILES입니다.\n\n"
            f"`{raw[:80]}{'…' if len(raw) > 80 else ''}`\n\n"
            "구조식을 다시 그리거나 SMILES 형식·원자가를 확인해 주세요."
        )
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        return None, (
            f"**구조 정제(SanitizeMol) 실패** — 화학적으로 허용되지 않는 구조입니다.\n\n"
            f"오류: {exc}"
        )
    return mol, None


# 가수분해로 생성된 산물은 방향족 수산화 검사에서 면제한다.
#
# 배경: ShelfLifeChildFilter 는 "부모에 없던 방향족 -OH 가 생성물에 나타나면
# 생체 대사(방향족 수산화) 산물로 보고 배제"한다. 보관 조건에서는 일어나지 않는
# 반응이므로 타당한 규칙이다.
#
# 그러나 아릴 에스터가 가수분해되면 원래 분자에 있던 산소가 -OH 로 드러날 뿐인데,
# 이 역시 "새로 생긴 방향족 OH"로 판정되어 함께 배제된다.
#
#   아스피린  CH3-C(=O)-O-c1ccccc1-COOH   (산소가 에스터로 가려져 있음)
#        ↓ 에스터 가수분해
#   살리실산  HO-c1ccccc1-COOH            (같은 산소가 -OH 로 노출)
#
# 수산화(c-H → c-OH)가 아니라 가수분해이므로 배제 대상이 아니다.
# 따라서 아래 기전 코드로 생성된 산물은 이 검사를 건너뛴다.
_HYDROLYSIS_MECHANISM_CODES = frozenset({
    "EST_HYD", "AMD_HYD", "CAR_HYD", "LAC_HYD",
    "SMIRKS_EST", "SMIRKS_AMD", "SMIRKS_CARB",
    "FB_EST", "FB_AMD",
})


def _cull_shelf_life_ineligible_hits(hits, parent_mol, parent_canonical: str):
    """유연물질 리스트 — 라디칼·방향족 수산화(`[OX2H]-[c;r6]`) Child 배제.

    단, 가수분해 기전으로 생성된 산물은 면제한다(위 주석 참조).
    """
    if not hits:
        return hits
    kept = []
    for hit in hits:
        # 가수분해 산물은 방향족 수산화 검사 면제
        if getattr(hit, "mechanism_code", "") in _HYDROLYSIS_MECHANISM_CODES:
            kept.append(hit)
            continue

        child_smi = _canonical_smiles(hit.smiles) or hit.smiles
        child_mol = Chem.MolFromSmiles(child_smi) if child_smi else None
        try:
            if child_mol is not None:
                Chem.SanitizeMol(child_mol)
        except Exception:
            child_mol = None
        try:
            if ShelfLifeChildFilter.evaluate_mols(
                parent_mol,
                child_mol,
                parent_smiles=parent_canonical,
                child_smiles=child_smi,
            ).eligible:
                kept.append(hit)
        except Exception:
            kept.append(hit)
    return kept


# 2026-09-03 추가: ENGINE_CACHE_VERSION이 import만 되고 어디에도 쓰이지 않아
# "엔진 로직 변경 시 캐시 자동 무효화"라는 p02 주석의 약속이 지켜지지 않았습니다.
# 예측 결과 캐시를 두고 그 키의 첫 요소로 버전을 넣어, 버전을 올리면 이전 캐시가
# 자동으로 무효화되도록 합니다.
_PREDICTION_CACHE: dict[tuple[str, str, str, float], DegradationPredictionResult] = {}
_PREDICTION_CACHE_MAX = 256


def degradation_cache_key(
    input_smiles: str,
    stress_condition: str,
    temperature: float,
) -> tuple[str, str, str, float]:
    """엔진 캐시 키 — 첫 요소가 ENGINE_CACHE_VERSION입니다."""
    return (
        ENGINE_CACHE_VERSION,
        (input_smiles or "").strip(),
        stress_condition or "",
        round(float(temperature), 2),
    )


def clear_degradation_cache() -> None:
    """예측 캐시 비우기 (테스트·엔진 교체용)."""
    _PREDICTION_CACHE.clear()


def _copy_prediction_result(src: DegradationPredictionResult) -> DegradationPredictionResult:
    """캐시 원본이 호출부에서 변형되지 않도록 얕은 복사본을 돌려줍니다."""
    return DegradationPredictionResult(
        products=list(src.products),
        error_message=src.error_message,
        parent_smiles=src.parent_smiles,
        parent_canonical_smiles=src.parent_canonical_smiles,
        parent_exact_mw=src.parent_exact_mw,
        parent_formula=src.parent_formula,
        stress_condition=src.stress_condition,
        temperature_c=src.temperature_c,
    )


def predict_advanced_degradants(
    input_smiles: str,
    stress_condition: str,
    temperature: float = 25.0,
) -> DegradationPredictionResult:
    """
    고급 분해산물 예측 엔진 — RDKit 물성 연산 + 가혹 조건 메커니즘 + DB RAG 신뢰도.

    동일 입력은 ENGINE_CACHE_VERSION이 포함된 키로 캐시됩니다.

    Parameters
    ----------
    input_smiles : Parent 화합물 SMILES
    stress_condition : 가혹 분해 환경 (STRESS_ENV_CHOICES 중 하나)
    temperature : 시험 온도 (°C) — 높을수록 DB 매칭 신뢰도 가중 상승

    Returns
    -------
    DegradationPredictionResult
        products: 예측 분해산물 목록 (ICH M7·Confidence 포함)
        error_message: SMILES 오류 등 예측 불가 시 한국어 안내
    """
    result = DegradationPredictionResult(
        parent_smiles=(input_smiles or "").strip(),
        stress_condition=stress_condition,
        temperature_c=float(temperature),
    )

    if stress_condition == "선택 안함":
        return result

    # 2026-09-03 수정: ENGINE_CACHE_VERSION을 실제 캐시 키에 반영.
    cache_key = degradation_cache_key(input_smiles, stress_condition, temperature)
    cached = _PREDICTION_CACHE.get(cache_key)
    if cached is not None:
        return _copy_prediction_result(cached)

    parent_mol, analysis = parse_parent_structure(input_smiles)
    if parent_mol is None or analysis.parse_error:
        result.error_message = analysis.parse_error or "Parent SMILES 파싱에 실패했습니다."
        return result

    result.parent_canonical_smiles = analysis.canonical_smiles
    result.parent_exact_mw = analysis.exact_molecular_weight
    result.parent_formula = analysis.molecular_formula

    hits = filter_valid_mechanism_hits(
        deduplicate_hits(
            generate_mechanism_hits(parent_mol, analysis, stress_condition),
            parent_canonical=analysis.canonical_smiles,
        ),
        analysis.canonical_smiles,
    )
    hits = _cull_shelf_life_ineligible_hits(
        hits, parent_mol, analysis.canonical_smiles
    )

    products = [
        _hit_to_product(
            h,
            stress_environment=stress_condition,
            temperature_c=float(temperature),
            parent_analysis=analysis,
        )
        for h in hits
    ]

    products.sort(key=lambda p: (-p.molecular_weight, -p.confidence_score))
    if not products and parent_mol is not None:
        prefix_map = {
            "산성 가혹 가수분해 (Acidic Hydrolysis)": "SD-AH",
            "염기성 가혹 가수분해 (Alkaline Hydrolysis)": "SD-AL",
            "수분 노출 (Neutral Hydrolysis)": "SD-NH",
            "과산화 산화 조건 (Oxidative Stress)": "SD-OX",
            "광분해 조건 (Photolytic Degradation)": "SD-PH",
        }
        fb_hits = filter_valid_mechanism_hits(
            deduplicate_hits(
                _generate_diverse_fallback_hits(
                    parent_mol,
                    analysis,
                    stress_condition,
                    id_prefix=prefix_map.get(stress_condition, "SD-XX"),
                ),
                parent_canonical=analysis.canonical_smiles,
            ),
            analysis.canonical_smiles,
        )
        fb_hits = _cull_shelf_life_ineligible_hits(
            fb_hits, parent_mol, analysis.canonical_smiles
        )
        products = [
            _hit_to_product(
                h,
                stress_environment=stress_condition,
                temperature_c=float(temperature),
                parent_analysis=analysis,
            )
            for h in fb_hits
        ]
    result.products = _finalize_degradant_products(
        sanitize_degradant_products(products, analysis.canonical_smiles),
        analysis.canonical_smiles,
        analysis.exact_molecular_weight,
    )
    # 2026-09-03 수정: 결과를 버전이 포함된 키로 캐시(무한 증식 방지 상한 적용).
    if len(_PREDICTION_CACHE) >= _PREDICTION_CACHE_MAX:
        _PREDICTION_CACHE.clear()
    _PREDICTION_CACHE[cache_key] = _copy_prediction_result(result)
    return result


def predict_degradation_products(
    parent_smiles: str,
    stress_environment: str,
    temperature_c: float = 25.0,
) -> list[DegradantProduct]:
    """
    스마트 조건부 분해 엔진 — predict_advanced_degradants 래퍼 (하위 호환).

    SMILES 오류 시 빈 목록 반환. UI에서는 DegradationPredictionResult를 직접 사용 권장.
    """
    outcome = predict_advanced_degradants(parent_smiles, stress_environment, temperature_c)
    return outcome.products
