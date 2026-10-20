"""
DP(Drug Product) 부형제–API 상호작용 유연물질 시뮬레이션 엔진.

구성
----
    반응 규칙        p18_excipient_reaction_rules — 14종
                     (Maillard 이민·엔아민 · 포르밀화 · 아마이드 축합 ·
                      에스터화 · 에스터교환 · 에스터 아미놀리시스 ·
                      Michael 부가 · 아황산부가물 · 과산화물 산화 · 헤미아세탈)
    부형제 라이브러리  p70 + p19 확장 = 143종
    니트로소화        p30_reaction_engine
    화학 타당성       p74_nitrosamine_chemical_validity
    아질산염 프로파일  이 파일의 EXCIPIENT_NITRITE_MATRIX (7종)

클래스
------
    ExcipientImpurityHit           결과 한 건
    AdvancedExcipientMatrixEngine  본체
    ExcipientInteractionEngine     구버전 호환 래퍼

변경 이력
---------
2026-08-25  반응 규칙을 내장 2종에서 p18의 14종으로 교체.
            라이브러리에 p19 확장 적용 (76종 → 143종).
2026-08-26  p73_nitrosamine_excipient_engine 을 흡수하고 폐지했습니다.
            니트로소화 엔진이 두 벌이었고 서로 다른 기준으로 판단해
            «한쪽이 놓치는 것을 다른 쪽이 잡는» 상태였습니다 —
            p73은 API 게이트가 3급 아민을 막고, 이쪽은 부형제 게이트가
            크로스포비돈·MCC를 막았습니다. 이제 경로가 하나입니다.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from p70_excipient_master_library import (
    COMPLETE_EXCIPIENT_MASTER as _LOCAL_MASTER,
    EXCIPIENT_CATEGORY_ORDER,
    LEGACY_EXCIPIENT_ALIASES,
)
from p71_excipient_opensource_registry import (
    get_full_excipient_library,
    opensource_catalog_stats,
    probe_smiles_for_entry,
)
from p30_reaction_engine import (
    NITRITE_SMILES,
    detect_amine_sites,
    is_nitrite_source_smiles,
    predict_nitrosation_with_reagent,
)

_LEGACY_EXCIPIENT_ALIASES = LEGACY_EXCIPIENT_ALIASES

EXCIPIENT_LIBRARY: dict[str, dict[str, Any]] = get_full_excipient_library()

# 2026-08-25 — 구조 보강(25종) + 신규 등재.
#   기존 항목 중 다수가 smiles·reactive_probes 없이 등록되어 있어 반응
#   예측에 전혀 참여하지 못했습니다. p19가 구조와 반응 조각을 채우고
#   경구고형제 상용 부형제를 추가합니다. 확장 모듈이 없어도 동작합니다.
try:
    from p19_excipient_library_ext import apply_extensions as _apply_excipient_extensions

    EXCIPIENT_LIBRARY = _apply_excipient_extensions(EXCIPIENT_LIBRARY)
except Exception:  # pragma: no cover
    pass
COMPLETE_EXCIPIENT_MASTER: dict[str, dict[str, Any]] = EXCIPIENT_LIBRARY
EXCIPIENT_CHOICES: tuple[str, ...] = tuple(EXCIPIENT_LIBRARY.keys())
EXCIPIENT_CATEGORIES: tuple[str, ...] = EXCIPIENT_CATEGORY_ORDER

NDMA_CANONICAL_SMILES = "CN(C)N=O"
NITROSO_SMARTS = "[#7]-[#7]=[#8]"

# ---------------------------------------------------------------------------
# 부형제별 아질산염(Nitrite) 잔류 프로파일 — QbD 매트릭스
#
# 2026-08-26 — p73_nitrosamine_excipient_engine.py 에서 옮겨 왔습니다.
#   니트로소화 엔진이 p72·p73 두 벌로 나뉘어 있었고, 서로 다른 기준으로
#   판단해 «한쪽이 놓치는 것을 다른 쪽이 잡는» 상태였습니다.
#   p73의 엔진은 API 게이트가 3급 아민을 막았고, p72는 부형제 게이트가
#   크로스포비돈·MCC를 막았습니다. 하나로 합치면서 데이터를 이리로 옮깁니다.
#
# has_nitrite_risk 는 «평가한 결과» 입니다. False 는 자료가 없다는 뜻이
# 아니라 무시할 수준이라고 판단했다는 뜻입니다.
#
# 2026-09-03 — 다만 **False 를 «완전 차단»으로 쓰던 것을 «감쇄»로 바꿨습니다.**
#   락토스의 프로파일은 `nitrite_level_ppm: "< 0.01 (Low)"` 로 **0이 아닌
#   미량 잔류를 스스로 인정**하면서 경로는 완전히 막고 있었습니다.
#   EMA/FDA 니트로사민 가이던스는 부형제 아질산염이 ppb 수준이어도
#   NDSRI 생성 가능성을 «평가»하도록 요구하며, 락토스는 아민 API와의
#   상호작용이 문서화된 부형제입니다. 완전 차단은 위음성을 구조적으로
#   고정시킵니다.
#
#   이제 —
#     nitrite_level_ppm 이 "Negligible"        → 차단 (스테아르산마그네슘)
#     그 외 (미량이라도 수치가 기록됨)          → 반응은 수행하되
#         confidence="screening-only, below reporting threshold" 를 붙여
#         보고 임계 미만임을 결과에 명시
# ---------------------------------------------------------------------------
EXCIPIENT_NITRITE_MATRIX: dict[str, dict[str, Any]] = {
    "Crospovidone": {
        "matrix_id": "Crospovidone_PVP",
        "smiles": "CCN1CCCC1=O",
        "has_nitrite_risk": True,
        "nitrite_level_ppm": "0.1 - 5.0 (High)",
        "desc": "🔥 포비돈 중합·보관 중 미량 Nitrite 및 Peroxide 동적 축적률이 높음",
    },
    "Polyethylene_Glycol_4000": {
        "matrix_id": "Polyethylene_Glycol_PEG",
        "smiles": "OCCOCCOCCO",
        "has_nitrite_risk": True,
        "nitrite_level_ppm": "0.1 - 2.0 (Medium-High)",
        "desc": "🔥 에테르 사슬 분해 시 산화물·미량 Nitrite 동시 용출 경고",
    },
    "Polyethylene_Glycol_6000": {
        "matrix_id": "Polyethylene_Glycol_PEG",
        "smiles": "OCCOCCOCCO",
        "has_nitrite_risk": True,
        "nitrite_level_ppm": "0.1 - 2.0 (Medium-High)",
        "desc": "🔥 PEG 계열 — 에테르 산화·Nitrite leaching 위험",
    },
    "Microcrystalline_Cellulose": {
        "matrix_id": "Microcrystalline_Cellulose_MCC",
        "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "has_nitrite_risk": True,
        "nitrite_level_ppm": "0.05 - 1.0 (Medium)",
        "desc": "천연 펄프 가공 기원 미량 무기 아질산염 잔류 가능",
    },
    "Lactose_Monohydrate": {
        "matrix_id": "Lactose",
        "smiles": "OCC1OC(OC2C(O)C(O)OC(CO)C2O)C(O)C(O)C1O",
        "has_nitrite_risk": False,
        "nitrite_level_ppm": "< 0.01 (Low)",
        "desc": "Nitrite 잔류 극히 낮음 — 환원당 Maillard 갈변 리스크는 별도",
    },
    "Lactose_Anhydrous": {
        "matrix_id": "Lactose",
        "smiles": "OCC1OC(OC2C(O)C(O)OC(CO)C2O)C(O)C(O)C1O",
        "has_nitrite_risk": False,
        "nitrite_level_ppm": "< 0.01 (Low)",
        "desc": "Nitrite 잔류 극히 낮음 — 환원당 Maillard 갈변 리스크는 별도",
    },
    "Magnesium_Stearate": {
        "matrix_id": "Magnesium_Stearate",
        "smiles": "CCCCCCCCCCCCCCCCCC(=O)O",
        "has_nitrite_risk": False,
        "nitrite_level_ppm": "Negligible",
        "desc": "아질산염 리스크 낮음 — 금속 이온 산화 촉매 주의",
    },
}

# 레거시·별칭 → 마스터 키
#
# 2026-09-03 수정 — **p70.LEGACY_EXCIPIENT_ALIASES 를 단일 출처로 삼습니다.**
#   전에는 여기에 같은 매핑을 6건 그대로 베껴 두고 있었습니다
#   (Crospovidone_PVP · PVP_Crospovidone · Polyethylene_Glycol_PEG ·
#    PEG_Polyethylene_Glycol · Microcrystalline_Cellulose_MCC · Lactose).
#   한쪽만 고치면 부형제 해석과 아질산염 프로파일 조회가 서로 다른 키를
#   가리켜 조용히 어긋납니다. 이제 p70을 import 해 병합하고, 이 모듈
#   고유의 별칭만 아래 `_MATRIX_ONLY_ALIASES` 에 둡니다.
#   새 별칭은 원칙적으로 **p70에** 추가하세요.

# 아질산염 매트릭스에서만 쓰는 별칭 (p70에 없는 것)
_MATRIX_ONLY_ALIASES: dict[str, str] = {}

_NITRITE_MATRIX_ALIASES: dict[str, str] = {
    **LEGACY_EXCIPIENT_ALIASES,
    **_MATRIX_ONLY_ALIASES,
}

# 아질산(HNO2). 프로파일에 잔류가 기록된 부형제는 이 시약을 씁니다.
NITROUS_ACID_SMILES = "O=NO"

_NDSRI_ALERT_KEYWORDS = (
    "ndsri",
    "nitros",
    "nitroso",
    "아질산",
    "formic",
    "formaldehyde",
    "포름",
    "nitrite",
)


def _normalize_excipient_entry(data: dict[str, Any]) -> dict[str, Any]:
    out = dict(data)
    out["risk_alerts"] = out.get("alerts") or out.get("risk_alerts") or ""
    out.setdefault("reactive_probes", [])
    return out


def resolve_excipient_name(name: str) -> str:
    key = (name or "").strip()
    return _LEGACY_EXCIPIENT_ALIASES.get(key, key)

# 결과 한 건을 담는 그릇
#
# 2026-09-03 — 부형제 귀속을 **두 축으로 분리**했습니다.
#
#   3단계(산물 재니트로소화)는 «모든 부형제»의 매트릭스 산물을 기질로
#   받는데, 결과 hit의 `excipient_name`에는 **아질산염을 공급한 부형제**만
#   기록됐습니다. 반응 파트너를 제공한 부형제 정보가 통째로 사라져,
#   MCC에 귀속된 건 안에 락토스 유래 Maillard 엔아민 골격과 스테아르산
#   마그네슘 유래 아미드 골격이 섞여 들어갔습니다. `DEFAULT_EXCIPIENTS`
#   순서와 `seen_smiles` 공유가 겹치면서 **논문 표의 «부형제별» 순위가
#   인공물**이 됐습니다 (아질산염이 가장 높은 크로스포비돈에 NDSRI 0건).
#
#   `excipient_name` 은 하위 호환을 위해 그대로 두되(= 아질산염 공급자),
#   집계가 두 축으로 나뉠 수 있도록 다음 두 필드를 새로 둡니다.
#
#     nitrite_source_excipient    아질산염(니트로싱 환경)을 공급한 부형제
#     substrate_origin_excipient  반응 기질(골격)을 제공한 부형제
#
#   매트릭스 상호작용(1단계)에서는 둘 다 같은 부형제이고,
#   3단계에서는 서로 다를 수 있습니다.
@dataclass
class ExcipientImpurityHit:
    smiles: str
    mechanism_origin: str
    alert_context: str
    excipient_name: str
    rule_name: str
    category: str = ""
    risk_class: str = "Matrix"
    nitrosamine_label: str = ""
    nitrite_exposure: str = ""
    impurity_type: str = ""
    ich_m7_class_hint: str = ""
    alpha_hydrogen_present: str = ""
    steric_hindrance: str = ""
    kinetics_comment: str = ""
    chemical_validity_notes: str = ""
    nitrite_exposure_raw: str = ""
    # 2026-09-03 추가 — 부형제 귀속 두 축 + 신뢰도
    nitrite_source_excipient: str = ""
    substrate_origin_excipient: str = ""
    confidence: str = ""

    def to_dict(self) -> dict[str, str]:
        return {
            "smiles": self.smiles,
            "mechanism_origin": self.mechanism_origin,
            "alert_context": self.alert_context,
            "excipient_name": self.excipient_name,
            "rule_name": self.rule_name,
            "category": self.category,
            "risk_class": self.risk_class,
            "nitrosamine_label": self.nitrosamine_label,
            "nitrite_exposure": self.nitrite_exposure,
            "impurity_type": self.impurity_type,
            "ich_m7_class_hint": self.ich_m7_class_hint,
            "alpha_hydrogen_present": self.alpha_hydrogen_present,
            "steric_hindrance": self.steric_hindrance,
            "kinetics_comment": self.kinetics_comment,
            "chemical_validity_notes": self.chemical_validity_notes,
            "nitrite_exposure_raw": self.nitrite_exposure_raw,
            "nitrite_source_excipient": self.nitrite_source_excipient,
            "substrate_origin_excipient": self.substrate_origin_excipient,
            "confidence": self.confidence,
        }


def _mol_to_canonical_smiles(mol: Chem.Mol | None) -> str | None:
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return smi if Chem.MolFromSmiles(smi) else None
    except Exception:
        return None


def _product_smiles_from_tuple(prod_tuple: tuple[Chem.Mol, ...]) -> str | None:
    candidates: list[str] = []
    for prod in prod_tuple:
        try:
            work = Chem.RemoveHs(prod)
            smi = _mol_to_canonical_smiles(work)
            if smi:
                candidates.append(smi)
        except Exception:
            continue
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    return max(
        candidates,
        key=lambda s: (
            Descriptors.MolWt(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else 0.0
        ),
    )


def _classify_nitrosamine_product(smiles: str) -> tuple[str, str]:
    """Return (risk_class, nitrosamine_label) or ('', '')."""
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return "", ""
    try:
        canon = Chem.MolToSmiles(mol, canonical=True)
        ndma_mol = Chem.MolFromSmiles(NDMA_CANONICAL_SMILES)
        if ndma_mol is not None and canon == Chem.MolToSmiles(ndma_mol, canonical=True):
            return "NDMA", "NDMA (N-nitrosodimethylamine)"
    except Exception:
        pass
    try:
        pat = Chem.MolFromSmarts(NITROSO_SMARTS)
        if pat and mol.HasSubstructMatch(pat):
            return "NDSRI", "N-nitroso (NDSRI)"
    except Exception:
        pass
    return "", ""


def nitrite_profile_for_excipient(excipient_key: str) -> dict[str, Any] | None:
    """부형제 키 → 아질산염 잔류 프로파일. 등재되지 않았으면 None."""
    resolved = resolve_excipient_name(excipient_key)
    master_key = _NITRITE_MATRIX_ALIASES.get(resolved, resolved)
    return EXCIPIENT_NITRITE_MATRIX.get(master_key)


def nitrite_matrix_catalog_rows() -> list[dict[str, str]]:
    """화면·보고서용 아질산염 매트릭스 표."""
    return [
        {
            "부형제키": key,
            "Nitrite_리스크": "Yes" if data.get("has_nitrite_risk") else "No",
            "잔류_ppm": str(data.get("nitrite_level_ppm") or ""),
            "가이드": str(data.get("desc") or ""),
        }
        for key, data in EXCIPIENT_NITRITE_MATRIX.items()
    ]


# 2026-09-03 추가 — «무시할 수준»으로 못박힌 등급만 완전 차단합니다.
_NEGLIGIBLE_NITRITE_TOKENS = ("negligible", "none", "not detected", "n/a")

# 감쇄 경로로 나온 결과에 붙는 신뢰도 문구
SCREENING_ONLY_CONFIDENCE = "screening-only, below reporting threshold"


def _is_negligible_nitrite_level(level: Any) -> bool:
    """아질산염 잔류 등급이 «무시할 수준»으로 명시됐는지."""
    text = str(level or "").strip().lower()
    if not text:
        # 등급 자체가 비어 있으면 판단 근거가 없으므로 막지 않습니다
        return False
    return any(token in text for token in _NEGLIGIBLE_NITRITE_TOKENS)


def _excipient_supports_nitrosation(exc_data: dict[str, Any]) -> bool:
    """경고문 키워드 폴백. 프로파일이 없는 부형제에만 씁니다."""
    blob = f"{exc_data.get('risk_alerts', '')} {exc_data.get('alerts', '')}".lower()
    if any(k in blob for k in _NDSRI_ALERT_KEYWORDS):
        return True
    for smi in exc_data.get("reactive_probes") or []:
        if is_nitrite_source_smiles(str(smi)):
            return True
    return False


def nitrite_reagents_for_excipient(
    excipient_key: str, exc_data: dict[str, Any]
) -> tuple[list[str], str]:
    """
    이 부형제로 니트로소화를 돌릴지, 돌린다면 어떤 시약으로 할지 정합니다.

    판단 순서 —
      ① 아질산염 프로파일이 등재되어 있으면 그것만 따릅니다.
           has_nitrite_risk=True  → 부형제 자체 반응 조각 + 아질산(HNO2)
           has_nitrite_risk=False → **감쇄**. 2026-09-03 이전에는 여기서
                                    완전히 차단했습니다. 프로파일이 미량
                                    잔류를 인정(예: "< 0.01 (Low)")하는데도
                                    경로를 막아 위음성이 고정됐습니다.
                                    이제 "Negligible" 등급만 차단하고,
                                    수치가 기록된 것은 돌리되 결과에
                                    «보고 임계 미만» 신뢰도를 붙입니다.
      ② 프로파일이 없으면 반응 조각에 아질산 공급원이 있는지 봅니다.
      ③ 그것도 없으면 경고문 키워드로 마지막 판단을 하되,
         근거가 약하므로 «trace 가정» 으로 표시합니다.

    Returns: (시약 SMILES 목록, 근거 등급)
             근거 등급 — "profile" · "profile-attenuated" · "probe" ·
                         "assumed" · "" (차단)
    """
    probes = [str(s) for s in (exc_data.get("reactive_probes") or [])]
    nitrite_probes = [p for p in probes if is_nitrite_source_smiles(p)]

    profile = nitrite_profile_for_excipient(excipient_key)
    if profile is not None:
        reagents = list(nitrite_probes)
        if NITROUS_ACID_SMILES not in reagents:
            reagents.append(NITROUS_ACID_SMILES)

        if profile.get("has_nitrite_risk"):
            return reagents, "profile"

        # 완전 차단은 «무시할 수준»이라고 등급까지 못박은 경우에만
        if _is_negligible_nitrite_level(profile.get("nitrite_level_ppm")):
            return [], ""

        # 미량이지만 0이 아님 → 돌리되 감쇄 표시
        return reagents, "profile-attenuated"

    if nitrite_probes:
        return nitrite_probes, "probe"

    if _excipient_supports_nitrosation(exc_data):
        return [NITRITE_SMILES], "assumed"

    return [], ""


# 근거 등급을 화면 문구로
_NITRITE_BASIS_NOTE = {
    "profile": "부형제 아질산염 잔류 프로파일 등재",
    # 2026-09-03 추가 — has_nitrite_risk=False 이지만 잔류가 0은 아닌 경우
    "profile-attenuated": (
        "프로파일상 보고 임계 미만 — 스크리닝 목적으로만 평가 (완전 배제 아님)"
    ),
    "probe": "부형제 반응 조각에 아질산 공급원 존재",
    "assumed": "⚠ 실측 프로파일 없음 — 경고문 근거 trace 가정",
}

# 근거 등급 → 결과에 실리는 신뢰도
_NITRITE_BASIS_CONFIDENCE = {
    "profile": "profile-based",
    "profile-attenuated": SCREENING_ONLY_CONFIDENCE,
    "probe": "probe-based",
    "assumed": "assumed (alert-keyword only)",
}

# 실제 계산 (evalaute_formulation_risk)
class AdvancedExcipientMatrixEngine:
    # 2026-08-25 — 규칙 확장.
    #   기존에는 Maillard 이민과 포르밀화 두 개뿐이라 부형제를 늘려도
    #   잡히는 것이 없었습니다. p18이 아실화·Michael 부가·중아황산 부가·
    #   과산화물 산화 등을 추가합니다.
    #
    #   기존 Maillard 규칙 `[NX3;H2,H1]`은 두 가지가 틀렸습니다.
    #     · 1급 아마이드 질소에도 이민을 붙였습니다 (친핵성 부족)
    #     · 2급 아민에 중성 이민을 만들려 했습니다 (질소 원자가 초과)
    #   p18은 1급 아민 이민과 2급 아민 엔아민으로 나눴습니다.
    _FALLBACK_RULES = (
        ("Maillard_Imine_Bond", "[CH1:1](=O).[NX3;H2,H1:2]>>[CH1:1]=[NX3:2]"),
        ("Formylation_Artifact", "C(=O)O.[NX3;H2,H1:1]>>[NX3:1]C=O"),
    )

    def __init__(self) -> None:
        self.matrix_rules: dict[str, Chem.rdChemReactions.ChemicalReaction] = {}
        try:
            from p18_excipient_reaction_rules import build_reaction_rules

            self.matrix_rules = build_reaction_rules()
        except Exception:
            self.matrix_rules = {}

        if not self.matrix_rules:
            for name, smirks in self._FALLBACK_RULES:
                try:
                    rxn = AllChem.ReactionFromSmarts(smirks)
                    if rxn is not None:
                        self.matrix_rules[name] = rxn
                except Exception:
                    continue

    def _evaluate_matrix_interactions(
        self,
        api_mol: Chem.Mol,
        api_canon: str,
        selected_excipients_list: list[str],
        seen_smiles: set[str],
    ) -> list[ExcipientImpurityHit]:
        detected: list[ExcipientImpurityHit] = []
        for raw_name in selected_excipients_list:
            exc_name = resolve_excipient_name(raw_name)
            if exc_name not in EXCIPIENT_LIBRARY:
                continue

            exc_data = _normalize_excipient_entry(EXCIPIENT_LIBRARY[exc_name])
            probe_smiles = probe_smiles_for_entry(exc_data)

            exc_mols = [
                m
                for smi in probe_smiles
                if (m := Chem.MolFromSmiles(smi)) is not None
            ]
            if not exc_mols:
                continue

            for exc_mol in exc_mols:
                for rule_name, rxn in self.matrix_rules.items():
                    try:
                        product_sets = rxn.RunReactants((exc_mol, api_mol)) or rxn.RunReactants(
                            (api_mol, exc_mol)
                        )
                    except Exception:
                        continue
                    if not product_sets:
                        continue
                    for prod_tuple in product_sets:
                        prod_smiles = _product_smiles_from_tuple(prod_tuple)
                        if not prod_smiles or prod_smiles in seen_smiles:
                            continue
                        if prod_smiles == api_canon:
                            continue
                        seen_smiles.add(prod_smiles)
                        detected.append(
                            ExcipientImpurityHit(
                                smiles=prod_smiles,
                                mechanism_origin=f"API-{exc_name} 상호작용 ({rule_name})",
                                alert_context=str(exc_data.get("risk_alerts") or ""),
                                excipient_name=exc_name,
                                rule_name=rule_name,
                                category=str(exc_data.get("category") or ""),
                                risk_class="Matrix",
                                # 2026-09-03: 매트릭스 반응은 골격 제공자와
                                #   반응 상대가 같은 부형제입니다.
                                #   아질산염은 관여하지 않으므로 공급자는 없습니다.
                                nitrite_source_excipient="",
                                substrate_origin_excipient=exc_name,
                                confidence="matrix-reaction",
                            )
                        )
        return detected

    @staticmethod
    def _excipient_hit_from_validated_dict(h: dict[str, Any]) -> ExcipientImpurityHit:
        return ExcipientImpurityHit(
            smiles=str(h.get("smiles") or ""),
            mechanism_origin=str(h.get("mechanism") or h.get("mechanism_origin") or ""),
            alert_context=str(h.get("regulatory_alert") or h.get("alert_context") or ""),
            excipient_name=str(h.get("excipient_name") or ""),
            rule_name=str(h.get("rule_name") or "Direct_NDSRI_Formation"),
            category=str(h.get("category") or "Nitrosamine"),
            risk_class=str(h.get("risk_class") or "NDSRI"),
            nitrosamine_label=str(h.get("nitrosamine_label") or h.get("impurity_type") or ""),
            nitrite_exposure=str(h.get("nitrite_exposure") or ""),
            nitrite_exposure_raw=str(h.get("nitrite_exposure_raw") or h.get("nitrite_exposure") or ""),
            impurity_type=str(h.get("impurity_type") or ""),
            ich_m7_class_hint=str(h.get("ich_m7_class_hint") or ""),
            alpha_hydrogen_present=str(h.get("alpha_hydrogen_present", "")),
            steric_hindrance=str(h.get("steric_hindrance", "")),
            kinetics_comment=str(h.get("kinetics_comment") or ""),
            chemical_validity_notes=str(h.get("chemical_validity_notes") or ""),
            # 2026-09-03 추가 — 부형제 귀속 두 축 + 신뢰도
            nitrite_source_excipient=str(
                h.get("nitrite_source_excipient") or h.get("excipient_name") or ""
            ),
            substrate_origin_excipient=str(h.get("substrate_origin_excipient") or ""),
            confidence=str(h.get("confidence") or ""),
        )

    def _evaluate_nitrosamine_risk(
        self,
        api_smiles: str,
        selected_excipients_list: list[str],
        seen_smiles: set[str],
        *,
        substrates: list[tuple[str, ...]] | None = None,
    ) -> list[ExcipientImpurityHit]:
        """
        API 또는 상호작용 생성물 + 부형제 니트로싱 환경 → NDSRI/NDMA.

        2026-08-26 — 이 메서드가 니트로소화의 **유일한 경로**입니다.
        전에는 p73의 QbD 엔진(_evaluate_qbd_nitrosamine_risk)이 같은 일을
        다른 기준으로 한 번 더 했습니다. 둘이 `seen_smiles`를 공유해
        먼저 도는 쪽이 산물을 가져갔고, 어느 쪽이 잡느냐에 따라 표시되는
        메타데이터가 달라졌습니다. 하나로 합쳤습니다.

        2026-09-03 — `substrates` 항목이 **기질을 만든 부형제 키**를 함께
        나릅니다. 3단계에서 «아질산염을 준 부형제»와 «골격을 준 부형제»가
        서로 다른데 결과에는 전자만 남아 귀속이 소실되고 있었습니다.

            substrates = [(smiles, label), ...]                 # 예전 형식 (계속 허용)
            substrates = [(smiles, label, origin_excipient), ...]  # 새 형식

        예전 2-튜플이 오면 `origin_excipient`는 빈 문자열이 됩니다.
        """
        from p74_nitrosamine_chemical_validity import apply_chemical_validity_layer

        amine_info = detect_amine_sites(api_smiles)
        if not amine_info.has_nitrosatable_amine and not substrates:
            return []

        # (smiles, label, origin_excipient) 로 정규화 — 2-튜플 하위 호환
        raw_subs: list[tuple[str, ...]] = list(
            substrates if substrates is not None else [(api_smiles.strip(), "API")]
        )
        sub_list: list[tuple[str, str, str]] = []
        for item in raw_subs:
            smi = str(item[0])
            label = str(item[1]) if len(item) > 1 else ""
            origin = str(item[2]) if len(item) > 2 else ""
            sub_list.append((smi, label, origin))

        detected: list[ExcipientImpurityHit] = []

        for raw_name in selected_excipients_list:
            exc_name = resolve_excipient_name(raw_name)
            if exc_name not in EXCIPIENT_LIBRARY:
                continue
            exc_data = _normalize_excipient_entry(EXCIPIENT_LIBRARY[exc_name])
            reagents, basis = nitrite_reagents_for_excipient(exc_name, exc_data)
            if not reagents:
                continue

            for sub_smi, sub_label, sub_origin in sub_list:
                sub_amine = detect_amine_sites(sub_smi)
                if not sub_amine.has_nitrosatable_amine:
                    continue
                for reagent_smi in reagents:
                    try:
                        _, mols, meta = predict_nitrosation_with_reagent(sub_smi, reagent_smi)
                    except Exception:
                        continue
                    nitrite_note = _NITRITE_BASIS_NOTE.get(basis, basis)
                    for mol, (rule_id, rule_label, _) in zip(mols, meta):
                        try:
                            prod_smiles = Chem.MolToSmiles(mol, canonical=True)
                        except Exception:
                            continue
                        if not prod_smiles or prod_smiles in seen_smiles:
                            continue
                        risk_class, nitro_label = _classify_nitrosamine_product(prod_smiles)
                        if not risk_class:
                            continue
                        seen_smiles.add(prod_smiles)
                        profile = nitrite_profile_for_excipient(exc_name)
                        nitrite_ppm = str(
                            (profile or {}).get("nitrite_level_ppm") or "trace (assumed)"
                        )
                        base_alert = str(
                            (profile or {}).get("desc")
                            or exc_data.get("risk_alerts")
                            or ""
                        )
                        if basis in ("assumed", "profile-attenuated") and base_alert:
                            base_alert = f"{base_alert} · {_NITRITE_BASIS_NOTE[basis]}"

                        # 2026-09-03: 기질 골격을 제공한 부형제.
                        #   API에서 직접 온 기질이면 부형제가 아니라 API입니다.
                        substrate_origin = sub_origin or (
                            "" if sub_label == "API" else sub_origin
                        )
                        origin_phrase = (
                            f"{sub_label} [골격 제공: {substrate_origin}]"
                            if substrate_origin and substrate_origin != exc_name
                            else sub_label
                        )
                        raw_hit = {
                            "smiles": prod_smiles,
                            # 하위 호환 — 예전부터 «아질산염 공급 부형제»였습니다
                            "excipient_name": exc_name,
                            "nitrite_source_excipient": exc_name,
                            "substrate_origin_excipient": substrate_origin,
                            "confidence": _NITRITE_BASIS_CONFIDENCE.get(basis, basis),
                            "impurity_type": (
                                "NDMA (고유 유전독성 불순물)"
                                if risk_class == "NDMA"
                                else "API 유래 NDSRI (니트로사민 불순물)"
                            ),
                            "mechanism_origin": (
                                f"{origin_phrase} + {exc_name} → {rule_label} ({nitrite_note})"
                            ),
                            "nitrite_exposure": nitrite_ppm,
                            "nitrite_exposure_raw": nitrite_ppm,
                            "regulatory_alert": base_alert,
                            "alert_context": base_alert,
                            "rule_name": rule_id,
                            "risk_class": risk_class,
                            "nitrosamine_label": nitro_label,
                            "category": str(exc_data.get("category") or "Nitrosamine"),
                        }
                        validated = apply_chemical_validity_layer(raw_hit)
                        detected.append(self._excipient_hit_from_validated_dict(validated))
        return detected

    def evaluate_formulation_risk(
        self,
        api_smiles: str,
        selected_excipients_list: list[str],
    ) -> list[ExcipientImpurityHit]:
        """
        세 단계를 순서대로 돌리고 결과를 이어 붙입니다.

            ① 매트릭스 상호작용   Maillard · 에스터화 · 산화 등 14종 (p18)
            ② 니트로소화          API + 부형제 니트로싱 환경
            ③ 산물 재처리         ①의 생성물을 다시 니트로소화

        `seen_smiles`를 공유해 같은 불순물이 두 줄로 나오지 않게 합니다.
        (부형제 **사이**의 공유는 위음성을 만들므로, 화면은 부형제마다
         엔진을 새로 만드는 이원 1:1 모드를 씁니다 — p17)

        2026-08-26 — 전에는 ② 자리에 니트로소화가 두 벌 있었습니다.
        p73의 QbD 엔진과 이 모듈의 것이 서로 다른 기준으로 판단했고,
        p73은 API 게이트가 3급 아민을 막고 이쪽은 부형제 게이트가
        크로스포비돈·MCC를 막아 «서로가 놓치는 것을 상대가 잡는» 상태였습니다.
        하나로 합치면서 p73은 폐지했습니다.
        """
        api_mol = Chem.MolFromSmiles((api_smiles or "").strip())
        if api_mol is None:
            return []
        try:
            Chem.SanitizeMol(api_mol)
        except Exception:
            return []

        api_canon = Chem.MolToSmiles(api_mol, canonical=True)
        seen_smiles: set[str] = set()

        matrix_hits = self._evaluate_matrix_interactions(
            api_mol, api_canon, selected_excipients_list, seen_smiles
        )
        nitro_hits = self._evaluate_nitrosamine_risk(
            api_smiles, selected_excipients_list, seen_smiles
        )

        # 2026-09-03: 각 매트릭스 산물이 **어느 부형제에서 나온 골격인지**를
        #   함께 넘깁니다. 전에는 (smiles, label) 2-튜플이라 골격 제공
        #   부형제가 소실되고, 결과에는 아질산염 공급 부형제만 남았습니다.
        matrix_substrates = [
            (
                h.smiles,
                f"Matrix product ({h.rule_name})",
                h.substrate_origin_excipient or h.excipient_name,
            )
            for h in matrix_hits
        ]
        product_nitro = self._evaluate_nitrosamine_risk(
            api_smiles,
            selected_excipients_list,
            seen_smiles,
            substrates=matrix_substrates,
        )

        return matrix_hits + nitro_hits + product_nitro

    def evaluate_formulation_risk_legacy_matrix_only(
        self,
        api_smiles: str,
        selected_excipients_list: list[str],
    ) -> list[ExcipientImpurityHit]:
        """Maillard/Formylation only (internal regression helper)."""
        api_mol = Chem.MolFromSmiles((api_smiles or "").strip())
        if api_mol is None:
            return []
        api_canon = Chem.MolToSmiles(api_mol, canonical=True)
        return self._evaluate_matrix_interactions(
            api_mol, api_canon, selected_excipients_list, set()
        )

# 구버전 호환 래퍼
class ExcipientInteractionEngine(AdvancedExcipientMatrixEngine):
    @property
    def interaction_rules(self) -> dict[str, Chem.rdChemReactions.ChemicalReaction]:
        return self.matrix_rules

    def simulate_excipient_interaction(
        self, api_smiles: str, selected_excipient_names: list[str]
    ) -> list[ExcipientImpurityHit]:
        return self.evaluate_formulation_risk(api_smiles, selected_excipient_names)


def excipient_labels_for_ui() -> dict[str, str]:
    labels: dict[str, str] = {}
    for key, data in EXCIPIENT_LIBRARY.items():
        entry = _normalize_excipient_entry(data)
        cat = entry.get("category", "")
        alert = str(entry.get("risk_alerts") or "")
        short = f"{alert[:28]}…" if len(alert) > 28 else alert
        labels[key] = f"{entry.get('label_ko', key)} [{cat}] — {short}"
    return labels


def excipients_by_category() -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {cat: [] for cat in EXCIPIENT_CATEGORY_ORDER}
    for key, data in EXCIPIENT_LIBRARY.items():
        cat = str(data.get("category") or "Other")
        grouped.setdefault(cat, []).append(key)
    return {k: v for k, v in grouped.items() if v}


def filter_excipient_keys(
    *,
    categories: list[str] | None = None,
    search_query: str = "",
) -> list[str]:
    """UI 필터 — 카테고리·검색어로 부형제 키 목록 반환."""
    q = (search_query or "").strip().lower()
    labels = excipient_labels_for_ui()
    cat_set = set(categories or []) - {"전체"}
    out: list[str] = []
    for key in EXCIPIENT_CHOICES:
        entry = _normalize_excipient_entry(EXCIPIENT_LIBRARY[key])
        if cat_set and entry.get("category") not in cat_set:
            continue
        if q:
            blob = f"{key} {entry.get('label_ko', '')} {entry.get('category', '')}".lower()
            if q not in blob:
                continue
        out.append(key)
    return out


def excipient_catalog_rows() -> list[dict[str, str]]:
    """전체 부형제 목록 테이블용."""
    rows: list[dict[str, str]] = []
    for key in EXCIPIENT_CHOICES:
        e = _normalize_excipient_entry(EXCIPIENT_LIBRARY[key])
        rows.append(
            {
                "부형제": str(e.get("label_ko") or key),
                "카테고리": str(e.get("category") or ""),
                "출처": str(e.get("source") or "PRIME-DS"),
                "주의보": str(e.get("risk_alerts") or ""),
                "key": key,
            }
        )
    return rows


if __name__ == "__main__":
    engine = AdvancedExcipientMatrixEngine()
    sample_api = "Nc1ccccc1"
    user_choices = ["Lactose_Monohydrate", "Polyethylene_Glycol_4000", "Crospovidone"]
    print("=" * 80)
    print(f" [PRIME-DS] 부형제 마스터 풀: {len(COMPLETE_EXCIPIENT_MASTER)}종")
    print("=" * 80)
    for imp in engine.evaluate_formulation_risk(sample_api, user_choices):
        print(f"▶ {imp.smiles} | {imp.excipient_name} | {imp.rule_name}")
