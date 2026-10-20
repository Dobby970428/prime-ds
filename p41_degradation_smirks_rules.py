"""
Universal SMIRKS degradation rule catalog — hydrolysis + autoxidation.

Multi-generation degradation tree with optional Layer C GNN / DeepChem
toxicity and ICH M7 hybrid ensemble per node.
"""

# 규칙 목록

# 2026-09-03 수정: docstring이 실제 구현과 어긋나 있어 바로잡았습니다.
#   (구) "규칙 6종 / 최대 4세대" → (실제) DEGRADATION_RULES 7종 · 기본 max_depth=3

# DEGRADATION_RULES 7종
# 가수분해 3종:   에스터 · 아마이드 · 카바메이트
# 산화    3종:   티오에터 · 3차아민 · 에터 탈알킬화
#                 └ 에터 탈알킬화는 CYP Phase I 대사이므로 metabolic=True 로 표시되어
#                   shelf_life / tree_engine 규칙 집합에서 제외됩니다(실행 규칙은 6종).
# NDSRI  1종:    니트로사민 탈니트로소화

# 다세대 추적 기능
# 아스피린 (G0)
#    ├─ 살리실산 (G1)
#    │     └─ 더 분해되나? (G2)
#    └─ 아세트산 (G1)
#          └─ 더 분해되나? (G2)
# 기본 max_depth=3 → 세대 0·1·2 (3개 층)까지 생성
# 같은 구조가 다시 나오면 무한루프 막기 위해 기록 후 건너뜀(_visited_smiles)
# 화면의 "분해 트리"

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem

from p16_ich_m7_terms import class_definition, class_label
from p23_chem_search_engine import safe_canonical_smiles
from p44_degradation_shelf_life_filter import ShelfLifeChildFilter

# -------------------------------------------------------------------------
# 분해 반응 SMIRKS — ICH Q3A/B·M7 in-vitro 보관 안정성 전용
# 허용: 가수분해 · 비라디칼 대기 자동산화 · NDSRI(니트로사민 화학 분해)
# 금지: Phase I 대사(방향족 수산화, 에폭시화, 벤질릭 수산화 등)
# -------------------------------------------------------------------------
SHELF_LIFE_RULE_CATEGORIES = frozenset({"Hydrolysis", "Oxidation", "NDSRI"})
TREE_ENGINE_RULE_CATEGORIES = frozenset({"Hydrolysis", "Oxidation", "NDSRI"})

# DegradationTreeEngine — 이름·SMIRKS 이중 차단 (생체 대사 규칙 잔존 방지)
_METABOLIC_RULE_NAME_BLOCKLIST = frozenset(
    {
        "Aromatic_Hydroxylation",
        "Benzylic_Hydroxylation",
        "Aromatic_Demethylation",
        "Epoxidation",
        "N_Dealkylation",
        "Aryl_Dehalogenation",
        "Nitro_Reduction",
        "Phenol_Hydroxylation",
        "Phase1_Aromatic_Oxidation",
    }
)

# Phase I / 라디칼 유도 반응 SMIRKS fragment (부분 일치 시 로드 거부)
# ⚠ 여기 있는 것은 **실행할 반응 템플릿이 아니다.**
#   `_smirks_is_metabolic_or_radical_prone()` 가 문자열 **부분일치**로
#   «이 SMIRKS 가 대사·라디칼 계열인가» 를 판정할 때 쓰는 차단 목록이다.
#   그래서 "[#7;H0;D3" 처럼 일부러 잘린 조각이 들어 있고, RDKit 반응으로
#   파싱될 필요가 없다. 규칙 수를 셀 때 이 튜플을 포함시키지 말 것
#   (s19_rule_inventory.py 가 제외한다).
_METABOLIC_SMIRKS_FRAGMENTS: tuple[str, ...] = (
    "[c:1][cH:2]>>[c:1][c:2][OH]",  # 방향족 고리 C-H → C-OH
    "[c:1][CH2:2]>>[c:1][CH:2][OH]",  # 벤질릭 α-수산화
    "[c:1][CH3:2]>>[c:1][H:3]",  # 방향족 탈메틸화 (CYP)
    "[C:1]=[C:2]>>[C:1]1[C:2]1",  # 벤젠/알켄 에폭시화
    "[c:1][F,Cl,Br,I:2]>>[c:1][H:3]",  # 방향족 탈할로겐화
    # 2026-09-15 수정: "[N:1+]" 는 잘못된 표기다. 원자 매핑 뒤에 전하를 쓸 수
    #   없으므로 실제 규칙은 "[N+:1]" 로 쓰인다. 이 항목은 어떤 SMIRKS 와도
    #   부분일치하지 않아 **아무것도 걸러내지 못하고 있었다**. 지금은
    #   DEGRADATION_RULES 에 니트로 환원 규칙이 없어 무해하지만, 추가되면
    #   차단돼야 할 대사 규칙이 조용히 통과한다. 두 표기를 모두 둔다.
    "[N+:1](=[O:2])[O-:3]>>[N:1]",  # 니트로 환원(대사)
    "[N:1+](=[O:2])[O-:3]>>[N:1]",  # 〃 (옛 표기 — 기존 데이터 호환)
    "[#7;H0;D3",  # N-탈알킬화 계열
    "[nH:1]>>[n+:1][O-]",  # 방향족 N-산화 (대사형)
    "[OX2H]-[c;r6]",  # 벤젠 고리 직접 수산화 (Phase I)
    ">>[c:2][OH]",  # 방향족 C-H → C-OH 산물
)

DEGRADATION_RULES: dict[str, Any] = {
    "Hydrolysis_Ester": {
        "smirks": "[CX3:1](=[O:2])[OX2:3][#6:4]>>[CX3:1](=[O:2])[OX2H].[#6:4][OX2H]",
        "category": "Hydrolysis",
    },
    "Hydrolysis_Amide": {
        # 2026-09-03 수정: 산물 템플릿의 [NX3H:3]가 수소를 1개로 못 박아, 파라세타몰에서
        #   [NH]c1ccc(O)cc1 라디칼이 생성되고 라디칼 필터에 걸려 조용히 삭제되었습니다.
        #   H 고정을 제거해 RDKit이 수소를 재계산하도록 했습니다 → 4-아미노페놀이 나옵니다.
        "smirks": "[CX3:1](=[O:2])[NX3:3][c,C:4]>>[CX3:1](=[O:2])[OX2H].[NX3:3][c,C:4]",
        "category": "Hydrolysis",
    },
    "Hydrolysis_Carbamate": {
        # 2026-09-03 수정: 같은 모듈 안에서 Hydrolysis_Ester는 중성인데 여기만 음이온([O-])
        #   이었습니다. DB canonical SMILES 대조가 영구 실패하므로 중성형으로 통일.
        # 2026-09-29 수정: 산물 템플릿이 질소를 떼어낸 뒤 카본에 [OX2H]만 붙여
        #   놓아, RDKit 이 남는 원자가를 수소로 채우면서 **폼산(HCOOH)** 이
        #   나왔습니다. 카바메이트 가수분해는 아민 + CO2 + 알코올을 줍니다
        #   (카바믹산을 거쳐 즉시 탈카복실화). 카바메이트 탄소를 CO2 로
        #   매핑해 교정합니다.
        #     검증: CNC(=O)Oc1ccccc1 → CN · O=C=O · Oc1ccccc1
        #   ※ 이 규칙은 사례연구 세 건 어디에서도 발화하지 않았으므로
        #     (해당 API 에 카바메이트 없음) 보고된 수치는 바뀌지 않습니다.
        "smirks": "[#7:1][#6:2](=O)O[#6:3]>>[#7:1].O=[#6:2]=O.[#6:3][OH]",
        "category": "Hydrolysis",
    },
    "Oxidation_Thioether": {
        "smirks": "[#6:1][SX2:2][#6:3]>>[#6:1][SX3:2](=[O])[#6:3]",
        "category": "Oxidation",
    },
    "Oxidation_Tertiary_Amine": {
        "smirks": "[#6:1][NX3;!$(N=O);!$(NC=O):2]([#6:3])[#6:4]>>[#6:1][NX3+:2]([#6:3])([#6:4])[O-]",
        "category": "Oxidation",
    },
    "Oxidation_Ether_Dealkylation": {
        # 2026-09-03 수정: O-탈메틸화는 전형적인 CYP Phase I 대사(아니솔 → 페놀 + 포름알데하이드)로,
        #   이 모듈 docstring이 "금지: Phase I 대사"라고 명시한 대상입니다. 그런데 이름 기반
        #   블록리스트에는 존재하지 않는 이름(Aromatic_Demethylation)만 있었고 SMIRKS 조각
        #   블록리스트에도 없어 두 방어선을 모두 통과했습니다. 오타 하나로 뚫리는 이름 매칭 대신
        #   규칙 dict에 명시 플래그를 두어 차단합니다.
        "smirks": "[c,C:1][OX2:2][CH3:3]>>[c,C:1][OX2H].[CH2:3]=O",
        "category": "Oxidation",
        "metabolic": True,
    },
    "Denitrosation_NDSRI": {
        # 2026-09-03 수정: 구 규칙 Nitrosamination_NDSRIs 는 존재하지 않는 알킬 이동을
        #   모델링했습니다 — NDMA → CN + CN=O(니트로소메탄, C-nitroso). 니트로사민의
        #   탈니트로소화는 R2NH + HNO2를 주고 알킬기가 니트로소 질소로 옮겨가지 않습니다.
        #   게다가 산물 C-N=O가 니트로소기를 갖는 탓에 하위 모듈에서 cohort of concern으로
        #   재판정되어 등급이 부풀려졌습니다. 이름도 실제 동작(분해)에 맞춰 변경했습니다
        #   (_RULE_MECHANISM_CODE는 이미 SMIRKS_NDSRI_DENITRO로 적혀 있었습니다).
        "smirks": "[NX3:1]([#6:2])([#6:3])[NX2:4]=[OX1:5]>>[NX3;H1:1]([#6:2])[#6:3]",
        "category": "NDSRI",
    },
}

_RULE_MECHANISM_CODE: dict[str, str] = {
    "Hydrolysis_Ester": "SMIRKS_EST",
    "Hydrolysis_Amide": "SMIRKS_AMD",
    "Hydrolysis_Carbamate": "SMIRKS_CARB",
    "Oxidation_Thioether": "SMIRKS_S_OX",
    "Oxidation_Tertiary_Amine": "SMIRKS_N_OX",
    "Oxidation_Ether_Dealkylation": "SMIRKS_ETH_OX",
    # 2026-09-03 수정: 규칙 이름 변경에 맞춰 키 갱신 (구: Nitrosamination_NDSRIs)
    "Denitrosation_NDSRI": "SMIRKS_NDSRI_DENITRO",
}

# 2026-09-03 추가: 가수분해 기전 코드 — 방향족 수산화 컬링 면제 판정에 사용합니다.
#   아릴 에스터/아마이드가 가수분해되면 원래 분자에 있던 산소가 -OH로 "드러날" 뿐인데,
#   ShelfLifeChildFilter는 이를 부모에 없던 방향족 OH로 보고 대사성 수산화로 오판해
#   컬링합니다(아스피린 → 살리실산이 통째로 사라지던 원인). p43에는 같은 면제가 이미
#   있었으나 p41에는 적용되어 있지 않았습니다 — 여기서도 rule_name 기준으로 면제합니다.
_HYDROLYSIS_MECHANISM_CODES = frozenset({"SMIRKS_EST", "SMIRKS_AMD", "SMIRKS_CARB"})


def _is_hydrolysis_rule(rule_name: str | None) -> bool:
    """rule_name이 가수분해 기전인지 (_RULE_MECHANISM_CODE 기준)."""
    if not rule_name:
        return False
    return _RULE_MECHANISM_CODE.get(rule_name, "") in _HYDROLYSIS_MECHANISM_CODES


def passes_shelf_life_child_filter(
    parent_canon: str,
    child_smiles: str,
    rule_name: str | None = None,
) -> bool:
    """보관 안정성 Child 필터 — 가수분해 산물은 방향족 수산화 검사만 면제.

    라디칼·Parent 동일 등 다른 배제 사유는 그대로 적용됩니다.
    같은 규칙을 재적용하는 다른 모듈(p32)도 이 함수를 쓰세요.
    """
    verdict = ShelfLifeChildFilter.evaluate(parent_canon, child_smiles)
    if verdict.eligible:
        return True
    if verdict.reason == "aromatic_r6_hydroxylation" and _is_hydrolysis_rule(rule_name):
        return True
    return False


# 모듈 내부 호출용 별칭 (하위 호환)
_passes_shelf_life_child_filter = passes_shelf_life_child_filter


@dataclass(frozen=True)
class SmirksDegradantHit:
    rule_name: str
    category: str
    degradant_smiles: str

    def to_dict(self) -> dict[str, str]:
        return {
            "Rule_Name": self.rule_name,
            "Degradant_SMILES": self.degradant_smiles,
            "Category": self.category,
        }


def degradation_smirks_by_rule() -> dict[str, str]:
    return {name: info["smirks"] for name, info in tree_engine_degradation_rules().items()}


def _smirks_is_metabolic_or_radical_prone(smirks: str) -> bool:
    """Phase I 대사·라디칼 유도 SMIRKS 여부."""
    normalized = (smirks or "").replace(" ", "")
    if not normalized:
        return True
    return any(frag in normalized for frag in _METABOLIC_SMIRKS_FRAGMENTS)


def shelf_life_degradation_rules() -> dict[str, dict[str, str]]:
    """가수분해·자동산화·NDSRI 분해 규칙 (생체 대사 SMIRKS 제외)."""
    out: dict[str, dict[str, str]] = {}
    for name, info in DEGRADATION_RULES.items():
        if info.get("category") not in SHELF_LIFE_RULE_CATEGORIES:
            continue
        # 2026-09-03 추가: 명시 플래그 우선 — 이름 오타로 뚫리지 않는 1차 방어선
        if info.get("metabolic"):
            continue
        if name in _METABOLIC_RULE_NAME_BLOCKLIST:
            continue
        smirks = info.get("smirks", "")
        if _smirks_is_metabolic_or_radical_prone(smirks):
            continue
        out[name] = info
    return out


def tree_engine_degradation_rules() -> dict[str, dict[str, str]]:
    """
    DegradationTreeEngine 전용 규칙 — 가수분해 + 비라디칼 자동산화 + NDSRI.

    `self.rules`는 이 함수 결과만 사용합니다.
    """
    out: dict[str, dict[str, str]] = {}
    for name, info in DEGRADATION_RULES.items():
        if info.get("category") not in TREE_ENGINE_RULE_CATEGORIES:
            continue
        # 2026-09-03 추가: 명시 플래그 우선 — 이름 오타로 뚫리지 않는 1차 방어선
        if info.get("metabolic"):
            continue
        if name in _METABOLIC_RULE_NAME_BLOCKLIST:
            continue
        smirks = info.get("smirks", "")
        if _smirks_is_metabolic_or_radical_prone(smirks):
            continue
        out[name] = info
    return out


def _normalize_smirks(smirks: str) -> str:
    return (smirks or "").replace(" ", "").strip()


def _canonicalize_smiles(smiles: str) -> str | None:
    canon, _ = safe_canonical_smiles((smiles or "").strip())
    if not canon:
        return None
    return canon if Chem.MolFromSmiles(canon) else None


def _mol_to_canonical_smiles(mol: Chem.Mol) -> str | None:
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    try:
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return smi if Chem.MolFromSmiles(smi) else None
    except Exception:
        return None


def _run_single_step_rules(smiles: str) -> list[dict[str, str]]:
    parent_canon = _canonicalize_smiles(smiles)
    if not parent_canon:
        return []

    mol = Chem.MolFromSmiles(parent_canon)
    if mol is None:
        return []

    step_results: list[dict[str, str]] = []
    seen: set[str] = {parent_canon}

    for rule_name, rule_info in shelf_life_degradation_rules().items():
        smirks = _normalize_smirks(rule_info["smirks"])
        try:
            rxn = AllChem.ReactionFromSmarts(smirks)
        except Exception:
            continue
        try:
            product_sets = rxn.RunReactants((mol,))
        except Exception:
            continue

        for product_set in product_sets:
            for prod in product_set:
                try:
                    prod_smiles = _mol_to_canonical_smiles(prod)
                    if not prod_smiles or prod_smiles in seen or prod_smiles == parent_canon:
                        continue
                    # 2026-09-03 수정: 가수분해 산물은 방향족 수산화 컬링에서 면제
                    #   (아스피린 → 살리실산이 여기서 통째로 사라지고 있었습니다).
                    if not _passes_shelf_life_child_filter(
                        parent_canon, prod_smiles, rule_name
                    ):
                        continue
                    seen.add(prod_smiles)
                    step_results.append({"rule": rule_name, "smiles": prod_smiles})
                except Exception:
                    continue

    return step_results


def _predict_node_toxicity(smiles: str) -> dict[str, Any]:
    """Layer C GNN / DeepChem RF — Ames mutagenicity probability per tree node."""
    gnn_model = ""
    prob: float | None = None
    fallback_reason = ""

    try:
        from p64_gnn_mutagenicity_predictor import predict_mutagenicity_probability

        gnn = predict_mutagenicity_probability(smiles)
        prob = gnn.p_mutagenic
        gnn_model = gnn.model_name
        if gnn.fallback_reason:
            fallback_reason = gnn.fallback_reason
    except Exception as exc:
        fallback_reason = str(exc)[:120]

    if prob is None:
        try:
            from p63_deepchem_predictor import predict_advanced_properties

            row = predict_advanced_properties([smiles]).get(smiles, {})
            raw = row.get("Ames_Mutagenicity_Prob")
            prob = float(raw) if raw is not None else None
            if not gnn_model:
                gnn_model = "ECFP-2048 + RandomForest (Tox21 SR-p53)"
        except Exception:
            pass

    if prob is None:
        return {
            "Ames_Mutagenicity_Prob": None,
            "Ames_AI_Probability": "N/A",
            "Prediction_Class": "Unknown",
            "GNN_Model": gnn_model or "unavailable",
            "Fallback_Log": fallback_reason,
        }

    pred_class = "Positive" if prob >= 0.5 else "Negative"
    return {
        "Ames_Mutagenicity_Prob": round(prob, 4),
        "Ames_AI_Probability": f"{round(prob * 100, 1)}%",
        "Prediction_Class": pred_class,
        "GNN_Model": gnn_model,
        "Fallback_Log": fallback_reason,
    }


# 2026-09-03 수정: 대체 M7 판정이 ICH M7 정의와 어긋나 있었습니다.
#   ① prob is None(모델 없음/실패)에 Class 5를 주면, 모델 파일이 없는 배포본에서 모든
#      분해산물이 트리 JSON·보고서에 "Class 5 — Non-mutagenic impurity"로 찍힙니다.
#      Class 5는 "알러트 없음 또는 비변이원성을 입증할 충분한 데이터"이지 모델이 답을
#      못 낸 상태가 아니므로, 규제 문서에 들어가면 안 되는 fail-open이었습니다.
#      → 미판정 + 전문가 검토 필요로 반환합니다.
#   ② "Class 3 (Mutagenic impurity)"는 Class 3의 실제 정의(알러트 구조·변이원성 데이터
#      없음)가 아니고, "Class 4 (Alerting structure — confirmatory Ames recommended)"도
#      Class 4의 실제 조치(비변이원성으로 취급)와 반대였습니다.
#      → 등급 문자열을 직접 적지 않고 p16_ich_m7_terms의 공식 용어를 가져다 씁니다.
_QSAR_UNAVAILABLE_STATUS = "Not classified — (Q)SAR unavailable, expert review required"


def _evaluate_ich_m7_from_prob(prob: float | None, rule_name: str | None) -> str:
    """Conservative fallback when full ensemble unavailable.

    등급 문자열은 p16_ich_m7_terms(ICH M7(R2) 원문 단일 출처)에서 가져옵니다.
    """
    if prob is None:
        return _QSAR_UNAVAILABLE_STATUS
    if prob >= 0.70:
        # 알러트 구조 + 변이원성 시험 데이터 없음 → Class 3
        label = class_label(3)
        if not rule_name:
            label = f"{label} [AI Inferred Liability]"
        return f"{label} — P(Mutagenic)={prob:.2f}"
    if prob >= 0.50:
        # 경계 구간 — 확정 등급을 배정하지 않고 Ames 확인시험 권고
        return (
            f"Not classified — borderline (Q)SAR signal P={prob:.2f}; "
            f"bacterial mutagenicity assay recommended ({class_definition(3)})"
        )
    # 알러트 없음 / 비변이원성 신호 → Class 5
    return f"{class_label(5)} — P(Mutagenic)={prob:.2f}"


def _ich_m7_node_analysis(smiles: str, ames_prob: float | None, rule_name: str | None) -> dict[str, Any]:
    """Hybrid ensemble ICH M7 per node; falls back to probability heuristic."""
    try:
        from p66_ich_m7 import classify_ich_m7

        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            m7 = classify_ich_m7(mol)
            status = m7.assigned_class_label or f"Class {m7.class_number}"
            if m7.ai_inferred_liability:
                status = "Class 3-Alert (AI Inferred Liability)"
            return {
                "ICH_M7_Class": m7.class_number,
                "ICH_M7_Status": status,
                "ICH_M7_Rationale": (m7.rationale or "")[:160],
                "Ensemble_P_Mutagenic": m7.gnn_mutagenic_probability,
            }
    except Exception:
        pass

    return {
        "ICH_M7_Class": None,
        "ICH_M7_Status": _evaluate_ich_m7_from_prob(ames_prob, rule_name),
        "ICH_M7_Rationale": "",
        "Ensemble_P_Mutagenic": ames_prob,
    }


def _build_toxicity_analysis(smiles: str, rule_name: str | None = None) -> dict[str, Any]:
    tox = _predict_node_toxicity(smiles)
    prob = tox.get("Ames_Mutagenicity_Prob")
    ich = _ich_m7_node_analysis(smiles, prob if isinstance(prob, (int, float)) else None, rule_name)
    return {**tox, **ich}


def _build_tree_node(
    parent_canon: str,
    generation: int,
    *,
    include_toxicity: bool,
    rule_name: str | None = None,
) -> dict[str, Any]:
    node: dict[str, Any] = {
        "name": parent_canon,
        "generation": generation,
        "children": [],
    }
    if include_toxicity:
        node["toxicity_analysis"] = _build_toxicity_analysis(parent_canon, rule_name)
    if rule_name:
        node["triggered_by"] = rule_name
    return node


def generate_degradation_tree(
    current_smiles: str,
    max_depth: int = 3,
    *,
    _current_depth: int = 1,
    _visited_smiles: set[str] | None = None,
    include_toxicity: bool = False,
) -> dict[str, Any] | None:
    """Recursive SMIRKS tree — optional toxicity / ICH M7 per node."""
    if max_depth < 1:
        max_depth = 1

    if _visited_smiles is None:
        _visited_smiles = set()

    parent_canon = _canonicalize_smiles(current_smiles)
    if not parent_canon or parent_canon in _visited_smiles:
        return None

    node = _build_tree_node(
        parent_canon,
        _current_depth - 1,
        include_toxicity=include_toxicity,
    )
    _visited_smiles.add(parent_canon)

    if _current_depth >= max_depth:
        return node

    local_seen: set[str] = set()
    for step in _run_single_step_rules(parent_canon):
        next_smiles = step["smiles"]
        rule_triggered = step["rule"]
        if next_smiles in local_seen:
            continue
        local_seen.add(next_smiles)

        child_tree = generate_degradation_tree(
            next_smiles,
            max_depth=max_depth,
            _current_depth=_current_depth + 1,
            _visited_smiles=_visited_smiles.copy(),
            include_toxicity=include_toxicity,
        )
        if child_tree is None:
            child_tree = _build_tree_node(
                next_smiles,
                _current_depth,
                include_toxicity=include_toxicity,
                rule_name=rule_triggered,
            )
        else:
            child_tree["triggered_by"] = rule_triggered
            if include_toxicity and "toxicity_analysis" not in child_tree:
                child_tree["toxicity_analysis"] = _build_toxicity_analysis(
                    next_smiles, rule_triggered
                )
        node["children"].append(child_tree)

    return node


def generate_toxic_degradation_tree(
    current_smiles: str,
    max_depth: int = 3,
    *,
    _current_depth: int = 1,
    _visited_smiles: set[str] | None = None,
) -> dict[str, Any] | None:
    """
    AI toxicity + ICH M7 hybrid ensemble integrated degradation tree.
    """
    return generate_degradation_tree(
        current_smiles,
        max_depth=max_depth,
        _current_depth=_current_depth,
        _visited_smiles=_visited_smiles,
        include_toxicity=True,
    )


def degradation_tree_to_json(tree: dict[str, Any] | None, *, indent: int = 2) -> str:
    return json.dumps(tree or {}, indent=indent, ensure_ascii=False)


def degradation_tree_to_cytoscape(
    tree: dict[str, Any] | None,
) -> dict[str, list[dict[str, dict[str, Any]]]]:
    """
    Cytoscape.js 등 그래프 시각화용 elements 포맷.

    제안하신 NetworkX DiGraph + to_cytoscape_json() 과 동일한 용도이며,
    기존 nested JSON 트리(generate_degradation_tree)를 그대로 변환합니다.
    """
    nodes: list[dict[str, dict[str, Any]]] = []
    edges: list[dict[str, dict[str, Any]]] = []
    node_ids: set[str] = set()
    edge_ids: set[str] = set()

    def _walk(node: dict[str, Any] | None, parent_id: str | None) -> None:
        if not node:
            return
        nid = str(node.get("name") or "").strip()
        if not nid:
            return

        if nid not in node_ids:
            node_ids.add(nid)
            gen = int(node.get("generation", 0))
            tox = node.get("toxicity_analysis") or {}
            data: dict[str, Any] = {
                "id": nid,
                "label": nid,
                "generation": gen,
                "type": "Parent" if gen == 0 else "Impurity",
            }
            if node.get("triggered_by"):
                data["triggered_by"] = node["triggered_by"]
            if tox:
                data["ames_prob"] = tox.get("Ames_AI_Probability")
                data["ich_m7_class"] = tox.get("ICH_M7_Class")
                data["ich_m7_status"] = tox.get("ICH_M7_Status")
            nodes.append({"data": data})

        if parent_id and parent_id != nid:
            rule = str(node.get("triggered_by") or "degradation")
            eid = f"{parent_id}__{nid}__{rule}"
            if eid not in edge_ids:
                edge_ids.add(eid)
                edges.append(
                    {
                        "data": {
                            "id": eid,
                            "source": parent_id,
                            "target": nid,
                            "reaction_type": rule,
                        }
                    }
                )

        for child in node.get("children") or []:
            _walk(child, nid)

    _walk(tree, None)
    return {"nodes": nodes, "edges": edges}


def degradation_tree_cytoscape_json(tree: dict[str, Any] | None, *, indent: int = 2) -> str:
    return json.dumps(
        degradation_tree_to_cytoscape(tree),
        indent=indent,
        ensure_ascii=False,
    )


def flatten_degradation_tree(tree: dict[str, Any] | None) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def _walk(node: dict[str, Any] | None, parent: str | None) -> None:
        if not node:
            return
        tox = node.get("toxicity_analysis") or {}
        rows.append(
            {
                "SMILES": node.get("name", ""),
                "Generation": node.get("generation", 0),
                "Triggered_By": node.get("triggered_by", ""),
                "Parent_SMILES": parent or "",
                "Ames_AI_Probability": tox.get("Ames_AI_Probability", ""),
                "Prediction_Class": tox.get("Prediction_Class", ""),
                "ICH_M7_Class": tox.get("ICH_M7_Class", ""),
                "ICH_M7_Status": tox.get("ICH_M7_Status", ""),
            }
        )
        for child in node.get("children") or []:
            _walk(child, str(node.get("name", "")))

    _walk(tree, None)
    return rows


def predict_all_degradants(smiles: str) -> dict[str, list[dict[str, str]]]:
    parent_canon = _canonicalize_smiles(smiles)
    empty: dict[str, list[dict[str, str]]] = {
        "Hydrolysis": [],
        "Oxidation": [],
        "NDSRI": [],
    }
    if not parent_canon:
        return empty

    parent_mol = Chem.MolFromSmiles(parent_canon)
    if parent_mol is None:
        return empty

    try:
        Chem.SanitizeMol(parent_mol)
    except Exception:
        return empty

    results: dict[str, list[dict[str, str]]] = {k: [] for k in empty}
    seen: set[str] = {parent_canon}

    for rule_name, rule_info in shelf_life_degradation_rules().items():
        category = rule_info["category"]
        smirks = _normalize_smirks(rule_info["smirks"])
        try:
            rxn = AllChem.ReactionFromSmarts(smirks)
        except Exception:
            continue
        try:
            product_sets = rxn.RunReactants((parent_mol,))
        except Exception:
            continue

        for product_set in product_sets:
            for prod in product_set:
                prod_smiles = _mol_to_canonical_smiles(prod)
                if not prod_smiles or prod_smiles in seen:
                    continue
                # 2026-09-03 수정: 가수분해 산물은 방향족 수산화 컬링에서 면제
                #   (아스피린 → 살리실산이 여기서 통째로 사라지고 있었습니다).
                if not _passes_shelf_life_child_filter(parent_canon, prod_smiles, rule_name):
                    continue
                seen.add(prod_smiles)
                results.setdefault(category, []).append(
                    SmirksDegradantHit(rule_name, category, prod_smiles).to_dict()
                )

    return results


def smirks_degradant_smiles_by_category(
    smiles: str,
    *,
    categories: frozenset[str] | None = None,
) -> list[tuple[str, str, str]]:
    allowed = categories or frozenset({"Hydrolysis", "Oxidation"})
    out: list[tuple[str, str, str]] = []
    grouped = predict_all_degradants(smiles)
    for cat in allowed:
        for row in grouped.get(cat, []):
            smi = row.get("Degradant_SMILES", "")
            rule = row.get("Rule_Name", "")
            if smi and rule:
                out.append((rule, cat, smi))
    return out


def rule_mechanism_code(rule_name: str) -> str:
    return _RULE_MECHANISM_CODE.get(rule_name, "SMIRKS")


if __name__ == "__main__":
    sample_api = "CCOC(=O)CCSCCc1ccccc1"
    print("[PRIME-DS] AI toxicity + ICH M7 integrated degradation tree")
    toxic_tree = generate_toxic_degradation_tree(sample_api, max_depth=3)
    print(degradation_tree_to_json(toxic_tree))
