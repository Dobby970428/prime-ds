"""
PRIME-DS 스마트 조건부 분해 엔진 — RDKit 물리화학 연산 + 1만 건 오리지널 DB RAG 교차 보정.

Mock 템플릿 없이 Parent SMILES·가혹 조건에 따라 메커니즘별 분해산물을 동적 생성하고,
compounds / drug_primary 인덱스와 Tanimoto 유사도로 Confidence Score(0~100%)를 산출합니다.
"""

# 실제 생성
# 2-1. 분자를 10개 항목으로 진단
# has_ester    True    ← 아스피린은 여기 걸림
# has_amide    False
# has_lactam   False
# has_tertiary_amine  False
# has_thioether       False
# ...

# 2-2. 조건에 맞는 규칙만 골라 적용
# 산성 가수분해 고른 경우,
# 에스터 절단 규칙 발동
# [C:1](=[O:2])[O:3][#6:4]  >>  [C:1](=[O:2])O . [O:3][#6:4]
# 결과:
# CC(=O)Oc1ccccc1C(=O)O
#         ↓
# O=C(O)c1ccccc1O    (살리실산)
# CC(=O)O            (아세트산)

# 2-3. 중복·잘못된 것 제거
# deduplicate_hits          중복 3종 → 2종
# filter_valid_mechanism_hits  부모와 같은 것, 분자량이 안 바뀐 것 제거

# 2-4. 신뢰도 점수-"RAG(DB 유사도 조회, LLM 검색증강생성 아님)"
# 만든 구조를 DB 13,000개와 비교
# p23_chem_search_engine 에게 물어봄
#    ↓
# DB에 똑같은 게 있음   → 96점
# 비슷한 게 있음        → Tanimoto 유사도 %
# 아무것도 없음         → 기본점수 × 0.58 (최대 48점)
#    ↓
# 60°C 초과면 가산점
# 아스피린 결과는 47.6점

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from p23_chem_search_engine import (
    exact_match_search,
    safe_canonical_smiles,
    search_index_status,
    similarity_search_top_n,
)
from p41_degradation_smirks_rules import rule_mechanism_code, smirks_degradant_smiles_by_category

_LOG = logging.getLogger(__name__)

# ── SMARTS: 가혹 조건별 유기화학 필터 (ICH Q1A/Q1B 메커니즘 매칭) ─────────────
_ESTER_SMARTS = "[C:1](=[O:2])[O:3][#6:4]"
_AMIDE_SMARTS = "[C:1](=[O:2])[N:3]"
_ESTER_SUBSTRUCT = "[CX3](=O)[OX2H0]"
_AMIDE_SUBSTRUCT = "[CX3](=O)[NX3]"
# 2026-09-03 수정: p40의 복제본과 같은 결함 — 기존 "[#6;R]1[#6;R][#6;R][#6;R][#7;R]1=O"는
#   고리 질소에 =O(N=O)를 요구해 2-피롤리디논·δ-락탐·NMP·페니실린 G 어느 것도 잡지 못했고,
#   그 결과 아래 락탐 개환 경로(_generate_hydrolysis_hits)가 통째로 죽어 있었습니다.
#   고리 아미드(N-C=O)를 직접 기술하도록 교체했습니다(p40과 동일 표기).
_LACTAM_SMARTS = "[NX3;R][CX3;R](=O)"
_CARBAMATE_SMARTS = "[#7][#6](=O)O[#6]"
_AROMATIC_RING = "[c]1[c][c][c][c][c]1"
_TERTIARY_AMINE = "[NX3;H0;!$(NC=O)]([#6])([#6])[#6]"
_SECONDARY_AMINE = "[NX3;H1;!$(NC=O)][#6]"
_ALKENE = "[CX3]=[CX3]"
_THIOETHER = "[#16;!$(S=O)]"
_NITRO = "[NX3+](=O)[O-]"
_PHENOL = "[c][OX2H1]"
_ALDEHYDE = "[CX3H1](=O)"
_ARYL_HALIDE_SMARTS = "[c][F,Cl,Br,I]"
_BENZYLIC_BETA_SMARTS = "[c][CH2][C](=O)"

# ── RDKit 반응 SMARTS (보관 안정성: 가수분해 · 자동산화만) ─────────────────────
# 가수분해 — 사용자 지정 결합 절단식 (RunReactants → fragment SMILES)
_RXN_ESTER_HYDROLYSIS = "[C:1](=[O:2])[O:3][#6:4]>>[C:1](=[O:2])O.[O:3][#6:4]"
# CX3/CX4 typed ester cleavage → neutral carboxylic acid + alcohol fragments
_RXN_ESTER_HYDROLYSIS_NEUTRAL = (
    "[CX3:1](=[O:2])[OX2:3][#6:4]>>[CX3:1](=[O:2])[OX2H].[#6:4][OX2H]"
)
_RXN_AMIDE_HYDROLYSIS = "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])O.[N:3]"
# 2026-09-03 수정: 구 락탐 SMIRKS는 ① 반응물 SMARTS가 아무 락탐도 잡지 못하고
#   ② 반응물에 원자 맵이 하나도 없는데 산물이 :1,:2를 참조해, 매치되더라도 기질과 무관한
#   조각만 나왔으며 ③ 산물이 카복실레이트 음이온이라 DB 대조도 실패했습니다.
#   고리 N-C(=O) 결합만 끊고 나머지 고리는 이어진 채로 남는 개환(하나의 산물 성분)으로
#   재작성했습니다 — 2-피롤리디논 → GABA, 페니실린 G → 페니실로산.
_RXN_LACTAM_HYDROLYSIS = (
    "[NX3;R:1][CX3;R:2](=[OX1:3])>>([NX3:1].[CX3:2](=[OX1:3])[OX2H])"
)
# 자동산화 — S / 3차·2차 아민 N-oxide, 알데하이드 산화
_RXN_S_SIMPLE_OXIDATION = "[S:1]>>[S:1]=O"
# 2026-09-03 수정: 원자 매핑 :1은 대괄호 표현식 맨 뒤에 와야 합니다. 기존 표기는
#   ChemicalReactionParserException을 냈고, 그 예외가 run_reaction_product_matrix의
#   `except Exception: return []`에 삼켜져 로그 없이 빈 결과가 되었습니다
#   (N_OX·S_OX_ALT·FB_N_OX·FB_S_OX가 영구히 0건, 아민만 있는 원료는 광분해 조건에서도
#   산화 산물 0건). 위치를 교정했습니다.
_RXN_N_SIMPLE_OXIDATION = "[N;!$(N=O);!$(NC=O):1]>>[N:1]=O"
_RXN_S_OXIDATION = "[S;!$(S=O);!$(S(=O)=O):1]>>[S+:1][O-]"
_RXN_ALDEHYDE_OXIDATION = "[CX3H1:1](=O)>>[CX3:1](=O)[OH]"
_RXN_PHOTO_NORRISH = "[C:1](=[O:2])[O:3][#6:4]>>[C:1](=[O:2])[O:5].[#6:4][O:6]"

_HYDROLYTIC_STRESS = frozenset(
    {
        "산성 가혹 가수분해 (Acidic Hydrolysis)",
        "염기성 가혹 가수분해 (Alkaline Hydrolysis)",
        "수분 노출 (Neutral Hydrolysis)",
    }
)

# Morgan Fingerprint — 1만 건+ 내부 DB RAG 매트릭스 (radius=2, 2048-bit)
MORGAN_RADIUS = 2
MORGAN_BITS = 2048

# 2026-09-03 추가: 이 모듈이 쓰는 반응식 전체 목록 — validate()가 임포트 시점에
#   한 번 파싱 검증합니다(p18_excipient_reaction_rules.validate()와 같은 방식).
#   SMARTS 오타가 조용한 0건으로 묻히지 않게 하는 안전장치입니다.
ALL_REACTION_SMARTS: tuple[tuple[str, str], ...] = (
    ("_RXN_ESTER_HYDROLYSIS", _RXN_ESTER_HYDROLYSIS),
    ("_RXN_ESTER_HYDROLYSIS_NEUTRAL", _RXN_ESTER_HYDROLYSIS_NEUTRAL),
    ("_RXN_AMIDE_HYDROLYSIS", _RXN_AMIDE_HYDROLYSIS),
    ("_RXN_LACTAM_HYDROLYSIS", _RXN_LACTAM_HYDROLYSIS),
    ("_RXN_S_SIMPLE_OXIDATION", _RXN_S_SIMPLE_OXIDATION),
    ("_RXN_N_SIMPLE_OXIDATION", _RXN_N_SIMPLE_OXIDATION),
    ("_RXN_S_OXIDATION", _RXN_S_OXIDATION),
    ("_RXN_ALDEHYDE_OXIDATION", _RXN_ALDEHYDE_OXIDATION),
    ("_RXN_PHOTO_NORRISH", _RXN_PHOTO_NORRISH),
)

# 구조 판정용 SMARTS — 반응식과 함께 파싱 검증합니다.
ALL_STRUCTURE_SMARTS: tuple[tuple[str, str], ...] = (
    ("_ESTER_SMARTS", _ESTER_SMARTS),
    ("_AMIDE_SMARTS", _AMIDE_SMARTS),
    ("_ESTER_SUBSTRUCT", _ESTER_SUBSTRUCT),
    ("_AMIDE_SUBSTRUCT", _AMIDE_SUBSTRUCT),
    ("_LACTAM_SMARTS", _LACTAM_SMARTS),
    ("_CARBAMATE_SMARTS", _CARBAMATE_SMARTS),
    ("_AROMATIC_RING", _AROMATIC_RING),
    ("_TERTIARY_AMINE", _TERTIARY_AMINE),
    ("_SECONDARY_AMINE", _SECONDARY_AMINE),
    ("_ALKENE", _ALKENE),
    ("_THIOETHER", _THIOETHER),
    ("_NITRO", _NITRO),
    ("_PHENOL", _PHENOL),
    ("_ALDEHYDE", _ALDEHYDE),
    ("_ARYL_HALIDE_SMARTS", _ARYL_HALIDE_SMARTS),
    ("_BENZYLIC_BETA_SMARTS", _BENZYLIC_BETA_SMARTS),
)


def validate() -> dict[str, Any]:
    """반응식·구조 SMARTS가 전부 파싱되는지 확인 — 검증용.

    Returns
    -------
    {"total": int, "ok": [name, ...], "bad": [(name, reason), ...]}
    """
    ok: list[str] = []
    bad: list[tuple[str, str]] = []

    for name, smarts in ALL_REACTION_SMARTS:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is None:
                bad.append((name, "ReactionFromSmarts returned None"))
            else:
                ok.append(name)
        except Exception as exc:
            bad.append((name, str(exc)))

    for name, smarts in ALL_STRUCTURE_SMARTS:
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat is None:
                bad.append((name, "MolFromSmarts returned None"))
            else:
                ok.append(name)
        except Exception as exc:
            bad.append((name, str(exc)))

    return {
        "total": len(ALL_REACTION_SMARTS) + len(ALL_STRUCTURE_SMARTS),
        "ok": ok,
        "bad": bad,
    }


# 임포트 시점 1회 검증 — 실패한 식은 로그로 남깁니다(조용한 0건 방지).
_VALIDATION_REPORT = validate()
if _VALIDATION_REPORT["bad"]:
    for _name, _reason in _VALIDATION_REPORT["bad"]:
        _LOG.error("p42 SMARTS 파싱 실패: %s — %s", _name, _reason)


@dataclass
class ParentStructureAnalysis:
    """Parent SMILES RDKit 구조·물성 분석 결과."""

    smiles: str
    canonical_smiles: str
    parse_error: str | None = None
    molecular_weight: float = 0.0
    exact_molecular_weight: float = 0.0
    molecular_formula: str = ""
    has_ester: bool = False
    has_amide: bool = False
    has_lactam: bool = False
    has_carbamate: bool = False
    has_aromatic: bool = False
    has_tertiary_amine: bool = False
    has_secondary_amine: bool = False
    has_alkene: bool = False
    has_thioether: bool = False
    has_nitro: bool = False
    has_phenol: bool = False
    has_aldehyde: bool = False
    has_aryl_halide: bool = False
    has_benzylic_beta: bool = False
    labile_bond_count: int = 0


@dataclass
class MechanismHit:
    """단일 분해 경로 후보 (RDKit 연산 직후)."""

    degradant_id: str
    pathway_name: str
    name: str
    smiles: str
    mechanism_code: str
    mechanism_detail: str
    mechanism_base_score: float
    relative_abundance: str


def parse_parent_structure(smiles: str) -> tuple[Chem.Mol | None, ParentStructureAnalysis]:
    """
    Chem.MolFromSmiles 파싱 + ExactMolWt·SMARTS 기반 가혹 조건 필터 스캔.
    파싱 실패 시 parse_error에 메시지를 담아 반환합니다.
    """
    raw = (smiles or "").strip()
    empty = ParentStructureAnalysis(smiles=raw, canonical_smiles="", parse_error="SMILES가 비어 있습니다.")
    if not raw:
        return None, empty

    canon, err = safe_canonical_smiles(raw)
    if err or not canon:
        return None, ParentStructureAnalysis(
            smiles=raw,
            canonical_smiles="",
            parse_error=err or "RDKit SMILES 파싱 실패",
        )

    try:
        mol = Chem.MolFromSmiles(canon)
        if mol is None:
            return None, ParentStructureAnalysis(
                smiles=raw,
                canonical_smiles=canon,
                parse_error="Chem.MolFromSmiles가 None을 반환했습니다.",
            )
        Chem.SanitizeMol(mol)
    except Exception as exc:
        return None, ParentStructureAnalysis(
            smiles=raw,
            canonical_smiles=canon,
            parse_error=f"구조 정제(Sanitize) 실패: {exc}",
        )

    def _match(smarts: str) -> bool:
        pat = Chem.MolFromSmarts(smarts)
        return bool(pat and mol.HasSubstructMatch(pat))

    flags = {
        "has_ester": _match(_ESTER_SMARTS) or _match(_ESTER_SUBSTRUCT),
        "has_amide": _match(_AMIDE_SMARTS) or _match(_AMIDE_SUBSTRUCT),
        "has_lactam": _match(_LACTAM_SMARTS),
        "has_carbamate": _match(_CARBAMATE_SMARTS),
        "has_aromatic": _match(_AROMATIC_RING),
        "has_tertiary_amine": _match(_TERTIARY_AMINE),
        "has_secondary_amine": _match(_SECONDARY_AMINE),
        "has_alkene": _match(_ALKENE),
        "has_thioether": _match(_THIOETHER),
        "has_nitro": _match(_NITRO),
        "has_phenol": _match(_PHENOL),
        "has_aldehyde": _match(_ALDEHYDE),
        "has_aryl_halide": _match(_ARYL_HALIDE_SMARTS),
        "has_benzylic_beta": _match(_BENZYLIC_BETA_SMARTS),
    }
    labile = sum(
        1
        for k in (
            "has_ester",
            "has_amide",
            "has_lactam",
            "has_carbamate",
            "has_tertiary_amine",
            "has_alkene",
            "has_thioether",
            "has_nitro",
            "has_phenol",
            "has_aldehyde",
            "has_aryl_halide",
        )
        if flags[k]
    )

    try:
        mw = round(Descriptors.MolWt(mol), 2)
        exact_mw = round(Descriptors.ExactMolWt(mol), 4)
        formula = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        mw, exact_mw, formula = 0.0, 0.0, ""

    analysis = ParentStructureAnalysis(
        smiles=raw,
        canonical_smiles=canon,
        molecular_weight=mw,
        exact_molecular_weight=exact_mw,
        molecular_formula=formula,
        labile_bond_count=labile,
        **flags,
    )
    return mol, analysis


def exact_mw_from_smiles(smiles: str) -> float:
    """분해산물 SMILES → RDKit Descriptors.ExactMolWt (g/mol)."""
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return 0.0
    try:
        Chem.SanitizeMol(mol)
        return round(Descriptors.ExactMolWt(mol), 4)
    except Exception:
        return round(Descriptors.MolWt(mol), 2)


def morgan_fingerprint(mol: Chem.Mol) -> DataStructs.ExplicitBitVect | None:
    """Morgan bit-vector fingerprint (radius=2, nBits=2048)."""
    if mol is None:
        return None
    try:
        return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_BITS)
    except Exception:
        return None


def _is_valid_degradant_smiles(smiles: str) -> bool:
    """더미 map 번호 잔여([1H], [10H] 등), 파싱 불가, 라디칼 중간체 SMILES 제외."""
    from p44_degradation_shelf_life_filter import is_spin_stabilized_neutral_smiles

    if not smiles or not isinstance(smiles, str):
        return False
    if re.search(r"\[\d+H\]", smiles):
        return False
    smi = smiles.strip()
    if Chem.MolFromSmiles(smi) is None:
        return False
    return is_spin_stabilized_neutral_smiles(smi)


def mol_to_degradant_smiles(prod_mol: Chem.Mol) -> str | None:
    """
    RunReactants 생성물 1분자 → SanitizeMol → canonical SMILES.
    Product Matrix의 각 prod_mol을 독립 분해산물 SMILES로 변환합니다.
    """
    from p44_degradation_shelf_life_filter import is_spin_stabilized_neutral_mol

    if prod_mol is None or prod_mol.GetNumAtoms() == 0:
        return None
    if not is_spin_stabilized_neutral_mol(prod_mol):
        return None
    try:
        rw = Chem.RWMol(prod_mol)
        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() == 0:
                # 2026-09-03 수정: 더미를 H로 바꿀 때 동위원소 번호를 지우지 않아
                #   [3H]·[1H]가 남고 _is_valid_degradant_smiles의 \[\d+H\] 정규식에
                #   전부 걸렸습니다. 동위원소 라벨도 함께 제거합니다.
                atom.SetAtomicNum(1)
                atom.SetIsotope(0)
        Chem.SanitizeMol(rw)
        smi = Chem.MolToSmiles(rw, canonical=True)
        return smi if _is_valid_degradant_smiles(smi) else None
    except Exception:
        pass
    try:
        Chem.SanitizeMol(prod_mol)
        smi = Chem.MolToSmiles(prod_mol, canonical=True)
        return smi if _is_valid_degradant_smiles(smi) else None
    except Exception:
        return None


# 하위 호환 별칭
sanitize_reaction_fragment = mol_to_degradant_smiles


def exclude_parent_smiles(
    smiles_list: list[str],
    parent_canon: str,
) -> list[str]:
    """입력 Parent와 동일한 canonical SMILES 제거 (자기 자신 필터)."""
    if not parent_canon:
        return list(dict.fromkeys(smiles_list))
    out: list[str] = []
    seen: set[str] = set()
    for smi in smiles_list:
        if not smi or smi in seen:
            continue
        if not _is_valid_degradant_smiles(smi):
            continue
        if _is_parent_smiles(smi, parent_canon):
            continue
        seen.add(smi)
        out.append(smi)
    return out


def run_reaction_product_matrix(
    parent_mol: Chem.Mol,
    rxn_smarts: str,
    *,
    parent_canon: str = "",
    max_unique: int = 32,
) -> list[str]:
    """
    AllChem.RunReactants 2차원 Product Matrix 완전 순회.

    반환 형태: ((product_1, product_2), (product_3, product_4), ...)
    각 prod_set 내부의 **모든** prod_mol을 추출 — 소분자 파편 + 거대 모핵 잔기 모두 포함.
    input_smiles(Parent)와 동일한 SMILES는 제외합니다.
    """
    if parent_mol is None:
        return []

    # 2026-09-03 수정: 기존에는 ReactionFromSmarts 실패(SMARTS 오타 등)까지 통째로
    #   `except Exception: return []`에 삼켜져 로그 한 줄 없이 영구 0건이 되었습니다.
    #   파싱 실패는 삼키지 않고 로그를 남깁니다(RunReactants 실패와 구분).
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    except Exception as exc:
        _LOG.error("반응식 파싱 실패 — SMARTS=%r: %s", rxn_smarts, exc)
        return []
    if rxn is None:
        _LOG.error("반응식 파싱 실패 — ReactionFromSmarts가 None 반환: %r", rxn_smarts)
        return []

    try:
        product_matrix = rxn.RunReactants((parent_mol,))
    except Exception as exc:
        _LOG.debug("RunReactants 실패 — SMARTS=%r: %s", rxn_smarts, exc)
        return []

    if not product_matrix:
        return []

    collected: list[str] = []
    for prod_set in product_matrix:
        if prod_set is None:
            continue
        for prod_mol in prod_set:
            smi = mol_to_degradant_smiles(prod_mol)
            if smi:
                collected.append(smi)

    unique = exclude_parent_smiles(collected, parent_canon)
    unique.sort(key=exact_mw_from_smiles, reverse=True)
    if max_unique and len(unique) > max_unique:
        unique = unique[:max_unique]
    return unique


def _mol_has_ester(mol: Chem.Mol) -> bool:
    pat = Chem.MolFromSmarts(_ESTER_SMARTS)
    return bool(pat and mol.HasSubstructMatch(pat))


def _mol_has_amide(mol: Chem.Mol) -> bool:
    pat = Chem.MolFromSmarts(_AMIDE_SMARTS)
    return bool(pat and mol.HasSubstructMatch(pat))


def temperature_linear_boost_celsius(temp_c: float) -> float:
    """가혹 온도 > 60°C 시 신뢰도 점수 선형 가중 (최대 +15%)."""
    if float(temp_c) <= 60.0:
        return 0.0
    return min(15.0, (float(temp_c) - 60.0) * 0.3)


def _is_parent_smiles(smiles: str, parent_canon: str) -> bool:
    if not parent_canon:
        return False
    return _canonical_smiles(smiles) == _canonical_smiles(parent_canon)


def _run_reaction(
    parent_mol: Chem.Mol,
    rxn_smarts: str,
    *,
    max_products: int = 32,
    parent_canon: str = "",
) -> list[str]:
    """RunReactants Product Matrix — run_reaction_product_matrix 래퍼."""
    return run_reaction_product_matrix(
        parent_mol,
        rxn_smarts,
        parent_canon=parent_canon,
        max_unique=max_products,
    )


def _mol_weight(smiles: str) -> float:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0
    return float(Descriptors.MolWt(mol))


# ---------------------------------------------------------------------------
# 2026-09-03 삭제: 결합 절단(FragmentOnBonds) 기반 가수분해 fallback
#
#   `_bond_indices_from_smarts` / `_fragments_from_bond_cleavage` /
#   `_ester_bond_indices` / `_amide_bond_indices` 를 제거했습니다.
#
#   왜 지웠나 — 이 경로는 100% 무효였습니다.
#     ① FragmentOnBonds(addDummies=True)는 더미 원자에 결합 양끝 원자 인덱스를
#        **동위원소 번호**로 붙입니다. 여기에 SetAtomicNum(1)을 하면 [3H]·[1H]가 되고,
#        `_is_valid_degradant_smiles`의 `\[\d+H\]` 정규식이 이를 전부 거부합니다.
#        실측: 아스피린·파라세타몰 모두 결과가 [] 였습니다.
#     ② 설령 통과했더라도 아실 조각 [3H]C(C)=O 는 아세트알데하이드이지 아세트산이
#        아니므로 화학적으로도 틀린 산물입니다(가수분해는 아실 측에 OH가 붙습니다).
#
#   가수분해는 `_RXN_ESTER_HYDROLYSIS_NEUTRAL`(중성 카복실산 + 알코올) 및
#   `_RXN_AMIDE_HYDROLYSIS` SMIRKS로 일원화했습니다 — 실제로도 그것만 동작했습니다.
# ---------------------------------------------------------------------------


def _fragment_label(smiles: str) -> str:
    """분해 fragment SMILES → 간단 한국어 명칭."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "분해 fragment"
    carboxyl = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    alcohol = Chem.MolFromSmarts("[OX2H1]")
    amine = Chem.MolFromSmarts("[NX3;!$(NC=O)]")
    if carboxyl and mol.HasSubstructMatch(carboxyl):
        return "카르복실산 잔기"
    if amine and mol.HasSubstructMatch(amine):
        return "아민/아마인 잔기"
    if alcohol and mol.HasSubstructMatch(alcohol):
        return "알코올/페놀릭 잔기"
    return "유기 분해 fragment"


def _degradant_display_name(smiles: str, parent_mw: float) -> str:
    """fragment SMILES → 주형 모핵 잔기 vs 소분자 파편 명칭 (ExactMolWt 기준)."""
    fw = exact_mw_from_smiles(smiles)
    role = _fragment_label(smiles)
    if parent_mw > 0 and fw >= parent_mw * 0.45:
        return f"주형 모핵 잔기 — {role}"
    if fw <= 90:
        return f"소분자 파편 — {role}"
    return f"분해산물 잔기 — {role}"


def predict_hydrolysis_degradants(smiles: str) -> list[str]:
    """
    SMILES 내 에스테르(Ester) 결합을 탐색하여 가수분해 분해산물 SMILES 목록을 반환합니다.

    CX3/CX4 typed SMIRKS — 카르복실산 + 알코올 중성 fragment (RunReactants).
    """
    raw = (smiles or "").strip()
    if not raw:
        return []
    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        return []
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return []
    parent_canon, _ = safe_canonical_smiles(raw)
    return run_reaction_product_matrix(
        mol,
        _RXN_ESTER_HYDROLYSIS_NEUTRAL,
        parent_canon=parent_canon or raw,
    )


def _hydrolysis_fragments(mol: Chem.Mol, *, mode: str, parent_canon: str = "") -> list[str]:
    """
    가수분해 — SMIRKS Product Matrix 전체 추출 (아실 잔기 + 알코올/아민 잔기).

    2026-09-03 수정: 무효였던 결합 절단(FragmentOnBonds) fallback을 제거하고
    SMIRKS 경로로 일원화했습니다. 에스터는 중성 카복실산을 내는
    `_RXN_ESTER_HYDROLYSIS_NEUTRAL`을 함께 적용해 DB 대조가 가능한 형태로 만듭니다.
    """
    rxn = _RXN_ESTER_HYDROLYSIS if mode == "ester" else _RXN_AMIDE_HYDROLYSIS
    rxn_frags = run_reaction_product_matrix(mol, rxn, parent_canon=parent_canon)

    neutral_ester: list[str] = []
    if mode == "ester" and parent_canon:
        neutral_ester = predict_hydrolysis_degradants(parent_canon)

    merged = exclude_parent_smiles(rxn_frags + neutral_ester, parent_canon)

    # 다중 에스테르: 1차 RunReactants 모핵 잔기에 잔존 ester/amide가 있으면 2차 절단 (1단계만)
    cascade: list[str] = []
    has_check = _mol_has_ester if mode == "ester" else _mol_has_amide
    for seed in rxn_frags:
        if not _is_valid_degradant_smiles(seed):
            continue
        sub = Chem.MolFromSmiles(seed)
        if sub is None or not has_check(sub):
            continue
        cascade.extend(
            run_reaction_product_matrix(sub, rxn, parent_canon=parent_canon, max_unique=8)
        )

    merged = exclude_parent_smiles(merged + cascade, parent_canon)
    merged.sort(key=exact_mw_from_smiles, reverse=True)
    return merged


def _abundance_for_rank(rank: int, stress: str) -> str:
    if rank == 0:
        return "주요 (Major, ~40–55%)"
    if rank == 1:
        return "부분 (Minor, ~20–35%)"
    if stress in _HYDROLYTIC_STRESS and "Neutral" in stress:
        return "미량 (Trace, ~5–15%)"
    return "미량 (Trace, <15%)"


def _halogen_count(smiles: str) -> int:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return 0
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (9, 17, 35, 53))


def _has_aryl_halide_smiles(smiles: str) -> bool:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return False
    pat = Chem.MolFromSmarts(_ARYL_HALIDE_SMARTS)
    return bool(pat and mol.HasSubstructMatch(pat))


def _canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return (smiles or "").strip()
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return (smiles or "").strip()


def _transform_smiles(smiles: str, rxn_smarts: str, parent_canon: str = "") -> list[str]:
    """단일 SMILES에 반응 SMARTS Product Matrix 적용."""
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return []
    return run_reaction_product_matrix(mol, rxn_smarts, parent_canon=parent_canon)


def _generate_hydrolysis_hits(
    parent_mol: Chem.Mol,
    analysis: ParentStructureAnalysis,
    stress: str,
    *,
    id_prefix: str,
) -> list[MechanismHit]:
    """산성/염기성/중성 가수분해 — 에스테르·아미드·락탐·카바메이트 동적 절단."""
    hits: list[MechanismHit] = []
    seq = 0
    base_score = 78.0 if "Acidic" in stress or "Alkaline" in stress else 68.0
    parent_canon = analysis.canonical_smiles

    if analysis.has_ester:
        frags = _hydrolysis_fragments(parent_mol, mode="ester", parent_canon=parent_canon)
        for smi in frags:
            hits.append(
                MechanismHit(
                    degradant_id=f"{id_prefix}-{seq:02d}",
                    pathway_name=f"에스테르 가수분해 ({stress.split('(')[-1].rstrip(')')})",
                    name=_degradant_display_name(smi, analysis.molecular_weight),
                    smiles=smi,
                    mechanism_code="EST_HYD",
                    mechanism_detail=(
                        "Ester C(=O)-O bond cleavage — acyl + alcohol/phenolic pair products"
                    ),
                    mechanism_base_score=base_score + 4.0,
                    relative_abundance=_abundance_for_rank(seq, stress),
                )
            )
            seq += 1

    if analysis.has_amide:
        frags = _hydrolysis_fragments(parent_mol, mode="amide", parent_canon=parent_canon)
        for smi in frags:
            hits.append(
                MechanismHit(
                    degradant_id=f"{id_prefix}-{seq:02d}",
                    pathway_name=f"아미드 가수분해 ({stress.split('(')[-1].rstrip(')')})",
                    name=_degradant_display_name(smi, analysis.molecular_weight),
                    smiles=smi,
                    mechanism_code="AMD_HYD",
                    mechanism_detail="Amide C(=O)-N bond cleavage — acyl + amine pair products",
                    mechanism_base_score=base_score + 2.0,
                    relative_abundance=_abundance_for_rank(seq, stress),
                )
            )
            seq += 1

    if analysis.has_lactam:
        frags = _run_reaction(
            parent_mol, _RXN_LACTAM_HYDROLYSIS, parent_canon=parent_canon
        )
        for smi in frags:
            hits.append(
                MechanismHit(
                    degradant_id=f"{id_prefix}-{seq:02d}",
                    pathway_name="락탐 고리 개환 가수분해 (Lactam Ring Opening)",
                    name=f"락탐 개환 — {_fragment_label(smi)}",
                    smiles=smi,
                    mechanism_code="LAC_HYD",
                    mechanism_detail="Lactam ring-opening hydrolysis",
                    mechanism_base_score=base_score + 1.0,
                    relative_abundance=_abundance_for_rank(seq, stress),
                )
            )
            seq += 1

    if analysis.has_carbamate and not analysis.has_ester:
        frags = _hydrolysis_fragments(parent_mol, mode="ester", parent_canon=parent_canon)
        for smi in frags:
            hits.append(
                MechanismHit(
                    degradant_id=f"{id_prefix}-{seq:02d}",
                    pathway_name="카바메이트 가수분해 (Carbamate Hydrolysis)",
                    name=_degradant_display_name(smi, analysis.molecular_weight),
                    smiles=smi,
                    mechanism_code="CAR_HYD",
                    mechanism_detail="Carbamate N-C(=O)-O cleavage (ester-type SMARTS)",
                    mechanism_base_score=base_score,
                    relative_abundance=_abundance_for_rank(seq, stress),
                )
            )
            seq += 1

    return hits


def _generate_oxidative_hits(
    parent_mol: Chem.Mol,
    analysis: ParentStructureAnalysis,
    *,
    id_prefix: str,
) -> list[MechanismHit]:
    """과산화·자동산화 — N-oxide, S=O, 알데하이드 산화 (보관 유연물질)."""
    hits: list[MechanismHit] = []
    seq = 0
    parent_canon = analysis.canonical_smiles
    rxn_specs: list[tuple[str, str, str, str, float, bool]] = [
        (
            "N_OX",
            _RXN_N_SIMPLE_OXIDATION,
            "N-자동산화 (N-Oxide Formation)",
            "N-옥사이드 자동산화체",
            76.0,
            analysis.has_tertiary_amine or analysis.has_secondary_amine,
        ),
        (
            "S_OX",
            _RXN_S_SIMPLE_OXIDATION,
            "S-자동산화 (S=O Formation)",
            "설폭시드/아울폭사이드형 산화체",
            74.0,
            analysis.has_thioether,
        ),
        (
            "S_OX_ALT",
            _RXN_S_OXIDATION,
            "S-자동산화 (Sulfoxide Formation)",
            "설폭시드 산화체",
            70.0,
            analysis.has_thioether,
        ),
        (
            "ALD_OX",
            _RXN_ALDEHYDE_OXIDATION,
            "알데하이드 자동산화 (Aldehyde → Carboxylic acid)",
            "카르복실산 산화 분해체",
            66.0,
            analysis.has_aldehyde,
        ),
    ]

    for code, rxn, pathway, name, score, enabled in rxn_specs:
        if not enabled:
            continue
        products = _run_reaction(parent_mol, rxn, parent_canon=parent_canon)
        for smi in products:
            hits.append(
                MechanismHit(
                    degradant_id=f"{id_prefix}-{seq:02d}",
                    pathway_name=pathway,
                    name=name,
                    smiles=smi,
                    mechanism_code=code,
                    mechanism_detail=f"Auto-oxidation — {code} RDKit reaction",
                    mechanism_base_score=score,
                    relative_abundance=_abundance_for_rank(seq, "과산화 산화 조건 (Oxidative Stress)"),
                )
            )
            seq += 1

    return hits


def _generate_photolytic_hits(
    parent_mol: Chem.Mol,
    analysis: ParentStructureAnalysis,
    *,
    id_prefix: str,
) -> list[MechanismHit]:
    """광분해 — Norrish-type 절단 + 광산화."""
    hits: list[MechanismHit] = []
    seq = 0

    if analysis.has_ester:
        for smi in _run_reaction(
            parent_mol, _RXN_PHOTO_NORRISH, parent_canon=analysis.canonical_smiles
        ):
            if _is_parent_smiles(smi, analysis.canonical_smiles):
                continue
            hits.append(
                MechanismHit(
                    degradant_id=f"{id_prefix}-{seq:02d}",
                    pathway_name="Norrish Type I 광화학 절단",
                    name=f"광분해 fragment — {_fragment_label(smi)}",
                    smiles=smi,
                    mechanism_code="PHOTO_NOR",
                    mechanism_detail="Photolytic ester bond homolysis/cleavage",
                    mechanism_base_score=72.0,
                    relative_abundance=_abundance_for_rank(seq, "광분해"),
                )
            )
            seq += 1

    if analysis.has_aromatic or analysis.has_aldehyde:
        ox = _generate_oxidative_hits(parent_mol, analysis, id_prefix=f"{id_prefix}-ox")
        for h in ox[:2]:
            h.degradant_id = f"{id_prefix}-{seq:02d}"
            h.pathway_name = "광산화 (Photo-oxidation) — " + h.pathway_name
            h.mechanism_base_score = max(55.0, h.mechanism_base_score - 8.0)
            hits.append(h)
            seq += 1

    return hits


def confidence_label_ko(score: float) -> str:
    """DB Morgan Tanimoto 매칭률(%) 기반 신뢰도 등급 (한국어)."""
    if score >= 85.0:
        return "🔥 실제 검출 이력과 구조 90% 이상 일치 (검증 완료)"
    if score >= 50.0:
        return "⚠️ 화학 메커니즘 및 유사 작용기 유효성 확인"
    return "💡 정보화학 모델 기반 가상 시뮬레이션 결과"


def morgan_tanimoto_db_match(
    degradant_smiles: str,
    *,
    top_n: int = 8,
) -> tuple[float, str]:
    """
    Morgan FP (2048-bit) · DataStructs.TanimotoSimilarity 기반 내부 DB Top-N 매칭.

    Returns
    -------
    (match_pct 0–100, db_detail_ko)
    """
    raw = (degradant_smiles or "").strip()
    if not raw:
        return 0.0, "SMILES 없음 — DB 매칭 생략"

    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        return 0.0, "RDKit 구조 해석 실패 — DB 매칭 생략"

    try:
        Chem.SanitizeMol(mol)
        _ = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_BITS)
    except Exception:
        return 0.0, "Morgan Fingerprint 생성 실패"

    try:
        exact = exact_match_search(raw)
        if exact:
            name = exact.get("name") or exact.get("drug_name") or "indexed compound"
            cid = exact.get("cid") or ""
            detail = f"Exact SMILES 일치 — {name}"
            if cid:
                detail += f" (CID {cid})"
            return 96.0, detail
    except Exception:
        pass

    try:
        sims, err = similarity_search_top_n(raw, top_n=top_n)
        if err:
            return 0.0, f"DB FP 매트릭스: {err}"
        if not sims:
            return 0.0, "DB FP 인덱스 미매칭 — 메커니즘 연산만 적용"
        best = sims[0]
        pct = float(best.get("tanimoto_pct", 0.0))
        ref = (
            best.get("drug_name")
            or best.get("compound_name")
            or best.get("drug_file")
            or "유사 화합물"
        )
        # 2026-09-03 수정: ① search_index_status()가 반환하는 키는
        #   drug_files / drug_primary / compounds_indexed / substruct_ready / db_path 이며
        #   "compounds"는 존재하지 않아 fallback이 항상 0이었습니다(화면에 `DB ~0건`).
        #   → compounds_indexed 로 교정.
        #   ② 최초 실행 시 이 호출이 `no such table` 예외를 내는데, 바깥 try에 묶여 있어
        #   이미 계산한 Tanimoto 결과까지 통째로 버려졌습니다. → 별도 try로 분리.
        db_n = 0
        try:
            idx_status = search_index_status()
            db_n = idx_status.get("drug_primary", 0) or idx_status.get("compounds_indexed", 0)
        except Exception as exc:
            _LOG.debug("search_index_status() 실패 — DB 건수 표기 생략: %s", exc)

        detail = f"Morgan Tanimoto {pct:.1f}% — {ref} (Top-1 / {len(sims)}"
        detail += f", DB ~{db_n:,}건 FP 매트릭스)" if db_n else ", DB 건수 미상)"
        return pct, detail
    except Exception as exc:
        return 0.0, f"DB RAG 연산 오류: {exc}"


def compute_rag_confidence(
    degradant_smiles: str,
    *,
    mechanism_base_score: float,
    stress_environment: str,
    temperature_c: float,
    mechanism_code: str,
) -> tuple[float, str, str]:
    """
    Morgan FP RAG + 가혹 조건 가중 — Confidence Score (0~100%) 산출.

    - DB Tanimoto 매칭률을 핵심 점수로 사용
    - 온도 > 60°C: 선형 가중 상승
    - DB 미구축 시 메커니즘 점수 기반 보수적 추정
    """
    _ = stress_environment, mechanism_code  # 향후 조건별 가중 확장용
    match_pct, db_detail = morgan_tanimoto_db_match(degradant_smiles, top_n=8)
    temp_boost = temperature_linear_boost_celsius(temperature_c)

    # 2026-09-03 수정: 점수가 역전되어 있었습니다.
    #   DB 매칭 0건      -> 45.2 + 12.0 = 57.2
    #   Tanimoto 5% 매칭 ->  5.0 + 12.0 = 17.0   ← DB 근거가 더 많은데 점수가 더 낮음
    #   메커니즘 점수를 하한으로 잡아, DB 매칭이 있어도 그 아래로 내려가지 않게 합니다.
    mech_base = min(48.0, float(mechanism_base_score) * 0.58)
    score = min(100.0, round(max(match_pct, mech_base) + temp_boost, 1))
    return score, confidence_label_ko(score), db_detail


# 2026-09-03 추가: fallback을 가혹 조건별로 분리합니다.
#   기존에는 조건과 무관하게 에스터/아미드 가수분해 + N/S 산화를 전부 돌려,
#   메티오닌 「산성 가혹 가수분해」에서 설폭사이드 산화체가, 아스피린 「과산화 산화 조건」에서
#   에스터 가수분해 산물이 나왔습니다(선택한 가혹 조건을 무시).
_FALLBACK_RXN_HYDROLYTIC: tuple[tuple[str, str, str, str, float], ...] = (
    ("FB_EST", _RXN_ESTER_HYDROLYSIS, "범용 에스테르 가수분해 (Universal Ester Cleavage)", "에스테르 절단 분해체", 52.0),
    ("FB_AMD", _RXN_AMIDE_HYDROLYSIS, "범용 아미드 가수분해 (Universal Amide Cleavage)", "아미드 절단 분해체", 50.0),
)

_FALLBACK_RXN_OXIDATIVE: tuple[tuple[str, str, str, str, float], ...] = (
    ("FB_N_OX", _RXN_N_SIMPLE_OXIDATION, "질소 자동산화 (N-Oxidation Fallback)", "N-옥사이드 자동산화체", 45.0),
    ("FB_S_OX", _RXN_S_SIMPLE_OXIDATION, "황 자동산화 (S-Oxidation Fallback)", "설폭시드형 산화체", 45.0),
)

NO_APPLICABLE_MECHANISM_MESSAGE = "해당 조건에서 적용 가능한 기전 없음"


def fallback_reactions_for_stress(
    stress: str,
) -> tuple[tuple[str, str, str, str, float], ...]:
    """가혹 조건에 화학적으로 대응하는 fallback 반응만 반환."""
    if stress in _HYDROLYTIC_STRESS:
        return _FALLBACK_RXN_HYDROLYTIC
    if stress == "과산화 산화 조건 (Oxidative Stress)":
        return _FALLBACK_RXN_OXIDATIVE
    if stress == "광분해 조건 (Photolytic Degradation)":
        # 광분해는 Norrish형 절단(에스터)과 광산화가 모두 성립합니다.
        return _FALLBACK_RXN_HYDROLYTIC + _FALLBACK_RXN_OXIDATIVE
    return ()


def no_applicable_mechanism_message(stress: str) -> str:
    """조건에 맞는 기전이 하나도 없을 때 화면·보고서에 쓸 한국어 안내."""
    label = (stress or "").strip() or "선택한 가혹 조건"
    return f"{label} — {NO_APPLICABLE_MECHANISM_MESSAGE}"


def _generate_diverse_fallback_hits(
    parent_mol: Chem.Mol,
    analysis: ParentStructureAnalysis,
    stress: str,
    *,
    id_prefix: str,
    notes: list[str] | None = None,
) -> list[MechanismHit]:
    """
    작용기 미검출·반응 실패 시의 fallback — **선택한 가혹 조건에 맞는 기전만** 시도합니다.

    조건에 맞는 기전이 하나도 없으면 빈 목록을 반환하고,
    `notes`가 주어지면 "해당 조건에서 적용 가능한 기전 없음" 메시지를 담아 줍니다.
    """
    parent_canon = analysis.canonical_smiles
    hits: list[MechanismHit] = []
    seq = 0
    seen: set[str] = set()

    def _add(code: str, pathway: str, name: str, smi: str, detail: str, score: float, abund: str) -> None:
        nonlocal seq
        canon = _canonical_smiles(smi)
        if not canon or _is_parent_smiles(canon, parent_canon) or canon in seen:
            return
        seen.add(canon)
        hits.append(
            MechanismHit(
                degradant_id=f"{id_prefix}-{seq:02d}",
                pathway_name=pathway,
                name=name,
                smiles=canon,
                mechanism_code=code,
                mechanism_detail=detail,
                mechanism_base_score=score,
                relative_abundance=abund,
            )
        )
        seq += 1

    fallback_rxns = fallback_reactions_for_stress(stress)
    if not fallback_rxns:
        msg = no_applicable_mechanism_message(stress)
        _LOG.info("%s", msg)
        if notes is not None:
            notes.append(msg)
        return []

    for code, rxn, pathway, name, score in fallback_rxns:
        for smi in _run_reaction(parent_mol, rxn, parent_canon=parent_canon):
            _add(
                code,
                pathway,
                f"{name} — {_fragment_label(smi)}",
                smi,
                f"Fallback RDKit RunReactants — {code}",
                score,
                "미량 (Trace, ~8–15%)",
            )

    if not hits:
        msg = no_applicable_mechanism_message(stress)
        _LOG.info("%s", msg)
        if notes is not None:
            notes.append(msg)

    return hits


def _generate_metabolic_fallback_hits(
    parent_mol: Chem.Mol,
    analysis: ParentStructureAnalysis,
    stress: str,
    *,
    id_prefix: str,
    notes: list[str] | None = None,
) -> list[MechanismHit]:
    """하위 호환 래퍼 — Parent 복사 없는 조건별 fallback."""
    return _generate_diverse_fallback_hits(
        parent_mol, analysis, stress, id_prefix=id_prefix, notes=notes
    )


def _generate_smirks_catalog_hits(
    analysis: ParentStructureAnalysis,
    stress: str,
    *,
    id_prefix: str,
    seq_start: int = 0,
    categories: frozenset[str],
) -> list[MechanismHit]:
    """Universal SMIRKS rule catalog — hydrolysis / autoxidation degradants."""
    hits: list[MechanismHit] = []
    seq = seq_start
    base_score = 74.0 if "Acidic" in stress or "Alkaline" in stress else 68.0
    ox_score = 76.0 if "Oxidative" in stress else 70.0

    _pathway_labels = {
        "Hydrolysis_Ester": "SMIRKS 에스테르 가수분해",
        "Hydrolysis_Amide": "SMIRKS 아미드 가수분해",
        "Hydrolysis_Carbamate": "SMIRKS 카바메이트 가수분해",
        "Oxidation_Thioether": "SMIRKS 티오에테르 자동산화 (Sulfoxide)",
        "Oxidation_Tertiary_Amine": "SMIRKS 3급 아민 N-Oxide",
        # 2026-09-03 수정: p41 규칙 이름 변경에 맞춰 키 갱신
        #   (구 Nitrosamination_NDSRIs — 실제 동작은 생성이 아니라 탈니트로소화였습니다).
        #   Oxidation_Ether_Dealkylation은 metabolic=True로 실행 규칙에서 빠졌으나,
        #   과거 캐시/보고서 재표시를 위해 라벨은 남겨 둡니다.
        "Oxidation_Ether_Dealkylation": "SMIRKS 에테르 탈메틸화 (Autoxidation)",
        "Denitrosation_NDSRI": "SMIRKS 니트로사민 탈니트로소화 (Denitrosation)",
    }

    for rule_name, category, smi in smirks_degradant_smiles_by_category(
        analysis.canonical_smiles,
        categories=categories,
    ):
        code = rule_mechanism_code(rule_name)
        score = ox_score if category == "Oxidation" else base_score
        hits.append(
            MechanismHit(
                degradant_id=f"{id_prefix}-{seq:02d}",
                pathway_name=_pathway_labels.get(rule_name, f"SMIRKS {rule_name}"),
                name=_degradant_display_name(smi, analysis.molecular_weight),
                smiles=smi,
                mechanism_code=code,
                mechanism_detail=f"Universal SMIRKS catalog — {rule_name} ({category})",
                mechanism_base_score=score,
                relative_abundance=_abundance_for_rank(seq, stress),
            )
        )
        seq += 1
    return hits


def generate_mechanism_hits(
    parent_mol: Chem.Mol,
    analysis: ParentStructureAnalysis,
    stress_environment: str,
    *,
    notes: list[str] | None = None,
) -> list[MechanismHit]:
    """가혹 환경별 메커니즘 매칭 → 분해산물 후보 목록."""
    prefix_map = {
        "산성 가혹 가수분해 (Acidic Hydrolysis)": "SD-AH",
        "염기성 가혹 가수분해 (Alkaline Hydrolysis)": "SD-AL",
        "수분 노출 (Neutral Hydrolysis)": "SD-NH",
        "과산화 산화 조건 (Oxidative Stress)": "SD-OX",
        "광분해 조건 (Photolytic Degradation)": "SD-PH",
    }
    prefix = prefix_map.get(stress_environment, "SD-XX")

    if stress_environment in _HYDROLYTIC_STRESS:
        hits = _generate_hydrolysis_hits(parent_mol, analysis, stress_environment, id_prefix=prefix)
        hits.extend(
            _generate_smirks_catalog_hits(
                analysis,
                stress_environment,
                id_prefix=prefix,
                seq_start=len(hits),
                categories=frozenset({"Hydrolysis"}),
            )
        )
    elif stress_environment == "과산화 산화 조건 (Oxidative Stress)":
        hits = _generate_oxidative_hits(parent_mol, analysis, id_prefix=prefix)
        hits.extend(
            _generate_smirks_catalog_hits(
                analysis,
                stress_environment,
                id_prefix=prefix,
                seq_start=len(hits),
                categories=frozenset({"Oxidation"}),
            )
        )
    elif stress_environment == "광분해 조건 (Photolytic Degradation)":
        hits = _generate_photolytic_hits(parent_mol, analysis, id_prefix=prefix)
        hits.extend(
            _generate_smirks_catalog_hits(
                analysis,
                stress_environment,
                id_prefix=prefix,
                seq_start=len(hits),
                categories=frozenset({"Oxidation"}),
            )
        )
    else:
        hits = []

    if not hits:
        # 2026-09-03 수정: fallback도 선택한 가혹 조건에 맞는 기전만 시도합니다.
        hits = _generate_diverse_fallback_hits(
            parent_mol, analysis, stress_environment, id_prefix=prefix, notes=notes
        )
    return hits


def deduplicate_hits(
    hits: list[MechanismHit],
    *,
    parent_canonical: str = "",
) -> list[MechanismHit]:
    """Canonical SMILES 기준 중복 제거 + Parent 동일 hit 전면 제거."""
    parent_canon = _canonical_smiles(parent_canonical) if parent_canonical else ""
    seen: set[str] = set()
    out: list[MechanismHit] = []
    for h in hits:
        mol = Chem.MolFromSmiles(h.smiles)
        if mol is None:
            continue
        canon = Chem.MolToSmiles(mol, canonical=True)
        if parent_canon and _is_parent_smiles(canon, parent_canon):
            continue
        if h.mechanism_code == "DEHAL" and parent_canon and _halogen_count(canon) >= _halogen_count(parent_canon):
            continue
        if canon in seen:
            continue
        seen.add(canon)
        out.append(h)
    return out


def filter_valid_mechanism_hits(
    hits: list[MechanismHit],
    parent_canonical: str,
) -> list[MechanismHit]:
    """Parent 동일·탈할로겐 미적용 등 잘못된 hit 제거 (Parent 복사 절대 불가)."""
    parent_canon = _canonical_smiles(parent_canonical)
    parent_h = _halogen_count(parent_canon)
    parent_mw = exact_mw_from_smiles(parent_canon) if parent_canon else 0.0
    out: list[MechanismHit] = []
    for h in hits:
        canon = _canonical_smiles(h.smiles)
        if not canon:
            continue
        if _is_parent_smiles(canon, parent_canon):
            continue
        if h.mechanism_code == "DEHAL" and _halogen_count(canon) >= parent_h:
            continue
        if "탈할로" in h.name and _has_aryl_halide_smiles(canon):
            continue
        frag_mw = exact_mw_from_smiles(canon)
        if parent_mw > 0 and abs(frag_mw - parent_mw) < 0.05 and h.mechanism_code not in ("FB_EST", "FB_AMD"):
            if h.mechanism_code in (
                "EST_HYD", "AMD_HYD", "CAR_HYD", "DEHAL", "N_OX", "S_OX", "N_DEALK",
            ):
                continue
        out.append(h)
    return out
