"""
PRIME-DS — Name Reaction Process Simulator (rdchiral + RDKit)

유명 Name Reaction SMARTS 템플릿으로 목적물·부반응 유연물질을 생성하고
내부 10K Morgan-Tanimoto RAG 매트릭스와 통합한다.
"""

# 경로 B
# Name Reaction 33종(ex) Suzuki, Grignard, Wittig)
# Suzuki 커플링
#   ├── main_smirks   의도한 생성물
#   └── side_smirks   부산물
# 의도한 것과 원치 않은 것을 구분해서 예측하는 게 핵심. 경로 A는 이 구분이 없음.

from __future__ import annotations

import json
from itertools import product
from pathlib import Path
from typing import Any, Callable

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRun, rdchiralRunText
# rdchiral: RDKit만 쓰면 반응 적용 시 입체 화학이 뒤집히는 오류가 생기는데, 이 라이브러리가 그걸 막아줌. 광학이성질체가 중요한 의약품에서는 필요한 기능

from engines.p03_branding import log_process_warning
from engines.p33_process_impurity_engine import (
    DB_TAG_MECHANISTIC,
    DB_TAG_MODEL,
    DB_TAG_VERIFIED,
    classify_ich_m7_genotoxicity,
    compute_db_verification_score,
    get_verification_fingerprint_db,
    mol_from_smiles,
    morgan_fp,
)
from engines.p05_safeguards import ProcessRunContext, log_invalid_smiles_bypass

# ── UI 선택지 · 반응 카탈로그 ─────────────────────────────────────────────────
NAME_RXN_NONE: str = "선택안함 (No Name Reaction)"
NAME_RXN_STANDARD_DEGRADATION: str = "표준 가수분해/산화"  # legacy · NAME_RXN_NONE 과 동일 skip

_CATALOG_DIR = Path(__file__).resolve().parent.parent / "data" / "reaction_catalogs"

# 베이스라인 Name Reaction / SMIRKS 템플릿 (USPTO·ORD JSON 병합 전 내장 카탈로그)
BASELINE_REACTION_CATALOG: dict[str, dict[str, Any]] = {
    "Suzuki-Miyaura Coupling": {
        "main_smarts": "[c:1]-[Cl,Br,I:2].[c:3]-B(O)O>>[c:1]-[c:3]",
        "byproduct_smarts": [
            ("Halide homo-coupling", "[c:1]-[Cl,Br,I:2].[c:3]-[Cl,Br,I:4]>>[c:1]-[c:3]"),
            ("Boronic homo-coupling", "[c:1]-B(O)O.[c:2]-B(O)O>>[c:1]-[c:2]"),
        ],
        "required_a": "[c]-[Cl,Br,I]",
        "required_b": "[c]-B(O)O",
        "hint_a": "Aryl/Vinyl halide",
        "hint_b": "Boronic acid",
        "source": "baseline",
    },
    "Friedel-Crafts Acylation": {
        "main_smarts": "[c:1]:[c:2]:[c:3].[C:4](=[O:5])Cl>>[c:1][C:4](=[O:5])[c:2]:[c:3]",
        "byproduct_smarts": [
            ("Di-acylation", "[c:1][C:2](=[O:3])[c:4]:[c:5]:[c:6].[C:7](=[O:8])Cl>>[c:1][C:2](=[O:3])[c:4]([C:7](=[O:8])[c:5]:[c:6])"),
        ],
        "required_a": "[c;R]",
        "required_b": "[C](=[O])Cl",
        "hint_a": "Aromatic ring",
        "hint_b": "Acyl chloride",
        "source": "baseline",
    },
    "Aldol Condensation": {
        "main_smarts": "[C:1]=[O:2].[C:3]=[O:4]>>[C:1]([C:3]=[O:4])[O:2]",
        "byproduct_smarts": [
            ("Self-condensation", "[CH:1]=[O:2].[CH:3]=[O:4]>>[CH:1]([CH:3]=[O:4])[O:2]"),
        ],
        "required_a": "[C]=[O]",
        "required_b": "[C]=[O]",
        "hint_a": "Carbonyl A",
        "hint_b": "Carbonyl B",
        "source": "baseline",
    },
    "Grignard Reaction": {
        "main_smarts": "[C:1]=[O:2].[C:3][Mg+0]>>[C:1]([O:2])[C:3]",
        "byproduct_smarts": [
            ("Over-addition", "[C:1]([O:2])[C:3].[C:4][Mg+0]>>[C:1]([O:2])([C:3])[C:4]"),
        ],
        "required_a": "[C]=[O]",
        "required_b": "[C][Mg+0]",
        "hint_a": "Carbonyl electrophile",
        "hint_b": "Grignard reagent",
        "source": "baseline",
    },
    "Amide Coupling": {
        "main_smarts": "[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4]",
        "byproduct_smarts": [
            ("Acid anhydride", "[C:1](=[O:2])[OH:3].[C:3](=[O:4])[OH:5]>>[C:1](=[O:2])OC(=O)[C:3]"),
        ],
        "required_a": "[C](=[O])O",
        "required_b": "[N;!$(N=*)]",
        "hint_a": "Carboxylic acid / acyl donor",
        "hint_b": "Amine nucleophile",
        "source": "baseline",
    },
    "Esterification (Fischer)": {
        "main_smarts": "[C:1](=[O:2])[OH:3].[OH:4]>>[C:1](=[O:2])O[O:4]",
        "byproduct_smarts": [
            ("Acid-catalyzed hydrolysis", "[C:1](=[O:2])O[O:3]>>[C:1](=[O:2])[OH].[OH:3]"),
        ],
        "required_a": "[C](=[O])O",
        "required_b": "[OH]",
        "hint_a": "Carboxylic acid",
        "hint_b": "Alcohol",
        "source": "baseline",
    },
    "Williamson Ether Synthesis": {
        "main_smarts": "[C:1][O-1:2].[C:3][Br,I,Cl:4]>>[C:1][O:2][C:3]",
        "byproduct_smarts": [],
        "required_a": "[C][O-]",
        "required_b": "[C][Cl,Br,I]",
        "hint_a": "Alkoxide",
        "hint_b": "Alkyl halide",
        "source": "baseline",
    },
    "SN2 N-Alkylation": {
        "main_smarts": "[N:1].[C:2][Cl,Br,I:3]>>[N:1][C:2]",
        "byproduct_smarts": [
            ("Over-alkylation", "[N:1][C:2].[C:3][Cl,Br,I:4]>>[N:1]([C:2])[C:3]"),
        ],
        "required_a": "[N;!$(N=*)]",
        "required_b": "[C][Cl,Br,I]",
        "hint_a": "Amine",
        "hint_b": "Alkyl halide",
        "source": "baseline",
    },
    "Mitsunobu Reaction": {
        "main_smarts": "[O:1][H].[C:2][OH:3]>>[O:1][C:2]",
        "byproduct_smarts": [],
        "required_a": "[O][H]",
        "required_b": "[C][OH]",
        "hint_a": "Alcohol (nucleophile)",
        "hint_b": "Alcohol (electrophile)",
        "source": "baseline",
    },
    "Heck Coupling": {
        "main_smarts": "[c:1]-[Cl,Br,I:2].[CH2:3]=[CH2:4]>>[c:1]-[CH2:3]-[CH2:4]",
        "byproduct_smarts": [
            ("Di-aryl coupling", "[c:1]-[Cl,Br,I:2].[c:3]-[Cl,Br,I:4]>>[c:1]-[c:3]"),
        ],
        "required_a": "[c]-[Cl,Br,I]",
        "required_b": "[CH2]=[CH2]",
        "hint_a": "Aryl halide",
        "hint_b": "Olefin",
        "source": "baseline",
    },
    "Sonogashira Coupling": {
        "main_smarts": "[c:1]-[Cl,Br,I:2].[C:3]#[CH:4]>>[c:1]-[C:3]#[CH:4]",
        "byproduct_smarts": [],
        "required_a": "[c]-[Cl,Br,I]",
        "required_b": "[C]#[CH]",
        "hint_a": "Aryl halide",
        "hint_b": "Terminal alkyne",
        "source": "baseline",
    },
    "Stille Coupling": {
        "main_smarts": "[c:1]-[Cl,Br,I:2].[c:3]-[Sn]>>[c:1]-[c:3]",
        "byproduct_smarts": [],
        "required_a": "[c]-[Cl,Br,I]",
        "required_b": "[c]-[Sn]",
        "hint_a": "Aryl halide",
        "hint_b": "Stannane",
        "source": "baseline",
    },
    "Buchwald-Hartwig Amination": {
        "main_smarts": "[c:1]-[Cl,Br,I:2].[N:3]>>[c:1]-[N:3]",
        "byproduct_smarts": [
            ("Bis-arylation", "[c:1]-[N:2].[c:3]-[Cl,Br,I:4]>>[c:1]-[N:2]-[c:3]"),
        ],
        "required_a": "[c]-[Cl,Br,I]",
        "required_b": "[N;!$(N=*)]",
        "hint_a": "Aryl halide",
        "hint_b": "Amine",
        "source": "baseline",
    },
    "Reductive Amination": {
        "main_smarts": "[C:1]=[O:2].[N:3]>>[C:1][N:3]",
        "byproduct_smarts": [
            ("Imine intermediate", "[C:1]=[O:2].[N:3]>>[C:1]=[N:3]"),
        ],
        "required_a": "[C]=[O]",
        "required_b": "[N;!$(N=*)]",
        "hint_a": "Aldehyde/Ketone",
        "hint_b": "Amine",
        "source": "baseline",
    },
    "NaBH4 Carbonyl Reduction": {
        "main_smarts": "[C:1]=[O:2]>>[C:1][OH:2]",
        "byproduct_smarts": [],
        "required_a": "[C]=[O]",
        "required_b": "[C]=[O]",
        "hint_a": "Carbonyl",
        "hint_b": "Hydride donor (implicit)",
        "source": "baseline",
    },
    "Swern Oxidation": {
        "main_smarts": "[C:1][CH2:2][OH:3]>>[C:1][CH2:2]=[O:3]",
        "byproduct_smarts": [],
        "required_a": "[CH2][OH]",
        "required_b": "[CH2][OH]",
        "hint_a": "Primary alcohol",
        "hint_b": "Oxidant (implicit)",
        "source": "baseline",
    },
    "Wittig Reaction": {
        "main_smarts": "[C:1]=[O:2].[CH2:3]=[P+0]>>[C:1]=[CH2:3]",
        "byproduct_smarts": [],
        "required_a": "[C]=[O]",
        "required_b": "[CH2]=[P]",
        "hint_a": "Carbonyl",
        "hint_b": "Ylide",
        "source": "baseline",
    },
    "Michael Addition": {
        "main_smarts": "[C:1]=[C:2]-[C:3]=[O:4].[N:5]>>[C:1]([N:5])-[C:2]-[C:3]=[O:4]",
        "byproduct_smarts": [],
        "required_a": "[C]=[C]-[C]=[O]",
        "required_b": "[N;!$(N=*)]",
        "hint_a": "Michael acceptor",
        "hint_b": "Nucleophile",
        "source": "baseline",
    },
    "Diels-Alder Reaction": {
        "main_smarts": "[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2][C:3][C:4][C:5][C:6]1",
        "byproduct_smarts": [],
        "required_a": "[C]=[C]-[C]=[C]",
        "required_b": "[C]=[C]",
        "hint_a": "Diene",
        "hint_b": "Dienophile",
        "source": "baseline",
    },
    "Boc Protection": {
        "main_smarts": "[N:1].[C:2](=[O:3])OC(=O)C(C)(C)C>>[N:1]C(=O)OC(C)(C)C",
        "byproduct_smarts": [],
        "required_a": "[N;!$(N=*)]",
        "required_b": "[C](=[O])OC(=O)C(C)(C)C",
        "hint_a": "Amine",
        "hint_b": "Boc anhydride",
        "source": "baseline",
    },
    "Fmoc Deprotection": {
        "main_smarts": "[N:1]C(=O)OCc1ccccc1>>[N:1]",
        "byproduct_smarts": [],
        "required_a": "[N]C(=O)OCc1ccccc1",
        "required_b": "[N]C(=O)OCc1ccccc1",
        "hint_a": "Fmoc-amine",
        "hint_b": "Base (implicit)",
        "source": "baseline",
    },
    "Aromatic Nitration": {
        "main_smarts": "[c:1]:[c:2]:[c:3]>>[c:1]([N+](=O)[O-])[c:2]:[c:3]",
        "byproduct_smarts": [
            ("Dinitration", "[c:1]([N+](=O)[O-])[c:2]:[c:3]>>[c:1]([N+](=O)[O-])[c:2]([N+](=O)[O-])[c:3]"),
        ],
        "required_a": "[c;R]",
        "required_b": "[c;R]",
        "hint_a": "Aromatic ring",
        "hint_b": "Nitronium (implicit)",
        "source": "baseline",
    },
    "Nitro Reduction to Amine": {
        "main_smarts": "[c:1][N+](=O)[O-:2]>>[c:1][NH2:2]",
        "byproduct_smarts": [
            ("Hydroxylamine intermediate", "[c:1][N+](=O)[O-]>>[c:1][NH]O"),
        ],
        "required_a": "[c][N+](=O)[O-]",
        "required_b": "[c][N+](=O)[O-]",
        "hint_a": "Nitroarene",
        "hint_b": "Reductant (implicit)",
        "source": "baseline",
    },
    "Halogenation (Aromatic)": {
        "main_smarts": "[c:1]:[c:2]:[c:3]>>[c:1]([Cl,Br,I:4])[c:2]:[c:3]",
        "byproduct_smarts": [],
        "required_a": "[c;R]",
        "required_b": "[c;R]",
        "hint_a": "Aromatic ring",
        "hint_b": "Halogen source",
        "source": "baseline",
    },
    "Hydrogenation (Alkene)": {
        "main_smarts": "[C:1]=[C:2]>>[C:1][C:2]",
        "byproduct_smarts": [],
        "required_a": "[C]=[C]",
        "required_b": "[C]=[C]",
        "hint_a": "Olefin",
        "hint_b": "H2 / catalyst (implicit)",
        "source": "baseline",
    },
    "Dehydration (Alcohol)": {
        "main_smarts": "[C:1][CH2:2][OH:3]>>[C:1]=[CH2:2]",
        "byproduct_smarts": [],
        "required_a": "[CH2][OH]",
        "required_b": "[CH2][OH]",
        "hint_a": "Alcohol",
        "hint_b": "Acid catalyst (implicit)",
        "source": "baseline",
    },
    "Claisen Rearrangement": {
        "main_smarts": "[C:1]=[C:2]-[O:3][C:4]=[C:5]>>[C:1]([C:4]=[O:3])-[C:2]-[C:5]",
        "byproduct_smarts": [],
        "required_a": "[C]=[C]-O[C]=[C]",
        "required_b": "[C]=[C]-O[C]=[C]",
        "hint_a": "Allyl vinyl ether",
        "hint_b": "Heat (implicit)",
        "source": "baseline",
    },
    "Cross Metathesis": {
        "main_smarts": "[CH2:1]=[CH:2][R:3].[CH2:4]=[CH:5][R:6]>>[CH2:1]=[CH:5][R:6]",
        "byproduct_smarts": [],
        "required_a": "[CH2]=[CH]",
        "required_b": "[CH2]=[CH]",
        "hint_a": "Olefin A",
        "hint_b": "Olefin B",
        "source": "baseline",
    },
    "Hydrolysis (Ester, Acidic)": {
        "main_smarts": "[C:1](=[O:2])O[C:3]>>[C:1](=[O:2])[OH].[OH:3]",
        "byproduct_smarts": [],
        "required_a": "[C](=[O])O[C]",
        "required_b": "[C](=[O])O[C]",
        "hint_a": "Ester",
        "hint_b": "Water / acid",
        "source": "baseline",
    },
    "Hydrolysis (Amide, Basic)": {
        "main_smarts": "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[O-].[N:3]",
        "byproduct_smarts": [],
        "required_a": "[C](=[O])[N]",
        "required_b": "[C](=[O])[N]",
        "hint_a": "Amide",
        "hint_b": "Base / water",
        "source": "baseline",
    },
    "Sandmeyer Reaction": {
        "main_smarts": "[c:1][N+]#N.[Cl-]>>[c:1][Cl]",
        "byproduct_smarts": [],
        "required_a": "[c][N+]#N",
        "required_b": "[Cl-]",
        "hint_a": "Diazonium salt",
        "hint_b": "Halide nucleophile",
        "source": "baseline",
    },
    "Knoevenagel Condensation": {
        "main_smarts": "[C:1]=[O:2].[C:3](=[O:4])C>>[C:1]=[C:3]([C:4]=[O])",
        "byproduct_smarts": [],
        "required_a": "[C]=[O]",
        "required_b": "[C](=[O])C",
        "hint_a": "Aldehyde/Ketone",
        "hint_b": "Active methylene",
        "source": "baseline",
    },
    "Esterification / Amide Bond Formation": {
        "main_smarts": "[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4]",
        "byproduct_smarts": [
            ("Esterification side", "[C:1](=[O:2])[OH:3].[OH:4]>>[C:1](=[O:2])O[O:4]"),
            ("Acid anhydride byproduct", "[C:1](=[O:2])[OH:3].[C:3](=[O:4])[OH:5]>>[C:1](=[O:2])OC(=O)[C:3]"),
        ],
        "required_a": "[C](=[O])O",
        "required_b": "[N;!$(N=*)]",
        "hint_a": "Carboxylic acid / acyl donor",
        "hint_b": "Amine nucleophile",
        "source": "baseline",
    },
}

_MERGED_CATALOG_CACHE: dict[str, dict[str, Any]] | None = None


def _normalize_reaction_entry(raw: dict[str, Any], *, source: str) -> dict[str, Any] | None:
    name = str(raw.get("name") or raw.get("reaction_name") or "").strip()
    main_smarts = str(raw.get("main_smarts") or raw.get("smirks") or raw.get("smarts") or "").strip()
    if not name or not main_smarts or ">>" not in main_smarts:
        return None
    byproducts = raw.get("byproduct_smarts") or raw.get("byproducts") or []
    normalized_bp: list[tuple[str, str]] = []
    for item in byproducts:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            normalized_bp.append((str(item[0]), str(item[1])))
        elif isinstance(item, dict):
            label = str(item.get("label") or item.get("name") or "Byproduct")
            smarts = str(item.get("smarts") or item.get("smirks") or "")
            if smarts:
                normalized_bp.append((label, smarts))
    return {
        "main_smarts": main_smarts,
        "byproduct_smarts": normalized_bp,
        "required_a": str(raw.get("required_a") or raw.get("reactant_a_smarts") or ""),
        "required_b": str(raw.get("required_b") or raw.get("reactant_b_smarts") or ""),
        "hint_a": str(raw.get("hint_a") or raw.get("reactant_a_hint") or "Reactant A"),
        "hint_b": str(raw.get("hint_b") or raw.get("reactant_b_hint") or "Reactant B"),
        "source": source,
    }


def load_reaction_catalog(path: str | Path) -> dict[str, dict[str, Any]]:
    """
    외부 JSON 반응 카탈로그 로드 (USPTO / ORD 템플릿 확장용).

    지원 형식
    ---------
    - {"reactions": [{name, main_smarts, ...}, ...]}
    - [{name, main_smarts, ...}, ...]
    - { "Reaction Name": {main_smarts: ...}, ... }
    """
    catalog: dict[str, dict[str, Any]] = {}
    file_path = Path(path)
    if not file_path.is_file():
        log_process_warning(f"Reaction catalog not found: {file_path}")
        return catalog
    try:
        with file_path.open(encoding="utf-8") as fh:
            payload = json.load(fh)
    except Exception as exc:
        log_process_warning(f"Reaction catalog JSON parse failed ({file_path}): {exc}")
        return catalog

    if isinstance(payload, dict):
        if "reactions" in payload and isinstance(payload["reactions"], list):
            items = payload["reactions"]
        else:
            items = []
            for key, value in payload.items():
                if isinstance(value, dict):
                    entry = dict(value)
                    entry.setdefault("name", key)
                    items.append(entry)
    elif isinstance(payload, list):
        items = payload
    else:
        return catalog

    for raw in items:
        if not isinstance(raw, dict):
            continue
        name = str(raw.get("name") or raw.get("reaction_name") or "").strip()
        normalized = _normalize_reaction_entry(raw, source=file_path.name)
        if normalized and name:
            catalog[name] = normalized
    return catalog


def load_reaction_catalogs_from_directory(directory: str | Path | None = None) -> dict[str, dict[str, Any]]:
    """디렉터리 내 모든 .json 카탈로그 병합 (후입 우선)."""
    merged: dict[str, dict[str, Any]] = {}
    dir_path = Path(directory) if directory else _CATALOG_DIR
    if not dir_path.is_dir():
        return merged
    for json_file in sorted(dir_path.glob("*.json")):
        merged.update(load_reaction_catalog(json_file))
    return merged


def get_reaction_catalog(*, reload: bool = False) -> dict[str, dict[str, Any]]:
    """베이스라인 + 외부 JSON 병합 카탈로그 (캐시)."""
    global _MERGED_CATALOG_CACHE
    if _MERGED_CATALOG_CACHE is not None and not reload:
        return _MERGED_CATALOG_CACHE
    merged = dict(BASELINE_REACTION_CATALOG)
    merged.update(load_reaction_catalogs_from_directory())
    _MERGED_CATALOG_CACHE = merged
    return merged


def get_reaction_ui_options() -> list[str]:
    """Name Reaction selectbox 옵션."""
    names = sorted(get_reaction_catalog().keys())
    return [NAME_RXN_NONE] + names


def get_reaction_spec(chosen_reaction: str) -> dict[str, Any] | None:
    return get_reaction_catalog().get(chosen_reaction)


NAME_REACTION_UI_OPTIONS: list[str] = get_reaction_ui_options()
NAME_REACTION_KB: dict[str, dict[str, Any]] = get_reaction_catalog()  # 레거시 별칭

_NAME_RXN_SKIP_VALUES: frozenset[str] = frozenset(
    {NAME_RXN_NONE, NAME_RXN_STANDARD_DEGRADATION, "선택안함", "선택 안함"}
)


def is_name_reaction_skipped(chosen_reaction: str | None) -> bool:
    """Name Reaction rdchiral 시뮬레이션 생략 여부."""
    return (chosen_reaction or "").strip() in _NAME_RXN_SKIP_VALUES

PATHWAY_PREFIX = "Name Reaction"


def rdchiralRunReactants(reaction_smarts: str, reactant_smiles: str | list[str]) -> list[str]:
    """
    rdchiral 1.1.0 호환 `rdchiralRunReactants` 실행 별칭.

    공식 패키지에는 `rdchiral.main.rdchiralRunReactants` 심볼이 없으며,
    2분자 Name Reaction은 `rdchiralReaction` + stereochemical `RunReactants`로 처리한다.
    1분자 분해/부반응은 `rdchiralRunText`를 우선 사용한다.
    """
    return rdchiral_run_reactants(reaction_smarts, reactant_smiles)


def rdchiral_run_reactants(reaction_smarts: str, reactant_smiles: str | list[str]) -> list[str]:
    """
    `rdchiralRunReactants` 호환 실행 래퍼.

    - 1분자 반응: rdchiralRunText
    - 2분자 이상: rdchiralReaction + rdchiralReactants → rdchiralRun
    """
    if isinstance(reactant_smiles, str):
        smiles_list = [s.strip() for s in reactant_smiles.split(".") if s.strip()]
    else:
        smiles_list = [str(s).strip() for s in reactant_smiles if str(s).strip()]

    if not smiles_list or ">>" not in reaction_smarts:
        return []

    reactant_templates = reaction_smarts.split(">>")[0].strip().split(".")
    n_templates = len(reactant_templates)

    if len(smiles_list) < n_templates:
        return []

    try:
        if n_templates == 1:
            return list(rdchiralRunText(reaction_smarts, smiles_list[0]) or [])
        reactant_dot = ".".join(smiles_list[:n_templates])
        rxn = rdchiralReaction(reaction_smarts)
        reactants = rdchiralReactants(reactant_dot)
        return list(rdchiralRun(rxn, reactants) or [])
    except Exception as exc:
        log_process_warning(f"rdchiralRun fallback → RDKit: {exc}")
        log_invalid_smiles_bypass(context="rdchiral_run_reactants")
        return _rdkit_run_reactants(reaction_smarts, smiles_list[:n_templates])


def _rdkit_run_reactants(reaction_smarts: str, smiles_list: list[str]) -> list[str]:
    """rdchiral 실패 시 RDKit AllChem fallback."""
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    except Exception:
        return []
    if rxn is None:
        return []
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    mols = [m for m in mols if m is not None]
    if not mols:
        return []
    try:
        outcomes = rxn.RunReactants(tuple(mols))
        return _outcomes_to_smiles(outcomes)
    except Exception:
        return []


def _outcomes_to_smiles(outcomes: tuple | list | None) -> list[str]:
    """RunReactants outcome → canonical SMILES 리스트."""
    if not outcomes:
        return []
    found: set[str] = set()
    for outcome in outcomes:
        if outcome is None:
            continue
        items = outcome if isinstance(outcome, (tuple, list)) else (outcome,)
        for prod in items:
            if prod is None:
                continue
            try:
                Chem.SanitizeMol(prod)
                found.add(Chem.MolToSmiles(prod, isomericSmiles=True))
            except Exception:
                log_invalid_smiles_bypass(context="rdchiral_outcome")
                continue
    return list(found)


def _has_smarts(smiles: str, smarts: str) -> bool:
    try:
        mol = mol_from_smiles(smiles, log_bypass=False)
        pat = Chem.MolFromSmarts(smarts)
        if mol is None or pat is None:
            return False
        return mol.HasSubstructMatch(pat)
    except Exception:
        log_invalid_smiles_bypass(context="has_smarts")
        return False


def _canonical(smiles: str) -> str:
    try:
        mol = mol_from_smiles(smiles, log_bypass=False)
        if mol is None:
            return smiles.strip()
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        log_invalid_smiles_bypass(context="name_reaction_canonical")
        return smiles.strip()


def _split_product_frags(combined_smiles: str) -> list[str]:
    """다분자 outcome SMILES (dot-separated) → 개별 fragment."""
    parts = [p.strip() for p in combined_smiles.split(".") if p.strip()]
    frags: list[str] = []
    for part in parts:
        try:
            canon = _canonical(part)
            if mol_from_smiles(canon, log_bypass=False) is not None:
                frags.append(canon)
        except Exception:
            log_invalid_smiles_bypass(context="split_product_frags")
            continue
    return frags


def _tanimoto_to_target(smiles: str, target_smiles: str) -> float:
    try:
        m1 = mol_from_smiles(smiles, log_bypass=False)
        m2 = mol_from_smiles(target_smiles, log_bypass=False)
        if m1 is None or m2 is None:
            return 0.0
        fp1 = morgan_fp(m1)
        fp2 = morgan_fp(m2)
        if fp1 is None or fp2 is None:
            return 0.0
        return float(DataStructs.TanimotoSimilarity(fp1, fp2))
    except Exception:
        log_invalid_smiles_bypass(context="tanimoto_to_target")
        return 0.0


def _classify_product_role(
    smiles: str,
    *,
    target_smiles: str | None,
    reaction_name: str,
    tag: str,
) -> str:
    """목적물 vs 부반응(유연물질) 역할 분류."""
    if tag != "main":
        return f"{PATHWAY_PREFIX}: {reaction_name} · Side-product / Byproduct (공정 유연물질)"
    if target_smiles and _tanimoto_to_target(smiles, target_smiles) >= 0.55:
        return f"{PATHWAY_PREFIX}: {reaction_name} · Main Product (목적물)"
    return f"{PATHWAY_PREFIX}: {reaction_name} · Main Product (목적물)"


def _boost_name_reaction_confidence(
    match_pct: float,
    reliability_pct: float,
    status_tag: str,
    *,
    is_byproduct: bool,
) -> tuple[float, str, str]:
    """
    Name Reaction 부산물이 10K DB와 고유사할 경우 Confidence 극대화 + 경보 태그.
    """
    display = f"{reliability_pct:.1f}% — {status_tag}"
    if not is_byproduct:
        return match_pct, status_tag, display

    if match_pct >= 90.0:
        boosted = min(100.0, reliability_pct + 12.0)
        tag = f"🚨 VERIFIED IMPURITY ALERT — {DB_TAG_VERIFIED}"
        display = f"{boosted:.1f}% — {tag}"
        return match_pct, tag, display
    if match_pct >= 75.0:
        boosted = min(100.0, reliability_pct + 8.0)
        tag = f"⚠️ HIGH-RISK BYPRODUCT — {DB_TAG_MECHANISTIC}"
        display = f"{boosted:.1f}% — {tag}"
        return match_pct, tag, display
    if match_pct >= 50.0 and is_byproduct:
        tag = f"⚠️ Name Rxn Byproduct — {DB_TAG_MECHANISTIC}"
        display = f"{reliability_pct:.1f}% — {tag}"
    return match_pct, status_tag, display


def _impurity_record_merge_key(rec: dict[str, Any]) -> tuple[str, ...]:
    """SMILES + 메커니즘 규칙 기준 병합 키 — 동일 SMILES라도 Suzuki 규칙별 유지."""
    smi = str(rec.get("smiles") or "")
    mechanism = str(
        rec.get("suzuki_rule")
        or rec.get("reaction_name")
        or rec.get("name_reaction")
        or ""
    ).strip()
    if mechanism:
        return (smi, mechanism)
    return (smi,)


def simulate_name_reaction_impurities(
    raw_a_smiles: str,
    raw_b_smiles: str,
    chosen_reaction: str,
    *,
    target_smiles: str | None = None,
    temperature_c: float = 25.0,
    db_fps: list | None = None,
    db_smiles: list[str] | None = None,
    reagent_context: dict[str, Any] | None = None,
    ctx: ProcessRunContext | None = None,
) -> list[dict[str, Any]]:
    """
    rdchiral 기반 Name Reaction 시뮬레이션.

    Parameters
    ----------
    raw_a_smiles, raw_b_smiles :
        Ketcher 주원료 A/B SMILES.
    chosen_reaction :
        NAME_REACTION_UI_OPTIONS 중 하나.
    target_smiles :
        목적물(API) SMILES — Main product 분류·유사도 가중.
    """
    if is_name_reaction_skipped(chosen_reaction):
        return []

    try:
        spec = get_reaction_spec(chosen_reaction)
        if spec is None:
            return []

        smi_a = (raw_a_smiles or "").strip()
        smi_b = (raw_b_smiles or "").strip()
        if not smi_a or not smi_b:
            return []

        req_a = str(spec.get("required_a") or "").strip()
        req_b = str(spec.get("required_b") or "").strip()
        if req_a and not _has_smarts(smi_a, req_a):
            return []
        if req_b and not _has_smarts(smi_b, req_b):
            return []

        if db_fps is None or db_smiles is None:
            db_fps, db_smiles = get_verification_fingerprint_db()

        reactant_dot = f"{smi_a}.{smi_b}"
        seen: set[str] = set()
        records: list[dict[str, Any]] = []

        def _process_outcome(combined: str, *, tag: str, sub_label: str) -> None:
            if ctx is not None and not ctx.check_continue():
                return
            try:
                frags = _split_product_frags(combined)
                if not frags:
                    frags = [_canonical(combined)]
                for smi in frags:
                    if ctx is not None and not ctx.check_continue():
                        return
                    if smi in seen:
                        continue
                    mol = mol_from_smiles(smi, log_bypass=False)
                    if mol is None:
                        log_invalid_smiles_bypass(context="name_reaction_product")
                        continue
                    seen.add(smi)
                    is_byproduct = tag != "main"
                    pathway = _classify_product_role(
                        smi,
                        target_smiles=target_smiles,
                        reaction_name=chosen_reaction,
                        tag=tag,
                    )
                    if sub_label and is_byproduct:
                        pathway = f"{pathway} · {sub_label}"

                    exact_mw = round(Descriptors.ExactMolWt(mol), 2)
                    ich_class, ich_detail = classify_ich_m7_genotoxicity(mol)
                    match_pct, status_tag, reliability_pct = compute_db_verification_score(
                        mol, db_fps, db_smiles, temperature_c=temperature_c
                    )
                    _, status_tag, confidence_display = _boost_name_reaction_confidence(
                        match_pct,
                        reliability_pct,
                        status_tag,
                        is_byproduct=is_byproduct,
                    )

                    records.append(
                        {
                            "pathway": pathway,
                            "smiles": smi,
                            "exact_mw": exact_mw,
                            "ich_m7": ich_class,
                            "ich_m7_detail": ich_detail,
                            "db_match_pct": match_pct,
                            "db_status_tag": status_tag,
                            "reliability_pct": reliability_pct,
                            "confidence_display": confidence_display,
                            "source_reactant": f"Name Rxn · {spec['hint_a']} + {spec['hint_b']}",
                            "name_reaction": chosen_reaction,
                            "is_name_reaction_byproduct": is_byproduct,
                        }
                    )
                    if ctx:
                        ctx.register_product()
            except Exception:
                log_invalid_smiles_bypass(context="name_reaction_process_outcome")
                return

        for outcome in rdchiralRunReactants(spec["main_smarts"], reactant_dot):
            if ctx is not None and not ctx.check_continue():
                break
            _process_outcome(outcome, tag="main", sub_label="")

        for sub_label, bp_smarts in spec.get("byproduct_smarts", []):
            if ctx is not None and not ctx.check_continue():
                break
            for outcome in rdchiralRunReactants(bp_smarts, reactant_dot):
                if ctx is not None and not ctx.check_continue():
                    break
                _process_outcome(outcome, tag="byproduct", sub_label=sub_label)

        records.sort(
            key=lambda r: (r.get("is_name_reaction_byproduct", False), -r["db_match_pct"])
        )

        return ctx.truncate_records(records) if ctx else records
    except Exception as exc:
        log_process_warning(f"simulate_name_reaction_impurities aborted: {exc}")
        return []


def simulate_suzuki_competitive_screening(
    raw_a_smiles: str,
    raw_b_smiles: str,
    *,
    chosen_reaction: str | None = None,
    temperature_c: float = 25.0,
    target_smiles: str | None = None,
    reagent_context: dict[str, Any] | None = None,
    ctx: ProcessRunContext | None = None,
) -> list[dict[str, Any]]:
    """
    Suzuki-Miyaura 메커니즘 기반 경쟁 부반응 5종 스크리닝.

    Main · Halide/Boronate homo-coupling · Des-bromo · Protodeboronation
    """
    from engines.p36_suzuki_impurity_rules import (
        annotate_impurity_records,
        get_suzuki_predictor,
        is_suzuki_coupling_context,
        resolve_suzuki_reactants,
    )

    if not is_suzuki_coupling_context(chosen_reaction, reagent_context):
        return []
    smi_a = (raw_a_smiles or "").strip()
    smi_b = (raw_b_smiles or "").strip()
    if not smi_a or not smi_b:
        return []
    if ctx is not None and not ctx.check_continue():
        return []

    pair = resolve_suzuki_reactants(smi_a, smi_b)
    if not pair:
        return []

    label = f"Suzuki competitive · {chosen_reaction or 'Pd catalysis'}"
    halide, boronic = pair
    bundle = get_suzuki_predictor().predict_all_scenarios_marshaled(halide, boronic)
    preds = bundle.get("predicted_impurities") or []
    records = get_suzuki_predictor().predictions_to_impurity_records(
        preds,
        source_label=label,
        temperature_c=temperature_c,
    )
    if target_smiles:
        target_canon = _canonical(target_smiles)
        for rec in records:
            if rec.get("smiles") == target_canon and rec.get("impurity_role") == "Main":
                rec["impurity_classification"] = "Main Product (Target)"
                rec["matches_target_api"] = True
    return annotate_impurity_records(records)


def run_universal_side_reaction_prediction(
    entries: list[tuple[str, str, dict[str, Any]]],
    *,
    target_smiles: str | None = None,
    selected_environment: dict[str, Any] | None = None,
    temperature_c: float = 25.0,
    name_reaction: str | None = None,
    ctx: ProcessRunContext | None = None,
) -> list[dict[str, Any]]:
    """
    Universal Side-Reaction Predictive Engine — Process Mode 메인 라우터.

  출발물질 배열 + 목적물 + 4단 Registry [Solvents, Acids, Bases, Catalysts] → 부반응 예측.
    """
    from engines.p35_universal_reaction_impurity_engine import (
        extract_reactant_smiles_from_entries,
        run_universal_side_reaction_pipeline,
    )

    if ctx is not None and not ctx.check_continue():
        return []

    reactant_pool = extract_reactant_smiles_from_entries(
        entries,
        exclude_target=target_smiles,
    )
    if not reactant_pool:
        return []

    bundle = run_universal_side_reaction_pipeline(
        reactant_pool,
        target_smiles=target_smiles or "",
        selected_environment=selected_environment,
        temperature_c=temperature_c,
        source_label="Universal Side-Reaction · full reactant pool",
        name_reaction=name_reaction,
        ctx=ctx,
    )
    records = bundle.get("records") or []
    if records:
        meta = {
            "ai_bridge_status": bundle.get("ai_bridge_status"),
            "ai_bridge_backend": bundle.get("ai_bridge_backend"),
            "ai_scenario_count": bundle.get("ai_scenario_count", 0),
        }
        records[0] = {**records[0], "_pipeline_meta": meta}
    return records


def merge_name_reaction_into_batch(
    base_records: list[dict[str, Any]],
    name_records: list[dict[str, Any]],
    *,
    max_products: int | None = None,
) -> list[dict[str, Any]]:
    """Step 3 batch 결과 + Name Reaction 레코드 병합 (메커니즘 규칙별 키)."""
    cap = max_products or len(base_records) + len(name_records)
    merged = list(base_records)
    seen = {_impurity_record_merge_key(r) for r in merged if r.get("smiles")}
    for rec in name_records:
        if len(merged) >= cap:
            break
        smi = rec.get("smiles")
        if not smi:
            continue
        key = _impurity_record_merge_key(rec)
        if key in seen:
            continue
        seen.add(key)
        merged.append(rec)
    return merged


# ── Phase 6: 공정 컨텍스트 기반 반응 카탈로그 엔진 ─────────────────────────────
REGISTRY_ENVIRONMENT_CATEGORIES: tuple[str, ...] = ("Solvents", "Acids", "Bases", "Catalysts")

# USPTO/ORD 스타일 내장 매트릭스 — 외부 JSON 없이도 가동 (required_environment 필수)
DEFAULT_PROCESS_REACTION_CATALOG: dict[str, dict[str, Any]] = {
    "Amide_Coupling": {
        "smirks": "[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4]",
        "required_environment": ["Bases"],
        "description": "아미드 커플링 — 염기성 조건 부반응",
    },
    "Esterification": {
        "smirks": "[C:1](=[O:2])[OH:3].[OH:4]>>[C:1](=[O:2])O[O:4]",
        "required_environment": ["Acids"],
        "description": "Fischer 에스테르화 — 산 촉매",
    },
    "Suzuki_Coupling": {
        "smirks": "[c:1]-[Cl,Br,I:2].[c:3]-B(O)O>>[c:1]-[c:3]",
        "required_environment": ["Catalysts"],
        "required_catalyst_metal": "Pd",
        "description": "스즈키 교차결합 — Pd 촉매 필수",
    },
    "Friedel_Crafts_Acylation": {
        "smirks": "[c:1]:[c:2]:[c:3].[C:4](=[O:5])Cl>>[c:1][C:4](=[O:5])[c:2]:[c:3]",
        "required_environment": ["Acids", "Catalysts"],
        "description": "Friedel-Crafts 아실화 부반응",
    },
    "Aldol_Condensation": {
        "smirks": "[C:1]=[O:2].[C:3]=[O:4]>>[C:1]([C:3]=[O:4])[O:2]",
        "required_environment": ["Bases"],
        "description": "Aldol 축합 부반응",
    },
    "Nitro_Reduction": {
        "smirks": "[c:1][N+1](=O)[O-1]>>[c:1][NH2]",
        "required_environment": ["Catalysts"],
        "description": "니트로기 환원 부산물",
    },
    "Nitrosamine_Formation": {
        "smirks": "[NH1:1][C:2].[O:3]=[N:4]=[O:5]>>[N:1]([C:2])[N:4]=[O:5]",
        "required_environment": ["Acids"],
        "description": "NDSRI 니트로사민 형성 — 산성 아질산 조건",
    },
    "Williamson_Ether": {
        "smirks": "[C:1][O-1:2].[C:3][Br,I,Cl:4]>>[C:1][O:2][C:3]",
        "required_environment": ["Bases"],
        "description": "Williamson 에테르 합성",
    },
    "SN2_N_Alkylation": {
        "smirks": "[N:1].[C:2][Cl,Br,I:3]>>[N:1][C:2]",
        "required_environment": ["Bases"],
        "description": "N-알킬화 / 과알킬화 부반응",
    },
    "Heck_Coupling": {
        "smirks": "[c:1]-[Cl,Br,I:2].[CH2:3]=[CH2:4]>>[c:1]-[CH2:3]-[CH2:4]",
        "required_environment": ["Catalysts", "Bases"],
        "required_catalyst_metal": "Pd",
        "description": "Heck 결합 — Pd 촉매 + 염기",
    },
    "Ester_Hydrolysis_Acidic": {
        "smirks": "[C:1](=[O:2])O[O:3]>>[C:1](=[O:2])[OH].[OH:3]",
        "required_environment": ["Acids"],
        "description": "산성 에스터 가수분해",
    },
    "Ester_Hydrolysis_Basic": {
        "smirks": "[C:1](=[O:2])O[O:3]>>[C:1](=[O:2])[O-1].[O:3]",
        "required_environment": ["Bases"],
        "description": "염기성 에스터 가수분해",
    },
    "Sonogashira_Coupling": {
        "smirks": "[c:1]-[Cl,Br,I:2].[C:3]#[CH:4]>>[c:1]-[C:3]#[CH:4]",
        "required_environment": ["Catalysts"],
        "required_catalyst_metal": "Pd",
        "description": "Sonogashira 결합 — Pd/Cu 촉매 환경",
    },
}

# 하위 호환 별칭
CORE_PROCESS_REACTION_CATALOG = DEFAULT_PROCESS_REACTION_CATALOG


def _normalize_engine_catalog_entry(name: str, raw: dict[str, Any]) -> dict[str, Any] | None:
    smirks = str(
        raw.get("smirks")
        or raw.get("main_smarts")
        or raw.get("smarts")
        or ""
    ).strip()
    if not smirks or ">>" not in smirks:
        return None

    req_env = raw.get("required_environment")
    if isinstance(req_env, str):
        req_env = [req_env]
    elif not isinstance(req_env, list):
        req_env = []
    req_env = [str(c) for c in req_env if c and str(c) in REGISTRY_ENVIRONMENT_CATEGORIES]

    entry: dict[str, Any] = {
        "smirks": smirks,
        "description": str(
            raw.get("description")
            or raw.get("hint_a")
            or name.replace("_", " ")
        ),
        "required_environment": req_env,
    }
    metal = raw.get("required_catalyst_metal")
    if metal:
        entry["required_catalyst_metal"] = str(metal)
    return entry


def _build_active_environment(reagent_context: dict[str, Any] | None) -> dict[str, bool]:
    """Registry UI 선택 → 활성 카테고리 플래그."""
    ctx = reagent_context or {}
    return {
        "Solvents": bool(ctx.get("solvents") or ctx.get("solvent_ids")),
        "Acids": bool(ctx.get("acids") or ctx.get("acid_ids")),
        "Bases": bool(ctx.get("bases") or ctx.get("base_ids")),
        "Catalysts": bool(ctx.get("catalysts") or ctx.get("catalyst_ids")),
    }


def _selected_catalyst_metals(reagent_context: dict[str, Any] | None) -> set[str]:
    ctx = reagent_context or {}
    metals: set[str] = set()
    for entry in ctx.get("registry_entries") or []:
        if str(entry.get("category")) != "Catalysts":
            continue
        metal = str(entry.get("metal") or "").strip()
        if metal and metal != "None":
            metals.add(metal)
    return metals


def _reaction_environment_satisfied(
    rxn_info: dict[str, Any],
    active_environment: dict[str, bool] | None,
    reagent_context: dict[str, Any] | None,
) -> bool:
    """required_environment 및 선택적 Pd/Ni/Cu 촉매 금속 조건 검사."""
    required = list(rxn_info.get("required_environment") or [])
    if required:
        active = active_environment or {}
        for category in required:
            if not active.get(category, False):
                return False

    metal_req = str(rxn_info.get("required_catalyst_metal") or "").strip()
    if metal_req:
        selected_metals = _selected_catalyst_metals(reagent_context)
        if metal_req not in selected_metals:
            return False
    return True


def _sanitize_product_mol(prod: Chem.Mol | None) -> str | None:
    """원자가 무결성 검증 — Valency 깨진 가짜 양성 구조 제거."""
    if prod is None:
        return None
    try:
        mol = Chem.Mol(prod)
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol, canonical=True)
        if not smi:
            return None
        verify = Chem.MolFromSmiles(smi)
        if verify is None:
            return None
        verify.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(verify)
        return Chem.MolToSmiles(verify, canonical=True)
    except Exception:
        log_invalid_smiles_bypass(context="engine_process_product")
        return None


# ── 다성분(Raw Material) 부반응 SMIRKS — 용매·산·염기·촉매 직접 참여 ─────────
RAW_MATERIAL_MULTI_COMPONENT_CATALOG: dict[str, dict[str, Any]] = {
    "Acid_Catalyzed_Solvent_Esterification": {
        "smirks": "[CX3:1](=[O:2])[OH:3].[OX2H:4][CX4:5]>>[CX3:1](=[O:2])[OX2:4][CX4:5].O",
        "requires": ["Acids", "Solvents"],
        "modes": ["parent_env"],
        "reactant_slots": ["parent", "Solvents"],
        "description": "산성 조건 하 메탄올/에탄올 용매 상호작용에 의한 에스테르 유연물질 형성",
    },
    "Base_Mediated_Solvent_Alkylation": {
        "smirks": "[NX3;H2,H1:1].[Cl:2][CX4:3][Cl:4]>>[NX3:1][CX4:3][Cl:4].[Cl-]",
        "requires": ["Bases", "Solvents"],
        "modes": ["parent_env", "env_only"],
        "reactant_slots": ["parent", "Solvents"],
        "env_only_slots": ["Bases", "Solvents"],
        "description": "염기성 조건 하 Halogenated 용매(DCM)의 분해 및 아민 잔류 알킬화 불순물 형성",
    },
    "Catalyst_Induced_Aromatic_Homocoupling": {
        "smirks": "[c:1]-[I,Br,Cl].[c:2]-[I,Br,Cl]>>[c:1]-[c:2]",
        "requires": ["Catalysts"],
        "modes": ["parent_env"],
        "description": "금속 촉매(Pd/Ni) 환경 하 할로겐화 아릴의 자체 호모커플링 이량체 유연물질 형성",
    },
    "Acid_Nitrite_Nitrosamination": {
        "smirks": "[NX3;H1:1].[N:2](=[O:3])[O-]>>[NX3:1][N:2]=[O:3]",
        "requires": ["Acids"],
        "modes": ["parent_env", "env_only"],
        "reactant_slots": ["parent", "Acids"],
        "env_only_slots": ["Bases", "Acids"],
        "description": "산성 조건 하 아질산염 혼입에 따른 N-니트로사민 발암성 불순물 형성 리스크",
    },
}

REACTION_MODE_PARENT_ENV: str = "parent_env"
REACTION_MODE_ENV_ONLY: str = "env_only"


def _build_selected_materials_from_context(
    reagent_context: dict[str, Any] | None,
) -> dict[str, list[str]]:
    """
    Registry UI 컨텍스트 → 카테고리별 SMILES 맵.

    {"Solvents": [...], "Acids": [...], "Bases": [...], "Catalysts": [...]}
    """
    from engines.p13_chemical_registry import build_selected_materials_payload

    return build_selected_materials_payload(reagent_context)


def _normalize_selected_materials(
    selected_materials: dict[str, list[str]] | None,
    reagent_context: dict[str, Any] | None = None,
) -> dict[str, list[str]]:
    """Capitalized 카테고리 키 + SMILES 배열로 정규화."""
    materials = _build_selected_materials_from_context(reagent_context)
    if not selected_materials:
        return materials

    normalized: dict[str, list[str]] = {
        cat: [] for cat in REGISTRY_ENVIRONMENT_CATEGORIES
    }
    for cat in REGISTRY_ENVIRONMENT_CATEGORIES:
        candidates = (
            selected_materials.get(cat)
            or selected_materials.get(cat.lower())
            or []
        )
        seen: set[str] = set()
        for smi in candidates:
            smi = str(smi).strip()
            if smi and smi not in seen:
                seen.add(smi)
                normalized[cat].append(smi)

    if any(normalized.values()):
        return normalized
    return materials


def _pool_env_mols_by_category(
    selected_materials: dict[str, list[str]],
) -> tuple[dict[str, list[Chem.Mol]], set[str]]:
    """카테고리별 SMILES → RDKit Mol 맵 + 활성 카테고리."""
    env_by_cat: dict[str, list[Chem.Mol]] = {
        cat: [] for cat in REGISTRY_ENVIRONMENT_CATEGORIES
    }
    active_categories: set[str] = set()
    seen_smiles: set[str] = set()

    for category in REGISTRY_ENVIRONMENT_CATEGORIES:
        for smi in selected_materials.get(category) or []:
            smi = str(smi).strip()
            if not smi or smi in seen_smiles:
                continue
            mol = mol_from_smiles(smi, log_bypass=False)
            if mol is None:
                continue
            seen_smiles.add(smi)
            env_by_cat[category].append(mol)
            active_categories.add(category)

    return env_by_cat, active_categories


def _flatten_env_mols(env_by_cat: dict[str, list[Chem.Mol]]) -> list[Chem.Mol]:
    flat: list[Chem.Mol] = []
    seen: set[str] = set()
    for category in REGISTRY_ENVIRONMENT_CATEGORIES:
        for mol in env_by_cat.get(category) or []:
            smi = Chem.MolToSmiles(mol, canonical=True)
            if smi not in seen:
                seen.add(smi)
                flat.append(mol)
    return flat


def _pool_env_mols_from_materials(
    selected_materials: dict[str, list[str]],
) -> tuple[list[Chem.Mol], set[str]]:
    """카테고리별 SMILES → flat RDKit Mol 풀 + 활성 카테고리 집합."""
    env_by_cat, active = _pool_env_mols_by_category(selected_materials)
    return _flatten_env_mols(env_by_cat), active


def _resolve_slot_mol_pool(
    slot: str,
    *,
    parent_mol: Chem.Mol | None,
    env_by_cat: dict[str, list[Chem.Mol]],
) -> list[Chem.Mol]:
    """반응 템플릿 슬롯 → 후보 Mol 리스트 (parent | 카테고리명)."""
    if slot == "parent":
        return [parent_mol] if parent_mol is not None else []
    if slot in env_by_cat:
        return list(env_by_cat.get(slot) or [])
    return []


def _enumerate_categorized_reactant_sets(
    slot_specs: list[str],
    *,
    parent_mol: Chem.Mol | None,
    env_by_cat: dict[str, list[Chem.Mol]],
) -> list[tuple[Chem.Mol, ...]]:
    """슬롯별 카테고리/parent 지정 → RunReactants 입력 Cartesian product."""
    if not slot_specs:
        return []
    pools = [
        _resolve_slot_mol_pool(slot, parent_mol=parent_mol, env_by_cat=env_by_cat)
        for slot in slot_specs
    ]
    if any(not pool for pool in pools):
        return []
    return [tuple(combo) for combo in product(*pools)]


def _enumerate_env_only_pairs(env_mols: list[Chem.Mol]) -> list[tuple[Chem.Mol, ...]]:
    """시약 풀 내 이성분 조합 (주성분 제외, 서로 다른 분자 우선)."""
    pairs: list[tuple[Chem.Mol, ...]] = []
    for i, m1 in enumerate(env_mols):
        for j, m2 in enumerate(env_mols):
            if i == j:
                continue
            pairs.append((m1, m2))
    return pairs


def _enumerate_reactant_sets(
    parent_mol: Chem.Mol,
    env_mols: list[Chem.Mol],
    n_templates: int,
) -> list[tuple[Chem.Mol, ...]]:
    """
    RDKit RunReactants 입력 튜플 생성 (레거시 flat 풀).

    - 1-template: (parent_mol,)
    - 2+-template: slot 0 = parent_mol 고정, slot 1+ = env_mols 교차 대입
    """
    if n_templates <= 0:
        return []
    if n_templates == 1:
        return [(parent_mol,)]
    if not env_mols:
        return []
    if n_templates == 2:
        return [(parent_mol, env_mol) for env_mol in env_mols]
    return [
        (parent_mol,) + combo
        for combo in product(env_mols, repeat=n_templates - 1)
    ]


def _enumerate_reactant_sets_for_rule(
    rxn_info: dict[str, Any],
    n_templates: int,
    *,
    parent_mol: Chem.Mol | None,
    env_by_cat: dict[str, list[Chem.Mol]],
    env_mols_flat: list[Chem.Mol],
) -> list[tuple[tuple[Chem.Mol, ...], str]]:
    """
    규칙별 반응 모드(parent_env / env_only)에 따른 RunReactants 입력 생성.

    Returns list of (reactant_tuple, simulation_subtype).
    """
    modes = list(rxn_info.get("modes") or [REACTION_MODE_PARENT_ENV])
    parent_slots = rxn_info.get("reactant_slots")
    env_slots = rxn_info.get("env_only_slots")

    if isinstance(parent_slots, str):
        parent_slots = [parent_slots]
    if isinstance(env_slots, str):
        env_slots = [env_slots]

    results: list[tuple[tuple[Chem.Mol, ...], str]] = []
    seen: set[tuple[str, ...]] = set()

    def _register(reactants: tuple[Chem.Mol, ...], subtype: str) -> None:
        if len(reactants) != n_templates:
            return
        key = tuple(Chem.MolToSmiles(m, canonical=True) for m in reactants)
        dedupe_key = (subtype, *key)
        if dedupe_key in seen:
            return
        seen.add(dedupe_key)
        results.append((reactants, subtype))

    if REACTION_MODE_PARENT_ENV in modes and parent_mol is not None:
        if isinstance(parent_slots, list) and len(parent_slots) == n_templates:
            for reactants in _enumerate_categorized_reactant_sets(
                parent_slots,
                parent_mol=parent_mol,
                env_by_cat=env_by_cat,
            ):
                _register(reactants, "parent_env_multi_component")
        else:
            for reactants in _enumerate_reactant_sets(
                parent_mol, env_mols_flat, n_templates
            ):
                _register(reactants, "parent_env_multi_component")

    if REACTION_MODE_ENV_ONLY in modes and env_mols_flat:
        if isinstance(env_slots, list) and len(env_slots) == n_templates:
            for reactants in _enumerate_categorized_reactant_sets(
                env_slots,
                parent_mol=None,
                env_by_cat=env_by_cat,
            ):
                _register(reactants, "env_only_multi_component")
        elif n_templates == 2:
            for reactants in _enumerate_env_only_pairs(env_mols_flat):
                _register(reactants, "env_only_multi_component")
        elif n_templates == 1:
            for mol in env_mols_flat:
                _register((mol,), "env_only_multi_component")
        else:
            for combo in product(env_mols_flat, repeat=n_templates):
                _register(tuple(combo), "env_only_multi_component")

    return results


def _append_products_from_outcomes(
    outcomes: Any,
    rxn_name: str,
    rxn_info: dict[str, Any],
    predicted: list[dict[str, Any]],
    seen_smiles: set[str],
    *,
    simulation_type: str = "multi_component",
) -> None:
    """RunReactants 결과 → canonical SMILES 예측 레코드."""
    for product_tuple in outcomes or []:
        if product_tuple is None:
            continue
        items = (
            product_tuple
            if isinstance(product_tuple, (tuple, list))
            else (product_tuple,)
        )
        for prod in items:
            prod_smiles = _sanitize_product_mol(prod)
            if not prod_smiles or prod_smiles in seen_smiles:
                continue
            seen_smiles.add(prod_smiles)
            predicted.append(
                {
                    "smiles": prod_smiles,
                    "origin_reaction": rxn_name,
                    "mechanism": rxn_info.get("description", rxn_name),
                    "simulation_type": simulation_type,
                }
            )


def _run_smirks_catalog_loop(
    parent_mol: Chem.Mol | None,
    env_by_cat: dict[str, list[Chem.Mol]],
    active_categories: set[str],
    catalog: dict[str, dict[str, Any]],
    *,
    requires_key: str = "requires",
    reaction_filter: Callable[[dict[str, Any], set[str]], bool] | None = None,
    simulation_type: str = "multi_component",
) -> list[dict[str, Any]]:
    """카탈로그 SMIRKS × parent_env / env_only reactant set 순회."""
    predicted: list[dict[str, Any]] = []
    seen_smiles: set[str] = set()
    env_mols_flat = _flatten_env_mols(env_by_cat)

    for rxn_name, rxn_info in catalog.items():
        required = set(rxn_info.get(requires_key) or [])
        if required and not required.issubset(active_categories):
            continue
        if reaction_filter is not None and not reaction_filter(rxn_info, active_categories):
            continue

        smirks = str(rxn_info.get("smirks") or "")
        if ">>" not in smirks:
            continue
        try:
            rxn = AllChem.ReactionFromSmarts(smirks)
        except Exception:
            continue
        if rxn is None:
            continue

        n_templates = rxn.GetNumReactantTemplates()
        reactant_entries = _enumerate_reactant_sets_for_rule(
            rxn_info,
            n_templates,
            parent_mol=parent_mol,
            env_by_cat=env_by_cat,
            env_mols_flat=env_mols_flat,
        )
        for reactants, subtype in reactant_entries:
            try:
                outcomes = rxn.RunReactants(reactants)
            except Exception:
                continue
            _append_products_from_outcomes(
                outcomes,
                rxn_name,
                rxn_info,
                predicted,
                seen_smiles,
                simulation_type=subtype or simulation_type,
            )

    return predicted


def run_reagent_pool_cross_simulation(
    selected_materials: dict[str, list[str]] | None,
    catalog: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """
    주성분 없이 시약 풀(env_mols) 상호작용만 독립 시뮬레이션.

    Base + Halogenated Solvent 알킬화 등 env_only 규칙 전용 batch hook.
    """
    materials = _normalize_selected_materials(selected_materials, None)
    env_by_cat, active_categories = _pool_env_mols_by_category(materials)
    if not _flatten_env_mols(env_by_cat):
        return []

    engine_catalog = catalog or RAW_MATERIAL_MULTI_COMPONENT_CATALOG
    env_catalog = {
        name: info
        for name, info in engine_catalog.items()
        if REACTION_MODE_ENV_ONLY in (info.get("modes") or [])
    }
    return _run_smirks_catalog_loop(
        None,
        env_by_cat,
        active_categories,
        env_catalog,
        requires_key="requires",
        simulation_type="env_only_multi_component",
    )


class AdvancedProcessImpurityEngine:
    """
    용매·산·염기·촉매가 구조적으로 결합/분해에 참여하는
    다성분(Multi-component) 부반응 SMIRKS 엔진.
    """

    def __init__(self) -> None:
        self.reaction_catalog: dict[str, dict[str, Any]] = {}
        self._initialize_raw_material_reactions()

    def _initialize_raw_material_reactions(self) -> None:
        """오픈소스 기반 다성분 부반응 SMIRKS 매트릭스 시드."""
        self.reaction_catalog = {
            name: dict(spec) for name, spec in RAW_MATERIAL_MULTI_COMPONENT_CATALOG.items()
        }

    def load_reaction_catalog_from_json(self, file_path: str) -> int:
        """requires 키를 갖는 raw-material JSON 병합."""
        path = Path(file_path)
        if not path.is_file():
            return 0
        added = 0
        try:
            with path.open(encoding="utf-8") as fh:
                payload = json.load(fh)
        except Exception:
            return 0
        items = (
            payload.items()
            if isinstance(payload, dict) and "reactions" not in payload
            else ((str(r.get("name", i)), r) for i, r in enumerate(payload.get("reactions", [])))
        )
        for name, raw in items:
            if not isinstance(raw, dict):
                continue
            smirks = str(raw.get("smirks") or raw.get("main_smarts") or "").strip()
            if ">>" not in smirks:
                continue
            requires = raw.get("requires") or raw.get("required_environment") or []
            if isinstance(requires, str):
                requires = [requires]
            key = str(name).replace(" ", "_")
            if key not in self.reaction_catalog:
                added += 1
            self.reaction_catalog[key] = {
                "smirks": smirks,
                "requires": [str(c) for c in requires if c],
                "description": str(raw.get("description") or key),
                "modes": raw.get("modes") or [REACTION_MODE_PARENT_ENV],
                "reactant_slots": raw.get("reactant_slots"),
                "env_only_slots": raw.get("env_only_slots"),
            }
        return added

    def run_raw_material_simulation(
        self,
        parent_smiles: str,
        selected_materials: dict[str, list[str]] | None,
        *,
        reagent_context: dict[str, Any] | None = None,
    ) -> list[dict[str, Any]]:
        """
        parent + 카테고리별 선택 시약 SMILES → 다성분 부반응 순회.

        selected_materials: {"Solvents": ["SMILES", ...], "Acids": [...], ...}
        """
        parent_mol = mol_from_smiles(parent_smiles, log_bypass=False)
        if parent_mol is None:
            return []

        materials = _normalize_selected_materials(selected_materials, reagent_context)
        env_by_cat, active_categories = _pool_env_mols_by_category(materials)
        if not _flatten_env_mols(env_by_cat) and not active_categories:
            return []

        return _run_smirks_catalog_loop(
            parent_mol,
            env_by_cat,
            active_categories,
            self.reaction_catalog,
            requires_key="requires",
            simulation_type="raw_material_multi_component",
        )

    def run_reagent_cross_simulation(
        self,
        selected_materials: dict[str, list[str]] | None,
        *,
        reagent_context: dict[str, Any] | None = None,
    ) -> list[dict[str, Any]]:
        """시약 풀 단독(env_only) 부반응 — 주성분 없이 RunReactants."""
        materials = _normalize_selected_materials(selected_materials, reagent_context)
        return run_reagent_pool_cross_simulation(materials, self.reaction_catalog)

    def _parse_products(
        self,
        products_set: Any,
        rxn_name: str,
        rxn_info: dict[str, Any],
        predicted_list: list[dict[str, Any]],
        seen_smiles: set[str],
    ) -> None:
        """생성 Mol → Sanitize → canonical SMILES 후 predicted_list 주입."""
        for product_tuple in products_set or []:
            if product_tuple is None:
                continue
            items = (
                product_tuple
                if isinstance(product_tuple, (tuple, list))
                else (product_tuple,)
            )
            for prod in items:
                prod_smiles = _sanitize_product_mol(prod)
                if not prod_smiles or prod_smiles in seen_smiles:
                    continue
                seen_smiles.add(prod_smiles)
                predicted_list.append(
                    {
                        "smiles": prod_smiles,
                        "origin_reaction": rxn_name,
                        "mechanism": rxn_info.get("description", rxn_name),
                        "simulation_type": "raw_material_multi_component",
                    }
                )

    def predictions_to_impurity_records(
        self,
        predictions: list[dict[str, Any]],
        *,
        source_label: str,
        temperature_c: float = 25.0,
    ) -> list[dict[str, Any]]:
        """raw material 예측 → ProcessImpurityEngine batch 레코드."""
        engine = NameReactionProcessEngine()
        return engine.predictions_to_impurity_records(
            predictions,
            source_label=source_label,
            temperature_c=temperature_c,
        )


_ADVANCED_ENGINE: AdvancedProcessImpurityEngine | None = None


def get_advanced_process_engine() -> AdvancedProcessImpurityEngine:
    global _ADVANCED_ENGINE
    if _ADVANCED_ENGINE is None:
        _ADVANCED_ENGINE = AdvancedProcessImpurityEngine()
    return _ADVANCED_ENGINE


class NameReactionProcessEngine:
    """
    Registry 시약 SMILES + 반응 카탈로그 SMIRKS 매칭 기반 가상 반응 트리 엔진.

    required_environment 조건을 만족하는 반응만 parent + 시약 환경에 대해 실행한다.
    """

    def __init__(self) -> None:
        self.reaction_catalog: dict[str, dict[str, Any]] = {}
        self._initialize_base_catalog()
        self.load_catalogs_from_directory()

    def _initialize_base_catalog(self) -> None:
        """제약 부반응 중심 내장 매트릭스 (10+ 규칙, required_environment 포함)."""
        self.reaction_catalog = {
            name: dict(info) for name, info in DEFAULT_PROCESS_REACTION_CATALOG.items()
        }

    def load_reaction_catalog_from_json(self, file_path: str) -> int:
        """
        대용량 오픈소스 반응 카탈로그 JSON 주입.

        지원 형식
        ---------
        - Phase 6 표준: {"Reaction_Name": {"smirks", "required_environment", "description"}}
        - USPTO/ORD: {"reactions": [{"name", "main_smarts", ...}, ...]}
        - 레거시 flat dict (main_smarts 키)
        """
        path = Path(file_path)
        if not path.is_file():
            log_process_warning(f"Reaction catalog JSON not found: {path}")
            return 0
        added = 0
        try:
            with path.open(encoding="utf-8") as fh:
                payload = json.load(fh)
        except Exception as exc:
            log_process_warning(f"Reaction catalog JSON parse failed ({path}): {exc}")
            return 0

        if isinstance(payload, dict) and "reactions" not in payload:
            items: list[tuple[str, dict[str, Any]]] = [
                (str(k), v) for k, v in payload.items() if isinstance(v, dict)
            ]
        else:
            loaded = load_reaction_catalog(path)
            items = [(name, spec) for name, spec in loaded.items()]

        for name, raw in items:
            normalized = _normalize_engine_catalog_entry(name, raw)
            if not normalized:
                continue
            key = name.replace(" ", "_").replace("/", "_")
            if key not in self.reaction_catalog:
                added += 1
            self.reaction_catalog[key] = normalized
        return added

    def load_external_uspto_catalog(self, json_path: str | Path) -> int:
        """USPTO/ORD JSON 로드 — load_reaction_catalog_from_json 별칭."""
        return self.load_reaction_catalog_from_json(str(json_path))

    def load_catalogs_from_directory(self, directory: str | Path | None = None) -> int:
        """data/reaction_catalogs/*.json 일괄 병합."""
        total = 0
        dir_path = Path(directory) if directory else _CATALOG_DIR
        if not dir_path.is_dir():
            return 0
        for json_file in sorted(dir_path.glob("*.json")):
            total += self.load_reaction_catalog_from_json(str(json_file))
        return total

    def run_process_simulation(
        self,
        parent_smiles: str,
        environment_smiles_list: list[str] | None,
        *,
        active_environment: dict[str, bool] | None = None,
        reagent_context: dict[str, Any] | None = None,
        selected_materials: dict[str, list[str]] | None = None,
    ) -> list[dict[str, Any]]:
        """
        parent + Registry 시약 SMILES 환경에 대해 카탈로그 반응을 순회 매칭.

        required_environment 조건이 충족된 반응만 RDKit RunReactants 실행.
        """
        parent_mol = mol_from_smiles(parent_smiles, log_bypass=False)
        if parent_mol is None:
            return []

        materials = _normalize_selected_materials(selected_materials, reagent_context)
        env_by_cat, active_categories = _pool_env_mols_by_category(materials)
        env_mols_flat = _flatten_env_mols(env_by_cat)

        # 레거시 flat SMILES 리스트 — 카테고리 풀에 미포함 구조만 추가
        pooled_smiles = {
            Chem.MolToSmiles(m, canonical=True)
            for m in env_mols_flat
            if m is not None
        }
        for smi in environment_smiles_list or []:
            smi = str(smi).strip()
            if not smi or smi in pooled_smiles:
                continue
            mol = mol_from_smiles(smi, log_bypass=False)
            if mol is not None:
                env_mols_flat.append(mol)
                env_by_cat.setdefault("Solvents", []).append(mol)
                active_categories.add("Solvents")
                pooled_smiles.add(Chem.MolToSmiles(mol, canonical=True))

        if active_environment:
            for category, is_active in active_environment.items():
                if is_active:
                    active_categories.add(category)

        def _process_filter(
            rxn_info: dict[str, Any],
            _active: set[str],
        ) -> bool:
            return _reaction_environment_satisfied(
                rxn_info, active_environment, reagent_context
            )

        return _run_smirks_catalog_loop(
            parent_mol,
            env_by_cat,
            active_categories,
            self.reaction_catalog,
            requires_key="required_environment",
            reaction_filter=_process_filter,
            simulation_type="process_catalog",
        )

    def predictions_to_impurity_records(
        self,
        predictions: list[dict[str, Any]],
        *,
        source_label: str,
        temperature_c: float = 25.0,
        db_fps: list | None = None,
        db_smiles: list[str] | None = None,
    ) -> list[dict[str, Any]]:
        """run_process_simulation 결과 → ProcessImpurityEngine batch 레코드."""
        if db_fps is None or db_smiles is None:
            db_fps, db_smiles = get_verification_fingerprint_db()

        records: list[dict[str, Any]] = []
        for pred in predictions:
            smi = str(pred.get("smiles") or "").strip()
            mol = mol_from_smiles(smi, log_bypass=False)
            if mol is None:
                continue
            origin = str(pred.get("origin_reaction", "Catalog"))
            mechanism = str(pred.get("mechanism", ""))
            pathway = f"{PATHWAY_PREFIX}: Catalog · {origin}"
            if mechanism:
                pathway = f"{pathway} · {mechanism}"

            exact_mw = round(Descriptors.ExactMolWt(mol), 2)
            ich_class, ich_detail = classify_ich_m7_genotoxicity(mol)
            match_pct, status_tag, reliability_pct = compute_db_verification_score(
                mol, db_fps, db_smiles, temperature_c=temperature_c
            )
            records.append(
                {
                    "pathway": pathway,
                    "smiles": smi,
                    "exact_mw": exact_mw,
                    "ich_m7": ich_class,
                    "ich_m7_detail": ich_detail,
                    "db_match_pct": match_pct,
                    "db_status_tag": status_tag,
                    "reliability_pct": reliability_pct,
                    "confidence_display": f"{reliability_pct:.1f}% — {status_tag}",
                    "source_reactant": source_label,
                    "name_reaction": origin,
                    "is_name_reaction_byproduct": True,
                    "catalog_mechanism": mechanism,
                }
            )
        return records

    def run_batch_environment_simulation(
        self,
        reactant_entries: list[tuple[str, str, dict[str, Any]]],
        environment_smiles_list: list[str],
        *,
        reagent_context: dict[str, Any] | None = None,
        selected_materials: dict[str, list[str]] | None = None,
        temperature_c: float = 25.0,
        ctx: ProcessRunContext | None = None,
    ) -> list[dict[str, Any]]:
        """다중 반응물 × Registry 시약 SMILES 다성분 카탈로그 순회."""
        all_records: list[dict[str, Any]] = []
        seen: set[str] = set()
        active_env = _build_active_environment(reagent_context)
        materials = _normalize_selected_materials(selected_materials, reagent_context)
        advanced = get_advanced_process_engine()

        for label, _key, props in reactant_entries:
            if ctx is not None and not ctx.check_continue():
                break
            if not props.get("valid") or not props.get("smiles"):
                continue
            parent_smi = str(props["smiles"]).strip()

            raw_preds = advanced.run_raw_material_simulation(
                parent_smi,
                materials,
                reagent_context=reagent_context,
            )
            proc_preds = self.run_process_simulation(
                parent_smi,
                environment_smiles_list,
                active_environment=active_env,
                reagent_context=reagent_context,
                selected_materials=materials,
            )

            preds: list[dict[str, Any]] = []
            seen_pred_smiles: set[str] = set()
            for pred in raw_preds + proc_preds:
                smi = pred.get("smiles")
                if smi and smi not in seen_pred_smiles:
                    seen_pred_smiles.add(smi)
                    preds.append(pred)

            batch = self.predictions_to_impurity_records(
                preds,
                source_label=f"{label} · Catalog/env",
                temperature_c=temperature_c,
            )
            for rec in batch:
                if ctx is not None and not ctx.check_continue():
                    break
                smi = rec.get("smiles")
                if not smi or smi in seen:
                    continue
                seen.add(smi)
                all_records.append(rec)

        if ctx is None or ctx.check_continue():
            env_preds = advanced.run_reagent_cross_simulation(
                materials,
                reagent_context=reagent_context,
            )
            if env_preds:
                env_batch = self.predictions_to_impurity_records(
                    env_preds,
                    source_label="Reagent pool · env_only",
                    temperature_c=temperature_c,
                )
                for rec in env_batch:
                    if ctx is not None and not ctx.check_continue():
                        break
                    smi = rec.get("smiles")
                    if not smi or smi in seen:
                        continue
                    seen.add(smi)
                    all_records.append(rec)

        return ctx.truncate_records(all_records) if ctx else all_records


# 싱글톤 — 앱 세션 간 카탈로그 캐시 재사용
_PROCESS_ENGINE: NameReactionProcessEngine | None = None


def get_process_reaction_engine(*, reload_catalog: bool = False) -> NameReactionProcessEngine:
    global _PROCESS_ENGINE
    if _PROCESS_ENGINE is None or reload_catalog:
        _PROCESS_ENGINE = NameReactionProcessEngine()
        _PROCESS_ENGINE.load_catalogs_from_directory()
    return _PROCESS_ENGINE


class NameReactionProcess:
    """PRIME-DS Name Reaction(rdchiral) 공정 유연물질 시뮬레이션 파사드."""

    STANDARD_DEGRADATION = NAME_RXN_STANDARD_DEGRADATION
    NONE = NAME_RXN_NONE
    NAME_RXN_NONE = NAME_RXN_NONE
    NAME_RXN_STANDARD_DEGRADATION = NAME_RXN_STANDARD_DEGRADATION
    NAME_REACTION_UI_OPTIONS = NAME_REACTION_UI_OPTIONS

    get_catalog = staticmethod(get_reaction_catalog)
    load_catalog = staticmethod(load_reaction_catalog)
    get_spec = staticmethod(get_reaction_spec)
    get_ui_options = staticmethod(get_reaction_ui_options)
    is_skipped = staticmethod(is_name_reaction_skipped)

    simulate_impurities = staticmethod(simulate_name_reaction_impurities)
    simulate_suzuki_competitive = staticmethod(simulate_suzuki_competitive_screening)
    run_universal_prediction = staticmethod(run_universal_side_reaction_prediction)
    merge_into_batch = staticmethod(merge_name_reaction_into_batch)
    run_rdchiral = staticmethod(rdchiral_run_reactants)
    get_process_engine = staticmethod(get_process_reaction_engine)
    get_advanced_engine = staticmethod(get_advanced_process_engine)
    ProcessEngine = NameReactionProcessEngine
    AdvancedProcessEngine = AdvancedProcessImpurityEngine
