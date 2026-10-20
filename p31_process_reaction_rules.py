"""
PRIME-DS process reaction SMIRKS — main product vs side impurity separation.

Binary reaction templates with regioselectivity-safe reactant pairing.
"""

# 경로 A

# 주/부생성물 구분(2반응)

from __future__ import annotations

import json
from itertools import combinations
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem

from p23_chem_search_engine import safe_canonical_smiles

# -------------------------------------------------------------------------
# 화학 공정 반응식 규칙 (주반응 + 부반응)
# -------------------------------------------------------------------------
PROCESS_REACTION_RULES: dict[str, dict[str, str]] = {
    "Amide_Coupling": {
        # 2026-09-03 수정: [NX3:3]에 배제 조건이 없어 아마이드 질소(→ 상온 공정에서 일어나지
        #   않는 이미드 생성)와 3급 아민(→ 원자가 초과로 조용히 폐기)까지 반응물로 받았습니다.
        #   p30의 대응 규칙들과 같은 기준(!$(NC=O) 등)으로 1·2급 아민만 받도록 좁혔습니다.
        "main_smirks": (
            "[CX3:1](=[O:2])[OX2H]"
            ".[NX3;H2,H1;!$(NC=O);!$(N=*);!$(NS(=O)=O):3]"
            ">>[CX3:1](=[O:2])[NX3:3]"
        ),
        "side_smirks": (
            "[CX3:1](=[O:2])[OX2H].[CX3:3](=[O:4])[OX2H]"
            ">>[CX3:1](=[O:2])[OX2][CX3:3](=[O:4])"
        ),
        "description": (
            "정상적인 아미드 결합 공정이지만, 과량의 산 조건 시 "
            "무수물(Anhydride) 불순물 리스크 존재."
        ),
    },
    "Nucleophilic_Substitution_with_Azide": {
        # 2026-09-03 수정: 기존 산물 [CX4:1][N-]=[N+]=N 은 음전하 질소가 3결합이 되어
        #   원자가 위반 → SanitizeMol 예외 → 산물 0건이었습니다. 아자이드 세 질소를 모두
        #   매핑하고 알킬 결합이 붙는 질소의 전하를 명시적으로 0으로 지정해 교정했습니다.
        "main_smirks": (
            "[CX4:1][Cl,Br,I].[N-:5]=[N+:6]=[N:7]"
            ">>[CX4:1][N+0:5]=[N+:6]=[N-:7]"
        ),
        "side_smirks": "[CX4:1][CX4:2][Cl,Br,I]>>[CX4:1]=[CX4:2]",
        "description": (
            "아지드 치환 주반응 진행 중, 경쟁적인 E2 제거반응으로 인한 "
            "알켄 불순물 혼입 리스크."
        ),
    },
}


def _normalize_smirks(smirks: str) -> str:
    return (smirks or "").replace(" ", "").strip()


def _canonical_product_smiles(mol: Chem.Mol) -> str | None:
    try:
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return smi if Chem.MolFromSmiles(smi) else None
    except Exception:
        return None


def _mols_from_smiles_list(reactant_smiles_list: list[str]) -> list[Chem.Mol]:
    mols: list[Chem.Mol] = []
    seen: set[str] = set()
    for raw in reactant_smiles_list:
        smi = (raw or "").strip()
        if not smi:
            continue
        canon, _ = safe_canonical_smiles(smi)
        canon = canon or smi
        if canon in seen:
            continue
        mol = Chem.MolFromSmiles(canon)
        if mol is None:
            continue
        seen.add(canon)
        mols.append(mol)
    return mols


def _run_binary_smirks(mol_a: Chem.Mol, mol_b: Chem.Mol, smirks: str) -> list[str]:
    smirks = _normalize_smirks(smirks)
    try:
        rxn = AllChem.ReactionFromSmarts(smirks)
    except Exception:
        return []
    if rxn is None:
        return []

    out: list[str] = []
    seen: set[str] = set()

    # 2026-09-03 수정: 반응물 템플릿이 1개인 SMIRKS(E2 제거 등)를 항상 2-튜플로 넘겨
    #   ChemicalParserException이 나고 except에서 조용히 버려졌습니다(부생성물 영구 0건).
    #   GetNumReactantTemplates()에 따라 1-튜플/2-튜플을 골라 호출합니다.
    if rxn.GetNumReactantTemplates() == 1:
        orders: tuple[tuple[Chem.Mol, ...], ...] = ((mol_a,), (mol_b,))
    else:
        orders = ((mol_a, mol_b), (mol_b, mol_a))

    for order in orders:
        try:
            product_sets = rxn.RunReactants(order)
        except Exception:
            continue
        for p_set in product_sets:
            for prod in p_set:
                smi = _canonical_product_smiles(prod)
                if smi and smi not in seen:
                    seen.add(smi)
                    out.append(smi)
    return out


def _reactant_pairs(mols: list[Chem.Mol]) -> list[tuple[Chem.Mol, Chem.Mol]]:
    if len(mols) < 2:
        return []
    return list(combinations(mols, 2))


def predict_process_reactions(reactant_smiles_list: list[str]) -> dict[str, Any]:
    """
    공정 투입물(반응물+시약) → 주반응 생성물 / 부반응 불순물 분류 리포트.
    """
    mols = _mols_from_smiles_list(reactant_smiles_list)
    if len(mols) < 2:
        return {
            "error": "반응을 진행하기 위해서는 최소 2개 이상의 화합물(반응물/시약)이 필요합니다.",
            "Main_Products": [],
            "Side_Products_Impurities": [],
            "Triggered_Mechanisms": [],
        }

    process_report: dict[str, Any] = {
        "Main_Products": [],
        "Side_Products_Impurities": [],
        "Triggered_Mechanisms": [],
    }
    main_seen: set[str] = set()
    side_seen: set[str] = set()
    pairs = _reactant_pairs(mols)

    for mech_name, smirks_info in PROCESS_REACTION_RULES.items():
        main_smirks = smirks_info["main_smirks"]
        side_smirks = smirks_info["side_smirks"]
        description = smirks_info["description"]
        mech_triggered = False

        for mol_a, mol_b in pairs:
            for p_smiles in _run_binary_smirks(mol_a, mol_b, main_smirks):
                if p_smiles in main_seen:
                    continue
                main_seen.add(p_smiles)
                process_report["Main_Products"].append(
                    {"SMILES": p_smiles, "Mechanism": mech_name}
                )
                mech_triggered = True

            for s_smiles in _run_binary_smirks(mol_a, mol_b, side_smirks):
                if s_smiles in main_seen or s_smiles in side_seen:
                    continue
                side_seen.add(s_smiles)
                process_report["Side_Products_Impurities"].append(
                    {
                        "SMILES": s_smiles,
                        "Mechanism": mech_name,
                        "Alert_Description": description,
                    }
                )
                mech_triggered = True

        if mech_triggered and mech_name not in process_report["Triggered_Mechanisms"]:
            process_report["Triggered_Mechanisms"].append(mech_name)

    return process_report


def process_reaction_report_to_json(report: dict[str, Any], *, indent: int = 2) -> str:
    return json.dumps(report, indent=indent, ensure_ascii=False)


if __name__ == "__main__":
    input_process_reagents = ["CC(=O)O", "CCN"]
    print("PRIME-DS process main/side reaction engine")
    result = predict_process_reactions(input_process_reagents)
    print(process_reaction_report_to_json(result))
