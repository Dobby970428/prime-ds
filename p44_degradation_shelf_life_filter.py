"""
ICH Q3A/B · M7 in-vitro 보관 안정성(Shelf-life) Child-node 여과 레이어.

DegradationTreeEngine · predict_advanced_degradants 공통 적용.
"""

# 검수- "보관 중에는 안 일어날 반응(대표적으로 방향족 수산화)"을 걸러냄

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

# 벤젠 고리(6원 방향족) 직접 수산화 — [OX2H]-[c;r6]
AROMATIC_R6_OH_SMARTS = "[OX2H]-[c;r6]"
AROMATIC_R6_OH_QUERY = Chem.MolFromSmarts(AROMATIC_R6_OH_SMARTS)


def _safe_mol(smiles: str) -> Chem.Mol | None:
    raw = (smiles or "").strip()
    if not raw:
        return None
    try:
        mol = Chem.MolFromSmiles(raw)
        if mol is None:
            return None
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def _safe_mol_obj(mol: Chem.Mol | None) -> Chem.Mol | None:
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def parent_has_aromatic_r6_oh(mol: Chem.Mol | None) -> bool:
    """`mol.HasSubstructMatch(aromatic_oh)` — Parent 벤젠고리-OH 존재 여부."""
    mol = _safe_mol_obj(mol)
    if mol is None or AROMATIC_R6_OH_QUERY is None:
        return False
    try:
        return bool(mol.HasSubstructMatch(AROMATIC_R6_OH_QUERY))
    except Exception:
        return False


def child_has_aromatic_r6_oh(mol: Chem.Mol | None) -> bool:
    """`prod_mol.HasSubstructMatch(aromatic_oh)` — Child 벤젠고리-OH 존재 여부."""
    return parent_has_aromatic_r6_oh(mol)


def should_cull_metabolic_aromatic_hydroxylation(
    parent_mol: Chem.Mol | None,
    child_mol: Chem.Mol | None,
) -> bool:
    """
    부모에 없던 벤젠고리-OH가 Child에 생겼으면 True (제외 대상).

    if child_has_oh and not parent_has_oh:
        continue
    """
    parent = _safe_mol_obj(parent_mol)
    child = _safe_mol_obj(child_mol)
    if parent is None or child is None or AROMATIC_R6_OH_QUERY is None:
        return False
    try:
        parent_has_oh = parent.HasSubstructMatch(AROMATIC_R6_OH_QUERY)
        child_has_oh = child.HasSubstructMatch(AROMATIC_R6_OH_QUERY)
        return bool(child_has_oh and not parent_has_oh)
    except Exception:
        return False


# 하위 호환 alias
is_de_novo_aromatic_hydroxylation_r6 = should_cull_metabolic_aromatic_hydroxylation


def is_de_novo_aromatic_hydroxylation_r6_smiles(
    parent_smiles: str,
    child_smiles: str,
) -> bool:
    return should_cull_metabolic_aromatic_hydroxylation(
        _safe_mol(parent_smiles),
        _safe_mol(child_smiles),
    )


def mol_total_radical_electrons(mol: Chem.Mol | None) -> int:
    if mol is None:
        return 0
    try:
        return sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    except Exception:
        return 0


def mol_has_radical_center(mol: Chem.Mol | None) -> bool:
    return mol_total_radical_electrons(mol) >= 1


def is_spin_stabilized_neutral_smiles(smiles: str) -> bool:
    return not mol_has_radical_center(_safe_mol(smiles))


def is_spin_stabilized_neutral_mol(mol: Chem.Mol | None) -> bool:
    return not mol_has_radical_center(_safe_mol_obj(mol))


@dataclass(frozen=True)
class ShelfLifeCullVerdict:
    eligible: bool
    reason: str = ""

    @property
    def culled(self) -> bool:
        return not self.eligible


class ShelfLifeChildFilter:
    """보관 유연물질 트리 / 리포트 — 라디칼 · 방향족 수산화 배제."""

    @staticmethod
    def evaluate_mols(
        parent_mol: Chem.Mol | None,
        child_mol: Chem.Mol | None,
        *,
        parent_smiles: str = "",
        child_smiles: str = "",
    ) -> ShelfLifeCullVerdict:
        parent_smi = (parent_smiles or "").strip()
        child_smi = (child_smiles or "").strip()

        if not child_smi and child_mol is None:
            return ShelfLifeCullVerdict(False, "empty_child")
        if parent_smi and child_smi and parent_smi == child_smi:
            return ShelfLifeCullVerdict(False, "identical_to_parent")

        if child_mol is not None and not is_spin_stabilized_neutral_mol(child_mol):
            return ShelfLifeCullVerdict(False, "radical_intermediate")
        if child_smi and not is_spin_stabilized_neutral_smiles(child_smi):
            return ShelfLifeCullVerdict(False, "radical_intermediate")

        pm = parent_mol if parent_mol is not None else _safe_mol(parent_smi)
        cm = child_mol if child_mol is not None else _safe_mol(child_smi)
        try:
            if should_cull_metabolic_aromatic_hydroxylation(pm, cm):
                return ShelfLifeCullVerdict(False, "aromatic_r6_hydroxylation")
        except Exception:
            pass

        return ShelfLifeCullVerdict(True, "shelf_life_eligible")

    @staticmethod
    def evaluate(parent_smiles: str, child_smiles: str) -> ShelfLifeCullVerdict:
        return ShelfLifeChildFilter.evaluate_mols(
            _safe_mol(parent_smiles),
            _safe_mol(child_smiles),
            parent_smiles=parent_smiles,
            child_smiles=child_smiles,
        )

    @staticmethod
    def passes(parent_smiles: str, child_smiles: str) -> bool:
        return ShelfLifeChildFilter.evaluate(parent_smiles, child_smiles).eligible

    @staticmethod
    def passes_mols(parent_mol: Chem.Mol | None, child_mol: Chem.Mol | None) -> bool:
        return ShelfLifeChildFilter.evaluate_mols(parent_mol, child_mol).eligible


is_shelf_life_eligible_child = ShelfLifeChildFilter.passes


if __name__ == "__main__":
    from rdkit import Chem as C

    aromatic_oh = Chem.MolFromSmarts("[OX2H]-[c;r6]")
    cases = [
        ("c1ccccc1", "Oc1ccccc1", True),
        ("CC(=O)Oc1ccccc1", "Oc1ccccc1", True),
        ("CC(=O)Oc1ccc(C=N)c(O)c1", "CC(=O)Oc1ccc(C=N)c(O)c1", False),
    ]
    for p_smi, c_smi, expect_cull in cases:
        pm = C.MolFromSmiles(p_smi)
        cm = C.MolFromSmiles(c_smi)
        parent_has = pm.HasSubstructMatch(aromatic_oh)
        child_has = cm.HasSubstructMatch(aromatic_oh)
        got = child_has and not parent_has
        print(
            f"[{'OK' if got == expect_cull else 'FAIL'}] "
            f"parent_oh={parent_has} child_oh={child_has} cull={got}"
        )
