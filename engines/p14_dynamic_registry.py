"""
PRIME-DS Dynamic Chemical Registry — PubChem PUG-REST + RDKit SMARTS 자동 분류.

마스터 DB 42종 외 임의 시약/용매를 이름·CAS·SMILES로 실시간 등록한다.
"""

# 목록에 없는 시약 (목록 42종에 갇히지 않게 해주는 장치)
# ① PubChem에 물어본다
# ② 구조를 받아온다
# ③ RDKit SMARTS로 "이건 산인가 염기인가 촉매인가" 자동 분류
# ④ 해당 목록에 추가

from __future__ import annotations

import re
from typing import Any
from urllib.parse import quote

import requests
from rdkit import Chem

_CATALYST_METALS: tuple[str, ...] = ("Pd", "Ni", "Cu", "Pt", "Ru", "Rh")


class DynamicChemicalRegistry:
    """PubChem 조회 및 RDKit 하위구조 기반 카테고리 자동 추론(Tagger)."""

    def __init__(self) -> None:
        self.category_patterns: dict[str, list[Chem.Mol | None]] = {
            "Acids": [
                Chem.MolFromSmarts("C(=O)[OH]"),
                Chem.MolFromSmarts("S(=O)(=O)[OH]"),
                Chem.MolFromSmarts("S(=O)(=O)Cl"),
                Chem.MolFromSmarts("[Cl,Br,I,F]"),
            ],
            "Bases": [
                Chem.MolFromSmarts("[NX3;H2,H1,H0]"),
                Chem.MolFromSmarts("[OH-]"),
                Chem.MolFromSmarts("[O-]"),
            ],
            "Catalysts": [
                Chem.MolFromSmarts("[Pd,Ni,Cu,Pt,Ru,Rh]"),
            ],
        }

    def search_pubchem_structure(self, identifier: str) -> str | None:
        """PubChem PUG-REST — 이름/CAS/상용명 → Canonical SMILES."""
        ident = str(identifier or "").strip()
        if not ident:
            return None

        encoded = quote(ident, safe="")
        endpoints = [
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/CanonicalSMILES/JSON",
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/SMILES/JSON",
        ]
        if re.fullmatch(r"\d{2,7}-\d{2}-\d", ident):
            cid_url = (
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                f"{encoded}/property/CanonicalSMILES/JSON"
            )
            endpoints.insert(0, cid_url)

        for url in endpoints:
            prop = "CanonicalSMILES" if "CanonicalSMILES" in url else "SMILES"
            try:
                response = requests.get(url, timeout=8)
                if response.status_code != 200:
                    continue
                data = response.json()
                props = data.get("PropertyTable", {}).get("Properties") or []
                if not props:
                    continue
                raw_smi = props[0].get(prop) or props[0].get("CanonicalSMILES") or props[0].get("SMILES")
                if not raw_smi:
                    continue
                mol = Chem.MolFromSmiles(str(raw_smi).strip())
                if mol is not None:
                    return Chem.MolToSmiles(mol, canonical=True)
            except Exception:
                continue
        return None

    def _infer_catalyst_metal(self, mol: Chem.Mol) -> str | None:
        for metal in _CATALYST_METALS:
            pattern = Chem.MolFromSmarts(f"[{metal}]")
            if pattern and mol.HasSubstructMatch(pattern):
                return metal
        return None

    def auto_classify_category(self, smiles: str) -> str:
        """RDKit 하위구조 매칭 → Solvents / Acids / Bases / Catalysts."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Solvents"

        for pattern in self.category_patterns["Catalysts"]:
            if pattern and mol.HasSubstructMatch(pattern):
                return "Catalysts"

        for pattern in self.category_patterns["Acids"]:
            if pattern and mol.HasSubstructMatch(pattern):
                return "Acids"

        for pattern in self.category_patterns["Bases"]:
            if pattern and mol.HasSubstructMatch(pattern):
                return "Bases"

        return "Solvents"

    def get_unlimited_entry(self, user_input: str) -> dict[str, Any]:
        """사용자 입력(이름·CAS·SMILES) → Registry 호환 entry dict."""
        raw = str(user_input or "").strip()
        if not raw:
            return {"success": False, "message": "입력값이 비어 있습니다."}

        smiles = self.search_pubchem_structure(raw)
        if not smiles:
            test_mol = Chem.MolFromSmiles(raw)
            if test_mol is not None:
                smiles = Chem.MolToSmiles(test_mol, canonical=True)
            else:
                return {"success": False, "message": "물질을 찾을 수 없습니다."}

        mol = Chem.MolFromSmiles(smiles)
        category = self.auto_classify_category(smiles)
        safe_id = "".join(ch if ch.isalnum() else "_" for ch in raw.lower()).strip("_")[:48]

        entry: dict[str, Any] = {
            "success": True,
            "id": f"dyn_{safe_id or 'entry'}",
            "name": raw,
            "smiles": smiles,
            "category": category,
            "role": category.rstrip("s") if category.endswith("s") else category,
            "source": "dynamic_pubchem",
        }
        if mol is not None and category == "Catalysts":
            metal = self._infer_catalyst_metal(mol)
            if metal:
                entry["metal"] = metal
        return entry


_DYNAMIC_REGISTRY: DynamicChemicalRegistry | None = None


def get_dynamic_registry() -> DynamicChemicalRegistry:
    global _DYNAMIC_REGISTRY
    if _DYNAMIC_REGISTRY is None:
        _DYNAMIC_REGISTRY = DynamicChemicalRegistry()
    return _DYNAMIC_REGISTRY
