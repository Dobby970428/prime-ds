"""
PRIME-DS Commercial Chemical Registry — 용매 · 산 · 염기 · 촉매 마스터 DB.

정보화학 표준 속성: id, name, smiles, role (category), risk / strength / metal / class.
"""

# 용매 13종   DCM, DMF, MeOH, THF, EtOAc, DMSO, 톨루엔...
# 산   8종    HCl, H2SO4, TFA, AcOH...
# 염기 11종   TEA, DIPEA, 피리딘, NaOH, K2CO3, KOtBu...
# 촉매 10종   Pd(PPh3)4, Pd(dppf)Cl2, Pd/C, CuI, DMAP, AIBN...

from __future__ import annotations

from typing import Any

REGISTRY_CATEGORIES: tuple[str, ...] = ("Solvents", "Acids", "Bases", "Catalysts")
REGISTRY_NONE: str = "선택안함 (Not Selected)"

# 레거시 UI 세션 id → 현행 id (Phase 5~8 호환)
LEGACY_ID_ALIASES: dict[str, str] = {
    "sol_mecn": "sol_acn",
    "sol_h2o": "sol_water",
    "cat_pd_oac2": "cat_pd_ac",
    "cat_cui": "cat_cu_i",
    "bas_pip": "bas_piper",
}

# 제약 합성 공정에 특화된 카테고리별 마스터 DB (id 기반 · 확장 가능)
CHEMICAL_REGISTRY_DB: dict[str, list[dict[str, Any]]] = {
    "Solvents": [
        {"id": "sol_dcm", "name": "Dichloromethane (DCM)", "smiles": "ClCCl", "class": "Halogenated", "risk": "High", "role": "Solvent"},
        {"id": "sol_dmf", "name": "Dimethylformamide (DMF)", "smiles": "CN(C)C=O", "class": "Polar Aprotic", "risk": "High", "role": "Solvent"},
        {"id": "sol_meoh", "name": "Methanol (MeOH)", "smiles": "CO", "class": "Polar Protic", "risk": "Medium", "role": "Solvent"},
        {"id": "sol_etoh", "name": "Ethanol (EtOH)", "smiles": "CCO", "class": "Polar Protic", "risk": "Medium", "role": "Solvent"},
        {"id": "sol_ipa", "name": "Isopropanol (IPA)", "smiles": "CC(C)O", "class": "Polar Protic", "risk": "Medium", "role": "Solvent"},
        {"id": "sol_thf", "name": "Tetrahydrofuran (THF)", "smiles": "C1CCOC1", "class": "Ether", "risk": "Medium", "role": "Solvent"},
        {"id": "sol_ea", "name": "Ethyl Acetate (EtOAc)", "smiles": "CCOC(=O)C", "class": "Ester", "risk": "Low", "role": "Solvent"},
        {"id": "sol_ace", "name": "Acetone", "smiles": "CC(=O)C", "class": "Ketone", "risk": "Medium", "role": "Solvent"},
        {"id": "sol_acn", "name": "Acetonitrile (MeCN)", "smiles": "CC#N", "class": "Polar Aprotic", "risk": "Medium", "role": "Solvent"},
        {"id": "sol_dmso", "name": "Dimethyl Sulfoxide (DMSO)", "smiles": "CS(=O)C", "class": "Polar Aprotic", "risk": "Medium", "role": "Solvent"},
        {"id": "sol_tol", "name": "Toluene", "smiles": "Cc1ccccc1", "class": "Non-polar Aromatic", "risk": "Medium", "role": "Solvent"},
        {"id": "sol_hex", "name": "n-Hexane", "smiles": "CCCCCC", "class": "Non-polar Aliphatic", "risk": "Low", "role": "Solvent"},
        {"id": "sol_water", "name": "Water (H2O)", "smiles": "O", "class": "Polar Protic", "risk": "Low", "role": "Solvent"},
    ],
    "Acids": [
        {"id": "acd_hcl", "name": "Hydrochloric Acid (HCl)", "smiles": "Cl", "strength": "Strong", "type": "Inorganic", "risk": "High", "role": "Acid"},
        {"id": "acd_h2so4", "name": "Sulfuric Acid (H2SO4)", "smiles": "OS(=O)(=O)O", "strength": "Strong", "type": "Inorganic", "risk": "High", "role": "Acid"},
        {"id": "acd_hno3", "name": "Nitric Acid (HNO3)", "smiles": "O=N(=O)O", "strength": "Strong", "type": "Inorganic / Oxidizing", "risk": "High", "role": "Acid"},
        {"id": "acd_tfa", "name": "Trifluoroacetic acid (TFA)", "smiles": "FC(F)(F)C(=O)O", "strength": "Strong", "type": "Organic", "risk": "Medium", "role": "Acid"},
        {"id": "acd_tsoh", "name": "p-Toluenesulfonic acid (TsOH)", "smiles": "Cc1ccc(S(=O)(=O)O)cc1", "strength": "Strong", "type": "Organic", "risk": "Medium", "role": "Acid"},
        {"id": "acd_aa", "name": "Acetic Acid (AcOH)", "smiles": "CC(=O)O", "strength": "Weak", "type": "Organic", "risk": "Low", "role": "Acid"},
        {"id": "acd_fa", "name": "Formic Acid", "smiles": "C(=O)O", "strength": "Weak", "type": "Organic", "risk": "Low", "role": "Acid"},
        {"id": "acd_h3po4", "name": "Phosphoric Acid (H3PO4)", "smiles": "OP(=O)(O)O", "strength": "Medium", "type": "Inorganic", "risk": "Medium", "role": "Acid"},
    ],
    "Bases": [
        {"id": "bas_tea", "name": "Triethylamine (TEA)", "smiles": "CCN(CC)CC", "strength": "Weak", "type": "Organic Tertiary Amine", "risk": "Low", "role": "Base"},
        {"id": "bas_diea", "name": "N,N-Diisopropylethylamine (DIPEA/Hunig's Base)", "smiles": "CCN(C(C)C)C(C)C", "strength": "Weak", "type": "Organic Tertiary Amine", "risk": "Low", "role": "Base"},
        {"id": "bas_py", "name": "Pyridine", "smiles": "c1ccncc1", "strength": "Weak", "type": "Organic Aromatic", "risk": "Low", "role": "Base"},
        {"id": "bas_piper", "name": "Piperidine", "smiles": "C1CCNCC1", "strength": "Weak", "type": "Organic Secondary Amine", "risk": "Low", "role": "Base"},
        {"id": "bas_naoh", "name": "Sodium Hydroxide (NaOH)", "smiles": "[Na+].[OH-]", "strength": "Strong", "type": "Inorganic Strong", "risk": "High", "role": "Base"},
        {"id": "bas_koh", "name": "Potassium Hydroxide (KOH)", "smiles": "[K+].[OH-]", "strength": "Strong", "type": "Inorganic Strong", "risk": "High", "role": "Base"},
        {"id": "bas_k2co3", "name": "Potassium Carbonate (K2CO3)", "smiles": "O=C([O-])[O-].[K+].[K+]", "strength": "Moderate", "type": "Inorganic Mild", "risk": "Low", "role": "Base"},
        {"id": "bas_na2co3", "name": "Sodium Carbonate (Na2CO3)", "smiles": "O=C([O-])[O-].[Na+].[Na+]", "strength": "Moderate", "type": "Inorganic Mild", "risk": "Low", "role": "Base"},
        {"id": "bas_nahco3", "name": "Sodium Bicarbonate (NaHCO3)", "smiles": "C(=O)(O)[O-].[Na+]", "strength": "Weak", "type": "Inorganic Weak", "risk": "Low", "role": "Base"},
        {"id": "bas_naome", "name": "Sodium Methoxide (NaOMe)", "smiles": "C[O-].[Na+]", "strength": "Strong", "type": "Alkoxide Strong", "risk": "High", "role": "Base"},
        {"id": "bas_ktbub", "name": "Potassium tert-butoxide (KOtBu)", "smiles": "CC(C)(C)[O-].[K+]", "strength": "Strong", "type": "Alkoxide Strong", "risk": "High", "role": "Base"},
    ],
    "Catalysts": [
        {"id": "cat_pd_pph3", "name": "Pd(PPh3)4 (Tetrakis)", "smiles": "[Pd].P(c1ccccc1)(c2ccccc2)c3ccccc3", "class": "Pd(0) Cross-Coupling", "metal": "Pd", "type": "Cross-Coupling", "risk": "Medium", "role": "Catalyst"},
        {"id": "cat_pd_dcl", "name": "Pd(dppf)Cl2", "smiles": "Cl[Pd]Cl", "class": "Pd(II) Cross-Coupling", "metal": "Pd", "type": "Cross-Coupling", "risk": "Medium", "role": "Catalyst"},
        {"id": "cat_pd_ac", "name": "Palladium(II) Acetate (Pd(OAc)2)", "smiles": "CC(=O)[O-].CC(=O)[O-].[Pd+2]", "class": "Pd(II) Precatalyst", "metal": "Pd", "type": "Cross-Coupling", "risk": "Medium", "role": "Catalyst"},
        {"id": "cat_pd_c", "name": "Palladium on Carbon (Pd/C)", "smiles": "[Pd]", "class": "Heterogeneous Hydrogenation", "metal": "Pd", "type": "Cross-Coupling", "risk": "Medium", "role": "Catalyst"},
        {"id": "cat_ni_cl", "name": "Ni(dppp)Cl2", "smiles": "Cl[Ni]Cl", "class": "Nickel Cross-Coupling", "metal": "Ni", "type": "Cross-Coupling", "risk": "Medium", "role": "Catalyst"},
        {"id": "cat_nicl2", "name": "Nickel(II) Chloride (NiCl2)", "smiles": "Cl[Ni]Cl", "class": "Nickel Halide", "metal": "Ni", "type": "Cross-Coupling", "risk": "Medium", "role": "Catalyst"},
        {"id": "cat_cu_i", "name": "Copper(I) Iodide (CuI)", "smiles": "[Cu]I", "class": "Click Chemistry / Sonogashira", "metal": "Cu", "type": "Cross-Coupling", "risk": "Medium", "role": "Catalyst"},
        {"id": "cat_dmap", "name": "4-Dimethylaminopyridine (DMAP)", "smiles": "CN(C)c1ccncc1", "class": "Nucleophilic Organocatalyst", "metal": "None", "type": "Organocatalyst", "risk": "Low", "role": "Catalyst"},
        {"id": "cat_aibn", "name": "AIBN", "smiles": "CC(C)(C#N)/N=N/C(C)(C)C#N", "class": "Radical Initiator", "metal": "None", "type": "Radical", "risk": "Medium", "role": "Catalyst"},
        {"id": "cat_bzox", "name": "Benzoyl Peroxide (BPO)", "smiles": "O=C(OOC(=O)c1ccccc1)c2ccccc2", "class": "Radical Initiator", "metal": "None", "type": "Radical", "risk": "Medium", "role": "Catalyst"},
    ],
}

REGISTRY_MASTER_COUNT: int = sum(
    len(entries) for entries in CHEMICAL_REGISTRY_DB.values()
)

REAGENT_NAME_ALIASES: dict[str, str] = {
    "NaOH": "Sodium Hydroxide (NaOH)",
    "TEA": "Triethylamine (TEA)",
    "HCl": "Hydrochloric Acid (HCl)",
    "Hydrochloric Acid": "Hydrochloric Acid (HCl)",
    "Sodium Hydroxide": "Sodium Hydroxide (NaOH)",
    "DIPEA (Hünig's Base)": "N,N-Diisopropylethylamine (DIPEA/Hunig's Base)",
    "N,N-Diisopropylethylamine (DIPEA)": "N,N-Diisopropylethylamine (DIPEA/Hunig's Base)",
    "Potassium Carbonate": "Potassium Carbonate (K2CO3)",
    "Acetic Acid": "Acetic Acid (AcOH)",
    "Water": "Water (H2O)",
    "Pd(OAc)2": "Palladium(II) Acetate (Pd(OAc)2)",
    "Pd(PPh3)4": "Pd(PPh3)4 (Tetrakis)",
    "CuI": "Copper(I) Iodide (CuI)",
    "Piperidine": "Piperidine",
    "None": REGISTRY_NONE,
}

SOLVENT_NAME_ALIASES: dict[str, str] = {
    "Methanol": "Methanol (MeOH)",
    "Acetonitrile (MeCN)": "Acetonitrile (MeCN)",
    "Water": "Water (H2O)",
}

_HOMO_COUPLING_METALS: frozenset[str] = frozenset({"Pd", "Ni", "Cu"})

# 용매 물리화학 확장 (Registry Solvents와 name 매칭)
SOLVENT_PHYSICS: dict[str, dict[str, Any]] = {
    "Tetrahydrofuran (THF)": {
        "abbrev": "THF",
        "formula": "C4H8O",
        "bp_c": 66.0,
        "density_g_ml": 0.889,
        "dielectric": 7.6,
        "polarity": "Medium",
        "water_miscible": False,
    },
    "Methanol (MeOH)": {
        "abbrev": "MeOH",
        "formula": "CH4O",
        "bp_c": 64.7,
        "density_g_ml": 0.792,
        "dielectric": 32.7,
        "polarity": "High",
        "water_miscible": True,
    },
    "Acetonitrile (MeCN)": {
        "abbrev": "MeCN",
        "formula": "C2H3N",
        "bp_c": 81.6,
        "density_g_ml": 0.786,
        "dielectric": 37.5,
        "polarity": "High",
        "water_miscible": True,
    },
    "Dichloromethane (DCM)": {
        "abbrev": "DCM",
        "formula": "CH2Cl2",
        "bp_c": 39.6,
        "density_g_ml": 1.325,
        "dielectric": 8.9,
        "polarity": "Low",
        "water_miscible": False,
    },
    "Dimethylformamide (DMF)": {
        "abbrev": "DMF",
        "formula": "C3H7NO",
        "bp_c": 153.0,
        "density_g_ml": 0.944,
        "dielectric": 36.7,
        "polarity": "High",
        "water_miscible": True,
    },
    "Toluene": {
        "abbrev": "PhMe",
        "formula": "C7H8",
        "bp_c": 110.6,
        "density_g_ml": 0.867,
        "dielectric": 2.4,
        "polarity": "Low",
        "water_miscible": False,
    },
    "Water (H2O)": {
        "abbrev": "H2O",
        "formula": "H2O",
        "bp_c": 100.0,
        "density_g_ml": 1.000,
        "dielectric": 78.5,
        "polarity": "Very High",
        "water_miscible": True,
    },
    "Ethanol (EtOH)": {
        "abbrev": "EtOH",
        "formula": "C2H6O",
        "bp_c": 78.4,
        "density_g_ml": 0.789,
        "dielectric": 24.5,
        "polarity": "High",
        "water_miscible": True,
    },
    "Isopropanol (IPA)": {
        "abbrev": "IPA",
        "formula": "C3H8O",
        "bp_c": 82.6,
        "density_g_ml": 0.786,
        "dielectric": 18.3,
        "polarity": "High",
        "water_miscible": True,
    },
    "Ethyl Acetate (EtOAc)": {
        "abbrev": "EtOAc",
        "formula": "C4H8O2",
        "bp_c": 77.1,
        "density_g_ml": 0.902,
        "dielectric": 6.0,
        "polarity": "Medium",
        "water_miscible": False,
    },
    "Acetone": {
        "abbrev": "ACE",
        "formula": "C3H6O",
        "bp_c": 56.0,
        "density_g_ml": 0.791,
        "dielectric": 20.7,
        "polarity": "High",
        "water_miscible": True,
    },
    "Dimethyl Sulfoxide (DMSO)": {
        "abbrev": "DMSO",
        "formula": "C2H6OS",
        "bp_c": 189.0,
        "density_g_ml": 1.100,
        "dielectric": 46.7,
        "polarity": "High",
        "water_miscible": True,
    },
    "n-Hexane": {
        "abbrev": "nHex",
        "formula": "C6H14",
        "bp_c": 68.7,
        "density_g_ml": 0.659,
        "dielectric": 1.9,
        "polarity": "Low",
        "water_miscible": False,
    },
}


def _build_master_lookup() -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    by_name: dict[str, dict[str, Any]] = {}
    by_id: dict[str, dict[str, Any]] = {}
    for category in REGISTRY_CATEGORIES:
        for entry in CHEMICAL_REGISTRY_DB.get(category, []):
            payload = dict(entry)
            payload["category"] = category
            by_name[str(entry["name"])] = payload
            if entry.get("id"):
                by_id[str(entry["id"])] = payload
    for alias, canonical in REAGENT_NAME_ALIASES.items():
        if canonical in by_name:
            by_name[alias] = dict(by_name[canonical])
    for alias, canonical in SOLVENT_NAME_ALIASES.items():
        if canonical in by_name:
            by_name[alias] = dict(by_name[canonical])
    for legacy_id, current_id in LEGACY_ID_ALIASES.items():
        if current_id in by_id:
            by_id[legacy_id] = dict(by_id[current_id])
    return by_name, by_id


_MASTER_LOOKUP, _ID_LOOKUP = _build_master_lookup()


def get_chemical_by_id(category: str, chem_id: str) -> dict[str, Any]:
    """카테고리와 ID를 기반으로 단일 시약 정보(SMILES 포함)를 반환."""
    resolved_id = LEGACY_ID_ALIASES.get(chem_id, chem_id)
    for chem in CHEMICAL_REGISTRY_DB.get(category, []):
        if chem.get("id") == resolved_id:
            payload = dict(chem)
            payload["category"] = category
            return payload
    return {}


def get_chemical_by_id_any(chem_id: str) -> dict[str, Any]:
    """카테고리 무관 ID 조회."""
    entry = _ID_LOOKUP.get(chem_id)
    return dict(entry) if entry else {}


def get_category_ids(category: str) -> list[str]:
    """카테고리별 registry id 목록 (표시명 알파벳/가나다 순 정렬)."""
    if category not in REGISTRY_CATEGORIES:
        return []
    ids = [str(e["id"]) for e in CHEMICAL_REGISTRY_DB.get(category, []) if e.get("id")]
    ids.sort(key=lambda cid: format_registry_display_name(category, cid).casefold())
    return ids


def is_dynamic_registry_id(chem_id: str) -> bool:
    return str(chem_id or "").startswith("dyn_")


def get_dynamic_catalog_for_category(
    category: str,
    dynamic_catalog: list[dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    return [
        dict(e)
        for e in (dynamic_catalog or [])
        if str(e.get("category")) == category and e.get("id")
    ]


def get_merged_category_ids(
    category: str,
    dynamic_catalog: list[dict[str, Any]] | None = None,
) -> list[str]:
    """마스터 DB id + 동적 등록 id (표시명 정렬)."""
    master = get_category_ids(category)
    dynamic_ids = [
        str(e["id"])
        for e in get_dynamic_catalog_for_category(category, dynamic_catalog)
        if e.get("id") not in master
    ]
    dynamic_ids.sort(
        key=lambda cid: format_registry_display_name(
            category, cid, dynamic_catalog=dynamic_catalog
        ).casefold()
    )
    return master + dynamic_ids


def resolve_registry_token(
    category: str,
    token: str,
    dynamic_catalog: list[dict[str, Any]] | None = None,
) -> dict[str, Any] | None:
    """마스터 id/명칭 또는 동적 id → registry entry dict."""
    cleaned = str(token or "").strip()
    if not cleaned or is_registry_none(cleaned):
        return None

    entry = get_chemical_by_id(category, cleaned)
    if entry:
        return entry

    entry = lookup_registry_entry(cleaned, category=category)
    if entry:
        return dict(entry)

    for dyn in dynamic_catalog or []:
        if str(dyn.get("id")) != cleaned:
            continue
        if str(dyn.get("category")) != category:
            continue
        return {
            "id": dyn.get("id"),
            "name": dyn.get("name", cleaned),
            "smiles": dyn.get("smiles"),
            "category": category,
            "role": dyn.get("role", category),
            "metal": dyn.get("metal"),
            "source": dyn.get("source", "dynamic_pubchem"),
        }
    return None


def format_registry_display_name(
    category: str,
    chem_id: str,
    *,
    dynamic_catalog: list[dict[str, Any]] | None = None,
) -> str:
    """multiselect format_func — id → 표시명."""
    if is_dynamic_registry_id(chem_id):
        for dyn in dynamic_catalog or []:
            if str(dyn.get("id")) == chem_id:
                tag = str(dyn.get("category", category))[:3]
                return f"🌐 {dyn.get('name', chem_id)} [{tag}]"
        return chem_id

    entry = get_chemical_by_id(category, chem_id)
    return str(entry.get("name") or chem_id) if entry else chem_id


def get_smiles_for_ids(category: str, chem_ids: list[str]) -> list[str]:
    """선택 id 목록 → 해당 카테고리 SMILES 배열 (선택 순서 유지)."""
    smiles: list[str] = []
    for chem_id in chem_ids:
        entry = get_chemical_by_id(category, chem_id)
        smi = entry.get("smiles")
        if smi:
            smiles.append(str(smi))
    return smiles


def build_smiles_array_from_id_selection(
    selected: dict[str, list[str]],
) -> list[str]:
    """4카테고리 id 선택 → run_process_simulation용 통합 SMILES 배열."""
    return list(build_process_reagent_context_from_ids(selected).get("all_smiles") or [])


def is_registry_none(name: str | None) -> bool:
    cleaned = (name or "").strip()
    return not cleaned or cleaned in {REGISTRY_NONE, "None", "선택안함"}


def normalize_reagent_name(reagent: str | None) -> str:
    cleaned = (reagent or "").strip()
    if is_registry_none(cleaned):
        return REGISTRY_NONE
    return REAGENT_NAME_ALIASES.get(cleaned, cleaned)


def get_category_options(category: str) -> list[str]:
    """카테고리별 Registry 명칭 목록 (선택안함 제외)."""
    if category not in REGISTRY_CATEGORIES:
        return []
    return [str(e["name"]) for e in CHEMICAL_REGISTRY_DB.get(category, [])]


def lookup_registry_entry(name: str | None, *, category: str | None = None) -> dict[str, Any] | None:
    cleaned = normalize_reagent_name(name)
    if is_registry_none(cleaned):
        return None
    entry = _MASTER_LOOKUP.get(cleaned)
    if entry is None:
        return None
    if category and entry.get("category") != category:
        return None
    return dict(entry)


def lookup_reagent(reagent: str | None) -> dict[str, Any] | None:
    """레거시 호환 — Acids/Bases/Catalysts 조회."""
    entry = lookup_registry_entry(reagent)
    if entry and entry.get("category") in ("Acids", "Bases", "Catalysts"):
        return entry
    return entry


def get_reagent_smiles(reagent: str | None) -> str | None:
    entry = lookup_registry_entry(reagent)
    if not entry:
        return None
    smiles = entry.get("smiles")
    return str(smiles) if smiles else None


def get_reagent_ui_options() -> list[str]:
    """레거시 단일 selectbox 옵션."""
    options = [REGISTRY_NONE]
    for category in ("Acids", "Bases", "Catalysts"):
        options.extend(get_category_options(category))
    return options


def is_catalyst_reagent_name(reagent: str | None) -> bool:
    entry = lookup_registry_entry(reagent, category="Catalysts")
    if not entry:
        return False
    metal = str(entry.get("metal") or "").strip()
    return metal in _HOMO_COUPLING_METALS


def catalyst_metal_label(reagent: str | None) -> str:
    entry = lookup_registry_entry(reagent, category="Catalysts")
    if entry:
        metal = str(entry.get("metal") or "").strip()
        if metal in _HOMO_COUPLING_METALS:
            return metal
    return "Catalyst"


def build_solvent_catalog() -> dict[str, dict[str, Any]]:
    catalog: dict[str, dict[str, Any]] = {}
    for entry in CHEMICAL_REGISTRY_DB.get("Solvents", []):
        name = str(entry["name"])
        physics = dict(SOLVENT_PHYSICS.get(name, {}))
        catalog[name] = {
            "abbrev": physics.get("abbrev") or name[:6],
            "formula": physics.get("formula", "—"),
            "bp_c": physics.get("bp_c"),
            "density_g_ml": physics.get("density_g_ml"),
            "dielectric": physics.get("dielectric"),
            "polarity": physics.get("polarity") or entry.get("class", "—"),
            "water_miscible": physics.get("water_miscible", False),
            "smiles": entry.get("smiles"),
            "class": entry.get("class"),
            "risk": entry.get("risk"),
            "role": entry.get("role", "Solvent"),
        }
    return catalog


SOLVENT_CATALOG: dict[str, dict[str, Any]] = build_solvent_catalog()


def _resolve_category_selection(
    category: str,
    selected_tokens: list[str] | None,
    *,
    dynamic_catalog: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    resolved: list[dict[str, Any]] = []
    for raw in selected_tokens or []:
        entry = resolve_registry_token(category, str(raw), dynamic_catalog)
        if entry:
            resolved.append(entry)
    return resolved


def build_process_reagent_context_from_ids(
    selected: dict[str, list[str]],
    *,
    dynamic_catalog: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """
    id 기반 UI 선택 → Request Body.

    Keys: solvents, acids, bases, catalysts (각 id list).
    """
    return build_process_reagent_context(
        {
            "Solvents": list(selected.get("solvents") or []),
            "Acids": list(selected.get("acids") or []),
            "Bases": list(selected.get("bases") or []),
            "Catalysts": list(selected.get("catalysts") or []),
        },
        dynamic_catalog=dynamic_catalog,
    )


def build_process_reagent_context(
    selected: dict[str, list[str] | str | None],
    *,
    dynamic_catalog: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """
    4-카테고리 UI 선택 → 백엔드 Request Body.

    Parameters
    ----------
    selected : keys Solvents, Acids, Bases, Catalysts
        Solvents는 str(단일) 또는 list 허용.
    """
    solvents_raw = selected.get("Solvents")
    if isinstance(solvents_raw, str):
        solvent_names = [] if is_registry_none(solvents_raw) else [solvents_raw]
    else:
        solvent_names = list(solvents_raw or [])

    acid_names = list(selected.get("Acids") or [])
    base_names = list(selected.get("Bases") or [])
    catalyst_names = list(selected.get("Catalysts") or [])

    solvent_entries = _resolve_category_selection("Solvents", solvent_names, dynamic_catalog=dynamic_catalog)
    acid_entries = _resolve_category_selection("Acids", acid_names, dynamic_catalog=dynamic_catalog)
    base_entries = _resolve_category_selection("Bases", base_names, dynamic_catalog=dynamic_catalog)
    catalyst_entries = _resolve_category_selection("Catalysts", catalyst_names, dynamic_catalog=dynamic_catalog)

    registry_entries = solvent_entries + acid_entries + base_entries + catalyst_entries
    all_smiles = [e["smiles"] for e in registry_entries if e.get("smiles")]

    homo_catalysts = [
        e["name"]
        for e in catalyst_entries
        if is_catalyst_reagent_name(e["name"])
        or str(e.get("metal") or "").strip() in _HOMO_COUPLING_METALS
    ]
    primary_catalyst = homo_catalysts[0] if homo_catalysts else (
        catalyst_entries[0]["name"]
        if catalyst_entries
        else REGISTRY_NONE
    )

    primary_solvent = (
        solvent_entries[0]["name"]
        if solvent_entries
        else REGISTRY_NONE
    )
    solvent_props = SOLVENT_CATALOG.get(primary_solvent, {}) if solvent_entries else {}

    selected_materials: dict[str, list[str]] = {
        "Solvents": [e["smiles"] for e in solvent_entries if e.get("smiles")],
        "Acids": [e["smiles"] for e in acid_entries if e.get("smiles")],
        "Bases": [e["smiles"] for e in base_entries if e.get("smiles")],
        "Catalysts": [e["smiles"] for e in catalyst_entries if e.get("smiles")],
    }

    return {
        "solvents": [e["name"] for e in solvent_entries],
        "acids": [e["name"] for e in acid_entries],
        "bases": [e["name"] for e in base_entries],
        "catalysts": [e["name"] for e in catalyst_entries],
        "solvent_ids": [e.get("id") for e in solvent_entries if e.get("id")],
        "acid_ids": [e.get("id") for e in acid_entries if e.get("id")],
        "base_ids": [e.get("id") for e in base_entries if e.get("id")],
        "catalyst_ids": [e.get("id") for e in catalyst_entries if e.get("id")],
        "solvent_smiles": selected_materials["Solvents"],
        "acid_smiles": selected_materials["Acids"],
        "base_smiles": selected_materials["Bases"],
        "catalyst_smiles": selected_materials["Catalysts"],
        "selected_materials": selected_materials,
        "all_smiles": all_smiles,
        "registry_entries": registry_entries,
        "homo_coupling_catalysts": homo_catalysts,
        "primary_catalyst": primary_catalyst,
        "primary_solvent": primary_solvent,
        "solvent_props": solvent_props,
        "reagent": primary_catalyst if not is_registry_none(primary_catalyst) else REGISTRY_NONE,
        "dynamic_catalog": list(dynamic_catalog or []),
    }


def merge_dynamic_registry_entries(
    ctx: dict[str, Any],
    dynamic_entries: list[dict[str, Any]],
) -> dict[str, Any]:
    """PubChem 동적 등록 entry → reagent_context / selected_materials 병합."""
    if not dynamic_entries:
        return ctx

    merged = dict(ctx)
    materials = build_selected_materials_payload(merged)
    registry_entries = list(merged.get("registry_entries") or [])
    name_lists: dict[str, list[str]] = {
        "Solvents": list(merged.get("solvents") or []),
        "Acids": list(merged.get("acids") or []),
        "Bases": list(merged.get("bases") or []),
        "Catalysts": list(merged.get("catalysts") or []),
    }
    homo_catalysts = list(merged.get("homo_coupling_catalysts") or [])

    for entry in dynamic_entries:
        if not entry.get("success", True):
            continue
        smi = str(entry.get("smiles") or "").strip()
        category = str(entry.get("category") or "Solvents")
        if category not in REGISTRY_CATEGORIES or not smi:
            continue
        if smi in materials.get(category, []):
            continue

        materials.setdefault(category, []).append(smi)
        name = str(entry.get("name") or smi)
        name_lists.setdefault(category, []).append(name)
        registry_entries.append(
            {
                "id": entry.get("id"),
                "name": name,
                "smiles": smi,
                "category": category,
                "role": entry.get("role", category),
                "source": entry.get("source", "dynamic_pubchem"),
            }
        )
        if category == "Catalysts" and is_catalyst_reagent_name(name):
            homo_catalysts.append(name)

    merged["selected_materials"] = materials
    merged["registry_entries"] = registry_entries
    merged["solvent_smiles"] = materials["Solvents"]
    merged["acid_smiles"] = materials["Acids"]
    merged["base_smiles"] = materials["Bases"]
    merged["catalyst_smiles"] = materials["Catalysts"]
    merged["solvents"] = name_lists["Solvents"]
    merged["acids"] = name_lists["Acids"]
    merged["bases"] = name_lists["Bases"]
    merged["catalysts"] = name_lists["Catalysts"]
    merged["all_smiles"] = []
    for cat in REGISTRY_CATEGORIES:
        for smi in materials.get(cat, []):
            if smi and smi not in merged["all_smiles"]:
                merged["all_smiles"].append(smi)
    merged["dynamic_entries"] = list(dynamic_entries)
    merged["homo_coupling_catalysts"] = homo_catalysts
    if homo_catalysts:
        merged["primary_catalyst"] = homo_catalysts[0]
        merged["reagent"] = homo_catalysts[0]
    return merged


def build_selected_materials_payload(
    ctx: dict[str, Any] | None,
) -> dict[str, list[str]]:
    """
    Registry 컨텍스트 → API Request Body.

    {"Solvents": ["SMILES", ...], "Acids": [...], "Bases": [...], "Catalysts": [...]}
    """
    source = ctx or {}
    existing = source.get("selected_materials")
    if isinstance(existing, dict) and any(existing.values()):
        return {
            "Solvents": [str(s) for s in (existing.get("Solvents") or []) if s],
            "Acids": [str(s) for s in (existing.get("Acids") or []) if s],
            "Bases": [str(s) for s in (existing.get("Bases") or []) if s],
            "Catalysts": [str(s) for s in (existing.get("Catalysts") or []) if s],
        }
    return {
        "Solvents": [str(s) for s in (source.get("solvent_smiles") or []) if s],
        "Acids": [str(s) for s in (source.get("acid_smiles") or []) if s],
        "Bases": [str(s) for s in (source.get("base_smiles") or []) if s],
        "Catalysts": [str(s) for s in (source.get("catalyst_smiles") or []) if s],
    }


def format_reagent_context_summary(ctx: dict[str, Any]) -> str:
    """결과 캡션용 요약 문자열."""
    parts: list[str] = []
    if ctx.get("solvents"):
        parts.append(f"Solvent: {', '.join(ctx['solvents'])}")
    if ctx.get("acids"):
        parts.append(f"Acid: {', '.join(ctx['acids'])}")
    if ctx.get("bases"):
        parts.append(f"Base: {', '.join(ctx['bases'])}")
    if ctx.get("catalysts"):
        parts.append(f"Catalyst: {', '.join(ctx['catalysts'])}")
    return " · ".join(parts) if parts else "시약 미선택"
