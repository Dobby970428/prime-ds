"""
2-반응물 공정 예측 + 목적물 분해/미반응 스캔 엔진 (rdChemReactions).
"""

# 경로 A

# 2물질 조합 -> 4바구니 (출발물질 2개를 받아 네 가지 바구니로 나눔)
# 교차반응 생성물: 물질1+물질2가 반응해 생긴 것
# 시약 유래 부산물: 시약 때문에 생긴 것
# 미반응 출발물질: 반응 안 하고 남은 것
# 분해산물: 공정 중 깨진 것

# 미반응 출발 물질:
# 출발물질이 최종 목적물의 부분구조인가?
#   → 아니다  =  반응에 안 쓰였을 수 있다  =  잔류 가능

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import Descriptors, rdChemReactions, rdMolDescriptors

from p40_degradation_engine import EXTRA_DEGRADATION_SMIRKS
from p30_reaction_engine import (
    ImpurityClassification,
    ImpurityResult,
    ReagentProfile,
    _mol_from_smiles,
    canonical_smiles,
    classify_reagent,
    is_nitrite_source_smiles,
    predict_nitrosation_with_reagent,
    predict_reaction_impurities,
    predict_reagent_side_products,
)

# ── 2반응물 SMIRKS (출발물질1 + 출발물질2 → 부생성물) ──
BINARY_REACTIONS: tuple[tuple[str, str, str], ...] = (
    (
        "amide_formation",
        "카르복실산 + 2급 아민 → 아미드",
        "[#6:1](=O)[OH].[#7;H1;D2;!$([#7][#7])!$([#7]=*):2]([#6:3])[#6:4]>>[#6:1](=O)[#7:2]([#6:3])([#6:4])",
    ),
    (
        "reductive_amination_side",
        "알데하이드 + 2급 아민 → 환원아민화 부산물",
        "[#6:1](=O)[#6:2].[#7;H1;D2;!$([#7][#7])!$([#7]=*):3]([#6:4])[#6:5]>>[#7:3]([#6:4])([#6:5])[#6:2][#6:1]",
    ),
    (
        "nitrosation_cross",
        "2급 아민 + 아질산염 → N-니트로소amine",
        "[#7;H1;D2;!$([#7][#7])!$([#7]=*):1]([#6:2])[#6:3].[N+](=O)[O-]>>[#7;H0;D3:1]([#6:2])([#6:3])N=O",
    ),
)

# ── 목적물 분해(Degradant) SMIRKS — 보관 안정성 (가수분해) ──
# 2026-09-03 수정: 가수분해 산물을 카복실레이트 음이온([O-])으로 만들면 중성형으로 저장된
#   PubChem DB canonical SMILES와 절대 일치하지 않아, 실제 검출 물질도 항상 "이론적"으로
#   분류됩니다. 네 모듈(p30·p32·p40·p41)의 가수분해 산물을 모두 중성형으로 통일했습니다.
DEGRADATION_RULES: tuple[tuple[str, str, str], ...] = (
    ("ester_hydrolysis", "에스터 가수분해", "[#6:1](=O)O[#6:2]>>[#6:1](=O)[OX2H].[#6:2][OH]"),
    ("amide_hydrolysis", "아미드 가수분해", "[#6:1](=O)[#7:2]>>[#6:1](=O)[OX2H].[#7:2]"),
) + EXTRA_DEGRADATION_SMIRKS


@dataclass
class ProcessAnalysisResult:
    cross_reaction_products: list[ImpurityResult] = field(default_factory=list)
    reagent_reaction_products: list[ImpurityResult] = field(default_factory=list)
    unreacted_starting_materials: list[ImpurityResult] = field(default_factory=list)
    degradants: list[ImpurityResult] = field(default_factory=list)
    reagent_profile: ReagentProfile | None = None


def _run_binary_reaction(mol_a: Chem.Mol, mol_b: Chem.Mol, smirks: str) -> list[Chem.Mol]:
    try:
        rxn = rdChemReactions.ReactionFromSmarts(smirks)
    except Exception:
        return []
    if rxn is None:
        return []
    products: list[Chem.Mol] = []
    for order in ((mol_a, mol_b), (mol_b, mol_a)):
        for outcome in rxn.RunReactants(order):
            for prod in outcome:
                try:
                    Chem.SanitizeMol(prod)
                    products.append(prod)
                except Exception:
                    continue
    return products


def _run_unimolecular(mol: Chem.Mol, smirks: str) -> list[Chem.Mol]:
    try:
        rxn = rdChemReactions.ReactionFromSmarts(smirks)
    except Exception:
        return []
    if rxn is None:
        return []
    products: list[Chem.Mol] = []
    for outcome in rxn.RunReactants((mol,)):
        for prod in outcome:
            try:
                Chem.SanitizeMol(prod)
                products.append(prod)
            except Exception:
                continue
    return products


def _is_unreacted_starting_material(reactant_smiles: str, product_smiles: str) -> bool:
    """출발물질이 목적물에 구조적으로 포함되지 않으면 미반응 잔류 가능."""
    r_mol = _mol_from_smiles(reactant_smiles)
    p_mol = _mol_from_smiles(product_smiles)
    if r_mol is None or p_mol is None:
        return False
    r_can = canonical_smiles(reactant_smiles)
    p_can = canonical_smiles(product_smiles)
    if r_can == p_can:
        return False
    return not p_mol.HasSubstructMatch(r_mol)


def _safe_product_mol(mol: Chem.Mol) -> Chem.Mol | None:
    """입체화학 정보 보존 시도 후 Sanitize."""
    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        # 2026-09-03 수정: 기존에는 SanitizeMol(catchErrors=True)로 실패 코드만 받고 망가진
        #   분자를 그대로 통과시켰습니다(방향족 고리가 조용히 포화 고리로 바뀐 날조 구조에
        #   분자식·MW·M7 등급이 붙음). sanitize 실패 분자는 불순물로 보고하지 않습니다.
        #   호출부 `_make_impurity`는 None이면 곧바로 None을 반환하므로 안전합니다.
        return None


def _make_impurity(
    mol: Chem.Mol,
    rule_id: str,
    rule_label: str,
    classification_fn,
    alert_fn,
    category: str,
) -> ImpurityResult | None:
    try:
        mol = _safe_product_mol(mol)
        if mol is None:
            return None
        smi = Chem.MolToSmiles(mol, canonical=True)
        if not smi:
            return None
        try:
            cls_str, db_rec = classification_fn(smi)
        except Exception:
            cls_str, db_rec = "theoretical", None
        try:
            alerts, hl = alert_fn(mol)
        except Exception:
            alerts, hl = [], []
        try:
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mw = round(Descriptors.MolWt(mol), 2)
        except Exception:
            formula, mw = "", 0.0
        return ImpurityResult(
            smiles=smi,
            classification=(
                ImpurityClassification.VALIDATED
                if cls_str == "validated"
                else ImpurityClassification.THEORETICAL
            ),
            rule_id=f"{category}:{rule_id}",
            rule_label=rule_label,
            molecular_formula=formula,
            molecular_weight=mw,
            db_name=db_rec.get("name") if db_rec else None,
            db_cid=str(db_rec["cid"]) if db_rec and db_rec.get("cid") is not None else None,
            ashby_tennant_alerts=alerts,
            highlight_atoms=hl,
        )
    except Exception:
        return None


def _append_api_catalog_degradants(
    result: ProcessAnalysisResult,
    api_smiles: str,
    add_unique,
) -> None:
    """목적물(API) — Universal SMIRKS 카탈로그 + 보관 안정성 분해산물."""
    from p41_degradation_smirks_rules import (
        passes_shelf_life_child_filter,
        predict_all_degradants,
    )

    parent_mol = _mol_from_smiles(api_smiles)
    if parent_mol is None:
        return
    parent_canon = canonical_smiles(api_smiles)
    cat_labels = {
        "Hydrolysis": "API 가수분해 분해산물",
        "Oxidation": "API 자동산화 분해산물",
        "NDSRI": "API NDSRI 분해",
    }
    try:
        grouped = predict_all_degradants(api_smiles)
    except Exception:
        grouped = {}

    for category, rows in grouped.items():
        for row in rows:
            smi = str(row.get("Degradant_SMILES") or "").strip()
            rule_name = str(row.get("Rule_Name") or "SMIRKS")
            if not smi or smi == parent_canon:
                continue
            child_mol = _mol_from_smiles(smi)
            if child_mol is None:
                continue
            try:
                # 2026-09-03 수정: 여기서 필터를 재적용하면서 가수분해 기전 면제가 빠져 있어,
                #   p41이 이미 내보낸 살리실산(아스피린 에스터 가수분해)이 다시 컬링되었습니다.
                #   rule_name 기준 면제를 적용하는 p41의 공용 함수로 교체합니다.
                if not passes_shelf_life_child_filter(
                    parent_canon or api_smiles, smi, rule_name
                ):
                    continue
            except Exception:
                pass
            label = f"{cat_labels.get(category, category)} — {rule_name}"
            add_unique(result.degradants, child_mol, rule_name, label, "degradant")


def _append_api_multistep_tree_degradants(
    result: ProcessAnalysisResult,
    api_smiles: str,
    add_unique,
    *,
    max_depth: int = 3,
) -> None:
    """목적물(API) — SMIRKS 재귀 다단계 분해 트리 노드 병합."""
    from p41_degradation_smirks_rules import (
        flatten_degradation_tree,
        generate_degradation_tree,
        passes_shelf_life_child_filter,
    )

    parent_mol = _mol_from_smiles(api_smiles)
    parent_canon = canonical_smiles(api_smiles)

    try:
        tree = generate_degradation_tree(api_smiles.strip(), max_depth=max_depth)
    except Exception:
        return
    if not tree:
        return

    for row in flatten_degradation_tree(tree):
        gen = int(row.get("Generation", 0))
        if gen <= 0:
            continue
        smi = str(row.get("SMILES") or "").strip()
        triggered = str(row.get("Triggered_By") or "SMIRKS tree")
        if not smi:
            continue
        mol = _mol_from_smiles(smi)
        if mol is None:
            continue
        try:
            # 2026-09-03 수정: 트리 노드에도 가수분해 기전 면제를 적용합니다
            #   (Triggered_By가 곧 rule_name입니다). 면제가 없으면 아스피린 트리의
            #   살리실산 노드가 여기서 다시 걸러졌습니다.
            if parent_mol is not None and not passes_shelf_life_child_filter(
                parent_canon or api_smiles, smi, triggered
            ):
                continue
        except Exception:
            pass
        label = f"다단계 분해 (Gen {gen}) — {triggered}"
        add_unique(result.degradants, mol, f"tree:{triggered}", label, "degradant")


def analyze_process_reaction(
    reactant1: str,
    reactant2: str,
    product: str,
    *,
    reagent: str = "",
    classification_fn,
    alert_fn,
    nitrite_in_process: bool = False,
) -> ProcessAnalysisResult:
    """
    출발물질 1·2 + 반응 시약 + 목적물 — 교차반응 · 시약 부반응 · NDSRI · 분해/미반응 스캔.
    """
    try:
        return _analyze_process_reaction_impl(
            reactant1,
            reactant2,
            product,
            reagent=reagent,
            classification_fn=classification_fn,
            alert_fn=alert_fn,
            nitrite_in_process=nitrite_in_process,
        )
    except Exception:
        return ProcessAnalysisResult()


def _analyze_process_reaction_impl(
    reactant1: str,
    reactant2: str,
    product: str,
    *,
    reagent: str = "",
    classification_fn,
    alert_fn,
    nitrite_in_process: bool = False,
) -> ProcessAnalysisResult:
    result = ProcessAnalysisResult()
    seen: set[str] = set()

    def add_unique(
        items: list[ImpurityResult],
        mol: Chem.Mol,
        rule_id: str,
        label: str,
        cat: str,
    ) -> None:
        try:
            smi = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            return
        if not smi or smi in seen:
            return
        seen.add(smi)
        imp = _make_impurity(mol, rule_id, label, classification_fn, alert_fn, cat)
        if imp:
            items.append(imp)

    mol1 = _mol_from_smiles(reactant1)
    mol2 = _mol_from_smiles(reactant2)
    mol_p = _mol_from_smiles(product)
    reagent = (reagent or "").strip()

    nitrite_active = nitrite_in_process or is_nitrite_source_smiles(reagent)
    reagent_profile: ReagentProfile | None = None
    if reagent:
        try:
            reagent_profile = classify_reagent(reagent)
            result.reagent_profile = reagent_profile
            if reagent_profile.is_nitrite_source:
                nitrite_active = True
        except Exception:
            pass

    # ── 1) 교차반응 (R1 + R2) ──
    if mol1 and mol2:
        try:
            for rule_id, label, smirks in BINARY_REACTIONS:
                if "nitros" in rule_id and not nitrite_active:
                    continue
                for prod in _run_binary_reaction(mol1, mol2, smirks):
                    add_unique(result.cross_reaction_products, prod, rule_id, label, "cross")
        except Exception:
            pass

    # ── 2) 시약(Reagent) 개입 — NDSRI · 산화/환원/알킬화 부반응 ──
    if reagent:
        amine_sources: list[tuple[str, str]] = []
        for smi_src, src_label in (
            (reactant1, "출발물질 1"),
            (reactant2, "출발물질 2"),
            (product, "목적물"),
        ):
            if smi_src:
                amine_sources.append((smi_src, src_label))

        if nitrite_active:
            for smi_src, src_label in amine_sources:
                try:
                    if reagent:
                        _, mols, meta = predict_nitrosation_with_reagent(smi_src, reagent)
                    else:
                        _, mols, meta = predict_reaction_impurities(smi_src, nitrite_present=True)
                    for mol, (rid, rlabel, _) in zip(mols, meta):
                        add_unique(
                            result.reagent_reaction_products,
                            mol,
                            rid,
                            f"{rlabel} ({src_label})",
                            "reagent_ndsri",
                        )
                except Exception:
                    continue

        if reagent_profile:
            for smi_src, src_label in amine_sources:
                try:
                    side = predict_reagent_side_products(smi_src, reagent, reagent_profile)
                    for mol, rid, rlabel in side:
                        add_unique(
                            result.reagent_reaction_products,
                            mol,
                            rid,
                            f"{rlabel} ({src_label})",
                            "reagent_side",
                        )
                except Exception:
                    continue

        # 시약 + 출발물질 교차 (R + Reagent)
        mol_reg = _mol_from_smiles(reagent)
        if mol_reg:
            for smi_src, src_label, mol_src in (
                (reactant1, "출발물질 1", mol1),
                (reactant2, "출발물질 2", mol2),
            ):
                if not mol_src:
                    continue
                try:
                    for rule_id, label, smirks in BINARY_REACTIONS:
                        if "nitros" in rule_id and not nitrite_active:
                            continue
                        for prod in _run_binary_reaction(mol_src, mol_reg, smirks):
                            add_unique(
                                result.reagent_reaction_products,
                                prod,
                                f"reagent_cross:{rule_id}",
                                f"{label} — {src_label} + 시약",
                                "reagent_cross",
                            )
                except Exception:
                    continue

    # ── 3) NDSRI 단독 (아질산 체크박스만, 시약 미입력) ──
    elif nitrite_active:
        for smi_src in (reactant1, reactant2):
            try:
                _, mols, meta = predict_reaction_impurities(smi_src, nitrite_present=True)
                for mol, (rid, rlabel, _) in zip(mols, meta):
                    add_unique(result.cross_reaction_products, mol, rid, rlabel, "ndsri")
            except Exception:
                continue

    # ── 4) 미반응 출발물질 · 시약 ──
    if product:
        for smi_src, label in (
            (reactant1, "미반응 출발물질 1"),
            (reactant2, "미반응 출발물질 2"),
            (reagent, "미반응 반응 시약"),
        ):
            if smi_src and _is_unreacted_starting_material(smi_src, product):
                mol = _mol_from_smiles(smi_src)
                if mol:
                    add_unique(result.unreacted_starting_materials, mol, "unreacted", label, "process")

    # ── 5) 목적물(API) 분해산물 — SMIRKS 카탈로그 + 다단계 트리 ──
    if mol_p:
        try:
            _append_api_catalog_degradants(result, product, add_unique)
            _append_api_multistep_tree_degradants(result, product, add_unique, max_depth=3)
        except Exception:
            pass
        if not result.degradants:
            try:
                from p44_degradation_shelf_life_filter import is_spin_stabilized_neutral_mol

                for rule_id, label, smirks in DEGRADATION_RULES:
                    for prod in _run_unimolecular(mol_p, smirks):
                        if not is_spin_stabilized_neutral_mol(prod):
                            continue
                        p_can = canonical_smiles(Chem.MolToSmiles(prod, canonical=True))
                        prod_can = canonical_smiles(product)
                        if p_can and prod_can and p_can != prod_can:
                            add_unique(result.degradants, prod, rule_id, label, "degradant")
            except Exception:
                pass

    return result
