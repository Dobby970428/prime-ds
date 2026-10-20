"""
PRIME-DS — Universal Side-Reaction Predictive Engine.

특정 Name Reaction에 의존하지 않고 [출발물질 배열 + 목적물 + 시약/용매] 조건으로부터
유기화학적 부반응 매커니즘을 동적으로 유도한다.
"""

# 경로 B
# 범용 부반응 12유형(Name Reaction에 안 걸리는 일반적인 부반응들)
# homo-coupling      같은 물질끼리 붙음 (이량체)
# des-halo           할로겐이 떨어져 나감
# protodeboronation  보론산이 떨어져 나감
# 용매 반응          용매가 끼어듦
# 알킬화             시약이 알킬기를 붙임
# ...

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from engines.p34_name_reaction_process import is_name_reaction_skipped
from engines.p33_process_impurity_engine import (
    compute_db_verification_score,
    get_verification_fingerprint_db,
    mol_from_smiles,
)
from engines.p05_safeguards import ProcessRunContext

SCREENING_ENGINE_TAG = "Universal_Side_Reaction_Engine"
AI_SCREENING_ENGINE_TAG = "External_AI_Reaction_Bridge"
_PATHWAY_PREFIX = "Process: Universal Side-Reaction"
_AI_PATHWAY_PREFIX = "Process: AI Forward"

MECHANISM_DISPLAY: dict[str, dict[str, str]] = {
    "Cross-Coupling Product": {
        "impurity_classification": "Main Product (Target)",
        "mechanism_ui_text": "의도된 교차결합 주산물 (할로겐 + 보론산)",
        "impurity_role": "Main",
        "impurity_warn": False,
    },
    "Halide Homo-coupling Byproduct": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "Homo-coupling 불순물 — 할로겐 출발물질 자가/교차 짝지음",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Boronic Homo-coupling Byproduct": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "Homo-coupling 불순물 — 보론산 출발물질 자가/교차 짝지음",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Des-halogenated Impurity": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "Des-halogenated 불순물 — 탈할로겐화(Des-bromo/Des-chloro)",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Des-boronated Impurity": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "Des-boronated 불순물 — 탈보론화(Protodeboronation)",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Solvent-Interacted Impurity": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "Solvent 부반응 산물 — 용매 축합/에스테르화",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Reagent-Alkylating Impurity": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "시약/용매 알킬화 부반응 불순물",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Acid-Catalyzed Side Product": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "산 촉매 부반응 — 카르복실산/에스터 전환",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Pd-Mediated Halide Dimerization": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "Pd 촉매 부반응 — 할로겐화물 자가/교차 이량화 (교차결합 미설정)",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Pd-Mediated Dehalogenation": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "Pd 촉매 부반응 — 탈할로겐화 (교차결합 미설정)",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "Base-Promoted Dehalogenation": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "염기 조건 탈할로겐화 — SNAr/탈할로겐 경로",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
    "AI-Generated Side Product": {
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "AI 생성 불순물 — ASKCOS/AiZynth generative forward",
        "impurity_role": "Side_Impurity",
        "impurity_warn": True,
    },
}

_REGISTRY_CATEGORIES = ("Solvents", "Acids", "Bases", "Catalysts")

_CROSS_COUPLING_NAME_KEYWORDS = (
    "suzuki",
    "buchwald",
    "heck",
    "sonogashira",
    "stille",
    "negishi",
    "kumada",
    "hiyama",
    "cross-coupling",
    "cross coupling",
)

_HALIDE_SMARTS = Chem.MolFromSmarts("[c:1]-[Cl,Br,I]")
_BORONIC_SMARTS = Chem.MolFromSmarts("[c:1]-B(O)O")


@dataclass(frozen=True)
class SideReactionContext:
    """어떤 부반응 규칙 팩을 적용할지 — Suzuki 미선택 시 결합 전용 규칙 억제."""

    has_halide: bool
    has_boronic: bool
    has_pd_catalyst: bool
    explicit_cross_coupling: bool
    halide_boronic_pair: bool
    run_cross_coupling: bool
    run_halide_homo_coupling: bool
    run_boronic_homo_coupling: bool
    run_deboronation: bool
    run_pd_halide_side: bool
    run_generic_dehalogenation: bool
    context_label: str


def _mol_has_substruct(mol: Chem.Mol, query: Chem.Mol | None) -> bool:
    if mol is None or query is None:
        return False
    try:
        return mol.HasSubstructMatch(query)
    except Exception:
        return False


def _mol_has_aromatic_halide(mol: Chem.Mol) -> bool:
    return _mol_has_substruct(mol, _HALIDE_SMARTS)


def _mol_has_boronic(mol: Chem.Mol) -> bool:
    return _mol_has_substruct(mol, _BORONIC_SMARTS)


def _env_has_pd_catalyst(env: dict[str, list[str]]) -> bool:
    for cat in _REGISTRY_CATEGORIES:
        for smi in env.get(cat) or []:
            text = str(smi).casefold()
            if "pd" in text or "palladium" in text or "dppf" in text:
                return True
    return False


def _is_explicit_cross_coupling_name(name_reaction: str | None) -> bool:
    label = str(name_reaction or "").strip()
    if not label or is_name_reaction_skipped(label):
        return False
    lowered = label.casefold()
    return any(keyword in lowered for keyword in _CROSS_COUPLING_NAME_KEYWORDS)


def build_side_reaction_context(
    reactant_mols: list[Chem.Mol],
    env: dict[str, list[str]],
    *,
    name_reaction: str | None = None,
) -> SideReactionContext:
    has_halide = any(_mol_has_aromatic_halide(mol) for mol in reactant_mols)
    has_boronic = any(_mol_has_boronic(mol) for mol in reactant_mols)
    has_pd = _env_has_pd_catalyst(env)
    explicit = _is_explicit_cross_coupling_name(name_reaction)
    pair = has_halide and has_boronic

    cross_coupling = pair or (explicit and has_halide and has_boronic)
    pd_side = has_halide and has_pd and not cross_coupling and not explicit

    if cross_coupling:
        label = f"Cross-coupling chemistry ({name_reaction or 'halide + boronic'})"
    elif pd_side:
        label = "Pd-catalyzed side reactions (no cross-coupling name reaction selected)"
    elif has_halide:
        label = "Generic halide side reactions (base/solvent — not Suzuki-specific)"
    else:
        label = "General process side reactions (reagent/environment matrix)"

    return SideReactionContext(
        has_halide=has_halide,
        has_boronic=has_boronic,
        has_pd_catalyst=has_pd,
        explicit_cross_coupling=explicit,
        halide_boronic_pair=pair,
        run_cross_coupling=cross_coupling,
        run_halide_homo_coupling=cross_coupling and has_halide,
        run_boronic_homo_coupling=cross_coupling and has_boronic,
        run_deboronation=cross_coupling and has_boronic,
        run_pd_halide_side=pd_side,
        run_generic_dehalogenation=has_halide and not cross_coupling,
        context_label=label,
    )


def _canonical_smiles(smiles: str) -> str:
    mol = mol_from_smiles(str(smiles or "").strip(), log_bypass=False)
    if mol is None:
        return str(smiles or "").strip()
    return Chem.MolToSmiles(mol, canonical=True)


def _sanitize_product_mol(prod: Chem.Mol | None) -> str | None:
    if prod is None:
        return None
    try:
        mol = Chem.Mol(prod)
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol, canonical=True)
        if not smi or mol_from_smiles(smi, log_bypass=False) is None:
            return None
        return smi
    except Exception:
        return None


def extract_reactant_smiles_from_entries(
    entries: list[tuple[str, str, dict[str, Any]]],
    *,
    exclude_target: str | None = None,
) -> list[str]:
    """batch entries → 출발물질 SMILES 풀 (목적물 제외)."""
    target_canon = _canonical_smiles(exclude_target) if exclude_target else ""
    pool: list[str] = []
    seen: set[str] = set()
    skip_keys = frozenset({"target_api"})
    for _label, key, meta in entries:
        if key in skip_keys:
            continue
        smi = str(meta.get("smiles") or "").strip()
        if not smi:
            continue
        canon = _canonical_smiles(smi)
        if target_canon and canon == target_canon:
            continue
        if canon in seen:
            continue
        seen.add(canon)
        pool.append(smi)
    return pool


def normalize_selected_environment(
    selected_environment: dict[str, Any] | None,
) -> dict[str, list[str]]:
    """4단 Registry 패널 → {Solvents, Acids, Bases, Catalysts} SMILES 배열."""
    env: dict[str, list[str]] = {cat: [] for cat in _REGISTRY_CATEGORIES}
    if not selected_environment:
        return env
    for cat in _REGISTRY_CATEGORIES:
        raw = selected_environment.get(cat) or []
        env[cat] = [str(s).strip() for s in raw if str(s or "").strip()]
    return env


class UniversalReactionImpurityEngine:
    """
    Universal Side-Reaction Predictive Engine.

    (1) 출발물질 교차 homo-coupling · cross-coupling
    (2) 탈할로겐화 · 탈보론화
    (3) 용매/시약 축합·치환 매트릭스
    """

    def __init__(self) -> None:
        self.cross_coupling_rxn = AllChem.ReactionFromSmarts(
            "[c:1]-[Cl,Br,I].[c:2]-B(O)O>>[c:1]-[c:2]"
        )
        self.halide_homo_coupling_rxn = AllChem.ReactionFromSmarts(
            "[c:1]-[Cl,Br,I,F].[c:2]-[Cl,Br,I,F]>>[c:1]-[c:2]"
        )
        self.boronic_homo_coupling_rxn = AllChem.ReactionFromSmarts(
            "[c:1]-B(O)O.[c:2]-B(O)O>>[c:1]-[c:2]"
        )
        self.dehalogenation_rxn = AllChem.ReactionFromSmarts("[c:1]-[Cl,Br,I]>>[c:1]")
        self.deboronation_rxn = AllChem.ReactionFromSmarts("[c:1]-B(O)O>>[c:1]")
        self.solvent_esterification = AllChem.ReactionFromSmarts(
            "[CX3:1](=[O:2])[OH:3].[OX2H:4][CX4:5]>>[CX3:1](=[O:2])[OX2:4][CX4:5]"
        )
        self.amine_alkylation = AllChem.ReactionFromSmarts(
            "[NX3;H2,H1:1].[Cl,Br,I:2][CX4:3]>>[NX3:1][CX4:3]"
        )

    def predict_all_potential_impurities(
        self,
        reactants_smiles: list[str],
        target_smiles: str,
        selected_environment: dict[str, list[str]],
        *,
        name_reaction: str | None = None,
    ) -> tuple[list[dict[str, Any]], SideReactionContext]:
        reactant_mols: list[Chem.Mol] = []
        reactant_canon: set[str] = set()
        for raw in reactants_smiles:
            smi = str(raw or "").strip()
            if not smi:
                continue
            mol = mol_from_smiles(smi, log_bypass=False)
            if mol is None:
                continue
            canon = _canonical_smiles(smi)
            reactant_mols.append(mol)
            reactant_canon.add(canon)

        target_canon = _canonical_smiles(target_smiles) if target_smiles else ""
        env = normalize_selected_environment(selected_environment)
        rxn_ctx = build_side_reaction_context(
            reactant_mols,
            env,
            name_reaction=name_reaction,
        )

        env_mols: list[tuple[Chem.Mol, str]] = []
        for cat, s_list in env.items():
            for smi in s_list:
                mol = mol_from_smiles(smi, log_bypass=False)
                if mol is not None:
                    env_mols.append((mol, cat))

        predicted_impurities: list[dict[str, Any]] = []
        has_base = bool(env.get("Bases"))

        if rxn_ctx.run_cross_coupling:
            for mol1, mol2 in itertools.combinations(reactant_mols, 2):
                for reactants in ((mol1, mol2), (mol2, mol1)):
                    self._run_reaction(
                        predicted_impurities,
                        self.cross_coupling_rxn,
                        reactants,
                        reaction_name="Cross-Coupling Product",
                        mechanism="할로겐-보론산 교차결합 주산물 후보",
                        reactant_canon=reactant_canon,
                        target_canon=target_canon,
                        exclude_target=False,
                        is_main_candidate=True,
                    )

        if rxn_ctx.run_halide_homo_coupling:
            for mol1, mol2 in itertools.combinations_with_replacement(reactant_mols, 2):
                self._run_reaction(
                    predicted_impurities,
                    self.halide_homo_coupling_rxn,
                    (mol1, mol2),
                    reaction_name="Halide Homo-coupling Byproduct",
                    mechanism="할로겐 출발물질 간 원치 않는 자가/교차 짝지음(이량체)",
                    reactant_canon=reactant_canon,
                    target_canon=target_canon,
                    exclude_target=False,
                )

        if rxn_ctx.run_boronic_homo_coupling:
            for mol1, mol2 in itertools.combinations_with_replacement(reactant_mols, 2):
                self._run_reaction(
                    predicted_impurities,
                    self.boronic_homo_coupling_rxn,
                    (mol1, mol2),
                    reaction_name="Boronic Homo-coupling Byproduct",
                    mechanism="보론산 출발물질 간 원치 않는 자가/교차 짝지음(이량체)",
                    reactant_canon=reactant_canon,
                    target_canon=target_canon,
                    exclude_target=False,
                )

        if rxn_ctx.run_pd_halide_side:
            for mol1, mol2 in itertools.combinations_with_replacement(reactant_mols, 2):
                if not _mol_has_aromatic_halide(mol1) and not _mol_has_aromatic_halide(mol2):
                    continue
                self._run_reaction(
                    predicted_impurities,
                    self.halide_homo_coupling_rxn,
                    (mol1, mol2),
                    reaction_name="Pd-Mediated Halide Dimerization",
                    mechanism="Pd 촉매 하 할로겐화물 이량화 — Suzuki/교차결합 미선택",
                    reactant_canon=reactant_canon,
                    target_canon=target_canon,
                    exclude_target=False,
                )

        for mol in reactant_mols:
            if not _mol_has_aromatic_halide(mol):
                continue
            if rxn_ctx.run_cross_coupling or rxn_ctx.run_pd_halide_side:
                self._run_reaction(
                    predicted_impurities,
                    self.dehalogenation_rxn,
                    (mol,),
                    reaction_name=(
                        "Des-halogenated Impurity"
                        if rxn_ctx.run_cross_coupling
                        else "Pd-Mediated Dehalogenation"
                    ),
                    mechanism=(
                        "촉매 중간체 양성자 치환 — 탈할로겐화(Des-bromo/Des-chloro)"
                        if rxn_ctx.run_cross_coupling
                        else "Pd 촉매 조건 탈할로겐화 부반응"
                    ),
                    reactant_canon=reactant_canon,
                    target_canon=target_canon,
                    exclude_target=False,
                )
            elif rxn_ctx.run_generic_dehalogenation and has_base:
                self._run_reaction(
                    predicted_impurities,
                    self.dehalogenation_rxn,
                    (mol,),
                    reaction_name="Base-Promoted Dehalogenation",
                    mechanism="염기(K2CO3 등) 조건 탈할로겐화 — SNAr/탈할로겐",
                    reactant_canon=reactant_canon,
                    target_canon=target_canon,
                    exclude_target=False,
                )

            if rxn_ctx.run_deboronation and _mol_has_boronic(mol):
                self._run_reaction(
                    predicted_impurities,
                    self.deboronation_rxn,
                    (mol,),
                    reaction_name="Des-boronated Impurity",
                    mechanism="염기·수분 조건 — 탈보론화(Protodeboronation)",
                    reactant_canon=reactant_canon,
                    target_canon=target_canon,
                    exclude_target=False,
                )

        # (3) 용매/시약 축합·치환 매트릭스
        for mol in reactant_mols:
            for env_mol, cat_name in env_mols:
                if cat_name == "Solvents":
                    self._run_reaction(
                        predicted_impurities,
                        self.solvent_esterification,
                        (mol, env_mol),
                        reaction_name="Solvent-Interacted Impurity",
                        mechanism="투입 용매와 출발물질 축합/에스테르화 부반응",
                        reactant_canon=reactant_canon,
                        target_canon=target_canon,
                        exclude_target=True,
                    )
                if cat_name == "Acids":
                    self._run_reaction(
                        predicted_impurities,
                        self.solvent_esterification,
                        (env_mol, mol),
                        reaction_name="Acid-Catalyzed Side Product",
                        mechanism="투입 산과 출발물질의 산 촉매 부반응",
                        reactant_canon=reactant_canon,
                        target_canon=target_canon,
                        exclude_target=True,
                    )
                if cat_name in ("Bases", "Solvents"):
                    self._run_reaction(
                        predicted_impurities,
                        self.amine_alkylation,
                        (mol, env_mol),
                        reaction_name="Reagent-Alkylating Impurity",
                        mechanism="시약/용매 친핵 공격 — 알킬화 부반응",
                        reactant_canon=reactant_canon,
                        target_canon=target_canon,
                        exclude_target=True,
                    )

        return predicted_impurities, rxn_ctx

    def _run_reaction(
        self,
        target_list: list[dict[str, Any]],
        rxn: Any,
        reactants: tuple[Chem.Mol, ...],
        *,
        reaction_name: str,
        mechanism: str,
        reactant_canon: set[str],
        target_canon: str,
        exclude_target: bool,
        is_main_candidate: bool = False,
    ) -> None:
        if rxn is None:
            return
        try:
            products = rxn.RunReactants(reactants)
        except Exception:
            return
        rule_seen: set[str] = set()
        for prod_tuple in products or []:
            if prod_tuple is None:
                continue
            items = prod_tuple if isinstance(prod_tuple, (tuple, list)) else (prod_tuple,)
            for prod in items:
                psmiles = _sanitize_product_mol(prod)
                if not psmiles:
                    continue
                dedupe_key = (psmiles, reaction_name)
                if dedupe_key in rule_seen:
                    continue
                rule_seen.add(dedupe_key)
                if exclude_target and target_canon and psmiles == target_canon:
                    continue
                if not is_main_candidate and psmiles in reactant_canon:
                    continue
                is_main = bool(
                    is_main_candidate and target_canon and psmiles == target_canon
                )
                self._add_impurity(
                    target_list,
                    psmiles,
                    reaction_name,
                    mechanism,
                    is_main_product=is_main,
                )

    def _add_impurity(
        self,
        target_list: list[dict[str, Any]],
        smiles: str,
        reaction_name: str,
        mechanism: str,
        *,
        is_main_product: bool = False,
    ) -> None:
        key = (smiles, reaction_name)
        if any(
            (
                item.get("smiles"),
                item.get("origin_reaction") or item.get("reaction_name"),
            )
            == key
            for item in target_list
        ):
            return
        meta = MECHANISM_DISPLAY.get(reaction_name, {})
        target_list.append(
            {
                "smiles": smiles,
                "origin_reaction": reaction_name,
                "reaction_name": reaction_name,
                "mechanism": mechanism,
                "screening_engine": SCREENING_ENGINE_TAG,
                "impurity_role": "Main" if is_main_product else meta.get("impurity_role", "Side_Impurity"),
                "impurity_classification": (
                    "Main Product (Target)"
                    if is_main_product
                    else meta.get("impurity_classification", "Side Impurity")
                ),
                "mechanism_ui_text": meta.get("mechanism_ui_text", mechanism),
                "impurity_warn": False if is_main_product else bool(meta.get("impurity_warn", True)),
                "is_main_product": is_main_product,
            }
        )


_ENGINE: UniversalReactionImpurityEngine | None = None


def get_universal_reaction_engine() -> UniversalReactionImpurityEngine:
    global _ENGINE
    if _ENGINE is None:
        _ENGINE = UniversalReactionImpurityEngine()
    return _ENGINE


def _merge_prediction_lists(
    *lists: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """(SMILES, reaction_name) 기준 병합 — AI + RDKit 공존."""
    merged: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for preds in lists:
        for pred in preds or []:
            smi = str(pred.get("smiles") or "").strip()
            rxn = str(pred.get("reaction_name") or pred.get("origin_reaction") or "Universal")
            key = (smi, rxn)
            if not smi or key in seen:
                continue
            seen.add(key)
            merged.append(pred)
    return merged


def universal_predictions_to_impurity_records(
    predictions: list[dict[str, Any]],
    *,
    source_label: str = "Universal · reactant pool",
    temperature_c: float = 25.0,
    db_fps: list | None = None,
    db_smiles: list[str] | None = None,
) -> list[dict[str, Any]]:
    """범용 예측 → batch 레코드 (ICH M7 Class 1~5는 enrich_impurity_records에서 일괄 부여)."""
    if db_fps is None or db_smiles is None:
        db_fps, db_smiles = get_verification_fingerprint_db()

    records: list[dict[str, Any]] = []
    for pred in predictions:
        smi = str(pred.get("smiles") or "").strip()
        mol = mol_from_smiles(smi, log_bypass=False)
        if mol is None:
            continue
        rxn_name = str(pred.get("reaction_name") or pred.get("origin_reaction") or "Universal")
        mech = str(pred.get("mechanism") or "")
        is_ai = pred.get("screening_engine") == AI_SCREENING_ENGINE_TAG or pred.get("ai_confidence") is not None
        prefix = _AI_PATHWAY_PREFIX if is_ai else _PATHWAY_PREFIX
        pathway = f"{prefix} · {rxn_name}"
        if mech:
            pathway = f"{pathway} · {mech}"

        match_pct, status_tag, reliability_pct = compute_db_verification_score(
            mol, db_fps, db_smiles, temperature_c=temperature_c
        )
        ai_conf = pred.get("ai_confidence")
        if ai_conf is not None:
            try:
                conf_f = float(ai_conf)
                confidence_display = (
                    f"AI {conf_f * 100:.1f}% · DB {reliability_pct:.1f}% — {status_tag}"
                )
            except (TypeError, ValueError):
                confidence_display = f"{reliability_pct:.1f}% — {status_tag}"
        else:
            confidence_display = f"{reliability_pct:.1f}% — {status_tag}"
        is_main = bool(pred.get("is_main_product"))
        records.append(
            {
                "pathway": pathway,
                "smiles": smi,
                "exact_mw": round(Descriptors.ExactMolWt(mol), 2),
                "db_match_pct": match_pct,
                "db_status_tag": status_tag,
                "reliability_pct": reliability_pct,
                "confidence_display": confidence_display,
                "ai_confidence": ai_conf,
                "ai_backend": pred.get("ai_backend"),
                "source_reactant": source_label,
                "name_reaction": rxn_name,
                "reaction_name": rxn_name,
                "origin_reaction": rxn_name,
                "is_name_reaction_byproduct": not is_main,
                "catalog_mechanism": mech,
                "impurity_role": pred.get("impurity_role", "Side_Impurity"),
                "impurity_classification": pred.get(
                    "impurity_classification", "Side Impurity"
                ),
                "mechanism_ui_text": pred.get("mechanism_ui_text") or mech or rxn_name,
                "impurity_warn": pred.get("impurity_warn", not is_main),
                "screening_engine": pred.get("screening_engine", SCREENING_ENGINE_TAG),
            }
        )
    return records


def run_universal_side_reaction_pipeline(
    reactants_smiles: list[str],
    *,
    target_smiles: str = "",
    selected_environment: dict[str, list[str]] | None = None,
    temperature_c: float = 25.0,
    source_label: str | None = None,
    name_reaction: str | None = None,
    ctx: ProcessRunContext | None = None,
) -> dict[str, Any]:
    """
    Universal Side-Reaction Predictive Engine — 메인 API.

    Returns
    -------
    dict with keys:
        predicted_impurities — raw mechanism-tagged predictions
        records — ProcessImpurityEngine batch 레코드 (ICH enrich 전)
        scenario_count
    """
    if ctx is not None and not ctx.check_continue():
        return {"predicted_impurities": [], "records": [], "scenario_count": 0}

    cleaned = [str(s).strip() for s in reactants_smiles if str(s or "").strip()]
    if not cleaned:
        return {"predicted_impurities": [], "records": [], "scenario_count": 0}

    engine = get_universal_reaction_engine()
    env = normalize_selected_environment(selected_environment)
    local_impurities, rxn_ctx = engine.predict_all_potential_impurities(
        cleaned,
        target_smiles or "",
        env,
        name_reaction=name_reaction,
    )

    from engines.p38_external_ai_reaction_bridge import get_external_ai_reaction_bridge

    bridge = get_external_ai_reaction_bridge()
    ai_bundle = bridge.predict_side_products(
        cleaned,
        target_smiles=target_smiles or "",
        selected_environment=env,
    )
    ai_impurities = ai_bundle.get("predicted_impurities") or []

    predicted_impurities = _merge_prediction_lists(local_impurities, ai_impurities)
    label = source_label or "Universal Side-Reaction · reactant + environment"
    records = universal_predictions_to_impurity_records(
        predicted_impurities,
        source_label=label,
        temperature_c=temperature_c,
    )
    if records:
        records[0] = {
            **records[0],
            "_pipeline_meta": {
                "ai_bridge_status": ai_bundle.get("status"),
                "ai_bridge_backend": ai_bundle.get("backend_used"),
                "ai_scenario_count": len(ai_impurities),
                "reaction_context": rxn_ctx.context_label,
                "cross_coupling_mode": rxn_ctx.run_cross_coupling,
            },
        }
    return {
        "predicted_impurities": predicted_impurities,
        "records": records,
        "scenario_count": len(predicted_impurities),
        "screening_engine": SCREENING_ENGINE_TAG,
        "ai_bridge_status": ai_bundle.get("status"),
        "ai_bridge_backend": ai_bundle.get("backend_used"),
        "ai_scenario_count": len(ai_impurities),
        "local_scenario_count": len(local_impurities),
        "reaction_context": rxn_ctx.context_label,
        "cross_coupling_mode": rxn_ctx.run_cross_coupling,
    }


def run_universal_impurity_simulation(
    reactants_smiles: list[str],
    *,
    target_smiles: str = "",
    selected_environment: dict[str, list[str]] | None = None,
    temperature_c: float = 25.0,
    source_label: str | None = None,
    ctx: ProcessRunContext | None = None,
) -> list[dict[str, Any]]:
    """하위 호환 — records 리스트만 반환."""
    bundle = run_universal_side_reaction_pipeline(
        reactants_smiles,
        target_smiles=target_smiles,
        selected_environment=selected_environment,
        temperature_c=temperature_c,
        source_label=source_label,
        ctx=ctx,
    )
    return bundle.get("records") or []
