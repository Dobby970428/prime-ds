"""
PRIME-DS — Suzuki-Miyaura 종합 부반응 SMIRKS 규칙 엔진.

주반응 · 할로겐 homo-coupling · 보론산 homo-coupling · 탈할로겐화 · 탈보론화를
RunReactants 기반으로 일괄 예측한다.
"""

# 경로 B
# Suzuki커플링이 의약품 합성에서 워낙 많이 쓰여서 만든 전용 파일
# 규칙 5개:
# ① 주반응 (원하는 것)
# ② 할라이드 homo-coupling
# ③ 보론산 homo-coupling
# ④ des-halo (탈할로겐)
# ⑤ protodeboronation

from __future__ import annotations

from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from engines.p33_process_impurity_engine import (
    classify_ich_m7_genotoxicity,
    compute_db_verification_score,
    get_verification_fingerprint_db,
    mol_from_smiles,
)

_HALIDE_QUERY = Chem.MolFromSmarts("[c;!$(c-B)]-[F,Cl,Br,I]")
_BORONIC_QUERY = Chem.MolFromSmarts("[c]-B(O)O")
_BORONIC_ESTER_QUERY = Chem.MolFromSmarts("[c]-B(-O)-O")

SUZUKI_COMPREHENSIVE_RULES: dict[str, dict[str, str]] = {
    "Suzuki_Main_Reaction": {
        "smirks": "[c:1]-[F,Cl,Br,I].[c:2]-B(O)O>>[c:1]-[c:2]",
        "type": "Main",
        "description": "의도된 주성분 Suzuki Cross-Coupling 결합 구조",
    },
    "Halide_Homo_Coupling_Impurity": {
        "smirks": "[c:1]-[Cl,Br,I].[c:2]-[Cl,Br,I]>>[c:1]-[c:2]",
        "type": "Side_Impurity",
        "description": "할로겐화 출발물질 간의 자체 Pd 촉매 호모커플링(이량체 유연물질)",
    },
    "Boronate_Homo_Coupling_Impurity": {
        "smirks": "[c:1]-B(O)O.[c:2]-B(O)O>>[c:1]-[c:2]",
        "type": "Side_Impurity",
        "description": "보론산 출발물질 간의 산화적 자가 짝지음(이량체 유연물질)",
    },
    "Hydro_Dehalogenation_Impurity": {
        "smirks": "[c:1]-[Cl,Br,I]>>[c:1]",
        "type": "Side_Impurity",
        "description": "산화적 첨가중간체의 양성자 치환에 의한 탈할로겐화(Des-bromo/Des-iodo 유연물질)",
    },
    "Protodeboronation_Impurity": {
        "smirks": "[c:1]-B(O)O>>[c:1]",
        "type": "Side_Impurity",
        "description": "염기 및 수분 조건 하 보론산 기능기 탈락 분해산물",
    },
}

_PATHWAY_PREFIX = "Process: Suzuki"
SCREENING_ENGINE_TAG = "Mechanism_Competitive_Suzuki"

SUZUKI_REACTION_DISPLAY_NAMES: dict[str, str] = {
    "Suzuki_Main_Reaction": "Suzuki Cross-Coupling (Main Product)",
    "Halide_Homo_Coupling_Impurity": "Halide Homo-coupling",
    "Boronate_Homo_Coupling_Impurity": "Boronic Acid Homo-coupling",
    "Hydro_Dehalogenation_Impurity": "Hydro-dehalogenation (Des-bromo)",
    "Protodeboronation_Impurity": "Protodeboronation (Des-boro)",
}

# UI · 리포트용 메커니즘 디스크립션 매핑
SUZUKI_MECHANISM_UI: dict[str, dict[str, Any]] = {
    "Suzuki_Main_Reaction": {
        "impurity_role": "Main",
        "impurity_classification": "Main Product (Target)",
        "mechanism_ui_text": "의도된 Suzuki-Miyaura 교차결합 주산물",
        "impurity_warn": False,
    },
    "Halide_Homo_Coupling_Impurity": {
        "impurity_role": "Side_Impurity",
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "출발물질 자체 결합성 공정 불순물 (Aryl halide Homo-coupling)",
        "impurity_warn": True,
    },
    "Boronate_Homo_Coupling_Impurity": {
        "impurity_role": "Side_Impurity",
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "출발물질 자체 결합성 공정 불순물 (Boronic acid Homo-coupling)",
        "impurity_warn": True,
    },
    "Hydro_Dehalogenation_Impurity": {
        "impurity_role": "Side_Impurity",
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "촉매 중간체 프로톤 분해산물 (Hydro-dehalogenation / Des-bromo)",
        "impurity_warn": True,
    },
    "Protodeboronation_Impurity": {
        "impurity_role": "Side_Impurity",
        "impurity_classification": "Side Impurity",
        "mechanism_ui_text": "촉매 중간체 프로톤 분해산물 (Protodeboronation / Des-boro)",
        "impurity_warn": True,
    },
}


def annotate_impurity_classification(record: dict[str, Any]) -> dict[str, Any]:
    """레코드에 유연물질 분류·메커니즘 UI 꼬리표 부여."""
    enriched = dict(record)
    rule = str(
        enriched.get("suzuki_rule")
        or enriched.get("name_reaction")
        or ""
    )
    meta = SUZUKI_MECHANISM_UI.get(rule)
    if meta:
        enriched.update(meta)
        enriched["screening_engine"] = SCREENING_ENGINE_TAG
        return enriched

    if enriched.get("impurity_classification"):
        return enriched

    if enriched.get("is_name_reaction_byproduct"):
        enriched["impurity_role"] = "Side_Impurity"
        enriched["impurity_classification"] = "Side Impurity"
        enriched["mechanism_ui_text"] = str(
            enriched.get("catalog_mechanism")
            or enriched.get("pathway", "")
        )[:120]
        enriched["impurity_warn"] = True
    elif "main" in str(enriched.get("pathway", "")).casefold():
        enriched["impurity_role"] = "Main"
        enriched["impurity_classification"] = "Main Product (Target)"
        enriched["mechanism_ui_text"] = "Name Reaction 주산물"
        enriched["impurity_warn"] = False
    else:
        enriched["impurity_role"] = enriched.get("impurity_role") or "Process"
        enriched["impurity_classification"] = "Process Impurity"
        enriched["mechanism_ui_text"] = str(enriched.get("pathway", ""))[:120]
        enriched["impurity_warn"] = False
    return enriched


def annotate_impurity_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [annotate_impurity_classification(r) for r in records]


def _sanitize_product_mol(prod: Chem.Mol | None) -> str | None:
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
        Chem.SanitizeMol(verify)
        return Chem.MolToSmiles(verify, canonical=True)
    except Exception:
        return None


def _mol_has_halide(mol: Chem.Mol) -> bool:
    return bool(_HALIDE_QUERY and mol.HasSubstructMatch(_HALIDE_QUERY))


def _mol_has_boronic(mol: Chem.Mol) -> bool:
    if _BORONIC_QUERY and mol.HasSubstructMatch(_BORONIC_QUERY):
        return True
    return bool(_BORONIC_ESTER_QUERY and mol.HasSubstructMatch(_BORONIC_ESTER_QUERY))


def resolve_suzuki_reactants_from_pool(
    smiles_pool: list[str],
) -> tuple[str, str] | None:
    """출발물질·시약 풀에서 (aryl halide, boronic acid) 쌍 탐색."""
    halides: list[str] = []
    boronics: list[str] = []
    for raw in smiles_pool:
        smi = str(raw or "").strip()
        if not smi:
            continue
        mol = mol_from_smiles(smi, log_bypass=False)
        if mol is None:
            continue
        if _mol_has_halide(mol):
            halides.append(smi)
        if _mol_has_boronic(mol):
            boronics.append(smi)
    if halides and boronics:
        return halides[0], boronics[0]
    return None


def resolve_suzuki_reactants(
    smiles_a: str,
    smiles_b: str,
) -> tuple[str, str] | None:
    """(aryl halide SMILES, boronic acid SMILES) 또는 None."""
    smi_a = str(smiles_a or "").strip()
    smi_b = str(smiles_b or "").strip()
    mol_a = mol_from_smiles(smi_a, log_bypass=False)
    mol_b = mol_from_smiles(smi_b, log_bypass=False)
    if mol_a is None or mol_b is None:
        return None
    if _mol_has_halide(mol_a) and _mol_has_boronic(mol_b):
        return smi_a, smi_b
    if _mol_has_halide(mol_b) and _mol_has_boronic(mol_a):
        return smi_b, smi_a
    return None


def is_suzuki_coupling_context(
    name_reaction: str | None,
    reagent_context: dict[str, Any] | None = None,
) -> bool:
    """Name Reaction 또는 Pd 촉매 Registry 선택 시 Suzuki 엔진 가동."""
    label = str(name_reaction or "")
    if "suzuki" in label.casefold():
        return True
    ctx = reagent_context or {}
    catalysts = [str(c).casefold() for c in (ctx.get("catalysts") or [])]
    if any("pd" in c or "palladium" in c for c in catalysts):
        return True
    for entry in ctx.get("registry_entries") or []:
        if str(entry.get("category")) != "Catalysts":
            continue
        metal = str(entry.get("metal") or "").strip()
        if metal == "Pd":
            return True
        name = str(entry.get("name") or "").casefold()
        if "pd" in name or "palladium" in name:
            return True
    return False


class SuzukiImpurityPredictor:
    """Suzuki-Miyaura 5종 SMIRKS 시나리오 예측기."""

    def __init__(self) -> None:
        self.rules: dict[str, Any] = {}
        for rule_name, spec in SUZUKI_COMPREHENSIVE_RULES.items():
            smirks = str(spec.get("smirks") or "")
            if ">>" not in smirks:
                continue
            try:
                rxn = AllChem.ReactionFromSmarts(smirks)
                if rxn is not None:
                    self.rules[rule_name] = rxn
            except Exception:
                continue

    def predict_all_scenarios(
        self,
        halide_smiles: str,
        boronic_smiles: str,
    ) -> list[dict[str, Any]]:
        """
        5종 Suzuki 시나리오를 규칙별 RunReactants로 실행.

        동일 SMILES라도 메커니즘(rule)마다 별도 항목을 유지한다.
        """
        halide_mol = mol_from_smiles(halide_smiles, log_bypass=False)
        boronic_mol = mol_from_smiles(boronic_smiles, log_bypass=False)
        if halide_mol is None or boronic_mol is None:
            return []

        results: list[dict[str, Any]] = []
        rules_fired: set[str] = set()

        for rule_name, rxn in self.rules.items():
            rule_info = SUZUKI_COMPREHENSIVE_RULES[rule_name]
            reactant_sets: list[tuple[Chem.Mol, ...]] = []

            if rule_name == "Suzuki_Main_Reaction":
                reactant_sets = [(halide_mol, boronic_mol)]
            elif rule_name == "Halide_Homo_Coupling_Impurity":
                reactant_sets = [(halide_mol, halide_mol)]
            elif rule_name == "Boronate_Homo_Coupling_Impurity":
                reactant_sets = [(boronic_mol, boronic_mol)]
            elif rule_name == "Hydro_Dehalogenation_Impurity":
                reactant_sets = [(halide_mol,)]
            elif rule_name == "Protodeboronation_Impurity":
                reactant_sets = [(boronic_mol,)]
            else:
                continue

            rule_seen: set[str] = set()
            for reactants in reactant_sets:
                try:
                    products = rxn.RunReactants(reactants)
                except Exception:
                    continue
                for prod_tuple in products or []:
                    if prod_tuple is None:
                        continue
                    items = (
                        prod_tuple
                        if isinstance(prod_tuple, (tuple, list))
                        else (prod_tuple,)
                    )
                    for prod in items:
                        prod_smiles = _sanitize_product_mol(prod)
                        if not prod_smiles or prod_smiles in rule_seen:
                            continue
                        rule_seen.add(prod_smiles)
                        if rule_name in rules_fired:
                            continue
                        rules_fired.add(rule_name)
                        ui_meta = SUZUKI_MECHANISM_UI.get(rule_name, {})
                        display_name = SUZUKI_REACTION_DISPLAY_NAMES.get(rule_name, rule_name)
                        results.append(
                            {
                                "smiles": prod_smiles,
                                "reaction_name": display_name,
                                "suzuki_rule": rule_name,
                                "type": rule_info["type"],
                                "mechanism": rule_info["description"],
                                "impurity_role": ui_meta.get(
                                    "impurity_role",
                                    "Main" if rule_info["type"] == "Main" else "Side_Impurity",
                                ),
                                "impurity_classification": ui_meta.get(
                                    "impurity_classification",
                                    "Main Product (Target)"
                                    if rule_info["type"] == "Main"
                                    else "Side Impurity",
                                ),
                                "mechanism_ui_text": ui_meta.get(
                                    "mechanism_ui_text", rule_info["description"]
                                ),
                                "impurity_warn": ui_meta.get(
                                    "impurity_warn", rule_info["type"] != "Main"
                                ),
                                "screening_engine": SCREENING_ENGINE_TAG,
                            }
                        )
                        break
                    if rule_name in rules_fired:
                        break
                if rule_name in rules_fired:
                    break
        return results

    def predict_all_scenarios_marshaled(
        self,
        halide_smiles: str,
        boronic_smiles: str,
    ) -> dict[str, Any]:
        """5종 시나리오 → predicted_impurities 배열로 마샬링."""
        predicted = self.predict_all_scenarios(halide_smiles, boronic_smiles)
        return {
            "predicted_impurities": predicted,
            "scenario_count": len(predicted),
            "halide_smiles": halide_smiles,
            "boronic_smiles": boronic_smiles,
            "screening_engine": SCREENING_ENGINE_TAG,
        }

    def predictions_to_impurity_records(
        self,
        predictions: list[dict[str, Any]],
        *,
        source_label: str = "Suzuki · Aryl halide + Boronic acid",
        temperature_c: float = 25.0,
        db_fps: list | None = None,
        db_smiles: list[str] | None = None,
    ) -> list[dict[str, Any]]:
        """Suzuki 예측 → ProcessImpurityEngine batch 레코드."""
        if db_fps is None or db_smiles is None:
            db_fps, db_smiles = get_verification_fingerprint_db()

        records: list[dict[str, Any]] = []
        for pred in predictions:
            smi = str(pred.get("smiles") or "").strip()
            mol = mol_from_smiles(smi, log_bypass=False)
            if mol is None:
                continue
            rule = str(pred.get("suzuki_rule") or pred.get("reaction_name", "Suzuki"))
            display_rule = SUZUKI_REACTION_DISPLAY_NAMES.get(rule, str(pred.get("reaction_name", rule)))
            mech = str(pred.get("mechanism", ""))
            is_side = str(pred.get("type", "")) != "Main"
            ui_meta = SUZUKI_MECHANISM_UI.get(rule, {})
            tag = "Side Impurity" if is_side else "Main Product"
            pathway = f"{_PATHWAY_PREFIX} · {display_rule} · {tag}"
            if mech:
                pathway = f"{pathway} · {mech}"

            ich_class, ich_detail = classify_ich_m7_genotoxicity(mol)
            match_pct, status_tag, reliability_pct = compute_db_verification_score(
                mol, db_fps, db_smiles, temperature_c=temperature_c
            )
            record = {
                "pathway": pathway,
                "smiles": smi,
                "exact_mw": round(Descriptors.ExactMolWt(mol), 2),
                "ich_m7": ich_class,
                "ich_m7_detail": ich_detail,
                "db_match_pct": match_pct,
                "db_status_tag": status_tag,
                "reliability_pct": reliability_pct,
                "confidence_display": f"{reliability_pct:.1f}% — {status_tag}",
                "source_reactant": source_label,
                "name_reaction": display_rule,
                "reaction_name": display_rule,
                "suzuki_rule": rule,
                "is_name_reaction_byproduct": is_side,
                "catalog_mechanism": mech,
                "impurity_role": pred.get("impurity_role")
                or ui_meta.get("impurity_role", "Side_Impurity" if is_side else "Main"),
                "impurity_classification": pred.get("impurity_classification")
                or ui_meta.get(
                    "impurity_classification",
                    "Main Product (Target)" if not is_side else "Side Impurity",
                ),
                "mechanism_ui_text": pred.get("mechanism_ui_text")
                or ui_meta.get("mechanism_ui_text", mech),
                "impurity_warn": pred.get("impurity_warn", is_side),
                "screening_engine": SCREENING_ENGINE_TAG,
            }
            records.append(annotate_impurity_classification(record))
        return records


_PREDICTOR: SuzukiImpurityPredictor | None = None


def get_suzuki_predictor() -> SuzukiImpurityPredictor:
    global _PREDICTOR
    if _PREDICTOR is None:
        _PREDICTOR = SuzukiImpurityPredictor()
    return _PREDICTOR


def run_suzuki_impurity_simulation(
    smiles_a: str,
    smiles_b: str,
    *,
    temperature_c: float = 25.0,
    source_label: str | None = None,
) -> list[dict[str, Any]]:
    """할로겐화물 + 보론산 자동 매칭 후 5종 시나리오 예측."""
    pair = resolve_suzuki_reactants(smiles_a, smiles_b)
    if not pair:
        pool = [s for s in (smiles_a, smiles_b) if str(s or "").strip()]
        pair = resolve_suzuki_reactants_from_pool(pool)
    if not pair:
        return []
    halide, boronic = pair
    bundle = get_suzuki_predictor().predict_all_scenarios_marshaled(halide, boronic)
    preds = bundle.get("predicted_impurities") or []
    label = source_label or "Suzuki · Halide + Boronic acid"
    return get_suzuki_predictor().predictions_to_impurity_records(
        preds,
        source_label=label,
        temperature_c=temperature_c,
    )
