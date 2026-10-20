"""공정 유연물질 예측 파이프라인 — 화면과 무관한 순수 계산.

p12_process_conditions_ui.py에서 추출했습니다. 원본은 Streamlit UI 모듈
안에 계산 로직이 들어 있어 `import streamlit` 없이는 쓸 수 없었습니다.
로직은 바꾸지 않았고 화면 의존성만 제거했습니다.

추출한 함수
    map_rxn_stress_to_prime          공정 가혹조건 → 분해엔진 조건명
    build_process_reactant_entries   캔버스 + Registry → batch 입력
    run_process_impurity_prediction  경로 B 전체 파이프라인
    build_process_conditions         UI 입력 → process_conditions dict

2026-08-25 추출.
"""

from __future__ import annotations

from typing import Any

from rdkit import Chem

from engines.p03_branding import log_process_exception
from engines.p05_safeguards import (
    MAX_PROCESS_PRODUCTS,
    ProcessRunContext,
    log_invalid_smiles_bypass,
)
from engines.p13_chemical_registry import (
    REGISTRY_NONE,
    build_process_reagent_context_from_ids,
    build_selected_materials_payload,
)
from engines.p33_process_impurity_engine import ProcessImpurityEngine
from engines.p34_name_reaction_process import NameReactionProcess
from engines.p36_suzuki_impurity_rules import annotate_impurity_records
from engines.p67_ich_m7_classifier import enrich_impurity_records
from p43_degradation_predictor import STRESS_ENV_CHOICES

__all__ = [
    "map_rxn_stress_to_prime",
    "build_process_reactant_entries",
    "run_process_impurity_prediction",
    "build_process_conditions",
    "CATEGORY_KEYS",
]

# UI 카테고리 이름 ↔ Registry가 받는 소문자 키
CATEGORY_KEYS: dict[str, str] = {
    "Solvents": "solvents",
    "Acids": "acids",
    "Bases": "bases",
    "Catalysts": "catalysts",
}


def map_rxn_stress_to_prime(stress_condition: str) -> str:
    """공정 stress 라벨 → PRIME-DS degradation STRESS_ENV_CHOICES."""
    s = ProcessImpurityEngine.normalize_stress_condition(stress_condition or "")
    mapping = {
        ProcessImpurityEngine.STRESS_NONE: STRESS_ENV_CHOICES[0],
        ProcessImpurityEngine.STRESS_ACIDIC_HYDROLYSIS: STRESS_ENV_CHOICES[1],
        ProcessImpurityEngine.STRESS_ALKALINE_HYDROLYSIS: STRESS_ENV_CHOICES[2],
        ProcessImpurityEngine.STRESS_OXIDATIVE: STRESS_ENV_CHOICES[4],
    }
    return mapping.get(s, STRESS_ENV_CHOICES[0])


def _valid_entry(label: str, key: str, smiles: str | None) -> tuple[str, str, dict[str, Any]] | None:
    smi = (smiles or "").strip()
    if not smi or Chem.MolFromSmiles(smi) is None:
        return None
    return (label, key, {"valid": True, "smiles": smi})


def build_process_reactant_entries(
    *,
    reactant1: str | None,
    reactant2: str | None,
    product: str | None,
    reagent_smiles: str | None,
    reference_api: str | None,
    process_reagent_context: dict[str, Any] | None = None,
) -> list[tuple[str, str, dict[str, Any]]]:
    """캔버스 SMILES + Registry 시약 → 공정 유연물질 batch 입력."""
    entries: list[tuple[str, str, dict[str, Any]]] = []
    api = (reference_api or product or "").strip()
    row = _valid_entry("Target API (목적물/대조)", "target_api", api)
    if row:
        entries.append(row)
    for label, key, smi in (
        ("출발물질 1", "reactant_1", reactant1),
        ("출발물질 2", "reactant_2", reactant2),
        ("반응 시약 (구조)", "reagent_struct", reagent_smiles),
        ("목적물 (Product)", "product", product),
    ):
        row = _valid_entry(label, key, smi)
        if row and row[2]["smiles"] != api:
            entries.append(row)

    ctx = process_reagent_context or {}
    seen_smiles = {e[2]["smiles"] for e in entries}
    for reg in ctx.get("registry_entries", []):
        smi = str(reg.get("smiles") or "").strip()
        if not smi or smi in seen_smiles:
            continue
        category = str(reg.get("category", "Reagent"))
        name = str(reg.get("name", category))
        row = _valid_entry(f"Registry · {category} · {name}", f"registry_{category}_{name}", smi)
        if row:
            entries.append(row)
            seen_smiles.add(smi)
    return entries


def build_process_conditions(
    *,
    stress_condition: str = "",
    name_reaction: str = "",
    solvents: list[str] | None = None,
    acids: list[str] | None = None,
    bases: list[str] | None = None,
    catalysts: list[str] | None = None,
    temperature_c: float = 25.0,
    reaction_time_h: float = 1.0,
    ph: float = 7.0,
    dynamic_catalog: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """UI 입력 → 엔진이 읽는 process_conditions dict.

    build_process_reagent_context_from_ids는 **소문자** 키를 받습니다
    (solvents/acids/bases/catalysts). 대문자로 넘기면 선택이 조용히
    무시되어 시약이 반영되지 않습니다.
    """
    id_selection = {
        "solvents": list(solvents or []),
        "acids": list(acids or []),
        "bases": list(bases or []),
        "catalysts": list(catalysts or []),
    }

    ctx = build_process_reagent_context_from_ids(id_selection, dynamic_catalog=dynamic_catalog)
    ctx["id_selection"] = id_selection
    ctx["dynamic_catalog"] = list(dynamic_catalog or [])
    ctx["selected_materials"] = build_selected_materials_payload(ctx)

    # 컨텍스트가 쓰는 키 이름은 primary_solvent / solvent_props 입니다.
    solvent_props = ctx.get("solvent_props") or {}

    return {
        "stress_condition": stress_condition or ProcessImpurityEngine.STRESS_NONE,
        "name_reaction": name_reaction or NameReactionProcess.NAME_RXN_NONE,
        "reagent": ctx.get("reagent") or REGISTRY_NONE,
        "reagent_smiles": list(ctx.get("all_smiles") or []),
        "selected_materials": ctx.get("selected_materials") or {},
        "process_reagents": ctx,
        "solvents": id_selection["solvents"],
        "acids": id_selection["acids"],
        "bases": id_selection["bases"],
        "catalysts": id_selection["catalysts"],
        "solvent": ctx.get("primary_solvent") or "—",
        "solvent_abbrev": solvent_props.get("abbrev") or "",
        "solvent_bp_c": solvent_props.get("bp_c"),
        "solvent_density_g_ml": solvent_props.get("density_g_ml"),
        "primary_catalyst": ctx.get("primary_catalyst") or REGISTRY_NONE,
        "homo_coupling_catalysts": list(ctx.get("homo_coupling_catalysts") or []),
        "temperature_c": float(temperature_c),
        "reaction_time_h": float(reaction_time_h),
        "ph": float(ph),
    }


def run_process_impurity_prediction(
    process_conditions: dict[str, Any],
    *,
    reactant1: str | None,
    reactant2: str | None,
    product: str | None,
    reagent_smiles: str | None,
    reference_api: str | None,
) -> list[dict[str, Any]]:
    """Universal Side-Reaction + Name Reaction + Stress batch 예측."""
    try:
        entries = build_process_reactant_entries(
            reactant1=reactant1,
            reactant2=reactant2,
            product=product,
            reagent_smiles=reagent_smiles,
            reference_api=reference_api,
            process_reagent_context=process_conditions.get("process_reagents"),
        )
        if not entries:
            return []

        stress = ProcessImpurityEngine.normalize_stress_condition(
            process_conditions.get("stress_condition", ProcessImpurityEngine.STRESS_NONE)
        )
        temp_c = float(process_conditions.get("temperature_c", 25.0))
        reagent_c = str(process_conditions.get("reagent", REGISTRY_NONE))
        process_reagent_ctx = process_conditions.get("process_reagents") or {}
        selected_materials = process_conditions.get(
            "selected_materials"
        ) or build_selected_materials_payload(process_reagent_ctx)
        has_registry_materials = any(
            selected_materials.get(cat) for cat in ("Solvents", "Acids", "Bases", "Catalysts")
        )

        ctx = ProcessRunContext(autostart=False)
        records = ProcessImpurityEngine.run_batch_impurity_prediction(
            entries,
            stress,
            temp_c,
            reagent=reagent_c,
            process_reagents={**process_reagent_ctx, "selected_materials": selected_materials},
            ctx=ctx,
        )

        name_rxn = process_conditions.get("name_reaction", NameReactionProcess.NAME_RXN_NONE)
        api_smi = (reference_api or product or "").strip() or None

        # Universal Side-Reaction Predictive Engine
        if ctx.check_continue():
            universal_records = NameReactionProcess.run_universal_prediction(
                entries,
                target_smiles=api_smi,
                selected_environment=selected_materials,
                temperature_c=temp_c,
                name_reaction=name_rxn,
                ctx=ctx,
            )
            records = NameReactionProcess.merge_into_batch(
                records, universal_records, max_products=MAX_PROCESS_PRODUCTS
            )

        # 선택적 Name Reaction (rdchiral)
        r1 = (reactant1 or "").strip()
        r2 = (reactant2 or "").strip()
        if not NameReactionProcess.is_skipped(name_rxn) and r1 and r2 and ctx.check_continue():
            nr_records = NameReactionProcess.simulate_impurities(
                r1,
                r2,
                name_rxn,
                target_smiles=api_smi,
                temperature_c=temp_c,
                reagent_context=process_reagent_ctx,
                ctx=ctx,
            )
            records = NameReactionProcess.merge_into_batch(
                records, nr_records, max_products=MAX_PROCESS_PRODUCTS
            )

        env_smiles = list(process_reagent_ctx.get("all_smiles") or [])
        if has_registry_materials and ctx.check_continue():
            engine = NameReactionProcess.get_process_engine()
            catalog_records = engine.run_batch_environment_simulation(
                entries,
                env_smiles,
                reagent_context={
                    **process_reagent_ctx,
                    "selected_materials": selected_materials,
                    "selected_id_lists": process_reagent_ctx.get("id_selection")
                    or {
                        "solvents": process_reagent_ctx.get("solvent_ids") or [],
                        "acids": process_reagent_ctx.get("acid_ids") or [],
                        "bases": process_reagent_ctx.get("base_ids") or [],
                        "catalysts": process_reagent_ctx.get("catalyst_ids") or [],
                    },
                },
                selected_materials=selected_materials,
                temperature_c=temp_c,
                ctx=ctx,
            )
            records = NameReactionProcess.merge_into_batch(
                records, catalog_records, max_products=MAX_PROCESS_PRODUCTS
            )

        # 2026-09-08 — 생성 경로 정보(모 화합물 구조)를 등급 판정부로 전달합니다.
        #   ICH M7 Class 4 는 «유연물질의 구조 알럿이 원료의약품과 공유되고,
        #   그 원료의약품에 변이원성 음성 실측 자료가 있는 경우» 이므로,
        #   모 화합물 구조가 없으면 이 분기가 아예 열리지 않습니다.
        parent_api_smiles = (product or reference_api or reactant1 or "") or None
        records = enrich_impurity_records(records, parent_smiles=parent_api_smiles)
        records = annotate_impurity_records(records)
        return ctx.truncate_records(records)
    except Exception as exc:
        log_process_exception("공정 유연물질 예측 파이프라인 실패", exc)
        log_invalid_smiles_bypass(context="run_process_impurity_prediction")
        return []
