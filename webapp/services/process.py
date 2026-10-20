"""모드 2 — 공정 불순물 예측.

⚠️ 이 화면에는 계산 경로가 둘 있습니다.

    경로 A  p39_hybrid_process_engine   → 판정 p66  (상향 전용)
    경로 B  engines/p33·p34·p35·p36     → 판정 p67  (하향 허용)

서로 부르지 않는 독립된 시스템이라 같은 화합물이 다른 등급으로 나올 수
있습니다. 검증 스크립트(s5)와 논문 수치는 경로 A(p66) 기준입니다.
화면에서도 어느 경로 결과인지 명시합니다.
"""

from __future__ import annotations

from typing import Any


# ---------------------------------------------------------------------------
# 선택지
# ---------------------------------------------------------------------------


def options() -> dict[str, Any]:
    out: dict[str, Any] = {
        "stress": [],
        "name_reactions": [],
        "registry": {},
        "master_count": 0,
        "errors": [],
    }

    try:
        from engines.p33_process_impurity_engine import ProcessImpurityEngine

        out["stress"] = list(ProcessImpurityEngine.STRESS_CONDITION_OPTIONS or [])
    except Exception as exc:
        out["errors"].append(f"가혹 조건 목록: {exc}")

    try:
        from engines.p34_name_reaction_process import NameReactionProcess

        out["name_reactions"] = list(NameReactionProcess.get_ui_options() or [])
    except Exception as exc:
        out["errors"].append(f"Name Reaction 목록: {exc}")

    try:
        from engines.p13_chemical_registry import (
            REGISTRY_MASTER_COUNT,
            format_registry_display_name,
            get_merged_category_ids,
        )

        out["master_count"] = REGISTRY_MASTER_COUNT
        for category in ("Solvents", "Acids", "Bases", "Catalysts"):
            ids = list(get_merged_category_ids(category, None) or [])
            out["registry"][category] = [
                {"id": cid, "label": format_registry_display_name(category, cid, dynamic_catalog=None)}
                for cid in ids
            ]
    except Exception as exc:
        out["errors"].append(f"시약 목록: {exc}")

    return out


def build_conditions(form: dict[str, Any], dynamic_catalog: list | None = None) -> dict[str, Any]:
    """폼 입력을 엔진이 받는 process_conditions dict로 바꿉니다.

    engines/ 패키지는 rdchiral을 필요로 합니다. 없으면 경로 B 전체가
    동작하지 않으므로, 여기서 막히더라도 경로 A는 살려 둡니다.
    """
    try:
        from p15_process_pipeline import build_process_conditions

        conditions = build_process_conditions(
            stress_condition=form.get("stress_condition") or "",
            name_reaction=form.get("name_reaction") or "",
            solvents=form.get("solvents"),
            acids=form.get("acids"),
            bases=form.get("bases"),
            catalysts=form.get("catalysts"),
            temperature_c=float(form.get("temperature_c") or 25.0),
            reaction_time_h=float(form.get("reaction_time_h") or 1.0),
            ph=float(form.get("ph") or 7.0),
            dynamic_catalog=dynamic_catalog,
        )
        conditions["_registry_error"] = None
        return conditions
    except Exception as exc:
        return {
            "_registry_error": str(exc),
            "stress_condition": form.get("stress_condition") or "",
            "name_reaction": form.get("name_reaction") or "",
            "reagent": "선택안함 (Not Selected)",
            "reagent_smiles": [],
            "selected_materials": {},
            "process_reagents": {},
            "solvents": list(form.get("solvents") or []),
            "acids": list(form.get("acids") or []),
            "bases": list(form.get("bases") or []),
            "catalysts": list(form.get("catalysts") or []),
            "solvent": "—",
            "solvent_abbrev": "",
            "solvent_bp_c": None,
            "solvent_density_g_ml": None,
            "temperature_c": float(form.get("temperature_c") or 25.0),
            "reaction_time_h": float(form.get("reaction_time_h") or 1.0),
            "ph": float(form.get("ph") or 7.0),
        }


# ---------------------------------------------------------------------------
# 경로 A — 하이브리드 공정 엔진 (p39)
# ---------------------------------------------------------------------------


def run_hybrid(
    reactant1: str,
    reactant2: str,
    product: str,
    *,
    reagent: str = "",
    reference_drug: str = "",
    nitrite: bool = False,
) -> dict[str, Any]:
    out: dict[str, Any] = {
        "rows": [],
        "errors": [],
        "validated_count": 0,
        "theoretical_count": 0,
        "reagent_roles": [],
        "product_risk": None,
        "nitrite_active": bool(nitrite),
        "_raw": None,
    }

    # 시약 구조가 아질산염이면 체크박스와 무관하게 경로를 켭니다.
    try:
        from p30_reaction_engine import is_nitrite_source_smiles

        if reagent and is_nitrite_source_smiles(reagent):
            out["nitrite_active"] = True
    except Exception:
        pass

    try:
        from p24_db_smiles_lookup import resolve_local_df_from_structure

        local_df, _ = resolve_local_df_from_structure(reference_drug)
    except Exception:
        local_df = None

    try:
        from p39_hybrid_process_engine import run_hybrid_process_analysis

        result = run_hybrid_process_analysis(
            reactant1,
            reactant2,
            product or "",
            reagent=reagent,
            reference_drug_smiles=reference_drug,
            nitrite_in_process=out["nitrite_active"],
            local_df=local_df,
        )
    except TypeError as exc:
        if "reagent" not in str(exc):
            out["errors"].append(f"하이브리드 엔진 실패: {exc}")
            return out
        try:
            from p39_hybrid_process_engine import run_hybrid_process_analysis

            result = run_hybrid_process_analysis(
                reactant1,
                reactant2,
                product or "",
                reference_drug_smiles=reference_drug,
                nitrite_in_process=out["nitrite_active"],
                local_df=local_df,
            )
        except Exception as exc2:
            out["errors"].append(f"하이브리드 엔진 실패: {exc2}")
            return out
    except Exception as exc:
        out["errors"].append(f"하이브리드 엔진 실패: {exc}")
        return out

    out["_raw"] = result
    out["errors"].extend(list(getattr(result, "errors", []) or []))
    out["validated_count"] = getattr(result, "validated_count", 0)
    out["theoretical_count"] = getattr(result, "theoretical_count", 0)
    out["reagent_roles"] = list(getattr(result, "reagent_roles", []) or [])
    out["rows"] = [_row_view(r) for r in (getattr(result, "rows", []) or [])]

    risk = getattr(result, "product_risk", None)
    if risk is not None:
        out["product_risk"] = {
            "smiles": getattr(risk, "canon_smiles", ""),
            "formula": getattr(risk, "molecular_formula", ""),
            "mw": getattr(risk, "molecular_weight", 0.0),
            "logp": getattr(risk, "logp", None),
            "tpsa": getattr(risk, "tpsa", None),
            "rotatable_bonds": getattr(risk, "rotatable_bonds", 0),
            "labels": list(getattr(risk, "vulnerability_labels", []) or []),
            "summary": getattr(risk, "summary", ""),
        }

    return out


def _row_view(row: Any) -> dict[str, Any]:
    get = lambda n, d=None: getattr(row, n, d)  # noqa: E731
    return {
        "category": str(get("category", "")),
        "category_label": get("category_label", ""),
        "smiles": get("smiles", ""),
        "rule_label": get("rule_label", ""),
        "verdict": get("verdict", ""),
        "is_validated": bool(get("is_validated", False)),
        "db_name": get("db_name", ""),
        "db_cid": get("db_cid", ""),
        "pubchem_cid": get("pubchem_cid", ""),
        "pubchem_link": get("pubchem_link", ""),
        "pubchem_status": get("pubchem_status", ""),
        "genotox": get("genotox_summary", ""),
        "alert_status": get("alert_status", ""),
        "alert_emoji": get("alert_emoji", ""),
        "ashby_tennant": get("ashby_tennant", ""),
        "m7_class": get("ich_m7_class"),
        "m7_label": get("ich_m7_label", ""),
        "logp": get("logp", ""),
        "stability": get("stability", ""),
        "highlight_atoms": list(get("highlight_atoms", []) or []),
        "degradation_flags": list(get("degradation_flags", []) or []),
    }


# ---------------------------------------------------------------------------
# 경로 B — engines/ 공정 불순물 파이프라인
# ---------------------------------------------------------------------------


def run_impurities(
    conditions: dict[str, Any],
    *,
    reactant1: str,
    reactant2: str,
    product: str,
    reagent_smiles: str = "",
    reference_api: str = "",
) -> dict[str, Any]:
    out: dict[str, Any] = {"records": [], "columns": [], "summary": {}, "error": None, "_raw": []}

    try:
        from p15_process_pipeline import run_process_impurity_prediction

        records = list(
            run_process_impurity_prediction(
                conditions,
                reactant1=reactant1,
                reactant2=reactant2,
                product=product,
                reagent_smiles=reagent_smiles,
                reference_api=reference_api,
            )
            or []
        )
    except Exception as exc:
        out["error"] = f"공정 불순물 예측 실패: {exc}"
        return out

    out["_raw"] = records
    out["records"] = [dict(r) for r in records]

    try:
        from engines.p33_process_impurity_engine import ProcessImpurityEngine

        out["columns"] = list(ProcessImpurityEngine.IMPURITY_RESULT_COLUMNS or [])
    except Exception:
        out["columns"] = []

    verified = sum(1 for r in out["records"] if (r.get("db_match_pct") or 0) >= 85)
    ai_rows = sum(1 for r in out["records"] if r.get("ai_confidence") is not None)
    high_risk = sum(
        1 for r in out["records"] if r.get("ich_m7_class") is not None and r["ich_m7_class"] <= 3
    )
    out["summary"] = {
        "total": len(out["records"]),
        "verified": verified,
        "ai_rows": ai_rows,
        "high_risk": high_risk,
    }

    meta = (out["records"][0] or {}).get("_pipeline_meta") if out["records"] else None
    out["pipeline_meta"] = dict(meta) if isinstance(meta, dict) else None

    return out


def resolve_target(product: str, reference: str, reactant1: str, reactant2: str, records: list):
    """보고서 대상 구조를 정합니다."""
    try:
        from engines.p37_retro_synthesis_process import resolve_report_target_smiles

        return resolve_report_target_smiles(
            product or reference, reactant1, reactant2, records=records
        )
    except Exception:
        return (product or reference or "", {"valid": bool(product or reference)})


def map_stress_for_degradation(stress_condition: str) -> str:
    """공정 모드 가혹 조건을 분해 엔진 조건 이름으로 변환."""
    try:
        from p15_process_pipeline import map_rxn_stress_to_prime

        return map_rxn_stress_to_prime(stress_condition)
    except Exception:
        return "선택 안함"
