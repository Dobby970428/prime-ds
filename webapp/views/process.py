"""모드 2 — 공정 불순물 예측 화면."""

from __future__ import annotations

from flask import Blueprint, flash, redirect, render_template, request, url_for

from .. import state
from ..services import chem
from ..services import degradation as degradation_svc
from ..i18n import t
from ..services import process as svc

bp = Blueprint("process", __name__)

RUN_KEY = "process_run"
FORM_KEY = "process_form"


@bp.get("/")
def index():
    return render_template(
        "process.html",
        nav_title=t("process.title"),
        nav_sub=t("process.subtitle"),
        canvas=state.canvas(),
        options=svc.options(),
        dynamic=state.get("process_dynamic_catalog", []),
        form=state.get(FORM_KEY, {}),
        run=state.get(RUN_KEY),
        stress_none=degradation_svc.STRESS_NONE,
    )


@bp.post("/run")
def run():
    f = request.form

    r1 = (f.get("reactant_1") or "").strip()
    r2 = (f.get("reactant_2") or "").strip()
    reagent = (f.get("reagent") or "").strip()
    product = (f.get("product") or "").strip()
    reference = (f.get("reference_drug") or "").strip()

    for name, value in (
        ("reactant_1", r1), ("reactant_2", r2),
        ("reagent", reagent), ("product", product),
        ("reference_drug", reference),
    ):
        state.set_canvas(name, value)

    form_data = {
        "stress_condition": f.get("stress_condition") or "",
        "name_reaction": f.get("name_reaction") or "",
        "solvents": f.getlist("solvents"),
        "acids": f.getlist("acids"),
        "bases": f.getlist("bases"),
        "catalysts": f.getlist("catalysts"),
        "temperature_c": _float(f.get("temperature_c"), 25.0),
        "reaction_time_h": _float(f.get("reaction_time_h"), 1.0),
        "ph": _float(f.get("ph"), 7.0),
        "nitrite": f.get("nitrite") == "on",
    }
    state.put(FORM_KEY, form_data)

    if not r1 and not r2:
        flash(t("flash.need_starting_material"), "warn")
        return redirect(url_for("process.index"))

    for key, value in (("process.reactant1", r1), ("process.reactant2", r2),
                       ("process.reagent", reagent), ("process.product", product)):
        if value and chem.parse(value) is None:
            flash(f"{t(key)} — {t('flash.parse_failed_named')}", "danger")
            return redirect(url_for("process.index"))

    if not product and not reference:
        flash(t("flash.need_product"), "warn")

    result: dict = {
        "reactant1": r1, "reactant2": r2, "reagent_smiles": reagent,
        "product": product, "reference_api": reference,
    }

    # 경로 A — 하이브리드 공정 엔진
    result["hybrid"] = svc.run_hybrid(
        r1, r2, product, reagent=reagent,
        reference_drug=reference, nitrite=form_data["nitrite"],
    )

    # 경로 B — engines/ 파이프라인
    conditions = svc.build_conditions(
        form_data, dynamic_catalog=state.get("process_dynamic_catalog", [])
    )
    result["conditions"] = conditions
    result["impurities"] = svc.run_impurities(
        conditions, reactant1=r1, reactant2=r2, product=product,
        reagent_smiles=reagent, reference_api=reference,
    )

    target_smiles, target_props = svc.resolve_target(
        product, reference, r1, r2, result["impurities"].get("_raw") or []
    )
    result["target_smiles"] = target_smiles
    result["target_props"] = target_props
    result["_records_raw"] = result["impurities"].pop("_raw", [])

    # 분해산물 — 가혹 조건을 골랐을 때만
    prime_stress = svc.map_stress_for_degradation(conditions["stress_condition"])
    result["degradation_stress"] = prime_stress
    if prime_stress and prime_stress != degradation_svc.STRESS_NONE and target_smiles:
        result["degradation"] = degradation_svc.predict(
            target_smiles, prime_stress, conditions["temperature_c"]
        )
    else:
        result["degradation"] = None

    state.put(RUN_KEY, result)

    for msg in (result["hybrid"].get("errors") or [])[:3]:
        flash(msg, "warn")
    if result["impurities"].get("error"):
        flash(result["impurities"]["error"], "danger")

    return redirect(url_for("process.index") + "#result")


@bp.post("/clear")
def clear():
    state.clear_keys(RUN_KEY)
    flash(t("flash.results_cleared"), "ok")
    return redirect(url_for("process.index"))


def _float(raw, default):
    try:
        return float(str(raw).strip())
    except (TypeError, ValueError):
        return default
