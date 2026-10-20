"""분해산물 예측 화면."""

from __future__ import annotations

import json

from flask import Blueprint, flash, redirect, render_template, request, url_for

from .. import state
from ..i18n import t
from ..services import degradation as svc

bp = Blueprint("degradation", __name__)

RUN_KEY = "degradation_run"
TREE_KEY = "degradation_tree"
API_KEY = "degradation_api"


@bp.get("/")
def index():
    run = state.get(RUN_KEY)
    return render_template(
        "degradation.html",
        nav_title=t("degradation.title"),
        nav_sub=t("degradation.subtitle"),
        canvas=state.canvas(),
        stress_choices=svc.stress_choices(),
        run=run,
        summary=svc.summary(run) if run else None,
        tree=state.get(TREE_KEY),
        api=state.get(API_KEY),
        form=state.get("degradation_form", {}),
    )


@bp.post("/run")
def run():
    smiles = (request.form.get("degrad_parent") or "").strip()
    stress = (request.form.get("stress_env") or svc.STRESS_NONE).strip()
    temp = _float(request.form.get("stress_temp"), 25.0)

    state.set_canvas("degrad_parent", smiles)
    state.put("degradation_form", {"stress_env": stress, "stress_temp": temp})

    if not smiles:
        flash(t("flash.draw_structure"), "warn")
        return redirect(url_for("degradation.index"))

    if stress == svc.STRESS_NONE:
        flash(t("flash.select_stress"), "warn")
        return redirect(url_for("degradation.index"))

    result = svc.predict(smiles, stress, temp)
    state.put(RUN_KEY, result)
    state.pop(TREE_KEY)
    state.pop(API_KEY)

    for msg in result.get("errors", []):
        flash(msg, "danger")
    for msg in result.get("warnings", []):
        flash(msg, "warn")

    return redirect(url_for("degradation.index") + "#result")


@bp.post("/tree")
def tree():
    run_data = state.get(RUN_KEY)
    parent = (request.form.get("tree_parent") or "").strip()
    if not parent and run_data:
        parent = run_data.get("parent_smiles", "")
    if not parent:
        flash(t("flash.run_degradation_first"), "warn")
        return redirect(url_for("degradation.index"))

    depth = _int(request.form.get("tree_depth"), 3)
    include_ai = request.form.get("tree_include_ai") == "on"

    result = svc.build_tree(parent, max_depth=depth, include_ai=include_ai)
    state.put(TREE_KEY, result)

    if result.get("error"):
        flash(result["error"], "danger")

    return redirect(url_for("degradation.index") + "#tree")


@bp.post("/api")
def api():
    run_data = state.get(RUN_KEY)
    parent = (request.form.get("api_parent") or "").strip()
    if not parent and run_data:
        parent = run_data.get("parent_smiles", "")
    if not parent:
        flash(t("flash.structure_required"), "warn")
        return redirect(url_for("degradation.index"))

    result = svc.backend_core_api(
        parent,
        formulation_ph=_float(request.form.get("formulation_ph"), 7.4),
        max_depth=_int(request.form.get("api_depth"), 3),
        include_structure_b64=request.form.get("include_b64") == "on",
    )
    if result.get("ok"):
        try:
            result["json_text"] = json.dumps(
                result["payload"], ensure_ascii=False, indent=2, default=str
            )
        except Exception:
            result["json_text"] = "{}"
    else:
        flash(result.get("error") or t("api.failed"), "danger")

    state.put(API_KEY, result)
    return redirect(url_for("degradation.index") + "#api")


@bp.post("/clear")
def clear():
    state.clear_keys(RUN_KEY, TREE_KEY, API_KEY)
    flash(t("flash.results_cleared"), "ok")
    return redirect(url_for("degradation.index"))


def _float(raw, default):
    try:
        return float(str(raw).strip())
    except (TypeError, ValueError):
        return default


def _int(raw, default):
    try:
        return int(str(raw).strip())
    except (TypeError, ValueError):
        return default
