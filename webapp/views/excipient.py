"""부형제-API 상호작용 화면."""

from __future__ import annotations

from flask import Blueprint, flash, redirect, render_template, request, url_for

from .. import state
from ..i18n import current_lang, t
from ..services import excipient as svc

bp = Blueprint("excipient", __name__)

RUN_KEY = "excipient_run"
CTX_KEY = "excipient_report_ctx"
SEL_KEY = "excipient_selection"


@bp.get("/")
def index():
    lang = current_lang()
    cat = svc.catalog(lang)
    selected = state.get(SEL_KEY, [])
    query = request.args.get("q", "").strip()
    categories = request.args.getlist("cat")

    visible = svc.filter_keys(
        cat["choices"], categories=categories, query=query, lang=lang
    )
    grouped = svc.group_by_category(visible, cat)

    return render_template(
        "excipient.html",
        nav_title=t("excipient.title"),
        nav_sub=t("excipient.subtitle"),
        canvas=state.canvas(),
        catalog=cat,
        visible=visible,
        grouped=grouped,
        selected=selected,
        query=query,
        active_cats=categories,
        run=state.get(RUN_KEY),
        nitrite=svc.nitrite_profiles(lang),
        mode=state.get("excipient_mode", "binary"),
        ratio=state.get("excipient_ratio", [1.0, 1.0]),
    )


@bp.post("/select")
def select():
    action = (request.form.get("action") or "").strip()
    selected = list(state.get(SEL_KEY, []))

    if action == "add":
        for key in request.form.getlist("key"):
            if key and key not in selected:
                selected.append(key)
    elif action == "remove":
        drop = set(request.form.getlist("key"))
        selected = [k for k in selected if k not in drop]
    elif action == "clear":
        selected = []
    elif action == "set":
        selected = [k for k in request.form.getlist("key") if k]

    state.put(SEL_KEY, selected)
    return redirect(_back())


@bp.post("/run")
def run():
    api_smiles = (request.form.get("excipient_api") or "").strip()
    state.set_canvas("excipient_api", api_smiles)

    keys = request.form.getlist("key") or list(state.get(SEL_KEY, []))
    state.put(SEL_KEY, keys)

    if not api_smiles:
        flash(t("flash.draw_structure"), "warn")
        return redirect(url_for("excipient.index"))
    if not keys:
        flash(t("flash.select_excipient"), "warn")
        return redirect(url_for("excipient.index"))

    mode = (request.form.get("mode") or "binary").strip()
    ratio_api = _float(request.form.get("ratio_api"), 1.0)
    ratio_exc = _float(request.form.get("ratio_exc"), 1.0)
    state.put("excipient_mode", mode)
    state.put("excipient_ratio", [ratio_api, ratio_exc])

    if mode == "matrix":
        result = svc.run(api_smiles, keys, lang=current_lang())
    else:
        result = svc.run_binary(
            api_smiles, keys, ratio=(ratio_api, ratio_exc), lang=current_lang()
        )

    ctx = result.pop("_ctx", None)
    result.pop("_all_hits", None)
    state.put(RUN_KEY, result)
    state.put(CTX_KEY, ctx)

    for msg in result.get("errors", []):
        flash(msg, "warn")

    empty = (
        not result.get("hits")
        if mode == "matrix"
        else not (result.get("summary") or {}).get("total_hits")
    )
    if empty and not result.get("errors"):
        flash(t("flash.no_excipient_hits"), "info")

    return redirect(url_for("excipient.index") + "#result")


def _float(raw, default):
    try:
        return float(str(raw).strip())
    except (TypeError, ValueError):
        return default


@bp.post("/clear")
def clear():
    state.clear_keys(RUN_KEY, CTX_KEY)
    flash(t("flash.results_cleared"), "ok")
    return redirect(url_for("excipient.index"))


def _back() -> str:
    ref = request.form.get("next") or request.referrer
    return ref or url_for("excipient.index")
