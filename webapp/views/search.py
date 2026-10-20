"""검색 화면 — 의약품명 · 구조."""

from __future__ import annotations

from flask import Blueprint, flash, redirect, render_template, request, url_for

from .. import state
from ..i18n import t
from ..services import search as svc

bp = Blueprint("search", __name__)


@bp.get("/drug")
def drug():
    query = (request.args.get("q") or "").strip()
    pick = (request.args.get("pick") or "").strip()
    pick_name = (request.args.get("name") or "").strip()

    result = svc.drug_search(query) if query else None
    table = None
    if pick:
        table = svc.drug_table(pick, pick_name, max_show=state.display_prefs()["max_show"])
        if table.get("error"):
            flash(table["error"], "warn")
        # 리포트 다운로드가 같은 표를 다시 계산하지 않도록 세션에 둡니다.
        if table.get("rows"):
            state.put("search_drug_table", table)
            state.put("search_drug_query", query)

    return render_template(
        "search_drug.html",
        nav_title=t("search.drug_title"),
        nav_sub=t("search.drug_subtitle"),
        query=query,
        result=result,
        table=table,
        pick=pick,
        prefs=state.display_prefs(),
    )


@bp.get("/structure")
def structure():
    return render_template(
        "search_structure.html",
        nav_title=t("search.structure_title"),
        nav_sub=t("search.structure_subtitle"),
        canvas=state.canvas(),
        sub=state.get("search_substructure"),
        sim=state.get("search_similarity"),
        history=svc.history(),
        prefs=state.display_prefs(),
    )


@bp.post("/structure/substructure")
def run_substructure():
    smiles = (request.form.get("search_query") or "").strip()
    mode = "deep" if (request.form.get("mode") or "") == "deep" else "fast"
    state.set_canvas("search_query", smiles)

    if not smiles:
        flash(t("flash.draw_structure"), "warn")
        return redirect(url_for("search.structure"))

    result = svc.substructure(smiles, mode=mode)
    state.put("search_substructure", result)
    if result.get("error"):
        flash(result["error"], "warn")
    return redirect(url_for("search.structure") + "#sub")


@bp.post("/structure/similarity")
def run_similarity():
    smiles = (request.form.get("search_query") or "").strip()
    state.set_canvas("search_query", smiles)

    if not smiles:
        flash("먼저 구조를 그려 주세요.", "warn")
        return redirect(url_for("search.structure"))

    result = svc.similarity(smiles, top_n=3)
    state.put("search_similarity", result)
    if result.get("error"):
        flash(result["error"], "warn")
    return redirect(url_for("search.structure") + "#sim")


@bp.post("/structure/clear")
def clear_structure():
    state.clear_keys("search_substructure", "search_similarity")
    return redirect(url_for("search.structure"))


@bp.post("/history/clear")
def clear_history():
    svc.clear_history()
    flash(t("flash.history_cleared"), "ok")
    return redirect(request.referrer or url_for("search.structure"))
