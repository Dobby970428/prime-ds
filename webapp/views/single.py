"""모드 1 — 단일 화합물 독성 · 물성 스캔."""

from __future__ import annotations

from flask import Blueprint, flash, redirect, render_template, request, url_for

from .. import state
from ..i18n import t
from ..services import search as search_svc
from ..services import single as svc

bp = Blueprint("single", __name__)

SCAN_KEY = "single_scan"


@bp.get("/")
def index():
    scan = state.get(SCAN_KEY)
    return render_template(
        "single.html",
        nav_title=t("single.title"),
        nav_sub=t("single.subtitle"),
        canvas=state.canvas(),
        scan=scan,
        oecd=state.get("single_oecd"),
        prefs=state.display_prefs(),
    )


@bp.post("/scan")
def scan():
    smiles = (request.form.get("single_search") or "").strip()
    reference = (request.form.get("single_reference_drug") or "").strip()
    slow = request.form.get("allow_slow_scan") == "on"

    state.set_canvas("single_search", smiles)
    state.set_canvas("single_reference_drug", reference)

    if not smiles:
        flash(t("flash.draw_structure"), "warn")
        return redirect(url_for("single.index"))

    result = svc.run_scan(smiles, reference_smiles=reference, allow_slow_db_scan=slow)

    if not (result.get("identity") or {}).get("valid"):
        flash(
            (result.get("identity") or {}).get("error") or t("flash.parse_failed"),
            "danger",
        )
        return redirect(url_for("single.index"))

    state.put(SCAN_KEY, result)

    # 검색 이력 기록 — 실패해도 스캔 자체는 유효합니다.
    # p11_structure_search_history는 st.session_state를 직접 쓰므로
    # Flask에서는 쓸 수 없습니다. services.search가 대체합니다.
    try:
        search_svc.record_history(result["identity"]["canonical"], label="독성 스캔")
    except Exception:
        pass

    for err in result.get("errors", [])[:3]:
        flash(err, "warn")

    return redirect(url_for("single.index") + "#result")


@bp.post("/oecd")
def oecd():
    """OECD QSAR (AMBIT2) 조회 — 네트워크 요청이라 별도 액션으로 뒀습니다."""
    from ..services import chem

    scan = state.get(SCAN_KEY)
    smiles = (request.form.get("oecd_smiles") or "").strip()
    if not smiles and scan:
        smiles = (scan.get("identity") or {}).get("canonical", "")

    if not smiles:
        flash(t("flash.scan_first"), "warn")
        return redirect(url_for("single.index"))

    result = chem.oecd_assessment(smiles)
    state.put("single_oecd", result)
    if result.get("error"):
        flash(result["error"], "warn")
    return redirect(url_for("single.index") + "#oecd")


@bp.post("/clear")
def clear():
    state.clear_keys(SCAN_KEY, "single_oecd")
    flash(t("flash.results_cleared"), "ok")
    return redirect(url_for("single.index"))
