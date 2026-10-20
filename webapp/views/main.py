"""홈 · 설정 · 엔진 정보."""

from __future__ import annotations

from flask import Blueprint, flash, redirect, render_template, request, url_for

from .. import state
from ..i18n import t
from . import context

bp = Blueprint("main", __name__)


@bp.get("/")
def index():
    return redirect(url_for("single.index"))


@bp.get("/lang/<code>")
def set_language(code: str):
    """언어 전환. 보던 화면으로 되돌아갑니다."""
    from ..i18n import set_lang

    set_lang(code)
    return redirect(request.referrer or url_for("single.index"))


@bp.route("/settings", methods=["GET", "POST"])
def settings():
    if request.method == "POST":
        action = (request.form.get("action") or "").strip()

        if action == "display":
            state.update_display_prefs(
                max_show=_int(request.form.get("max_show"), 12),
                cols=_int(request.form.get("cols"), 3),
                img_size=(request.form.get("img_size") or "300x300").strip(),
                show_smiles=request.form.get("show_smiles") == "on",
            )
            flash(t("flash.display_saved"), "ok")
            return redirect(url_for("main.settings"))

        if action == "reindex":
            try:
                from p22_drug_data import DATA_DIR
                from p26_import_index import import_all

                ok, skip, fail = import_all(DATA_DIR, verbose=False)
                context.invalidate()
                flash(f"{t('flash.index_refreshed')} {t('flash.registered')} {ok:,} · {t('flash.skipped')} {skip:,} · {t('flash.failed')} {fail:,}", "ok")
            except Exception as exc:
                flash(f"{t('flash.index_failed')} {exc}", "danger")
            return redirect(url_for("main.settings"))

        if action == "build_primary":
            try:
                from p23_chem_search_engine import (
                    build_drug_primary_index,
                    clear_substruct_library_cache,
                )

                ins, total = build_drug_primary_index()
                clear_substruct_library_cache()
                context.invalidate()
                flash(f"{t('flash.fp_built')} {t('flash.new')} {ins:,} · {t('about.total')} {total:,}", "ok")
            except Exception as exc:
                flash(f"{t('flash.fp_failed')} {exc}", "danger")
            return redirect(url_for("main.settings"))

        if action == "build_substruct":
            try:
                from p23_chem_search_engine import build_substruct_library

                count, err = build_substruct_library()
                if err:
                    flash(f"{t('flash.substruct_error')} {err}", "warn")
                else:
                    flash(f"{t('flash.substruct_built')} {count:,}", "ok")
            except Exception as exc:
                flash(f"{t('flash.substruct_failed')} {exc}", "danger")
            return redirect(url_for("main.settings"))

        if action == "build_compound":
            try:
                from p23_chem_search_engine import build_compound_index

                ins, rows, files = build_compound_index()
                context.invalidate()
                flash(
                    f"{t('flash.exact_built')} {t('flash.new')} {ins:,} · "
                    f"{t('flash.scanned_rows')} {rows:,} · {t('flash.files')} {files:,}",
                    "ok",
                )
            except Exception as exc:
                flash(f"{t('flash.exact_failed')} {exc}", "danger")
            return redirect(url_for("main.settings"))

        if action == "reset_session":
            state.reset()
            flash(t("flash.session_reset"), "ok")
            return redirect(url_for("main.settings"))

    return render_template(
        "settings.html",
        nav_title=t("settings.title"),
        nav_sub=t("settings.subtitle"),
        prefs=state.display_prefs(),
        index_status=_index_status(),
    )


@bp.get("/about")
def about():
    return render_template(
        "about.html",
        nav_title=t("about.title"),
        nav_sub=t("about.subtitle"),
        engine=_engine_info(),
    )


# ---------------------------------------------------------------------------


def _int(raw: str | None, default: int) -> int:
    try:
        return int(str(raw).strip())
    except (TypeError, ValueError):
        return default


def _index_status() -> dict:
    info: dict = {}
    try:
        from p23_chem_search_engine import (
            count_compounds,
            count_drug_primary,
            search_index_status,
            substruct_library_ready,
        )

        info["primary"] = count_drug_primary()
        info["compounds"] = count_compounds()
        info["substruct_ready"] = substruct_library_ready()
        try:
            info["status"] = search_index_status()
        except Exception:
            info["status"] = None
    except Exception as exc:
        info["error"] = str(exc)

    try:
        from p20_db_index import count_indexed

        info["indexed"] = count_indexed()
    except Exception:
        info["indexed"] = 0

    return info


def _engine_info() -> dict:
    info: dict = {"layers": [], "notes": []}

    try:
        from p60_m7_empirical_db import empirical_db_stats

        info["empirical"] = empirical_db_stats()
    except Exception as exc:
        info["empirical"] = {"error": str(exc)}

    try:
        import os
        from pathlib import Path

        root = Path(__file__).resolve().parent.parent.parent
        ames = root / "models" / "ames_ecfp_rf.pkl"
        tox21 = root / "models" / "tox21_ecfp_rf_sklearn.pkl"
        forced = os.environ.get("PRIME_LAYER_C_MODEL", "").strip().lower()
        info["layer_c"] = {
            "ames_available": ames.exists(),
            "tox21_available": tox21.exists(),
            "forced": forced or None,
            "active": "Tox21 SR-p53" if forced in ("tox21", "tox") else "Ames",
        }
    except Exception as exc:
        info["layer_c"] = {"error": str(exc)}

    try:
        from p02_degradation_constants import ENGINE_CACHE_VERSION

        info["degradation_version"] = ENGINE_CACHE_VERSION
    except Exception:
        info["degradation_version"] = "?"

    try:
        from p71_excipient_opensource_registry import opensource_catalog_stats

        info["excipients"] = opensource_catalog_stats()
    except Exception as exc:
        info["excipients"] = {"error": str(exc)}

    return info
