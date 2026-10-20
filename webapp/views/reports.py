"""보고서 생성 · 다운로드.

Streamlit 판은 화면을 그릴 때마다 Word/PDF/Excel을 전부 만들어 두고
다운로드 버튼에 바이트를 물렸습니다. 여기서는 사용자가 실제로 누른
포맷만 그 시점에 생성합니다 — 화면이 훨씬 빨라집니다.
"""

from __future__ import annotations

import io
import re
from typing import Any

from flask import Blueprint, Response, abort, flash, redirect, render_template, send_file, url_for

from ..i18n import t as _t
from .. import state
from ..services import single as single_svc

bp = Blueprint("reports", __name__)

MIME = {
    "txt": "text/plain; charset=utf-8",
    "docx": "application/vnd.openxmlformats-officedocument.wordprocessingml.document",
    "pdf": "application/pdf",
    "xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    "json": "application/json",
    "csv": "text/csv; charset=utf-8",
    "html": "text/html; charset=utf-8",
}


def _safe_base(text: str, fallback: str = "PRIME_DS_report") -> str:
    try:
        from p81_report_export import safe_filename_base

        base = safe_filename_base(text)
        if base:
            return base
    except Exception:
        pass
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", (text or "").strip())[:60].strip("_")
    return cleaned or fallback


def _send(data: bytes | str, filename: str, fmt: str) -> Response:
    payload = data.encode("utf-8") if isinstance(data, str) else data
    return send_file(
        io.BytesIO(payload),
        mimetype=MIME.get(fmt, "application/octet-stream"),
        as_attachment=True,
        download_name=filename,
    )


@bp.get("/")
def index():
    return render_template(
        "reports.html",
        nav_title=_t("reports.title"),
        nav_sub=_t("reports.subtitle"),
        has_single=state.has("single_scan"),
        has_process=state.has("process_run"),
        has_degradation=state.has("degradation_run"),
        has_excipient=state.has("excipient_run"),
        has_drugdb=bool((state.get("search_drug_table") or {}).get("rows")),
    )


# ---------------------------------------------------------------------------
# 모드 1 — 단일 화합물
# ---------------------------------------------------------------------------


def _single_ctx_and_text() -> tuple[Any, str, str]:
    scan = state.get("single_scan")
    if not scan:
        abort(404, _t("flash.no_scan"))
    identity = scan.get("identity") or {}
    name = (
        ((scan.get("db") or {}).get("record") or {}).get("name")
        or (scan.get("pubchem") or {}).get("name")
        or identity.get("formula")
        or "compound"
    )
    ctx = single_svc.build_report_context(scan, name=name)
    text = single_svc.render_report_text(scan, name=name)
    return ctx, text, _safe_base(str(name))


@bp.get("/single/preview")
def single_preview():
    try:
        _, text, _ = _single_ctx_and_text()
    except Exception as exc:
        flash(str(exc), "warn")
        return redirect(url_for("single.index"))
    return render_template(
        "report_preview.html",
        nav_title=_t("reports.preview_title"),
        nav_sub=_t("reports.single"),
        report_text=text,
        back=url_for("single.index"),
        downloads=[
            ("docx", "📝 Word", url_for("reports.single_download", fmt="docx")),
            ("pdf", "📕 PDF", url_for("reports.single_download", fmt="pdf")),
            ("txt", "📄 TXT", url_for("reports.single_download", fmt="txt")),
        ],
    )


@bp.get("/single/<fmt>")
def single_download(fmt: str):
    fmt = fmt.lower().strip()
    if fmt not in ("txt", "docx", "pdf"):
        abort(400, "Unsupported format")

    ctx, text, base = _single_ctx_and_text()

    if fmt == "txt":
        return _send(text or "", f"{base}.txt", "txt")

    if ctx is None:
        flash(_t("flash.report_ctx_failed"), "warn")
        return _send(text or "", f"{base}.txt", "txt")

    try:
        if fmt == "docx":
            from p81_report_export import report_to_regulatory_docx

            return _send(report_to_regulatory_docx(ctx), f"{base}.docx", "docx")

        from p81_report_export import report_to_regulatory_pdf

        return _send(report_to_regulatory_pdf(ctx), f"{base}.pdf", "pdf")
    except Exception as exc:
        flash(f"{fmt.upper()} {_t('flash.report_failed')} {exc} {_t('flash.report_failed_txt')}", "danger")
        return redirect(url_for("single.index"))


# ---------------------------------------------------------------------------
# 분해산물
# ---------------------------------------------------------------------------


@bp.get("/degradation/<fmt>")
def degradation_download(fmt: str):
    run = state.get("degradation_run")
    if not run:
        abort(404, _t("flash.no_degradation"))

    fmt = fmt.lower().strip()
    parent = run.get("parent_smiles", "")
    base = _safe_base(parent, "degradation")
    try:
        from p86_degradation_report_export import (
            build_degradation_report_context,
            export_degradation_excel,
            export_degradation_excel_with_structures,
            export_degradation_pdf,
            export_degradation_txt,
            safe_degradation_filename_base,
        )

        base = safe_degradation_filename_base(parent) or base
        ctx = build_degradation_report_context(
            parent_smiles=parent,
            stress_environment=run.get("stress_environment", ""),
            temperature_c=run.get("temperature_c", 25.0),
            products=run.get("_products_raw") or [],
        )
        if fmt == "txt":
            return _send(export_degradation_txt(ctx), f"{base}.txt", "txt")
        if fmt == "pdf":
            return _send(export_degradation_pdf(ctx), f"{base}.pdf", "pdf")
        if fmt == "xlsx":
            return _send(export_degradation_excel(ctx), f"{base}.xlsx", "xlsx")
        if fmt == "xlsx_img":
            return _send(
                export_degradation_excel_with_structures(ctx), f"{base}_structures.xlsx", "xlsx"
            )
    except Exception as exc:
        flash(f"{_t('flash.report_failed')} {exc}", "danger")
        return redirect(url_for("degradation.index"))

    abort(400, "Unsupported format")


@bp.get("/degradation/tree/<fmt>")
def degradation_tree_download(fmt: str):
    tree = state.get("degradation_tree")
    if not tree:
        abort(404, _t("flash.no_tree"))

    fmt = fmt.lower().strip()
    base = _safe_base(tree.get("parent_smiles", ""), "degradation_tree")
    suffix = "_toxic_tree" if tree.get("include_ai") else "_degradation_tree"

    if fmt == "json":
        return _send(tree.get("json_text", "{}"), f"{base}{suffix}.json", "json")
    if fmt == "cytoscape":
        return _send(tree.get("cytoscape_text", "{}"), f"{base}{suffix}_cytoscape.json", "json")
    if fmt == "csv":
        return _send(tree.get("csv_text", ""), f"{base}{suffix}.csv", "csv")
    if fmt == "pdf":
        try:
            import json

            from p87_degradation_tree_report_export import export_degradation_tree_pdf

            data = export_degradation_tree_pdf(
                json.loads(tree.get("json_text") or "{}"),
                parent_smiles=tree.get("parent_smiles", ""),
                max_depth=tree.get("max_depth", 3),
                include_ai=tree.get("include_ai", True),
            )
            return _send(data, f"{base}{suffix}.pdf", "pdf")
        except Exception as exc:
            flash(f"PDF {_t('flash.report_failed')} {exc}", "danger")
            return redirect(url_for("degradation.index"))

    abort(400, "Unsupported format")


# ---------------------------------------------------------------------------
# 부형제
# ---------------------------------------------------------------------------


@bp.get("/excipient/<fmt>")
def excipient_download(fmt: str):
    run = state.get("excipient_run")
    ctx = state.get("excipient_report_ctx")
    if not run or ctx is None:
        abort(404, _t("flash.no_excipient"))

    fmt = fmt.lower().strip()
    base = _safe_base(run.get("api_smiles", ""), "excipient")
    try:
        from p88_excipient_report_export import (
            export_excipient_excel,
            export_excipient_json,
            export_excipient_pdf,
            export_excipient_txt,
            safe_excipient_filename_base,
        )

        base = safe_excipient_filename_base(run.get("api_smiles", "")) or base
        if fmt == "txt":
            return _send(export_excipient_txt(ctx), f"{base}.txt", "txt")
        if fmt == "pdf":
            return _send(export_excipient_pdf(ctx), f"{base}.pdf", "pdf")
        if fmt == "xlsx":
            return _send(export_excipient_excel(ctx), f"{base}.xlsx", "xlsx")
        if fmt == "json":
            return _send(export_excipient_json(ctx), f"{base}.json", "json")
    except Exception as exc:
        flash(f"{_t('flash.report_failed')} {exc}", "danger")
        return redirect(url_for("excipient.index"))

    abort(400, "Unsupported format")


# ---------------------------------------------------------------------------
# 공정 불순물
# ---------------------------------------------------------------------------


@bp.get("/process/<fmt>")
def process_download(fmt: str):
    run = state.get("process_run")
    if not run:
        abort(404, _t("flash.no_process"))

    fmt = fmt.lower().strip()
    records = run.get("_records_raw") or []
    conditions = run.get("conditions") or {}

    try:
        from engines.p90_report_generator import (
            EXCEL_FILE_NAME,
            PDF_FILE_NAME,
            build_process_mode_regulatory_context,
            generate_primeds_process_comparative_pdf,
            generate_primeds_process_excel,
        )

        if fmt == "xlsx":
            if not records:
                flash(_t("flash.no_records_excel"), "warn")
                return redirect(url_for("process.index"))
            data = generate_primeds_process_excel(records=records, process_conditions=conditions)
            if not data:
                flash(_t("flash.empty_excel"), "warn")
                return redirect(url_for("process.index"))
            return _send(data, EXCEL_FILE_NAME, "xlsx")

        target_props = run.get("target_props") or {}

        if fmt == "pdf":
            data = generate_primeds_process_comparative_pdf(
                target_props,
                process_conditions=conditions,
                records=records,
                reactant1=run.get("reactant1"),
                reactant2=run.get("reactant2"),
                reagent_smiles=run.get("reagent_smiles"),
                reference_api=run.get("reference_api"),
            )
            return _send(data, PDF_FILE_NAME, "pdf")

        if fmt == "docx":
            from p80_formal_assessment_report import REPORT_TITLE
            from p81_report_export import report_to_docx

            ctx = build_process_mode_regulatory_context(
                target_props,
                records,
                conditions,
                reactant1=run.get("reactant1"),
                reactant2=run.get("reactant2"),
                reagent_smiles=run.get("reagent_smiles"),
                reference_api=run.get("reference_api"),
            )
            data = report_to_docx("", title=f"{REPORT_TITLE} — Process Mode", ctx=ctx)
            return _send(data, PDF_FILE_NAME.replace(".pdf", ".docx"), "docx")

        if fmt == "json":
            import json

            payload = {
                "generated_at": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
                "inputs": {
                    "reactant1": run.get("reactant1"),
                    "reactant2": run.get("reactant2"),
                    "reagent": run.get("reagent_smiles"),
                    "product": run.get("product"),
                    "reference_api": run.get("reference_api"),
                    "target_smiles": run.get("target_smiles"),
                },
                "conditions": {
                    k: v for k, v in conditions.items()
                    if k not in ("process_reagents", "_registry_error")
                },
                "hybrid_rows": (run.get("hybrid") or {}).get("rows") or [],
                "impurity_records": records or (run.get("impurities") or {}).get("records") or [],
                "degradation": (run.get("degradation") or {}).get("products") or [],
            }
            return _send(
                json.dumps(payload, ensure_ascii=False, indent=2, default=str),
                "PRIME_DS_process_report.json",
                "json",
            )
    except Exception as exc:
        flash(f"{_t('flash.report_failed')} {exc}", "danger")
        return redirect(url_for("process.index"))

    abort(400, "Unsupported format")


# ---------------------------------------------------------------------------
# 의약품 DB 검색 리포트
#
# Streamlit 판에서는 p85_report_download_ui 가 p89 를 이어 주고 있었는데
# 그 모듈이 화면과 함께 사라지면서 p89 가 통째로 끊겨 있었습니다. 다시 잇습니다.
# ---------------------------------------------------------------------------


def _drugdb_table() -> dict[str, Any]:
    table = state.get("search_drug_table")
    if not table or not table.get("rows"):
        abort(404, _t("flash.no_drugdb"))
    return table


@bp.get("/drugdb/preview")
def drugdb_preview():
    try:
        table = _drugdb_table()
        from ..services import search as search_svc

        text = search_svc.drug_db_report_text(table, query=state.get("search_drug_query") or "")
    except Exception as exc:
        flash(str(exc), "warn")
        return redirect(url_for("search.drug"))

    return render_template(
        "report_preview.html",
        nav_title=_t("reports.preview_title"),
        nav_sub=_t("reports.drugdb"),
        report_text=text,
        back=url_for("search.drug"),
        downloads=[
            ("docx", "📝 Word", url_for("reports.drugdb_download", fmt="docx")),
            ("pdf", "📕 PDF", url_for("reports.drugdb_download", fmt="pdf")),
            ("xlsx", "📊 Excel", url_for("reports.drugdb_download", fmt="xlsx")),
            ("txt", "📄 TXT", url_for("reports.drugdb_download", fmt="txt")),
            ("json", "📋 JSON", url_for("reports.drugdb_download", fmt="json")),
        ],
    )


@bp.get("/drugdb/<fmt>")
def drugdb_download(fmt: str):
    table = _drugdb_table()
    fmt = fmt.lower().strip()
    if fmt not in ("txt", "docx", "pdf", "xlsx", "json"):
        abort(400, "Unsupported format")

    try:
        from p89_drug_db_report_export import (
            export_drug_db_docx,
            export_drug_db_excel,
            export_drug_db_json,
            export_drug_db_pdf,
            export_drug_db_txt,
            safe_drug_db_filename_base,
        )

        from ..services import search as search_svc

        ctx = search_svc.drug_db_report_context(
            table, query=state.get("search_drug_query") or ""
        )
        base = safe_drug_db_filename_base(ctx.drug_label) or _safe_base(ctx.drug_label, "drug_db")

        builder = {
            "txt": export_drug_db_txt,
            "docx": export_drug_db_docx,
            "pdf": export_drug_db_pdf,
            "xlsx": export_drug_db_excel,
            "json": export_drug_db_json,
        }[fmt]
        return _send(builder(ctx), f"{base}.{fmt}", fmt)
    except Exception as exc:
        flash(f"{_t('flash.report_failed')} {exc}", "danger")
        return redirect(url_for("search.drug"))
