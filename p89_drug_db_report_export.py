"""의약품명 DB 검색 — 관련 유연물질 리포트 (TXT · Word · PDF · Excel · JSON)."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Any

from fpdf import FPDF
from openpyxl import Workbook
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
# 2026-09-03 수정: Chem/Descriptors 직접 사용처(_exact_mw)가 p82 로 이동해 import 제거.

from p80_formal_assessment_report import DOCUMENT_VERSION, PLATFORM_FULL, PLATFORM_NAME
from p81_report_export import (
    _format_report_value,
    _pdf_add_structure_figure,
    _pdf_blank_lines,
    _pdf_kv_rows,
    _pdf_output_bytes,
    _pdf_sanitize,
    _register_pdf_fonts,
    _pdf_section,
    _pdf_write_wrapped,
    report_to_docx,
    safe_filename_base,
)
from p82_pdf_typesetting import (
    MOLECULAR_WEIGHT_LABEL,
    average_mw_text,
    report_control_header_rows,
    structure_png_bytes,
)
from p84_structure_images import DEFAULT_REPORT_IMAGE_SIZE

_DISCLAIMER = (
    "REGULATORY DISCLAIMER: In-house database impurity records support ICH M7 screening. "
    "Computed Ashby-Tennant alerts and ICH M7 class assignments are in silico predictions "
    "and do not replace validated analytical or toxicological data."
)


@dataclass
class DrugDbReportContext:
    drug_label: str
    search_query: str
    source_file: str
    total_records: int
    displayed_count: int
    rows: list[dict[str, Any]]
    report_id: str = field(default_factory=lambda: datetime.now().strftime("DRG-%Y%m%d-%H%M%S"))
    generated_at: str = field(default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC"))


# 2026-09-03 수정 (E-7 / F-4): p88 과 축자 동일했고, 또한 p86(MolWt, 평균)과
#   분자량 정의가 달라 라니티딘이 314.14 / 314.41 로 갈렸습니다.
#   규제 문서 관례인 평균 분자량(MolWt)으로 통일하고 p82 로 위임합니다.
def _exact_mw(smiles: str) -> str:
    return average_mw_text(smiles)


# 2026-09-03 수정 (F-4): p86/p87 판본과 갈려 있던 구조 렌더링을 p82 로 통합.
def _structure_png_bytes(smiles: str, *, size: str = DEFAULT_REPORT_IMAGE_SIZE) -> bytes | None:
    return structure_png_bytes(smiles, size=size)


def _pdf_line(text: str) -> str:
    return _pdf_sanitize(str(text or ""))


def build_drug_db_report_context(
    *,
    drug_label: str,
    search_query: str,
    source_file: str,
    total_records: int,
    rows: list[dict[str, Any]],
) -> DrugDbReportContext:
    return DrugDbReportContext(
        drug_label=(drug_label or "").strip() or "Unknown drug",
        search_query=(search_query or "").strip(),
        source_file=str(source_file or ""),
        total_records=int(total_records),
        displayed_count=len(rows),
        rows=list(rows),
    )


def _m7_class_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        key = str(row.get("ICH M7") or "—")
        counts[key] = counts.get(key, 0) + 1
    return counts


def _document_control_rows(ctx: DrugDbReportContext) -> list[tuple[str, str]]:
    m7_counts = _m7_class_counts(ctx.rows)
    m7_summary = ", ".join(f"{k}: {v}" for k, v in sorted(m7_counts.items())) or "—"
    # 2026-09-03 수정 (F-4): 앞 5행은 p88 과 축자 동일했습니다 → p82 공용 헬퍼.
    return [
        *report_control_header_rows(
            report_id=ctx.report_id,
            generated=ctx.generated_at,
            platform=f"{PLATFORM_NAME} — {PLATFORM_FULL}",
            document_version=DOCUMENT_VERSION,
            assessment_type="In-house drug impurity database screening report",
        ),
        ("Search query", ctx.search_query or "—"),
        ("Drug / product", ctx.drug_label),
        ("Source database file", Path(ctx.source_file).name if ctx.source_file else "—"),
        ("Total records in file", str(ctx.total_records)),
        ("Impurity entries in report", str(ctx.displayed_count)),
        ("ICH M7 class distribution", m7_summary),
    ]


def _impurity_evaluation_rows(row: dict[str, Any], idx: int) -> list[tuple[str, str]]:
    smi = str(row.get("SMILES") or "")
    alerts = str(row.get("Ashby-Tennant") or "—")
    return [
        ("Impurity ID", f"IMP-{idx:03d}"),
        ("Compound name", _format_report_value(row.get("Name"))),
        ("SMILES", smi or "N/A"),
        ("PubChem CID", _format_report_value(row.get("CID"))),
        # 2026-09-03 수정: 평균 분자량으로 통일 + 라벨 명기.
        ("Molecular weight", f"{_exact_mw(smi)} {MOLECULAR_WEIGHT_LABEL}" if smi else "N/A"),
        ("Structure depiction", f"Figure {idx + 1} — 2D structure"),
        ("ICH M7 hazard class", _format_report_value(row.get("ICH M7"))),
        ("LogP (computed)", _format_report_value(row.get("LogP"))),
        ("Process stability", _format_report_value(row.get("공정 안정성"))),
        ("Ashby-Tennant alerts", alerts),
        (
            "Screening outcome",
            "POSITIVE structural alert" if alerts and alerts != "—" else "No Ashby-Tennant alert",
        ),
        (
            "Recommended action",
            "Confirm identity by LC-MS/NMR; apply ICH M7 control strategy per hazard class",
        ),
    ]


def _executive_text_block(row: dict[str, Any], idx: int) -> list[str]:
    name = str(row.get("Name") or f"Impurity {idx}")
    lines = [
        "",
        "╔" + "═" * 86 + "╗",
        f"║  IMPURITY {idx:03d} — {name[:70]:<70}║",
        "╚" + "═" * 86 + "╝",
    ]
    for label, val in _impurity_evaluation_rows(row, idx):
        if label == "Structure depiction":
            continue
        lines.append(f"  ▶ {label:<32}: {val}")
    return lines


def generate_drug_db_report_text(ctx: DrugDbReportContext) -> str:
    lines: list[str] = [
        "=" * 88,
        f"{PLATFORM_NAME} — DRUG DATABASE IMPURITY SCREENING REPORT",
        PLATFORM_FULL,
        f"Document Version: {DOCUMENT_VERSION}",
        "=" * 88,
        "",
        "DOCUMENT CONTROL",
        "-" * 88,
    ]
    for label, val in _document_control_rows(ctx):
        lines.append(f"  {label:<32}: {val}")

    lines.extend(
        [
            "",
            f"1. DRUG PRODUCT — {ctx.drug_label}",
            "-" * 88,
            f"  Search query                    : {ctx.search_query}",
            f"  Database source                 : {Path(ctx.source_file).name if ctx.source_file else '—'}",
            f"  Total impurity records on file  : {ctx.total_records}",
            f"  Entries included in this report : {ctx.displayed_count}",
            "",
            f"2. IMPURITY APPENDIX ({ctx.displayed_count} entries)",
            "-" * 88,
        ]
    )
    if not ctx.rows:
        lines.append("  No parseable impurity records with valid SMILES.")
    else:
        for idx, row in enumerate(ctx.rows, start=1):
            lines.extend(_executive_text_block(row, idx))

    lines.extend(
        [
            "",
            "3. SUMMARY",
            f"   ICH M7 distribution: {_m7_class_counts(ctx.rows)}",
            "",
            "4. DISCLAIMER",
            f"   {_DISCLAIMER}",
            "",
        ]
    )
    return "\n".join(lines)


def export_drug_db_txt(ctx: DrugDbReportContext) -> bytes:
    return generate_drug_db_report_text(ctx).encode("utf-8")


def export_drug_db_docx(ctx: DrugDbReportContext) -> bytes:
    title = f"Drug Database Impurity Report — {ctx.drug_label}"
    return report_to_docx(generate_drug_db_report_text(ctx), title=title, ctx=None)


def export_drug_db_json(ctx: DrugDbReportContext) -> bytes:
    payload = {
        "report_id": ctx.report_id,
        "generated_at": ctx.generated_at,
        "document_version": DOCUMENT_VERSION,
        "drug_label": ctx.drug_label,
        "search_query": ctx.search_query,
        "source_file": ctx.source_file,
        "total_records": ctx.total_records,
        "displayed_count": ctx.displayed_count,
        "ich_m7_summary": _m7_class_counts(ctx.rows),
        "impurities": [
            {**row, "evaluation": dict(_impurity_evaluation_rows(row, i))}
            for i, row in enumerate(ctx.rows, start=1)
        ],
        "disclaimer": _DISCLAIMER,
    }
    return json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")


def export_drug_db_excel(ctx: DrugDbReportContext) -> bytes:
    wb = Workbook()
    header_fill = PatternFill("solid", fgColor="E3F2FD")
    bold = Font(bold=True)

    ws_sum = wb.active
    ws_sum.title = "Summary"
    for r, (k, v) in enumerate(_document_control_rows(ctx), start=1):
        ws_sum.cell(row=r, column=1, value=k).font = bold
        ws_sum.cell(row=r, column=1).fill = header_fill
        ws_sum.cell(row=r, column=2, value=v)
        ws_sum.cell(row=r, column=2).alignment = Alignment(wrap_text=True, vertical="top")
    ws_sum.column_dimensions["A"].width = 34
    ws_sum.column_dimensions["B"].width = 72

    ws = wb.create_sheet("Impurities")
    headers = [
        "No",
        "Name",
        "SMILES",
        "CID",
        "ICH M7",
        "LogP",
        "Process_Stability",
        "Ashby_Tennant",
        "MW_g_mol",
        "2D_Structure",
    ]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=1, column=c, value=h)
        cell.font = bold
        cell.fill = header_fill

    for row_idx, (idx, row) in enumerate(zip(range(1, len(ctx.rows) + 1), ctx.rows), start=2):
        smi = str(row.get("SMILES") or "")
        data = [
            idx,
            row.get("Name"),
            smi,
            row.get("CID"),
            row.get("ICH M7"),
            row.get("LogP"),
            row.get("공정 안정성"),
            row.get("Ashby-Tennant"),
            _exact_mw(smi) if smi else "N/A",
            "",
        ]
        for col_idx, val in enumerate(data, start=1):
            cell = ws.cell(row=row_idx, column=col_idx, value=val)
            cell.alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[row_idx].height = 100
        png = _structure_png_bytes(smi, size="220x165")
        if png:
            img_buf = BytesIO(png)
            img_buf.seek(0)
            img = XLImage(img_buf)
            img.width = 150
            img.height = 112
            ws.add_image(img, f"J{row_idx}")

    widths = [6, 28, 36, 12, 12, 10, 18, 28, 12, 22]
    for col_idx, width in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(col_idx)].width = width

    buf = BytesIO()
    wb.save(buf)
    return buf.getvalue()


def export_drug_db_pdf(ctx: DrugDbReportContext) -> bytes:
    pdf = FPDF()
    pdf.set_margins(18, 18, 18)
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.add_page()
    w = pdf.epw
    family, family_has_bold = _register_pdf_fonts(pdf)
    kv = lambda rows: _pdf_kv_rows(pdf, family, rows, w, family_has_bold=family_has_bold)

    _pdf_write_wrapped(
        pdf,
        family,
        _pdf_line(f"{PLATFORM_NAME} — DRUG DATABASE IMPURITY SCREENING REPORT"),
        w,
        bold=True,
        size=14,
    )
    _pdf_write_wrapped(
        pdf,
        family,
        _pdf_line(f"Product: {ctx.drug_label}"),
        w,
        size=10,
    )
    _pdf_blank_lines(pdf, 1)

    _pdf_section(pdf, family, "DOCUMENT CONTROL", w)
    kv(_document_control_rows(ctx))

    _pdf_section(pdf, family, f"1. DRUG PRODUCT — {ctx.drug_label}", w)
    kv(
        [
            ("Search query", ctx.search_query or "—"),
            ("Source file", Path(ctx.source_file).name if ctx.source_file else "—"),
            ("Total records", str(ctx.total_records)),
            ("Report entries", str(ctx.displayed_count)),
        ]
    )

    _pdf_section(pdf, family, f"2. IMPURITY APPENDIX ({ctx.displayed_count})", w)
    if not ctx.rows:
        _pdf_write_wrapped(pdf, family, _pdf_line("No valid impurity entries."), w)
    else:
        for idx, row in enumerate(ctx.rows, start=1):
            name = str(row.get("Name") or f"Impurity {idx}")
            _pdf_section(
                pdf,
                family,
                _pdf_line(f"IMPURITY {idx:03d} — {name}"),
                w,
                pre_gap_lines=2,
            )
            smi = str(row.get("SMILES") or "")
            png = _structure_png_bytes(smi, size="400x300")
            if png:
                _pdf_add_structure_figure(
                    pdf,
                    family,
                    png,
                    w,
                    caption=_pdf_line(f"Figure {idx}: {name}"),
                )
            kv(_impurity_evaluation_rows(row, idx))

    _pdf_section(pdf, family, "3. REGULATORY DISCLAIMER", w)
    _pdf_write_wrapped(pdf, family, _pdf_line(_DISCLAIMER), w)
    return _pdf_output_bytes(pdf)


def safe_drug_db_filename_base(drug_label: str) -> str:
    base = safe_filename_base(drug_label, fallback="drug_db")
    return f"PRIME_DS_DrugDB_{base}"

