"""부형제-API 상호작용 · NDSRI/NDMA 리포트 (TXT · Excel · PDF · JSON)."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO
from typing import Any

from fpdf import FPDF
from openpyxl import Workbook
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from rdkit import Chem
# 2026-09-03 수정: Descriptors 직접 사용처(_exact_mw)가 p82 로 이동해 import 제거.

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
)
from p82_pdf_typesetting import (
    MOLECULAR_WEIGHT_LABEL,
    average_mw_text,
    report_control_header_rows,
    structure_png_bytes,
)
from p84_structure_images import DEFAULT_REPORT_IMAGE_SIZE

_DISCLAIMER = (
    "REGULATORY DISCLAIMER: In silico excipient-API interaction and nitrosamine (NDSRI/NDMA) "
    "screening supports formulation risk assessment. Results do not replace validated "
    "analytical methods (LC-MS), nitrosamine confirmatory testing, or qualified toxicological review."
)

_RISK_OUTCOME = {
    "Matrix": "Matrix interaction — Maillard / Formylation / formulation artifact (ICH Q1A/Q8 context)",
    "NDSRI": "NDSRI nitrosamine risk — secondary/tertiary amine + nitrosating excipient environment",
    "NDMA": "NDMA-specific structural alert — N-nitrosodimethylamine topology match",
}


@dataclass
class ExcipientReportContext:
    api_smiles: str
    api_canonical: str
    excipient_keys: list[str]
    excipient_labels: list[str]
    hits: list[dict[str, Any]]
    amine_secondary: int = 0
    amine_tertiary: int = 0
    report_id: str = field(default_factory=lambda: datetime.now().strftime("EXC-%Y%m%d-%H%M%S"))
    generated_at: str = field(default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC"))


def _canon(smiles: str) -> str:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return (smiles or "").strip()
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return (smiles or "").strip()


# 2026-09-03 수정 (E-7 / F-4): 분자량 정의가 보고서마다 달랐습니다.
#   p88·p89 는 Descriptors.ExactMolWt(단동위원소), p86 은 Descriptors.MolWt(평균)
#   이라 라니티딘이 314.14 / 314.41 로 갈렸습니다. 규제 문서 관례인 **평균
#   분자량**(MolWt)으로 통일하고 라벨에 어느 쪽인지 명기합니다.
#   함수명 `_exact_mw` 는 p88·p89 내부에서 쓰이므로 유지하되 p82 로 위임합니다.
def _exact_mw(smiles: str) -> str:
    return average_mw_text(smiles)


# 2026-09-03 수정 (F-4): p86/p87 판본과 갈려 있던 구조 렌더링을 p82 로 통합.
def _structure_png_bytes(smiles: str, *, size: str = DEFAULT_REPORT_IMAGE_SIZE) -> bytes | None:
    return structure_png_bytes(smiles, size=size)


def _pdf_line(text: str) -> str:
    return _pdf_sanitize(str(text or "").replace("\u207a", "+").replace("\u207b", "-"))


def build_excipient_report_context(
    api_smiles: str,
    excipient_keys: list[str],
    hits: list[dict[str, Any]],
    *,
    excipient_label_map: dict[str, str] | None = None,
    amine_secondary: int = 0,
    amine_tertiary: int = 0,
) -> ExcipientReportContext:
    labels_map = excipient_label_map or {}
    labels = [labels_map.get(k, k) for k in excipient_keys]
    return ExcipientReportContext(
        api_smiles=api_smiles.strip(),
        api_canonical=_canon(api_smiles),
        excipient_keys=list(excipient_keys),
        excipient_labels=labels,
        hits=list(hits),
        amine_secondary=amine_secondary,
        amine_tertiary=amine_tertiary,
    )


def _summary_counts(ctx: ExcipientReportContext) -> tuple[int, int, int]:
    matrix = sum(1 for h in ctx.hits if str(h.get("risk_class")) == "Matrix")
    ndsri = sum(1 for h in ctx.hits if str(h.get("risk_class")) == "NDSRI")
    ndma = sum(1 for h in ctx.hits if str(h.get("risk_class")) == "NDMA")
    return matrix, ndsri, ndma


def _document_control_rows(ctx: ExcipientReportContext) -> list[tuple[str, str]]:
    matrix, ndsri, ndma = _summary_counts(ctx)
    # 2026-09-03 수정 (F-4): 앞 5행은 p89 와 축자 동일했습니다 → p82 공용 헬퍼.
    return [
        *report_control_header_rows(
            report_id=ctx.report_id,
            generated=ctx.generated_at,
            platform=f"{PLATFORM_NAME} — {PLATFORM_FULL}",
            document_version=DOCUMENT_VERSION,
            assessment_type="Excipient-API interaction & nitrosamine (NDSRI/NDMA) screening",
        ),
        ("Selected excipients", str(len(ctx.excipient_keys))),
        ("Total predicted hits", str(len(ctx.hits))),
        ("Matrix hits", str(matrix)),
        ("NDSRI hits", str(ndsri)),
        ("NDMA hits", str(ndma)),
    ]


def _api_assessment_rows(ctx: ExcipientReportContext) -> list[tuple[str, str]]:
    api_smi = ctx.api_canonical or ctx.api_smiles
    nitrosatable = "Yes" if (ctx.amine_secondary or ctx.amine_tertiary) else "No"
    return [
        ("Substance role", "Drug Substance (API) — primary evaluation subject"),
        ("Input SMILES", ctx.api_smiles),
        ("Canonical SMILES", api_smi),
        # 2026-09-03 수정: 평균 분자량으로 통일 + 라벨 명기.
        ("Molecular weight", f"{_exact_mw(api_smi)} {MOLECULAR_WEIGHT_LABEL}"),
        ("Structure depiction", "Figure 1 — 2D structure (see PDF / Excel)"),
        ("Secondary amine sites", str(ctx.amine_secondary)),
        ("Tertiary amine sites", str(ctx.amine_tertiary)),
        ("Nitrosatable amine present", nitrosatable),
        ("Formulation excipients screened", str(len(ctx.excipient_keys))),
    ]


def _excipient_profile_rows(ctx: ExcipientReportContext) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for i, (key, label) in enumerate(zip(ctx.excipient_keys, ctx.excipient_labels), start=1):
        rows.append((f"Excipient {i}", f"{label} [{key}]"))
    return rows or [("Excipients", "None selected")]


def _hit_evaluation_rows(hit: dict[str, Any], idx: int) -> list[tuple[str, str]]:
    risk = str(hit.get("risk_class") or "Matrix")
    smi = str(hit.get("smiles") or "")
    nitrite = str(hit.get("nitrite_exposure") or "—")
    nitrite_raw = str(hit.get("nitrite_exposure_raw") or "")
    if nitrite_raw and nitrite_raw != nitrite:
        nitrite = f"{nitrite} (QbD raw: {nitrite_raw})"
    validity = str(
        hit.get("chemical_validity_notes") or hit.get("kinetics_comment") or "—"
    )
    rows = [
        ("Hit ID", f"HIT-{idx:03d}"),
        ("Risk classification", risk),
        ("Executive outcome", _RISK_OUTCOME.get(risk, "Predicted excipient-API interaction hit")),
        ("Excipient", _format_report_value(hit.get("excipient_name"))),
        ("Excipient category", _format_report_value(hit.get("category"))),
        ("Mechanism origin", _format_report_value(hit.get("mechanism_origin"))),
        ("SMIRKS / rule", _format_report_value(hit.get("rule_name"))),
        ("Predicted product SMILES", smi or "N/A"),
        # 2026-09-03 수정: 평균 분자량으로 통일 + 라벨 명기.
        ("Molecular weight", f"{_exact_mw(smi)} {MOLECULAR_WEIGHT_LABEL}" if smi else "N/A"),
        ("Structure depiction", f"Figure {idx + 1} — predicted interaction / nitrosamine product"),
        ("Nitrosamine label", _format_report_value(hit.get("nitrosamine_label"))),
        ("Impurity type", _format_report_value(hit.get("impurity_type"))),
        ("Expected nitrite exposure", nitrite),
        ("ICH M7 class hint", _format_report_value(hit.get("ich_m7_class_hint"))),
        ("Alpha-hydrogen present", _format_report_value(hit.get("alpha_hydrogen_present"))),
        ("Steric hindrance factor", _format_report_value(hit.get("steric_hindrance"))),
        ("Chemical validity / kinetics", validity),
        ("Formulation alert context", _format_report_value(hit.get("alert_context"))),
        (
            "Regulatory action (screening)",
            "Confirm by LC-MS / nitrosamine testing; apply control strategy per ICH M7 / Q9 if positive",
        ),
    ]
    return rows


def _hit_rows(ctx: ExcipientReportContext) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for i, h in enumerate(ctx.hits, start=1):
        rows.append(
            {
                "No": str(i),
                "Risk_Class": str(h.get("risk_class") or h.get("rule_name") or ""),
                "Nitrosamine": str(h.get("nitrosamine_label") or h.get("impurity_type") or "—"),
                "Nitrite_ppm": str(h.get("nitrite_exposure") or "—"),
                "Impurity_Type": str(h.get("impurity_type") or "—"),
                "ICH_M7": str(h.get("ich_m7_class_hint") or "—"),
                "Alpha_H": str(h.get("alpha_hydrogen_present") or "—"),
                "Steric": str(h.get("steric_hindrance") or "—"),
                "Validity_Notes": str(
                    h.get("chemical_validity_notes") or h.get("kinetics_comment") or "—"
                ),
                "Excipient": str(h.get("excipient_name") or ""),
                "Mechanism": str(h.get("mechanism_origin") or ""),
                "Rule": str(h.get("rule_name") or ""),
                "SMILES": str(h.get("smiles") or ""),
                "Category": str(h.get("category") or ""),
                "Alert": str(h.get("alert_context") or ""),
                "MW": _exact_mw(str(h.get("smiles") or "")),
            }
        )
    return rows


def _executive_text_block(hit: dict[str, Any], idx: int) -> list[str]:
    risk = str(hit.get("risk_class") or "Matrix")
    title = f"  HIT {idx} EXECUTIVE OUTCOME — {risk}"
    lines = [
        "",
        "╔" + "═" * 86 + "╗",
        f"║{title:<86}║",
        "╚" + "═" * 86 + "╝",
    ]
    for label, val in _hit_evaluation_rows(hit, idx):
        if label in ("Structure depiction",):
            continue
        lines.append(f"  ▶ {label:<32}: {val}")
    return lines


def generate_excipient_report_text(ctx: ExcipientReportContext) -> str:
    matrix, ndsri, ndma = _summary_counts(ctx)
    lines: list[str] = [
        "=" * 88,
        f"{PLATFORM_NAME} — EXCIPIENT-API INTERACTION & NITROSAMINE REPORT",
        PLATFORM_FULL,
        f"Document Version: {DOCUMENT_VERSION}",
        "=" * 88,
        "",
        "DOCUMENT CONTROL",
        "-" * 88,
    ]
    for label, val in _document_control_rows(ctx):
        lines.append(f"  {label:<28}: {val}")

    lines.extend(["", "1. PRIMARY ASSESSMENT — DRUG SUBSTANCE (API)", "-" * 88])
    for label, val in _api_assessment_rows(ctx):
        lines.append(f"  {label:<28}: {val}")

    lines.extend(["", "2. FORMULATION EXCIPIENT PROFILE", "-" * 88])
    for label, val in _excipient_profile_rows(ctx):
        lines.append(f"  {label:<28}: {val}")

    lines.extend(
        [
            "",
            f"3. PREDICTED INTERACTION / NITROSAMINE HITS ({len(ctx.hits)})",
            "-" * 88,
        ]
    )
    if not ctx.hits:
        lines.append("  No SMIRKS-matched interaction or nitrosamine hits detected.")
    else:
        for idx, hit in enumerate(ctx.hits, start=1):
            lines.extend(_executive_text_block(hit, idx))

    lines.extend(
        [
            "",
            "4. SUMMARY",
            f"   Matrix interactions (Maillard/Formylation): {matrix}",
            f"   NDSRI nitrosamine hits: {ndsri}",
            f"   NDMA-specific hits: {ndma}",
            "",
            "5. DISCLAIMER",
            f"   {_DISCLAIMER}",
            "",
        ]
    )
    return "\n".join(lines)


def export_excipient_txt(ctx: ExcipientReportContext) -> bytes:
    return generate_excipient_report_text(ctx).encode("utf-8")


def export_excipient_json(ctx: ExcipientReportContext) -> bytes:
    matrix, ndsri, ndma = _summary_counts(ctx)
    payload = {
        "report_id": ctx.report_id,
        "generated_at": ctx.generated_at,
        "document_version": DOCUMENT_VERSION,
        "api_smiles": ctx.api_smiles,
        "api_canonical": ctx.api_canonical,
        "api_molecular_weight": _exact_mw(ctx.api_canonical or ctx.api_smiles),
        "excipient_keys": ctx.excipient_keys,
        "excipient_labels": ctx.excipient_labels,
        "amine_sites": {
            "secondary": ctx.amine_secondary,
            "tertiary": ctx.amine_tertiary,
        },
        "summary": {
            "total_hits": len(ctx.hits),
            "matrix": matrix,
            "ndsri": ndsri,
            "ndma": ndma,
        },
        "hits": [
            {**hit, "evaluation": dict(_hit_evaluation_rows(hit, i))}
            for i, hit in enumerate(ctx.hits, start=1)
        ],
        "disclaimer": _DISCLAIMER,
    }
    return json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")


def export_excipient_excel(ctx: ExcipientReportContext) -> bytes:
    wb = Workbook()
    header_fill = PatternFill("solid", fgColor="E8F5E9")
    bold = Font(bold=True)

    ws_sum = wb.active
    ws_sum.title = "Summary"
    summary_rows = (
        _document_control_rows(ctx)
        + [("", "")]
        + _api_assessment_rows(ctx)
    )
    for r, (k, v) in enumerate(summary_rows, start=1):
        ws_sum.cell(row=r, column=1, value=k).font = bold
        ws_sum.cell(row=r, column=1).fill = header_fill if k else header_fill
        ws_sum.cell(row=r, column=2, value=v)
        ws_sum.cell(row=r, column=2).alignment = Alignment(wrap_text=True, vertical="top")
    ws_sum.column_dimensions["A"].width = 34
    ws_sum.column_dimensions["B"].width = 78

    api_png = _structure_png_bytes(ctx.api_canonical or ctx.api_smiles, size="320x240")
    if api_png:
        ws_api = wb.create_sheet("API_Structure")
        ws_api.cell(row=1, column=1, value="Figure 1: Drug Substance (API)").font = bold
        ws_api.cell(row=2, column=1, value=ctx.api_canonical or ctx.api_smiles)
        img_buf = BytesIO(api_png)
        img_buf.seek(0)
        img = XLImage(img_buf)
        img.width = 280
        img.height = 210
        ws_api.add_image(img, "A4")
        ws_api.column_dimensions["A"].width = 48

    ws = wb.create_sheet("Interaction_Hits")
    headers = [
        "No",
        "Risk_Class",
        "Executive_Outcome",
        "Excipient",
        "Mechanism",
        "Rule",
        "SMILES",
        "MW_g_mol",
        "Nitrosamine",
        "Nitrite_ppm",
        "Impurity_Type",
        "ICH_M7",
        "Alpha_H",
        "Steric",
        "Validity_Notes",
        "Alert",
        "2D_Structure",
    ]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=1, column=c, value=h)
        cell.font = bold
        cell.fill = header_fill

    for row_idx, (hit, row) in enumerate(zip(ctx.hits, _hit_rows(ctx)), start=2):
        risk = str(hit.get("risk_class") or "Matrix")
        row_data = [
            row["No"],
            row["Risk_Class"],
            _RISK_OUTCOME.get(risk, ""),
            row["Excipient"],
            row["Mechanism"],
            row["Rule"],
            row["SMILES"],
            row["MW"],
            row["Nitrosamine"],
            row["Nitrite_ppm"],
            row["Impurity_Type"],
            row["ICH_M7"],
            row["Alpha_H"],
            row["Steric"],
            row["Validity_Notes"],
            row["Alert"],
            "",
        ]
        for col_idx, val in enumerate(row_data, start=1):
            cell = ws.cell(row=row_idx, column=col_idx, value=val)
            cell.alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[row_idx].height = 100
        png = _structure_png_bytes(row["SMILES"], size="220x165")
        if png:
            img_buf = BytesIO(png)
            img_buf.seek(0)
            img = XLImage(img_buf)
            img.width = 150
            img.height = 112
            col_letter = get_column_letter(len(headers))
            ws.add_image(img, f"{col_letter}{row_idx}")

    widths = [6, 10, 36, 18, 22, 16, 36, 10, 18, 14, 24, 14, 8, 8, 28, 28, 22]
    for col_idx, width in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(col_idx)].width = width

    buf = BytesIO()
    wb.save(buf)
    return buf.getvalue()


def export_excipient_pdf(ctx: ExcipientReportContext) -> bytes:
    pdf = FPDF()
    pdf.set_margins(18, 18, 18)
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.add_page()
    w = pdf.epw
    family, family_has_bold = _register_pdf_fonts(pdf)
    kv = lambda rows: _pdf_kv_rows(
        pdf, family, rows, w, family_has_bold=family_has_bold
    )

    _pdf_write_wrapped(
        pdf,
        family,
        _pdf_line(f"{PLATFORM_NAME} — EXCIPIENT-API INTERACTION & NITROSAMINE REPORT"),
        w,
        bold=True,
        size=14,
    )
    _pdf_write_wrapped(
        pdf,
        family,
        _pdf_line("Formulation Risk Screening — Maillard · Formylation · NDSRI/NDMA"),
        w,
        size=10,
    )
    _pdf_blank_lines(pdf, 1)

    _pdf_section(pdf, family, "DOCUMENT CONTROL", w)
    kv(_document_control_rows(ctx))

    api_smi = ctx.api_canonical or ctx.api_smiles
    api_png = _structure_png_bytes(api_smi, size="400x300")
    if api_png:
        _pdf_add_structure_figure(
            pdf,
            family,
            api_png,
            w,
            caption="Figure 1: Drug Substance (API) — Primary Assessment Subject",
        )

    _pdf_section(pdf, family, "1. PRIMARY ASSESSMENT — DRUG SUBSTANCE (API)", w)
    _pdf_write_wrapped(
        pdf,
        family,
        _pdf_line(
            "Executive characterization of the API precedes excipient matrix and "
            "nitrosamine interaction screening in this report."
        ),
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 1)
    kv(_api_assessment_rows(ctx))

    _pdf_section(pdf, family, "2. FORMULATION EXCIPIENT PROFILE", w)
    kv(_excipient_profile_rows(ctx))

    _pdf_section(
        pdf,
        family,
        f"3. PREDICTED INTERACTION / NITROSAMINE HITS ({len(ctx.hits)})",
        w,
    )
    matrix, ndsri, ndma = _summary_counts(ctx)
    kv(
        [
            ("Matrix hits", str(matrix)),
            ("NDSRI hits", str(ndsri)),
            ("NDMA hits", str(ndma)),
        ]
    )

    if not ctx.hits:
        _pdf_write_wrapped(
            pdf,
            family,
            _pdf_line("No interaction or nitrosamine hits detected for selected excipients."),
            w,
        )
    else:
        for idx, hit in enumerate(ctx.hits, start=1):
            risk = str(hit.get("risk_class") or "Matrix")
            exc = str(hit.get("excipient_name") or "")
            _pdf_section(
                pdf,
                family,
                _pdf_line(f"HIT {idx} — {risk} — {exc}"),
                w,
                pre_gap_lines=2,
            )
            hit_png = _structure_png_bytes(str(hit.get("smiles") or ""), size="400x300")
            if hit_png:
                _pdf_add_structure_figure(
                    pdf,
                    family,
                    hit_png,
                    w,
                    caption=_pdf_line(
                        f"Figure {idx + 1}: Predicted product — {risk} / {exc}"
                    ),
                )
            kv(_hit_evaluation_rows(hit, idx))

    _pdf_section(pdf, family, "4. REGULATORY DISCLAIMER", w)
    _pdf_write_wrapped(pdf, family, _pdf_line(_DISCLAIMER), w)
    return _pdf_output_bytes(pdf)


def safe_excipient_filename_base(api_smiles: str) -> str:
    base = _canon(api_smiles) or "excipient_matrix"
    base = re.sub(r"[^\w\-]+", "_", base)[:48]
    return f"PRIME_DS_Excipient_{base}"
