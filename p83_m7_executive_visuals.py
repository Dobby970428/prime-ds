"""
ICH M7(R2) Executive Outcome — shared visual emphasis for TXT / PDF / DOCX reports.

Highlights Classification, Definition, Primary M7 route, Classification basis,
and Recommended control with class-specific color themes.
"""

# 등급별 색상(보고서 맨 위의 ★ EXECUTIVE OUTCOME ★ 박스)
# Class 1  진한 빨강
# Class 2  주황
# Class 3  노랑
# Class 4  파랑
# Class 5  초록

from __future__ import annotations

from typing import TYPE_CHECKING

from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from fpdf import FPDF

from p82_pdf_typesetting import (
    BOX_PAD_MM,
    PdfTwoColumnRowLayout,
    pdf_blend_rgb,
    pdf_draw_two_column_row,
    pdf_draw_wrapped_text,
    pdf_measure_two_column_row,
    pdf_measure_wrapped_text,
    pdf_sanitize,
    pdf_set_font,
)

if TYPE_CHECKING:
    from docx.document import Document as DocxDocument

M7_CLASS_THEMES: dict[int, dict[str, object]] = {
    1: {
        "bg": (255, 235, 238),
        "border": (183, 28, 28),
        "accent": (183, 28, 28),
        "fill_hex": "FFEBEE",
        "accent_hex": "B71C1C",
        "label": "Highest concern",
    },
    2: {
        "bg": (255, 243, 224),
        "border": (230, 81, 0),
        "accent": (230, 81, 0),
        "fill_hex": "FFF3E0",
        "accent_hex": "E65100",
        "label": "Known mutagen — TTC control",
    },
    3: {
        "bg": (255, 248, 225),
        "border": (245, 127, 23),
        "accent": (245, 127, 23),
        "fill_hex": "FFF8E1",
        "accent_hex": "F57F17",
        "label": "Alerting structure — confirmatory testing",
    },
    4: {
        "bg": (227, 242, 253),
        "border": (21, 101, 192),
        "accent": (21, 101, 192),
        "fill_hex": "E3F2FD",
        "accent_hex": "1565C0",
        "label": "Drug-related alert — Q3A/Q3B limits",
    },
    5: {
        "bg": (232, 245, 233),
        "border": (46, 125, 50),
        "accent": (46, 125, 50),
        "fill_hex": "E8F5E9",
        "accent_hex": "2E7D32",
        "label": "No alert / sufficient negative data",
    },
}

_EXECUTIVE_EMPHASIS_LABELS = frozenset(
    {
        "FINAL HAZARD CLASS",
        "FINAL HAZARD CLASS (Parent API)",
        "Classification",
        "Classification basis",
    }
)


def m7_class_theme(class_num: int) -> dict[str, object]:
    try:
        return M7_CLASS_THEMES[int(class_num)]
    except (TypeError, ValueError, KeyError):
        return M7_CLASS_THEMES[5]


def build_m7_executive_rows(
    class_num: int,
    *,
    classification: str,
    definition: str,
    primary_route: str,
    classification_basis: str,
    recommended_control: str,
    hazard_label: str = "FINAL HAZARD CLASS",
) -> list[tuple[str, str]]:
    """Standard executive-outcome key-value rows."""
    return [
        (hazard_label, f"Class {class_num}"),
        ("Classification", classification),
        ("Definition", definition),
        ("Primary M7 route", primary_route),
        ("Classification basis", classification_basis),
        ("Recommended control", recommended_control),
    ]


def build_degradant_executive_rows(product) -> list[tuple[str, str]]:
    """Degradant appendix — same executive layout as main ICH M7 Section 5."""
    from p80_formal_assessment_report import ICH_M7_REGULATORY
    from p43_degradation_predictor import ich_m7_signal_label_en, translate_degradant_text

    m7 = ICH_M7_REGULATORY.get(product.ich_m7_class, ICH_M7_REGULATORY[5])
    alerts = ", ".join(product.alerts) if product.alerts else "None detected"
    basis = (
        f"ICH M7 hazard class assignment logic: Class {product.ich_m7_class} — "
        f"{m7['definition']} Ashby-Tennant alerts: {alerts}."
    )
    return build_m7_executive_rows(
        product.ich_m7_class,
        classification=m7["name"],
        definition=m7["definition"],
        primary_route="Option 4 — Expert rule-based assessment (PRIMARY)",
        classification_basis=basis,
        recommended_control=m7["control"],
        hazard_label="FINAL HAZARD CLASS",
    )


def format_m7_executive_text_block(
    class_num: int,
    rows: list[tuple[str, str]],
    *,
    width: int = 86,
    title: str | None = None,
) -> list[str]:
    """ASCII callout box for plain-text reports.

    2026-09-03 수정 (E-3): 테두리 줄과 본문 줄이 4칸 어긋나 있었습니다.
      - 테두리 = "╔" + "═" * width + "╗"  → width + 2 자
      - 본문   = "║" + prefix(2) + chunk + pad + "║" 에서 채움 폭이
                 inner_w = width - 4 라 width - 2 자 (86 폭이면 84 vs 88)
      inner_w 를 "두 ║ 사이의 컬럼 수" = width 로 재정의해 배너·본문·테두리가
      모두 width + 2 자로 일치하게 했습니다.
      또한 ★(U+2605)·▶(U+25B6)은 East Asian Wide 라 모노스페이스에서 2칸으로
      렌더돼 눈으로 봐도 어긋나므로 ASCII '*' · '>' 로 교체했습니다.
      (주의: CJK 값이 들어가면 글자 수는 일정하되 실제 렌더 폭은 넓어집니다.
       규제 보고서 본문은 영문이므로 글자 수 기준 정렬을 채택했습니다.)
    """
    theme = m7_class_theme(class_num)
    banner = title or f"EXECUTIVE OUTCOME — ICH M7 Class {class_num} ({theme['label']})"
    inner_w = width  # 두 ║ 사이의 컬럼 수
    top = "╔" + "═" * width + "╗"
    mid = "╠" + "═" * width + "╣"
    bot = "╚" + "═" * width + "╝"
    banner_text = f"* {banner} *"[:inner_w]
    lines = ["", top, f"║{banner_text:^{inner_w}}║", mid]
    for label, value in rows:
        text = f"{label}: {value}"
        emphasize = label in _EXECUTIVE_EMPHASIS_LABELS
        prefix = "> " if emphasize else "  "
        avail = inner_w - len(prefix)
        while True:
            chunk, text = text[:avail], text[avail:]
            lines.append(f"║{prefix}{chunk:<{avail}}║")
            if not text:
                break
    lines.append(bot)
    lines.append("")
    return lines


# 2026-09-03 수정: sanitize 단일 출처는 p82.pdf_sanitize 입니다.
#   이 이름은 하위 호환용 얇은 재export 로만 남깁니다.
def _pdf_sanitize(text: str, family: str | None = None) -> str:
    return pdf_sanitize(text, family)


def _pdf_set_font(pdf: FPDF, family: str, bold: bool, size: float) -> None:
    pdf_set_font(pdf, family, bold, size)


def pdf_add_m7_executive_callout(
    pdf: FPDF,
    family: str,
    w: float,
    class_num: int,
    rows: list[tuple[str, str]],
    *,
    family_has_bold: bool = True,
    title: str | None = None,
    label_frac: float = 0.38,
    line_h: float | None = None,
) -> None:
    """Colored banner + dynamic two-column callout; boxes sized from measured text."""
    _ = line_h  # legacy kwarg from degradation PDF export (layout is auto-measured)
    theme = m7_class_theme(class_num)
    bg: tuple[int, int, int] = theme["bg"]  # type: ignore[assignment]
    border: tuple[int, int, int] = theme["border"]  # type: ignore[assignment]
    accent: tuple[int, int, int] = theme["accent"]  # type: ignore[assignment]
    x = pdf.l_margin
    pad = BOX_PAD_MM
    inner_x = x + pad
    inner_w = w - 2 * pad
    label_col_w = inner_w * label_frac
    value_col_w = inner_w - label_col_w

    # 2026-09-03 수정: family 를 넘겨 Helvetica fallback 시 latin-1 강제 인코딩.
    banner_text = pdf_sanitize(
        title or f"* EXECUTIVE OUTCOME - ICH M7 Class {class_num} ({theme['label']}) *",
        family,
    )
    banner_block = pdf_measure_wrapped_text(
        pdf,
        family,
        bold=family_has_bold,
        width_mm=w - 2 * pad,
        text=banner_text,
        base_size_pt=10.0,
        min_size_pt=8.0,
        family_has_bold=family_has_bold,
        max_lines=3,
    )
    banner_h = banner_block.block_height + 4.0

    row_layouts: list[PdfTwoColumnRowLayout] = []
    for label, value in rows:
        emphasize = label in _EXECUTIVE_EMPHASIS_LABELS
        row_layouts.append(
            pdf_measure_two_column_row(
                pdf,
                family,
                label=label,
                value=value,
                label_width=label_col_w,
                value_width=value_col_w,
                label_size_pt=10.0,
                value_size_pt=11.0 if emphasize else 10.0,
                value_bold=emphasize,
                emphasize=emphasize,
                family_has_bold=family_has_bold,
                value_min_size_pt=7.0,
            )
        )

    body_h = sum(layout.row_height for layout in row_layouts) + pad * 2
    total_h = banner_h + 1.0 + body_h + 2.5
    if pdf.get_y() + total_h > pdf.page_break_trigger:
        pdf.add_page()

    y0 = pdf.get_y()

    pdf.set_fill_color(*accent)
    pdf.set_draw_color(*border)
    pdf.set_text_color(255, 255, 255)
    pdf.rect(x, y0, w, banner_h, style="F")
    banner_y = y0 + max(1.0, (banner_h - banner_block.block_height) / 2)
    pdf_draw_wrapped_text(
        pdf,
        family,
        bold=family_has_bold,
        x=x + pad,
        y=banner_y,
        width_mm=w - 2 * pad,
        text=banner_text,
        block=banner_block,
        family_has_bold=family_has_bold,
        align="C",
        text_color=(255, 255, 255),
    )

    y_body = y0 + banner_h + 1.0
    pdf.set_fill_color(*bg)
    pdf.set_draw_color(*border)
    pdf.set_line_width(0.35)
    pdf.rect(x, y_body, w, body_h, style="FD")

    y = y_body + pad
    value_row_fill = pdf_blend_rgb(accent, bg, weight_a=0.14)
    for (label, value), layout in zip(rows, row_layouts):
        emphasize = layout.emphasize
        label_color = accent if emphasize else (40, 40, 40)
        value_color = accent if emphasize else (40, 40, 40)
        pdf_draw_two_column_row(
            pdf,
            family,
            x=inner_x,
            y=y,
            label=label,
            value=value,
            label_width=label_col_w,
            value_width=value_col_w,
            layout=layout,
            family_has_bold=family_has_bold,
            label_color=label_color,
            value_color=value_color,
            value_fill=value_row_fill if emphasize else None,
        )
        y += layout.row_height

    pdf.set_text_color(0, 0, 0)
    pdf.set_y(y_body + body_h + 2.5)


def _shade_docx_cell(cell, fill_hex: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), fill_hex)
    tc_pr.append(shd)


def _set_docx_run_font(run, *, size: int = 10, bold: bool = False, color: RGBColor | None = None) -> None:
    run.font.name = "Times New Roman"
    run.font.size = Pt(size)
    run.font.bold = bold
    if color is not None:
        run.font.color.rgb = color
    r_pr = run._element.get_or_add_rPr()
    r_fonts = r_pr.get_or_add_rFonts()
    for attr in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"):
        r_fonts.set(qn(attr), "Times New Roman")


def docx_add_m7_executive_callout(
    doc: DocxDocument,
    class_num: int,
    rows: list[tuple[str, str]],
    *,
    title: str | None = None,
) -> None:
    """Shaded Word table callout for executive ICH M7 outcome."""
    theme = m7_class_theme(class_num)
    fill_hex = str(theme["fill_hex"])
    accent_hex = str(theme["accent_hex"])
    accent_rgb = RGBColor(int(accent_hex[0:2], 16), int(accent_hex[2:4], 16), int(accent_hex[4:6], 16))

    table = doc.add_table(rows=1 + len(rows), cols=2)
    table.style = "Table Grid"

    header = table.rows[0].cells
    merged = header[0].merge(header[1])
    _shade_docx_cell(merged, accent_hex)
    p = merged.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    banner = title or f"★ EXECUTIVE OUTCOME — ICH M7 Class {class_num} ({theme['label']}) ★"
    run = p.add_run(banner)
    _set_docx_run_font(run, size=11, bold=True, color=RGBColor(255, 255, 255))

    for idx, (label, value) in enumerate(rows, start=1):
        row = table.rows[idx]
        c0, c1 = row.cells[0], row.cells[1]
        c0.width = Inches(2.2)
        c1.width = Inches(4.3)
        _shade_docx_cell(c0, fill_hex)
        _shade_docx_cell(c1, fill_hex)
        emphasize = label in _EXECUTIVE_EMPHASIS_LABELS
        r0 = c0.paragraphs[0].add_run(label)
        _set_docx_run_font(r0, size=10, bold=True, color=accent_rgb if emphasize else None)
        r1 = c1.paragraphs[0].add_run(str(value))
        _set_docx_run_font(r1, size=11 if emphasize else 10, bold=emphasize, color=accent_rgb if emphasize else None)

    doc.add_paragraph()
