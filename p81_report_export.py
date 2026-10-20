"""Regulatory report export: structured Word (DOCX) and PDF."""

from __future__ import annotations

import logging
import os
import re
from dataclasses import replace
from io import BytesIO

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from fpdf import FPDF
# 2026-09-03 수정: fpdf2 deprecated `ln=1` 대체용 enum.
from fpdf.enums import XPos, YPos

# p80이 만든 내용을 Word와 PDF로 바꿈.
# 폰트 지정        Times New Roman (없으면 Calibri)
# 한글 폰트 탐색    Windows 폰트 폴더에서 자동 검색
# 표 그리기        셀 너비 · 음영 · 병합
# 구조식 삽입      p84 호출
# 경보 강조        구조식에서 위험 부위를 빨갛게
# 페이지 나눔      표가 잘리지 않게 미리 계산

from p80_formal_assessment_report import (
    DOCUMENT_VERSION,
    ICH_M7_REGULATORY,
    PLATFORM_FULL,
    PLATFORM_NAME,
    REPORT_TITLE,
    REPORT_TITLE_LINE1,
    REPORT_TITLE_LINE2,
    AssessmentContext,
    ProductAssessment,
    build_m7_regulatory_option_summary,
    hybrid_ensemble_report_rows,
    includes_process_route_section,
    indexed_db_count_label,
    m7_action_flags,
    substance_category_label,
)
from p83_m7_executive_visuals import (
    build_m7_executive_rows,
    docx_add_m7_executive_callout,
    pdf_add_m7_executive_callout,
)
# 2026-09-03 수정: PDF 공용 헬퍼 단일 출처는 p82 입니다.
#   `_pdf_sanitize` 가 p81·p82·p83 에 따로 있었고 p88·p89 는 p81 판을 import 해
#   실질 5벌이었습니다. `_pdf_output_bytes` 는 p81·p86·p87 에 축자 3벌,
#   `_pdf_line_height` 는 p81(4.2/0.45) 과 p82(4.4/0.46) 가 서로 달라
#   측정·그리기 경로가 어긋났습니다. 전부 p82 로 통합합니다.
from p82_pdf_typesetting import (
    pdf_line_height,
    pubchem_structure_fetch_enabled,
    pdf_measure_two_column_row,
    pdf_output_bytes,
    pdf_sanitize,
)
from p84_structure_images import (
    DEFAULT_REPORT_IMAGE_SIZE,
    rdkit_alert_highlight_png,
    rdkit_clean_structure_png,
    structure_png,
)
from p62_toxicity_alerts import scan_report_structural_alerts

logger = logging.getLogger(__name__)

REG_FONT = "Times New Roman"
REG_FONT_FALLBACK = "Calibri"
REPORT_TITLE_COLOR = RGBColor(0x1A, 0x36, 0x5D)
REPORT_TITLE_COLOR_PDF = (26, 54, 93)


def _windows_fonts() -> tuple[str | None, str | None]:
    """한글이 들어간 PDF를 위해 유니코드 TTF를 찾습니다.

    Windows에서는 맑은 고딕을 우선합니다. 폰트를 못 찾으면 fpdf가
    Helvetica로 떨어지는데, Helvetica는 latin-1만 지원하므로 보고서에
    한글이 한 글자라도 있으면 생성이 실패합니다.

    2026-08-25 — 리눅스·맥 경로를 추가했습니다. Windows 전용이던 탐색
    때문에 다른 환경에서 한글 PDF가 전부 실패했습니다.
    """
    windir = os.environ.get("WINDIR", r"C:\Windows")
    win_dir = os.path.join(windir, "Fonts")

    candidates: tuple[tuple[str, str, str | None], ...] = (
        # (폴더, 일반, 굵게)
        (win_dir, "malgun.ttf", "malgunbd.ttf"),
        (win_dir, "malgunsl.ttf", None),
        (win_dir, "times.ttf", "timesbd.ttf"),
        (win_dir, "arial.ttf", "arialbd.ttf"),
        (win_dir, "calibri.ttf", "calibrib.ttf"),
        # macOS
        ("/System/Library/Fonts", "AppleSDGothicNeo.ttc", None),
        ("/Library/Fonts", "AppleGothic.ttf", None),
        # Linux — 나눔고딕 · 노토 산스 CJK
        ("/usr/share/fonts/truetype/nanum", "NanumGothic.ttf", "NanumGothicBold.ttf"),
        ("/usr/share/fonts/truetype/nanum", "NanumBarunGothic.ttf", "NanumBarunGothicBold.ttf"),
        ("/usr/share/fonts/opentype/noto", "NotoSansCJK-Regular.ttc", "NotoSansCJK-Bold.ttc"),
        ("/usr/share/fonts/truetype/noto", "NotoSansCJK-Regular.ttc", None),
        ("/usr/share/fonts/truetype/dejavu", "DejaVuSans.ttf", "DejaVuSans-Bold.ttf"),
    )

    for folder, regular_name, bold_name in candidates:
        regular_path = os.path.join(folder, regular_name)
        if not os.path.isfile(regular_path):
            continue
        bold_path = None
        if bold_name:
            candidate = os.path.join(folder, bold_name)
            if os.path.isfile(candidate):
                bold_path = candidate
        return regular_path, bold_path
    return None, None


# 2026-09-03 \uC218\uC815: \uC790\uCCB4 sanitize \uAD6C\uD604\uC744 \uC0AD\uC81C\uD558\uACE0 p82.pdf_sanitize \uB85C \uC704\uC784\uD569\uB2C8\uB2E4.
#   \uAE30\uC874 \uAD6C\uD604\uC740 U+2600-27BF \uC640 U+10000 \uC774\uC0C1\uB9CC \uC81C\uAC70\uD574\uC11C U+2500(\u2500, p43 \uC5D0 32\uD68C)
#   U+2265(\u2265) U+2026(\u2026) U+2460-62(\u2460\u2461\u2462, p72) U+03B1/B2(\u03B1 \u03B2) \uADF8\uB9AC\uACE0 \uD55C\uAE00\uC774
#   \uC804\uBD80 \uD1B5\uACFC\uD588\uACE0, \uD3F0\uD2B8\uB97C \uBABB \uCC3E\uC544 Helvetica(latin-1)\uB85C \uB5A8\uC5B4\uC9C0\uBA74 fpdf2 \uAC00
#   FPDFUnicodeEncodingException \uC744 \uB358\uC9C0\uBA70 PDF \uC0DD\uC131\uC774 \uC8FD\uC5C8\uC2B5\uB2C8\uB2E4.
#   \uC774 \uC774\uB984\uC740 p88\u00B7p89 \uAC00 import \uD558\uACE0 \uC788\uC73C\uBBC0\uB85C \uD558\uC704 \uD638\uD658 \uC7ACexport \uB85C \uB0A8\uAE41\uB2C8\uB2E4.
def _pdf_sanitize(text: str, family: str | None = None) -> str:
    """Replace Unicode chars that break narrow PDF fonts (p82 \uB2E8\uC77C \uCD9C\uCC98 \uC704\uC784)."""
    return pdf_sanitize(text, family)


# 2026-09-03 \uC218\uC815: p86\u00B7p87 \uACFC \uCD95\uC790 \uB3D9\uC77C\uD588\uB358 \uAD6C\uD604 \u2192 p82.pdf_output_bytes \uB85C \uD1B5\uD569.
def _pdf_output_bytes(pdf: FPDF) -> bytes:
    return pdf_output_bytes(pdf)


def _register_pdf_fonts(pdf: FPDF) -> tuple[str, bool]:
    regular, bold = _windows_fonts()
    if regular:
        pdf.add_font("Reg", "", regular)
        if bold:
            pdf.add_font("Reg", "B", bold)
            return "Reg", True
        return "Reg", False
    # 2026-09-03 \uC218\uC815: \uC9C0\uAE08\uAE4C\uC9C0\uB294 \uC870\uC6A9\uD788 Helvetica(latin-1 \uCF54\uC5B4 \uD3F0\uD2B8)\uB85C \uB5A8\uC5B4\uC838
    #   \uD55C\uAE00\uC774 '?' \uB85C \uAE68\uC9C0\uAC70\uB098 \uC608\uC678\uB85C \uC8FD\uC5C8\uC2B5\uB2C8\uB2E4. \uBA85\uC2DC\uC801\uC73C\uB85C \uACBD\uACE0\uB97C \uB0A8\uAE41\uB2C8\uB2E4.
    #   \uC2E4\uC81C \uBC29\uC5B4(latin-1 \uAC15\uC81C \uC778\uCF54\uB529)\uB294 p82.pdf_sanitize \uAC00 \uC218\uD589\uD569\uB2C8\uB2E4.
    logger.warning(
        "\uC720\uB2C8\uCF54\uB4DC TTF \uB97C \uCC3E\uC9C0 \uBABB\uD574 latin-1 \uCF54\uC5B4 \uD3F0\uD2B8(Helvetica)\uB85C \uB300\uCCB4\uD569\uB2C8\uB2E4. "
        "\uD55C\uAE00\u00B7\uBC15\uC2A4\uBB38\uC790\u00B7\uAE30\uD638\uB294 '?' \uB85C \uCE58\uD658\uB418\uC5B4 \uCD9C\uB825\uB429\uB2C8\uB2E4. "
        "\uB9D1\uC740 \uACE0\uB515 / NanumGothic / NotoSansCJK / DejaVuSans \uC911 \uD558\uB098\uB97C \uC124\uCE58\uD558\uC138\uC694."
    )
    return "Helvetica", True


def _pdf_set_font(pdf: FPDF, family: str, bold: bool, size: float) -> None:
    if family == "Helvetica":
        pdf.set_font("Helvetica", "B" if bold else "", size)
    else:
        style = "B" if bold and family == "Reg" else ""
        pdf.set_font(family, style=style, size=size)


def _format_report_value(value: object) -> str:
    """Match DOCX table cells: None/blank → N/A; otherwise preserve str(value)."""
    if value is None:
        return "N/A"
    s = str(value).strip()
    return s if s else "N/A"


# 2026-09-03 수정: 정의만 있고 호출되지 않던 죽은 코드였습니다.
#   `_pdf_add_table` 의 값 셀 렌더링에서 실제로 사용하도록 연결했습니다.
def _pdf_safe_value(value: object, family: str | None = None) -> str:
    """PDF-safe cell text — same semantics as DOCX, plus narrow-font sanitization."""
    return _pdf_sanitize(_format_report_value(value), family)


_VERDICT_LABEL_KEYS = frozenset({
    "structural alert result",
    "alerts identified",
    "match status",
    "final hazard class",
    "classification",
    "ames prediction",
    "process stability",
    "acceptable intake",
    "regulatory action",
})

_NO_EMPHASIS_LABEL_KEYS = frozenset({
    "conclusion",
    "definition",
    "recommended control",
    "stability rationale",
    "exposure assessment",
    "reference",
    "formation mechanism",
    "engine",
    "external database",
})

_VERDICT_VALUE_MARKERS = (
    "negative",
    "positive",
    "none detected",
    "no alert",
    "no ashby",
    "validated",
    "theoretical",
    "not matched",
    "no match",
    "stable",
    "unstable",
    "moderately stable",
    "class 1",
    "class 2",
    "class 3",
    "class 4",
    "class 5",
    "mutagenic",
    "non-mutagenic",
    "not integrated",
    "see classification",
    "threshold of toxicological concern",
    "ich q3a",
    "ich q3b",
    "ich m7",
)


def _should_emphasize_verdict(label: str, value: object) -> bool:
    """True when a table value is a key assessment outcome (bold in DOCX/PDF)."""
    ll = label.lower().strip()
    if ll in _NO_EMPHASIS_LABEL_KEYS:
        return False
    val = _format_report_value(value) if value is not None else "N/A"
    if ll in _VERDICT_LABEL_KEYS:
        return True
    if any(k in ll for k in ("alert result", "hazard class", "match status", "ames prediction")):
        return True
    vl = val.lower()
    if len(vl) > 120:
        return False
    return any(m in vl for m in _VERDICT_VALUE_MARKERS)


# 2026-09-03 수정: p82.pdf_line_height(4.4 / 0.46) 와 값이 달라(4.2 / 0.45)
#   측정은 p82, 그리기는 p81 상수로 도는 혼합 경로 때문에 페이지 나눔 예측이
#   실제 높이와 최대 2.2% 어긋났습니다. 정의를 삭제하고 p82 를 재export 합니다.
_pdf_line_height = pdf_line_height


_PDF_BODY_LINE_MM = 5.0
_PDF_SECTION_GAP_LINES = 1


_PDF_FIGURE_IMG_MM = 72.0
_PDF_FIGURE_GAP_LINES = 2
_PDF_APPENDIX_IMG_MM = 50.0
_PDF_APPENDIX_IMAGE_SIZE = "400x300"


def _smiles_alert_highlight_info(
    smiles: str,
    alerts: list[str] | None = None,
) -> tuple[list[str], bool]:
    """SMILES + optional alert labels → (labels, highlight_available)."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles((smiles or "").strip()) if smiles else None
    scanned_names, atoms, _ = scan_report_structural_alerts(mol)
    labels = list(alerts or []) or scanned_names
    return labels, bool(atoms)


def _report_alert_labels(ctx: AssessmentContext, smiles: str) -> tuple[list[str], bool]:
    """리포트 Figure 2용 alert 라벨 · 하이라이트 가능 여부."""
    return _smiles_alert_highlight_info(smiles, ctx.alerts)


def _report_smiles(ctx: AssessmentContext) -> str:
    return (ctx.smiles or ctx.predicted_impurity_smiles or "").strip()


def _structure_png_for_smiles(
    smiles: str,
    *,
    alerts: list[str] | None = None,
    pubchem_cid: str | None = None,
    highlighted: bool = False,
    size: str = DEFAULT_REPORT_IMAGE_SIZE,
) -> bytes | None:
    """2D structure — clean or alert-highlighted RDKit rendering."""
    smi = (smiles or "").strip()
    if not smi:
        return None

    try:
        if highlighted:
            data = rdkit_alert_highlight_png(smi, alerts, size=size)
        else:
            data = rdkit_clean_structure_png(smi, size=size)
        if data:
            return data
        # 2026-09-03 수정: RDKit 실패 시 SMILES/CID 를 PubChem REST 로 보내던 경로를
        #   PRIME_ALLOW_PUBCHEM_STRUCTURE_FETCH 명시적 opt-in 으로 제한합니다.
        #   기본값에서는 조회하지 않고 구조 이미지를 생략합니다(오프라인 재현성 +
        #   예측된 불순물 구조의 외부 전송 방지).
        if pubchem_cid and pubchem_structure_fetch_enabled():
            data = structure_png(pubchem_cid, smi, size=size)
            if data:
                return data
        return rdkit_clean_structure_png(smi, size=size)
    except Exception:
        return None


def _report_structure_png(ctx: AssessmentContext, *, highlighted: bool = False) -> bytes | None:
    return _structure_png_for_smiles(
        _report_smiles(ctx),
        alerts=ctx.alerts,
        pubchem_cid=ctx.pubchem_cid,
        highlighted=highlighted,
    )


def _add_evaluation_structure_figures_docx(doc: Document, ctx: AssessmentContext) -> None:
    """Title 직후 — Figure 1 (clean) + Figure 2 (alert highlight, 감지 시만)."""
    smi = _report_smiles(ctx)
    clean_png = _structure_png_for_smiles(smi, pubchem_cid=ctx.pubchem_cid, highlighted=False)
    if not clean_png:
        return

    _add_structure_image_docx(
        doc,
        clean_png,
        title="",
        caption="Figure 1: Target Evaluation Structure",
    )

    alert_labels, has_highlight = _report_alert_labels(ctx, smi)
    if has_highlight:
        alert_png = _structure_png_for_smiles(
            smi,
            alerts=alert_labels,
            pubchem_cid=ctx.pubchem_cid,
            highlighted=True,
        )
        alert_text = ", ".join(alert_labels) if alert_labels else "Structural alert zone"
        _add_structure_image_docx(
            doc,
            alert_png or clean_png,
            title="",
            caption=(
                "Figure 2: Identified Structural Alert Zone (Highlighted in Red) — "
                f"{alert_text}"
            ),
        )
    else:
        _add_para(
            doc,
            "Figure 2: No Ashby-Tennant / ICH M7 structural alert zone detected "
            "(highlight figure omitted).",
            size=9,
            italic=True,
            align=WD_ALIGN_PARAGRAPH.CENTER,
            space_after=10,
        )


def _pdf_add_evaluation_structure_figures(pdf: FPDF, family: str, ctx: AssessmentContext, w: float) -> None:
    """Report title 직후 — 평가 대상 구조(Fig.1) + alert 하이라이트(Fig.2, 감지 시만)."""
    smi = _report_smiles(ctx)
    clean_png = _report_structure_png(ctx, highlighted=False)
    if not clean_png:
        return

    _pdf_add_structure_figure(
        pdf,
        family,
        clean_png,
        w,
        caption="Figure 1: Target Evaluation Structure",
    )

    alert_labels, has_highlight = _report_alert_labels(ctx, smi)
    if has_highlight:
        alert_png = _structure_png_for_smiles(
            smi,
            alerts=alert_labels,
            pubchem_cid=ctx.pubchem_cid,
            highlighted=True,
        ) or clean_png
        alert_text = ", ".join(alert_labels) if alert_labels else "Structural alert zone"
        _pdf_add_structure_figure(
            pdf,
            family,
            alert_png,
            w,
            caption=(
                "Figure 2: Identified Structural Alert Zone (Highlighted in Red) — "
                f"{alert_text}"
            ),
        )
    else:
        _pdf_write_centered(
            pdf,
            family,
            "Figure 2: No Ashby-Tennant / ICH M7 structural alert zone detected "
            "(highlight figure omitted).",
            w,
            size=9,
            italic=True,
        )
        pdf.ln(2)


def _add_structure_image_docx(
    doc: Document,
    png: bytes | None,
    *,
    title: str,
    caption: str,
) -> None:
    if not png:
        return
    if title.strip():
        _add_para(doc, title, size=11, bold=True, space_after=4)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    try:
        p.add_run().add_picture(BytesIO(png), width=Inches(2.75))
    except Exception:
        return
    _add_para(doc, caption, size=9, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER, space_after=10)


def _pdf_add_structure_figure(
    pdf: FPDF,
    family: str,
    png: bytes | None,
    w: float,
    *,
    caption: str,
) -> None:
    """Centered structure figure with caption beneath."""
    if not png:
        return
    img_w = min(_PDF_FIGURE_IMG_MM, w * 0.58)
    caption_h = _pdf_line_height(8) * 2
    _pdf_ensure_space(pdf, img_w + caption_h + 8)
    try:
        x = pdf.l_margin + max(0.0, (w - img_w) / 2)
        pdf.image(BytesIO(png), x=x, w=img_w)
    except Exception:
        return
    _pdf_blank_lines(pdf, 1)
    _pdf_write_wrapped(pdf, family, caption, w, size=8)
    _pdf_blank_lines(pdf, _PDF_FIGURE_GAP_LINES)


# 2026-09-03 수정: `_pdf_add_structure_image` 삭제 — 프로젝트 전체에서 호출자가
#   없었고, 제목 한 줄을 덧붙이는 것 외에는 `_pdf_add_structure_figure` 와 동일해
#   따로 유지할 이유가 없습니다. 실제 사용처(2.1절 포함)는 모두
#   `_pdf_add_structure_figure` 를 직접 호출합니다.


def _product_assessment_rows(pa: ProductAssessment) -> list[tuple[str, str]]:
    pa_m7 = pa.m7.class_number if pa.m7 else 5
    pa_info = ICH_M7_REGULATORY[pa_m7]
    pa_alert_status = "POSITIVE" if pa.alerts else "NEGATIVE"
    pa_alerts_text = ", ".join(pa.alerts) if pa.alerts else "None detected"
    pa_pub = pa.pubchem_status or (
        f"PubChem CID {pa.pubchem_cid}" if pa.pubchem_cid else "No PubChem lookup performed"
    )
    return [
        ("Substance name", pa.name),
        ("SMILES", pa.smiles),
        ("Molecular weight", f"{pa.mw} g/mol (RDKit)"),
        ("Structural alert result", pa_alert_status),
        ("Alerts identified", pa_alerts_text),
        ("ICH M7 hazard class", f"Class {pa_m7} — {pa_info['name']}"),
        ("Assessment conclusion", pa.expert_comment),
        ("Lipophilicity (Log P)", pa.logp),
        ("Process stability", pa.stability),
        ("Stability rationale", pa.stability_detail),
        ("External database", pa_pub),
        ("In-house DB match", pa.db_match or "N/A"),
    ]


def _add_product_assessment_docx(doc: Document, pa: ProductAssessment) -> None:
    _add_para(
        doc,
        "2.1 DRUG SUBSTANCE (FINAL PRODUCT) — STRUCTURE AND ASSESSMENT",
        size=12,
        bold=True,
        space_after=6,
    )
    t_prod = doc.add_table(rows=0, cols=2)
    t_prod.style = "Table Grid"
    id_rows = _product_assessment_rows(pa)[:3]
    eval_rows = _product_assessment_rows(pa)[3:]
    for label, val in id_rows:
        _add_table_row(t_prod, label, val, bold_label=True)
    clean_png = _structure_png_for_smiles(pa.smiles, pubchem_cid=pa.pubchem_cid, highlighted=False)
    if clean_png:
        _add_structure_image_docx(
            doc,
            clean_png,
            title="",
            caption="Figure 1: Target Evaluation Structure (final product / API)",
        )
    alert_labels, has_highlight = _smiles_alert_highlight_info(pa.smiles, pa.alerts)
    if has_highlight:
        alert_png = _structure_png_for_smiles(
            pa.smiles,
            alerts=alert_labels,
            pubchem_cid=pa.pubchem_cid,
            highlighted=True,
        )
        _add_structure_image_docx(
            doc,
            alert_png or clean_png,
            title="",
            caption=(
                "Figure 2: Identified Structural Alert Zone (Highlighted in Red) — "
                f"{', '.join(alert_labels)}"
            ),
        )
    t_eval = doc.add_table(rows=0, cols=2)
    t_eval.style = "Table Grid"
    for label, val in eval_rows:
        _add_table_row(
            t_eval,
            label,
            val,
            bold_label=label in ("Structural alert result", "ICH M7 hazard class"),
        )
    doc.add_paragraph()


def _pdf_add_product_assessment(
    pdf: FPDF,
    family: str,
    pa: ProductAssessment,
    w: float,
    *,
    family_has_bold: bool = True,
) -> None:
    """2.1 DRUG SUBSTANCE 절 — TXT(p80:699-730)와 동일 항목을 PDF 로."""
    _pdf_section(
        pdf,
        family,
        "2.1 DRUG SUBSTANCE (FINAL PRODUCT) — STRUCTURE AND ASSESSMENT",
        w,
    )
    # 2026-09-03 수정: family_has_bold 를 하드코딩(False)하지 않고 호출자 값을 씁니다.
    kv = lambda rows: _pdf_kv_rows(
        pdf,
        family,
        rows,
        w,
        family_has_bold=family_has_bold,
        gap_between_lines=_PDF_SECTION_GAP_LINES,
    )
    all_rows = _product_assessment_rows(pa)
    kv(all_rows[:3])
    clean_png = _structure_png_for_smiles(pa.smiles, pubchem_cid=pa.pubchem_cid, highlighted=False)
    if clean_png:
        _pdf_add_structure_figure(
            pdf,
            family,
            clean_png,
            w,
            caption="Figure 1: Target Evaluation Structure (final product / API)",
        )
    alert_labels, has_highlight = _smiles_alert_highlight_info(pa.smiles, pa.alerts)
    if has_highlight:
        alert_png = _structure_png_for_smiles(
            pa.smiles,
            alerts=alert_labels,
            pubchem_cid=pa.pubchem_cid,
            highlighted=True,
        )
        _pdf_add_structure_figure(
            pdf,
            family,
            alert_png or clean_png,
            w,
            caption=(
                "Figure 2: Identified Structural Alert Zone (Highlighted in Red) — "
                f"{', '.join(alert_labels)}"
            ),
        )
    kv(all_rows[3:])


# 2026-09-03 수정: `_pdf_add_structure` 삭제 — 호출자가 없는 1줄 위임 래퍼였고,
#   실제 경로는 `_pdf_add_evaluation_structure_figures` 를 직접 호출합니다.


def _pdf_blank_lines(pdf: FPDF, lines: int = _PDF_SECTION_GAP_LINES) -> None:
    """Insert blank vertical space equivalent to N text lines."""
    if lines > 0:
        pdf.ln(lines * _PDF_BODY_LINE_MM)


def _pdf_ensure_space(pdf: FPDF, needed_mm: float) -> None:
    """Add a page only when the remaining space is truly insufficient."""
    if pdf.get_y() + needed_mm > pdf.page_break_trigger:
        pdf.add_page()


def _pdf_write_wrapped(pdf: FPDF, family: str, text: str, w: float, *, bold: bool = False, size: float = 9) -> None:
    line_h = _pdf_line_height(size)
    # 2026-09-03 수정: family 를 넘겨 Helvetica fallback 시 latin-1 강제 인코딩.
    content = _pdf_sanitize(text, family) if str(text).strip() else "N/A"
    _pdf_ensure_space(pdf, line_h * 2)
    _pdf_set_font(pdf, family, bold, size)
    pdf.set_x(pdf.l_margin)
    pdf.multi_cell(w, line_h, content)


def _pdf_write_centered(
    pdf: FPDF,
    family: str,
    text: str,
    w: float,
    *,
    bold: bool = False,
    size: float = 9,
    italic: bool = False,
) -> None:
    line_h = _pdf_line_height(size)
    # 2026-09-03 수정: family 를 넘겨 Helvetica fallback 시 latin-1 강제 인코딩.
    content = _pdf_sanitize(text, family) if str(text).strip() else "N/A"
    _pdf_ensure_space(pdf, line_h * 2)
    _pdf_set_font(pdf, family, bold, size)
    if italic:
        try:
            pdf.set_font(family, "I", size)
        except Exception:
            pass
    pdf.set_x(pdf.l_margin)
    pdf.multi_cell(w, line_h, content, align="C")


def _pdf_write_report_title_block(pdf: FPDF, family: str, w: float) -> None:
    """Cover title — prominent centered banner."""
    _pdf_ensure_space(pdf, 28)
    pdf.set_text_color(*REPORT_TITLE_COLOR_PDF)
    _pdf_write_centered(pdf, family, REPORT_TITLE_LINE1, w, bold=True, size=18)
    pdf.ln(1)
    _pdf_write_centered(pdf, family, REPORT_TITLE_LINE2, w, bold=True, size=16)
    pdf.ln(2)
    pdf.set_text_color(0, 0, 0)
    _pdf_write_centered(
        pdf,
        family,
        "ICH M7(R2) — Complementary (Q)SAR and Expert Rule-Based Evaluation",
        w,
        size=10.5,
        italic=True,
    )
    pdf.ln(1)
    _pdf_write_centered(
        pdf,
        family,
        f"{PLATFORM_NAME} — {PLATFORM_FULL}",
        w,
        size=8.5,
        italic=True,
    )


def _pdf_section(
    pdf: FPDF,
    family: str,
    title: str,
    w: float,
    *,
    pre_gap_lines: int = _PDF_SECTION_GAP_LINES,
) -> None:
    gap_mm = pre_gap_lines * _PDF_BODY_LINE_MM
    _pdf_ensure_space(pdf, 14 + gap_mm)
    _pdf_blank_lines(pdf, pre_gap_lines)
    _pdf_write_wrapped(pdf, family, title, w, bold=True, size=11)


def _pdf_add_table(
    pdf: FPDF,
    family: str,
    rows: list[tuple[str, str]],
    w: float,
    *,
    family_has_bold: bool = True,
    bold_labels: bool = True,
    label_frac: float = 0.36,
) -> None:
    """Two-column bordered table with measured row heights (no clipped cells)."""
    if not rows:
        return

    label_w = w * label_frac
    value_w = w - label_w
    x0 = pdf.l_margin

    for label, value in rows:
        val = _format_report_value(value)
        emphasize = _should_emphasize_verdict(label, val)

        layout = pdf_measure_two_column_row(
            pdf,
            family,
            label=str(label),
            value=val,
            label_width=label_w,
            value_width=value_w,
            label_size_pt=9.0,
            value_size_pt=10.0 if emphasize else 9.0,
            value_bold=emphasize,
            emphasize=emphasize,
            family_has_bold=family_has_bold,
            value_min_size_pt=7.0,
        )

        _pdf_ensure_space(pdf, layout.row_height + 2)
        y0 = pdf.get_y()

        _pdf_set_font(pdf, family, bold_labels and family_has_bold, layout.label_block.font_size)
        pdf.set_xy(x0, y0)
        pdf.multi_cell(
            label_w,
            layout.label_block.line_height,
            # 2026-09-03 수정: family 를 넘겨 Helvetica fallback 시 latin-1 강제 인코딩.
            _pdf_sanitize(str(label), family),
            border="LTR",
            align="L",
        )
        y_after_label = pdf.get_y()

        _pdf_set_font(
            pdf,
            family,
            emphasize and family_has_bold,
            layout.value_block.font_size,
        )
        pdf.set_xy(x0 + label_w, y0)
        pdf.multi_cell(
            value_w,
            layout.value_block.line_height,
            # 2026-09-03 수정: 죽어 있던 `_pdf_safe_value` 를 실제 렌더 경로에 연결.
            _pdf_safe_value(val, family),
            border="TR",
            align="L",
        )
        y_after_value = pdf.get_y()

        y_bottom = max(y_after_label, y_after_value, y0 + layout.row_height)
        if y_after_label < y_bottom:
            pdf.set_xy(x0, y_after_label)
            pdf.cell(label_w, y_bottom - y_after_label, "", border="L")
        if y_after_value < y_bottom:
            pdf.set_xy(x0 + label_w, y_after_value)
            pdf.cell(value_w, y_bottom - y_after_value, "", border="R")
        pdf.set_xy(x0, y_bottom)
        # 2026-09-03 수정: fpdf2 deprecated `ln=1` → new_x/new_y.
        pdf.cell(w, 0, "", border="T", new_x=XPos.LMARGIN, new_y=YPos.NEXT)


def _pdf_kv_rows(
    pdf: FPDF,
    family: str,
    rows: list[tuple[str, str]],
    w: float,
    *,
    family_has_bold: bool = True,
    gap_between_lines: int = 0,
) -> None:
    """Render key-value rows as a DOCX-style two-column table."""
    _ = gap_between_lines  # table rows carry their own spacing
    _pdf_add_table(pdf, family, rows, w, family_has_bold=family_has_bold)


def _add_appendix_line_docx(doc: Document, line: str) -> None:
    """Appendix species row — bold class/classification; SMILES stays regular."""
    stripped = line.strip()
    if " | SMILES:" in stripped:
        head, tail = stripped.split(" | SMILES:", 1)
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(2)
        r0 = p.add_run(head)
        _set_run_font(r0, 9, bold=True)
        r1 = p.add_run(f" | SMILES:{tail}")
        _set_run_font(r1, 9)
    else:
        _add_para(doc, stripped, size=9, bold=True, space_after=2)


def _appendix_species_class_label(class_num: object) -> str:
    try:
        info = ICH_M7_REGULATORY[int(class_num)]
        return f"Class {class_num} — {info['name']}"
    except (TypeError, ValueError, KeyError):
        return f"Class {class_num}"


def _appendix_species_detail_rows(row: dict, index: int) -> list[tuple[str, str]]:
    alerts = row.get("alerts") or []
    alert_text = ", ".join(str(a) for a in alerts) if alerts else "None detected"
    gen = row.get("generation")
    parent = str(row.get("parent_smiles") or "").strip()
    triggered = str(row.get("triggered_by") or row.get("name") or "").strip()
    rows: list[tuple[str, str]] = [
        ("Node ID", f"N.{index:03d}"),
        ("SMILES", _format_report_value(row.get("smiles"))),
    ]
    if gen is not None:
        rows.append(("Generation", str(gen)))
    if parent:
        rows.append(("Parent SMILES", parent))
    if triggered:
        rows.append(("Transformation / mechanism", triggered))
    rows.extend(
        [
            ("Name / pathway", _format_report_value(row.get("name"))),
            ("ICH M7 hazard class", _appendix_species_class_label(row.get("ich_m7_class", "—"))),
            ("Classification", _format_report_value(row.get("classification"))),
            ("Structural alerts", alert_text),
            ("Log P", _format_report_value(row.get("logp"))),
            ("Process stability", _format_report_value(row.get("stability"))),
        ]
    )
    try:
        from p52_ambit_client import build_oecd_assessment

        oecd = build_oecd_assessment(str(row.get("smiles") or ""))
        rows.extend(
            [
                ("OECD QSAR status", oecd.get("oecd_connection_status", "N/A")),
                ("OECD hydrolysis", oecd.get("oecd_hydrolysis", "N/A")),
                ("OECD biodegradability", oecd.get("oecd_biodegradability", "N/A")),
                ("OECD mutagenicity alert", oecd.get("oecd_mutagenicity_alert", "N/A")),
            ]
        )
    except Exception:
        rows.append(("OECD QSAR status", "Unavailable"))
    return rows


def _oecd_qsar_report_rows(ctx: AssessmentContext) -> list[tuple[str, str]]:
    return [
        ("OECD QSAR engine", ctx.qsar_engine),
        ("AMBIT connection", ctx.oecd_connection_status),
        ("OECD hydrolysis stability", ctx.oecd_hydrolysis),
        ("OECD biodegradability", ctx.oecd_biodegradability),
        ("OECD mutagenicity alert", ctx.oecd_mutagenicity_alert),
        ("Ames / (Q)SAR summary", ctx.qsar_ames),
    ]


def _appendix_species_png(row: dict) -> bytes | None:
    smiles = str(row.get("smiles") or "").strip()
    if not smiles:
        return None
    alerts = row.get("alerts") or []
    if alerts:
        return _structure_png_for_smiles(
            smiles,
            alerts=alerts,
            highlighted=True,
            size=_PDF_APPENDIX_IMAGE_SIZE,
        )
    return _structure_png_for_smiles(
        smiles,
        highlighted=False,
        size=_PDF_APPENDIX_IMAGE_SIZE,
    )


def _pdf_add_appendix_species_section(
    pdf: FPDF,
    family: str,
    species_rows: list[dict],
    w: float,
    *,
    family_has_bold: bool,
    parent_api_smiles: str = "",
) -> None:
    """Appendix A — 다단계 분해 트리 보고서와 동일한 노드별 구조식 + 경로 메타 + ICH M7."""
    if not species_rows:
        return

    _pdf_ensure_space(pdf, 20)
    pdf.ln(4)
    _pdf_section(
        pdf,
        family,
        "APPENDIX A — PROCESS-RELATED SPECIES (NODE-BY-NODE)",
        w,
    )
    _pdf_write_wrapped(
        pdf,
        family,
        (
            "Multi-step pathway layout (aligned with SMIRKS degradation tree reports): "
            "each node includes a 2D structure depiction, parent linkage, transformation "
            "rule, and ICH M7(R2) hazard class."
        ),
        w,
        size=8,
    )
    _pdf_blank_lines(pdf, 1)

    _pdf_section(pdf, family, "PATHWAY OVERVIEW", w)
    _pdf_write_wrapped(
        pdf,
        family,
        f"Flat node index (Generation 0 = Parent API · {len(species_rows)} node(s)):",
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 1)
    for idx, row in enumerate(species_rows, start=1):
        gen = row.get("generation", 1)
        trig = str(row.get("triggered_by") or row.get("name") or "Side reaction")
        ich_cls = row.get("ich_m7_class", "—")
        smi = str(row.get("smiles") or "")
        smi_short = smi[:48]
        suffix = "…" if len(smi) > 48 else ""
        overview = f"  {idx:02d}. Gen {gen} | {trig[:40]} | Class {ich_cls} | {smi_short}{suffix}"
        _pdf_write_wrapped(pdf, family, overview, w, size=8)
    _pdf_blank_lines(pdf, 2)

    _pdf_section(pdf, family, "NODE-BY-NODE STRUCTURE APPENDIX", w)
    _pdf_write_wrapped(
        pdf,
        family,
        "Each node includes a centered 2D chemical structure depiction and pathway metadata.",
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 2)

    for index, row in enumerate(species_rows, start=1):
        gen = row.get("generation", 1)
        triggered = str(row.get("triggered_by") or row.get("name") or "Process side product")
        smiles = str(row.get("smiles") or "")
        figure_ref = f"Figure A.{index}"

        _pdf_ensure_space(pdf, 90)
        _pdf_write_wrapped(
            pdf,
            family,
            f"Node {index} — Generation {gen} — {triggered}",
            w,
            bold=family_has_bold,
            size=11,
        )
        _pdf_blank_lines(pdf, 1)

        png = _appendix_species_png(row)
        if png:
            alerts = row.get("alerts") or []
            cap = (
                f"{figure_ref}: Structural alert zone highlighted in red."
                if alerts
                else f"{figure_ref}: 2D structure depiction — {smiles[:60]}{'…' if len(smiles) > 60 else ''}"
            )
            _pdf_add_structure_figure(pdf, family, png, w, caption=cap)
        else:
            _pdf_write_wrapped(
                pdf,
                family,
                f"{figure_ref}: Structure depiction unavailable.",
                w,
                size=9,
            )

        _pdf_kv_rows(
            pdf,
            family,
            _appendix_species_detail_rows(row, index),
            w,
            family_has_bold=family_has_bold,
        )

        try:
            m7_cls = int(row.get("ich_m7_class"))
        except (TypeError, ValueError):
            m7_cls = None
        if m7_cls is not None and 1 <= m7_cls <= 5:
            m7_info = ICH_M7_REGULATORY[m7_cls]
            pdf_add_m7_executive_callout(
                pdf,
                family,
                w,
                m7_cls,
                build_m7_executive_rows(
                    m7_cls,
                    classification=m7_info["name"],
                    definition=m7_info["definition"],
                    primary_route="PRIME-DS Hybrid — Expert Rules + (Q)SAR + GNN",
                    classification_basis=_appendix_species_class_label(m7_cls),
                    recommended_control=m7_info["control"],
                    hazard_label=f"NODE {index} ICH M7 OUTCOME",
                ),
                family_has_bold=family_has_bold,
                title=f"NODE {index} — ICH M7 Class {m7_cls}",
            )

        _pdf_blank_lines(pdf, 2)

    if parent_api_smiles:
        _pdf_write_wrapped(
            pdf,
            family,
            f"Parent API (root): {parent_api_smiles}",
            w,
            size=8,
            italic=True,
        )


def _add_appendix_species_section_docx(doc: Document, species_rows: list[dict]) -> None:
    if not species_rows:
        return

    doc.add_page_break()
    _add_para(
        doc,
        "APPENDIX A — ALL PREDICTED PROCESS-RELATED SPECIES",
        size=12,
        bold=True,
        space_after=6,
    )
    _add_para(
        doc,
        (
            "Each entry includes a 2D structure depiction, ICH M7(R2) hazard class, "
            "structural alert status, and key physicochemical properties."
        ),
        size=9,
        italic=True,
        space_after=8,
    )

    for index, row in enumerate(species_rows, start=1):
        gen = row.get("generation", 1)
        triggered = str(row.get("triggered_by") or row.get("name") or "Process side product")
        _add_para(
            doc,
            f"Node {index} — Generation {gen} — {triggered}",
            size=11,
            bold=True,
            space_after=4,
        )

        png = _appendix_species_png(row)
        if png:
            alerts = row.get("alerts") or []
            cap = (
                "Structural alert zone highlighted in red."
                if alerts
                else "2D structure depiction (no alert highlight)."
            )
            _add_structure_image_docx(doc, png, title="", caption=cap)

        t = doc.add_table(rows=0, cols=2)
        t.style = "Table Grid"
        for label, val in _appendix_species_detail_rows(row, index):
            _add_table_row(t, label, val, bold_label=label in ("ICH M7 hazard class", "Classification"))
        doc.add_paragraph()


def _pdf_appendix_line(pdf: FPDF, family: str, line: str, w: float, *, family_has_bold: bool) -> None:
    stripped = line.strip()
    if " | SMILES:" in stripped:
        head, tail = stripped.split(" | SMILES:", 1)
        _pdf_write_wrapped(pdf, family, head, w, bold=family_has_bold, size=8)
        _pdf_write_wrapped(pdf, family, f" | SMILES:{tail}", w, bold=False, size=8)
    else:
        _pdf_write_wrapped(pdf, family, stripped, w, bold=family_has_bold, size=8)
    pdf.ln(0.3)


def _set_run_font(
    run,
    size_pt: int = 11,
    *,
    bold: bool = False,
    italic: bool = False,
    color: RGBColor | None = None,
) -> None:
    font = REG_FONT
    run.font.name = font
    run.font.size = Pt(size_pt)
    run.font.bold = bold
    run.font.italic = italic
    if color is not None:
        run.font.color.rgb = color
    r_pr = run._element.get_or_add_rPr()
    r_fonts = r_pr.get_or_add_rFonts()
    for attr in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"):
        r_fonts.set(qn(attr), font)


def _apply_document_styles(doc: Document) -> None:
    for style_name in ("Normal", "Heading 1", "Heading 2", "Heading 3", "Title"):
        try:
            style = doc.styles[style_name]
            style.font.name = REG_FONT
            r_pr = style.element.get_or_add_rPr()
            r_fonts = r_pr.get_or_add_rFonts()
            for attr in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"):
                r_fonts.set(qn(attr), REG_FONT)
        except KeyError:
            pass


def _add_report_title_docx(doc: Document) -> None:
    """Prominent regulatory report title block."""
    _add_para(
        doc,
        REPORT_TITLE_LINE1,
        size=20,
        bold=True,
        align=WD_ALIGN_PARAGRAPH.CENTER,
        space_after=2,
        color=REPORT_TITLE_COLOR,
    )
    _add_para(
        doc,
        REPORT_TITLE_LINE2,
        size=18,
        bold=True,
        align=WD_ALIGN_PARAGRAPH.CENTER,
        space_after=4,
        color=REPORT_TITLE_COLOR,
    )
    _add_para(
        doc,
        "ICH M7(R2) — Complementary (Q)SAR and Expert Rule-Based Evaluation",
        size=12,
        italic=True,
        align=WD_ALIGN_PARAGRAPH.CENTER,
        space_after=4,
    )
    _add_para(
        doc,
        f"{PLATFORM_NAME} — {PLATFORM_FULL}",
        size=10,
        italic=True,
        align=WD_ALIGN_PARAGRAPH.CENTER,
        space_after=10,
    )


def _add_para(
    doc: Document,
    text: str,
    *,
    size: int = 11,
    bold: bool = False,
    italic: bool = False,
    align=None,
    space_after: int = 6,
    color: RGBColor | None = None,
) -> None:
    p = doc.add_paragraph()
    if align is not None:
        p.alignment = align
    p.paragraph_format.space_after = Pt(space_after)
    run = p.add_run(text)
    _set_run_font(run, size, bold=bold, italic=italic, color=color)


def _add_table_row(
    table,
    label: str,
    value: str,
    *,
    bold_label: bool = True,
    bold_value: bool | None = None,
) -> None:
    row = table.add_row()
    c0, c1 = row.cells[0], row.cells[1]
    c0.width = Inches(2.2)
    c1.width = Inches(4.3)
    r0 = c0.paragraphs[0].add_run(label)
    _set_run_font(r0, 10, bold=bold_label)
    display = _format_report_value(value)
    emphasize = bold_value if bold_value is not None else _should_emphasize_verdict(label, display)
    r1 = c1.paragraphs[0].add_run(display)
    _set_run_font(r1, 11 if emphasize else 10, bold=emphasize)


def _process_route_rows(ctx: AssessmentContext) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = [
        ("Reactant 1 (SMILES)", ctx.reactant1 or "N/A"),
        ("Reactant 2 (SMILES)", ctx.reactant2 or "N/A"),
        ("Reagent (SMILES)", ctx.reagent or "N/A"),
        (
            "Drug substance / API (assessed in this report)",
            ctx.product or ctx.smiles or "N/A",
        ),
    ]
    if ctx.reference_structure:
        rows.append(
            (
                "Reference structure (Class 4 comparison only)",
                f"{ctx.reference_structure} — not assessed in Sections 1–7",
            )
        )
    if ctx.primary_impurity_smiles or ctx.primary_impurity_name:
        rows.append(
            (
                "Representative predicted impurity",
                f"{ctx.primary_impurity_name or 'N/A'} — {ctx.primary_impurity_smiles or 'N/A'} (Appendix A)",
            )
        )
    rows.extend(
        [
            ("Nitrite source present", "Yes" if ctx.nitrite else "No"),
            ("Acidic conditions", "Yes" if ctx.acidic_condition else "No"),
            ("Formation mechanism", ctx.mechanism or "N/A"),
            ("Predicted impurity (SMILES)", ctx.predicted_impurity_smiles or "N/A"),
        ]
    )
    return rows


def report_to_regulatory_docx(ctx: AssessmentContext, *, appendix_text: str = "") -> bytes:
    """Build a structured regulatory Word document from assessment context."""
    doc = Document()
    _apply_document_styles(doc)

    sections = doc.sections
    for section in sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)

    m7 = ctx.m7
    class_num = m7.class_number if m7 else 5
    m7_info = ICH_M7_REGULATORY.get(class_num, ICH_M7_REGULATORY[5])
    # 2026-09-03 수정: cohort of concern(니트로사민 등)은 ICH M7(R2) §7.5에서
    # TTC 적용 대상에서 제외됩니다. 이 인자를 넘기지 않아 DOCX/PDF 경로만
    # Class 1·CoC 불순물에 TTC 1.5 µg/day 를 계속 제시하고 있었습니다.
    m7_control, q3_control, ai_limit = m7_action_flags(
        class_num,
        cohort_of_concern=bool(getattr(m7, "cohort_of_concern", False)),
    )
    alert_status = "POSITIVE" if ctx.alerts else "NEGATIVE"
    alerts_text = ", ".join(ctx.alerts) if ctx.alerts else "None detected"
    pub = ctx.pubchem_status or (
        f"PubChem CID {ctx.pubchem_cid}"
        if ctx.pubchem_cid
        else "No PubChem lookup performed"
    )

    # Title block
    _add_report_title_docx(doc)
    _add_evaluation_structure_figures_docx(doc, ctx)

    # Document control table
    _add_para(doc, "DOCUMENT CONTROL", size=12, bold=True, space_after=4)
    dc = doc.add_table(rows=0, cols=2)
    dc.style = "Table Grid"
    for label, val in [
        ("Report Number", ctx.report_id),
        ("Document Version", DOCUMENT_VERSION),
        ("Issue Date", ctx.issue_date),
        ("Software Platform", f"{PLATFORM_NAME} — {PLATFORM_FULL}"),
        ("Regulatory Basis", "ICH M7(R2); ICH Q3A(R2); ICH Q3B(R2)"),
        ("Substance Name", ctx.name),
        ("OECD QSAR status", ctx.oecd_connection_status),
        ("OECD hydrolysis", ctx.oecd_hydrolysis),
        ("OECD biodegradability", ctx.oecd_biodegradability),
        ("OECD mutagenicity alert", ctx.oecd_mutagenicity_alert),
    ]:
        _add_table_row(dc, label, val)
    doc.add_paragraph()

    _add_para(
        doc,
        "CONFIDENTIALITY: This document supports regulatory submission and internal GMP use. "
        "Computational results require qualified expert review and do not replace confirmatory testing.",
        size=10,
        italic=True,
        space_after=12,
    )

    # Section 1
    _add_para(doc, "1. SUBSTANCE IDENTIFICATION", size=12, bold=True, space_after=6)
    t1 = doc.add_table(rows=0, cols=2)
    t1.style = "Table Grid"
    cat = substance_category_label(ctx)
    for label, val in [
        ("Substance category", cat),
        ("Substance name", ctx.name),
        ("SMILES", ctx.smiles),
        ("Molecular weight", f"{ctx.mw} g/mol (RDKit)"),
    ]:
        _add_table_row(t1, label, val, bold_label=True)

    # Section 2 — 공정 반응식·불순물 형성 (단일 화합물 스캔에서는 제외)
    if includes_process_route_section(ctx):
        # 2026-09-03 수정: 절 제목을 TXT 생성기(p80:672)와 일치시킵니다.
        _add_para(
            doc,
            "2. PROCESS ROUTE AND IMPURITY FORMATION PATHWAY",
            size=12,
            bold=True,
            space_after=6,
        )
        t2 = doc.add_table(rows=0, cols=2)
        t2.style = "Table Grid"
        for label, val in _process_route_rows(ctx):
            _add_table_row(t2, label, val)
        doc.add_paragraph()
        # 2026-09-03 수정 (E-2): `_add_product_assessment_docx` 가 정의만 되어 있고
        #   호출되지 않아 '2.1 DRUG SUBSTANCE' 절이 DOCX 에서 통째로 빠졌습니다.
        #   같은 평가의 TXT(p80:699)에는 있고 DOCX/PDF 에만 없던 상태입니다.
        #   조건을 p80:699 와 정확히 동일하게 맞춰 호출합니다.
        if ctx.product_assessment and not ctx.is_final_product:
            _add_product_assessment_docx(doc, ctx.product_assessment)

    # Section 3
    _add_para(doc, "3. COMPLEMENTARY HYBRID ASSESSMENT (ICH M7(R2))", size=12, bold=True, space_after=6)
    # 2026-09-03 수정: 3.1·3.2 소제목도 p80(737·743·744)에 맞춥니다.
    _add_para(doc, "3.1 Expert Rule-Based Assessment (Derek-style)", size=11, bold=True, space_after=4)
    t3a = doc.add_table(rows=0, cols=2)
    t3a.style = "Table Grid"
    for label, val in [
        ("Engine", ctx.expert_engine),
        ("Structural alert result", alert_status),
        ("Alerts identified", alerts_text),
        ("Conclusion", ctx.expert_comment),
    ]:
        _add_table_row(t3a, label, val)
    doc.add_paragraph()

    # 2026-09-03 수정: p80 은 3.2 와 3.2a 를 별도 소제목으로 찍습니다. 동일하게 분리.
    _add_para(
        doc,
        "3.2 Statistical / (Q)SAR-Based Assessment (Sarah-style)",
        size=11,
        bold=True,
        space_after=2,
    )
    _add_para(doc, "3.2a OECD QSAR Toolbox (AMBIT2)", size=11, bold=True, space_after=4)
    t3b = doc.add_table(rows=0, cols=2)
    t3b.style = "Table Grid"
    for label, val in _oecd_qsar_report_rows(ctx):
        _add_table_row(t3b, label, val)
    _add_table_row(t3b, "External database", pub)
    doc.add_paragraph()

    _add_para(
        doc,
        "3.2b Hybrid Ensemble Matrix (Layers A / B / C)",
        size=11,
        bold=True,
        space_after=4,
    )
    t3b2 = doc.add_table(rows=0, cols=2)
    t3b2.style = "Table Grid"
    for label, val in hybrid_ensemble_report_rows(m7):
        _add_table_row(t3b2, label, val)
    doc.add_paragraph()

    _add_para(
        doc,
        # 2026-09-03 수정: 3.3 소제목·건수 문구를 p80(762)과 일치시킵니다.
        f"3.3 In-House Master Database Cross-Validation ({indexed_db_count_label(ctx)})",
        size=11,
        bold=True,
        space_after=4,
    )
    t3c = doc.add_table(rows=0, cols=2)
    t3c.style = "Table Grid"
    for label, val in [
        ("Match status", ctx.db_match),
        ("Reference", ctx.db_ref or "N/A"),
    ]:
        _add_table_row(t3c, label, val)
    doc.add_paragraph()

    option_summary = build_m7_regulatory_option_summary(ctx)
    _add_para(
        doc,
        "3.4 ICH M7(R2) Regulatory Option Determination and Classification Rationale",
        size=11,
        bold=True,
        space_after=4,
    )
    t3d = doc.add_table(rows=0, cols=2)
    t3d.style = "Table Grid"
    for label, val in [
        ("Primary assessment route", option_summary.primary_option),
        ("Options applied in silico", "; ".join(option_summary.applied_options)),
        (
            "Options NOT applied",
            "; ".join(option_summary.not_applied_options[:4]),
        ),
        (
            "Supporting open data",
            "; ".join(option_summary.supporting_data) if option_summary.supporting_data else "N/A",
        ),
    ]:
        _add_table_row(t3d, label, val)
    doc.add_paragraph()
    _add_para(doc, "Evidence chain for Class assignment:", size=10, bold=True, space_after=4)
    for i, ev in enumerate(option_summary.evidence_chain, 1):
        _add_para(doc, f"({i}) {ev}", size=10, space_after=3)

    # Section 4
    # 2026-09-03 수정: 절 제목을 TXT 생성기(p80:782)와 일치시킵니다.
    _add_para(
        doc,
        "4. PHYSICOCHEMICAL PROPERTIES AND EXPOSURE CONSIDERATIONS",
        size=12,
        bold=True,
        space_after=6,
    )
    t4 = doc.add_table(rows=0, cols=2)
    t4.style = "Table Grid"
    for label, val in [
        ("Log P", ctx.logp),
        ("Log D / pKa / Boiling point", ctx.logd_pka_bp),
        ("OECD hydrolysis (AMBIT)", ctx.oecd_hydrolysis),
        ("OECD biodegradability", ctx.oecd_biodegradability),
        ("Process stability", ctx.stability),
        ("Stability rationale", ctx.stability_detail),
        ("Exposure assessment", ctx.exposure_comment),
    ]:
        _add_table_row(t4, label, val)
    doc.add_paragraph()

    # Section 5 — executive ICH M7 outcome (highlighted callout)
    # 2026-09-03 수정: 절 제목을 TXT 생성기(p80:793)와 일치시킵니다.
    _add_para(
        doc,
        "5. FINAL ICH M7(R2) HAZARD CLASSIFICATION",
        size=12,
        bold=True,
        space_after=6,
    )
    class_basis = (
        option_summary.evidence_chain[-1]
        if option_summary.evidence_chain
        else "See Section 3.4"
    )
    docx_add_m7_executive_callout(
        doc,
        class_num,
        build_m7_executive_rows(
            class_num,
            classification=m7_info["name"],
            definition=m7_info["definition"],
            primary_route=option_summary.primary_option,
            classification_basis=class_basis,
            recommended_control=m7_info["control"],
        ),
    )
    doc.add_paragraph()

    # Section 6
    _add_para(doc, "6. RECOMMENDED ACTIONS AND ACCEPTABLE INTAKE", size=12, bold=True, space_after=6)
    actions = []
    if m7_control:
        actions.append("Control per ICH M7 mutagenic impurity limits")
    if q3_control:
        actions.append("Manage per ICH Q3A/Q3B general impurity limits")
    t6 = doc.add_table(rows=0, cols=2)
    t6.style = "Table Grid"
    _add_table_row(t6, "Regulatory action", "; ".join(actions) or "See classification")
    _add_table_row(t6, "Acceptable intake", ai_limit)
    doc.add_paragraph()

    # Section 7
    _add_para(doc, "7. CONCLUSION AND LIMITATIONS", size=12, bold=True, space_after=6)
    _add_para(doc, option_summary.conclusion_paragraph, size=10, space_after=6)
    _add_para(
        doc,
        "Results remain preliminary until qualified toxicological review, literature search, "
        "and where warranted, confirmatory bacterial mutagenicity (Ames) testing are completed.",
        size=10,
        space_after=8,
    )
    _add_para(
        doc,
        f"Report status: In silico assessment complete — expert sign-off required. "
        f"Generated by {PLATFORM_NAME} v{DOCUMENT_VERSION}.",
        size=10,
        italic=True,
    )

    if ctx.appendix_species:
        _add_appendix_species_section_docx(doc, ctx.appendix_species)
    elif appendix_text.strip():
        doc.add_page_break()
        for line in appendix_text.splitlines():
            if line.strip().startswith("="):
                _add_para(doc, line.strip("= ").strip(), size=11, bold=True, space_after=4)
            elif line.strip().startswith("["):
                _add_appendix_line_docx(doc, line)
            elif line.strip():
                _add_para(doc, line.strip(), size=9, space_after=2)

    buf = BytesIO()
    doc.save(buf)
    return buf.getvalue()


def parse_appendix_species_from_text(report_text: str) -> list[dict]:
    """TXT 부록 블록 → appendix_species 구조 (PDF/DOCX 리치 레이아웃 복원용)."""
    if "APPENDIX A" not in report_text:
        return []

    section = report_text.split("APPENDIX A", 1)[1]
    species: list[dict] = []
    current: dict | None = None

    header_re = re.compile(
        r"^\[(\d+)\]\s*(.+?)\s*\|\s*ICH M7 Class\s*(\d+)\s*\|\s*(.+)$"
    )
    for line in section.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("="):
            continue
        if stripped.startswith("APPENDIX"):
            continue

        header_match = header_re.match(stripped)
        if header_match:
            if current:
                species.append(current)
            current = {
                "name": header_match.group(2).strip(),
                "ich_m7_class": int(header_match.group(3)),
                "classification": header_match.group(4).strip(),
                "alerts": [],
                "logp": "N/A",
                "stability": "N/A",
                "smiles": "",
            }
            continue

        if current is None:
            continue
        if stripped.startswith("Log P:"):
            tail = stripped.replace("Log P:", "", 1).strip()
            if "|" in tail:
                logp_part, rest = tail.split("|", 1)
                current["logp"] = logp_part.strip()
                if "Stability:" in rest:
                    current["stability"] = rest.split("Stability:", 1)[1].strip()
            else:
                current["logp"] = tail
        elif stripped.startswith("Alerts:"):
            alert_part = stripped.replace("Alerts:", "", 1).strip()
            if alert_part.lower() in ("none", "none detected", "—"):
                current["alerts"] = []
            else:
                current["alerts"] = [a.strip() for a in alert_part.split(",") if a.strip()]
        elif stripped.startswith("SMILES:"):
            current["smiles"] = stripped.replace("SMILES:", "", 1).strip()

    if current:
        species.append(current)
    return species


def _resolve_assessment_context(
    report_text: str,
    ctx: AssessmentContext,
) -> tuple[AssessmentContext, str]:
    """appendix_species가 비어 있으면 TXT 부록에서 복원."""
    appendix = ""
    if "APPENDIX A" in report_text:
        appendix = "APPENDIX A" + report_text.split("APPENDIX A", 1)[1]
    if ctx.appendix_species:
        return ctx, appendix
    parsed = parse_appendix_species_from_text(report_text)
    if parsed:
        return replace(ctx, appendix_species=parsed), appendix
    return ctx, appendix


def report_to_docx(report_text: str, title: str = "PRIME-DS Report", ctx: AssessmentContext | None = None) -> bytes:
    if ctx is not None:
        effective_ctx, appendix = _resolve_assessment_context(report_text, ctx)
        return report_to_regulatory_docx(effective_ctx, appendix_text=appendix)
    return _legacy_plain_docx(report_text, title)


def _legacy_plain_docx(report_text: str, title: str) -> bytes:
    doc = Document()
    _apply_document_styles(doc)
    _add_para(doc, title, size=16, bold=True, align=WD_ALIGN_PARAGRAPH.CENTER)
    doc.add_paragraph()
    for line in report_text.splitlines():
        if not line.strip():
            doc.add_paragraph()
            continue
        is_heading = bool(re.match(r"^=+$", line.strip())) or line.strip().startswith("APPENDIX")
        _add_para(doc, line, size=11, bold=is_heading)
    buf = BytesIO()
    doc.save(buf)
    return buf.getvalue()


def report_to_regulatory_pdf(ctx: AssessmentContext, *, appendix_text: str = "") -> bytes:
    """Structured regulatory PDF (matches DOCX layout; avoids clipped plain-text lines)."""
    pdf = FPDF()
    pdf.set_margins(18, 18, 18)
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.add_page()
    w = pdf.epw
    family, family_has_bold = _register_pdf_fonts(pdf)

    m7 = ctx.m7
    class_num = m7.class_number if m7 else 5
    m7_info = ICH_M7_REGULATORY.get(class_num, ICH_M7_REGULATORY[5])
    # 2026-09-03 수정: cohort of concern(니트로사민 등)은 ICH M7(R2) §7.5에서
    # TTC 적용 대상에서 제외됩니다. 이 인자를 넘기지 않아 DOCX/PDF 경로만
    # Class 1·CoC 불순물에 TTC 1.5 µg/day 를 계속 제시하고 있었습니다.
    m7_control, q3_control, ai_limit = m7_action_flags(
        class_num,
        cohort_of_concern=bool(getattr(m7, "cohort_of_concern", False)),
    )
    alert_status = "POSITIVE" if ctx.alerts else "NEGATIVE"
    alerts_text = ", ".join(ctx.alerts) if ctx.alerts else "None detected"
    pub = ctx.pubchem_status or (
        f"PubChem CID {ctx.pubchem_cid}" if ctx.pubchem_cid else "No PubChem lookup performed"
    )
    actions = []
    if m7_control:
        actions.append("Control per ICH M7 mutagenic impurity limits")
    if q3_control:
        actions.append("Manage per ICH Q3A/Q3B general impurity limits")
    option_summary = build_m7_regulatory_option_summary(ctx)
    class_basis = (
        option_summary.evidence_chain[-1]
        if option_summary.evidence_chain
        else "See Section 3.4"
    )

    kv = lambda rows: _pdf_kv_rows(
        pdf,
        family,
        rows,
        w,
        family_has_bold=family_has_bold,
        gap_between_lines=_PDF_SECTION_GAP_LINES,
    )

    _pdf_write_report_title_block(pdf, family, w)
    pdf.ln(3)
    _pdf_add_evaluation_structure_figures(pdf, family, ctx, w)
    pdf.ln(2)

    _pdf_section(pdf, family, "DOCUMENT CONTROL", w)
    kv(
        [
            ("Report Number", ctx.report_id),
            ("Document Version", DOCUMENT_VERSION),
            ("Issue Date", ctx.issue_date),
            ("Software Platform", f"{PLATFORM_NAME} — {PLATFORM_FULL}"),
            ("Regulatory Basis", "ICH M7(R2); ICH Q3A(R2); ICH Q3B(R2)"),
            ("Substance Name", ctx.name),
            ("OECD QSAR status", ctx.oecd_connection_status),
            ("OECD hydrolysis", ctx.oecd_hydrolysis),
            ("OECD biodegradability", ctx.oecd_biodegradability),
            ("OECD mutagenicity alert", ctx.oecd_mutagenicity_alert),
        ],
    )
    _pdf_write_wrapped(
        pdf,
        family,
        "CONFIDENTIALITY: This document supports regulatory submission and internal GMP use. "
        "Computational results require qualified expert review and do not replace confirmatory testing.",
        w,
        size=9,
    )

    cat = substance_category_label(ctx)
    _pdf_section(pdf, family, "1. SUBSTANCE IDENTIFICATION", w)
    kv(
        [
            ("Substance category", cat),
            ("Substance name", ctx.name),
            ("SMILES", ctx.smiles),
            ("Molecular weight", f"{ctx.mw} g/mol (RDKit)"),
        ],
    )

    if includes_process_route_section(ctx):
        # 2026-09-03 수정: 절 제목을 TXT 생성기(p80:672)와 일치시킵니다.
        _pdf_section(pdf, family, "2. PROCESS ROUTE AND IMPURITY FORMATION PATHWAY", w)
        kv(_process_route_rows(ctx))
        # 2026-09-03 수정 (E-2): `_pdf_add_product_assessment` 가 정의만 되어 있고
        #   호출되지 않아 '2.1 DRUG SUBSTANCE' 절이 PDF 에서 통째로 빠졌습니다.
        #   조건을 p80:699 와 정확히 동일하게 맞춰 호출합니다.
        if ctx.product_assessment and not ctx.is_final_product:
            _pdf_add_product_assessment(
                pdf, family, ctx.product_assessment, w, family_has_bold=family_has_bold
            )

    _pdf_section(pdf, family, "3. COMPLEMENTARY HYBRID ASSESSMENT (ICH M7(R2))", w)
    # 2026-09-03 수정: 3.1·3.2 소제목도 p80(737·743·744)에 맞춥니다.
    _pdf_write_wrapped(pdf, family, "3.1 Expert Rule-Based Assessment (Derek-style)", w, bold=True, size=10)
    kv(
        [
            ("Engine", ctx.expert_engine),
            ("Structural alert result", alert_status),
            ("Alerts identified", alerts_text),
            ("Conclusion", ctx.expert_comment),
        ],
    )
    _pdf_blank_lines(pdf, 1)
    # 2026-09-03 수정: p80 은 3.2 와 3.2a 를 별도 소제목으로 찍습니다. 동일하게 분리.
    _pdf_write_wrapped(
        pdf,
        family,
        "3.2 Statistical / (Q)SAR-Based Assessment (Sarah-style)",
        w,
        bold=True,
        size=10,
    )
    _pdf_write_wrapped(pdf, family, "3.2a OECD QSAR Toolbox (AMBIT2)", w, bold=True, size=10)
    kv(_oecd_qsar_report_rows(ctx))
    kv([("External database", pub)])
    _pdf_blank_lines(pdf, 1)
    _pdf_write_wrapped(
        pdf,
        family,
        "3.2b Hybrid Ensemble Matrix (Layers A / B / C)",
        w,
        bold=True,
        size=10,
    )
    kv(hybrid_ensemble_report_rows(m7))
    _pdf_blank_lines(pdf, 1)
    _pdf_write_wrapped(
        pdf,
        family,
        # 2026-09-03 수정: 3.3 소제목·건수 문구를 p80(762)과 일치시킵니다.
        f"3.3 In-House Master Database Cross-Validation ({indexed_db_count_label(ctx)})",
        w,
        bold=True,
        size=10,
    )
    kv(
        [
            ("Match status", ctx.db_match),
            ("Reference", ctx.db_ref or "N/A"),
        ],
    )
    _pdf_blank_lines(pdf, 1)
    _pdf_write_wrapped(
        pdf,
        family,
        "3.4 ICH M7(R2) Regulatory Option Determination and Classification Rationale",
        w,
        bold=True,
        size=10,
    )
    kv(
        [
            ("Primary assessment route", option_summary.primary_option),
            ("Options applied in silico", "; ".join(option_summary.applied_options)),
            ("Options NOT applied", "; ".join(option_summary.not_applied_options[:4])),
            (
                "Supporting open data",
                "; ".join(option_summary.supporting_data) if option_summary.supporting_data else "N/A",
            ),
        ],
    )
    _pdf_blank_lines(pdf, 1)
    _pdf_write_wrapped(pdf, family, "Evidence chain for Class assignment:", w, bold=True, size=9)
    for i, ev in enumerate(option_summary.evidence_chain, 1):
        _pdf_write_wrapped(pdf, family, f"({i}) {ev}", w, size=9)
        _pdf_blank_lines(pdf, 1)

    # 2026-09-03 수정: 절 제목을 TXT 생성기(p80:782)와 일치시킵니다.
    _pdf_section(pdf, family, "4. PHYSICOCHEMICAL PROPERTIES AND EXPOSURE CONSIDERATIONS", w)
    kv(
        [
            ("Log P", ctx.logp),
            ("Log D / pKa / Boiling point", ctx.logd_pka_bp),
            ("OECD hydrolysis (AMBIT)", ctx.oecd_hydrolysis),
            ("OECD biodegradability", ctx.oecd_biodegradability),
            ("Process stability", ctx.stability),
            ("Stability rationale", ctx.stability_detail),
            ("Exposure assessment", ctx.exposure_comment),
        ],
    )

    # 2026-09-03 수정: 절 제목을 TXT 생성기(p80:793)와 일치시킵니다.
    _pdf_section(pdf, family, "5. FINAL ICH M7(R2) HAZARD CLASSIFICATION", w)
    pdf_add_m7_executive_callout(
        pdf,
        family,
        w,
        class_num,
        build_m7_executive_rows(
            class_num,
            classification=m7_info["name"],
            definition=m7_info["definition"],
            primary_route=option_summary.primary_option,
            classification_basis=class_basis,
            recommended_control=m7_info["control"],
        ),
        family_has_bold=family_has_bold,
    )

    _pdf_section(pdf, family, "6. RECOMMENDED ACTIONS AND ACCEPTABLE INTAKE", w)
    kv(
        [
            ("Regulatory action", "; ".join(actions) or "See classification"),
            ("Acceptable intake", ai_limit),
        ],
    )

    _pdf_section(pdf, family, "7. CONCLUSION AND LIMITATIONS", w)
    _pdf_write_wrapped(pdf, family, option_summary.conclusion_paragraph, w, size=9)
    _pdf_blank_lines(pdf, 1)
    _pdf_write_wrapped(
        pdf,
        family,
        "Results remain preliminary until qualified toxicological review, literature search, "
        "and where warranted, confirmatory bacterial mutagenicity (Ames) testing are completed.",
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 1)
    _pdf_write_wrapped(
        pdf,
        family,
        f"Report status: In silico assessment complete — expert sign-off required. "
        f"Generated by {PLATFORM_NAME} v{DOCUMENT_VERSION}.",
        w,
        size=9,
    )

    if ctx.appendix_species:
        _pdf_add_appendix_species_section(
            pdf,
            family,
            ctx.appendix_species,
            w,
            family_has_bold=family_has_bold,
            parent_api_smiles=str(ctx.smiles or ctx.product or ""),
        )
    elif appendix_text.strip():
        _pdf_ensure_space(pdf, 20)
        pdf.ln(4)
        _pdf_section(pdf, family, "APPENDIX A — ALL PREDICTED PROCESS-RELATED SPECIES", w)
        for line in appendix_text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("="):
                continue
            if stripped.startswith("APPENDIX"):
                continue
            if stripped.startswith("["):
                _pdf_appendix_line(pdf, family, stripped, w, family_has_bold=family_has_bold)

    out = pdf.output()
    if isinstance(out, bytearray):
        return bytes(out)
    if isinstance(out, bytes):
        return out
    return out.encode("latin-1")


def _plain_text_to_pdf(report_text: str, title: str) -> bytes:
    """Fallback plain-text PDF with safe wrapping (no full-width separator lines)."""
    pdf = FPDF()
    pdf.set_margins(18, 18, 18)
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.add_page()
    w = pdf.epw
    family, _ = _register_pdf_fonts(pdf)

    _pdf_write_wrapped(pdf, family, title, w, bold=True, size=13)
    pdf.ln(3)

    for line in report_text.splitlines():
        stripped = line.strip()
        if not stripped:
            pdf.ln(2)
            continue
        if re.fullmatch(r"[=\-_]{4,}", stripped):
            pdf.ln(2)
            continue
        collapsed = re.sub(r"\s{2,}", " ", stripped)
        is_heading = bool(re.match(r"^[0-9]+\.", collapsed)) or collapsed.startswith(
            ("DOCUMENT CONTROL", "APPENDIX", "CONFIDENTIALITY")
        )
        _pdf_write_wrapped(pdf, family, collapsed, w, bold=is_heading, size=10 if is_heading else 9)

    out = pdf.output()
    if isinstance(out, bytearray):
        return bytes(out)
    if isinstance(out, bytes):
        return out
    return out.encode("latin-1")


def report_to_pdf(
    report_text: str,
    title: str = "PRIME-DS Report",
    ctx: AssessmentContext | None = None,
) -> bytes:
    if ctx is not None:
        effective_ctx, appendix = _resolve_assessment_context(report_text, ctx)
        return report_to_regulatory_pdf(effective_ctx, appendix_text=appendix)
    return _plain_text_to_pdf(report_text, title)


def safe_filename_base(name: str, fallback: str = "report") -> str:
    base = "".join(c if c.isalnum() or c in "-_" else "_" for c in (name or fallback)[:50])
    return base.strip("_") or fallback
