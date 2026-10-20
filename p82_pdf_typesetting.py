"""
FPDF2 dynamic typesetting helpers for PRIME-DS regulatory PDF reports.

Measures wrapped text before drawing boxes so highlight containers always
fit their content. Supports auto-shrink for long ICH M7 / toxicology strings.

2026-09-03 수정: p81·p86·p87·p88·p89 에 복제되어 있던 PDF 공용 헬퍼를 이
모듈로 통합했습니다 (pdf_sanitize / pdf_output_bytes / structure_png_bytes /
average_mw_text / report_control_header_rows). 나머지 모듈은 얇은 재export만
남깁니다. 순환 import 방지를 위해 이 모듈은 p80/p81/p83 을 module-level 로
import 하지 않습니다 (의존 방향이 p80 -> p83 -> p82 이므로).
"""

# 글자가 넘치지 않게 함(PDF는 Word와 달리 자동 줄바꿈이 없음)

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass

from fpdf import FPDF

logger = logging.getLogger(__name__)

CELL_PAD_MM = 1.8
BOX_PAD_MM = 3.5
MIN_FONT_PT = 7.0
MAX_WRAP_LINES = 12
SHRINK_STEP_PT = 0.5

# 2026-09-03 수정: 예측된 불순물 SMILES 를 외부(PubChem REST)로 보내는 조회는
#   규제 도구 기본값에서 꺼져 있어야 합니다. 명시적 opt-in 으로 바꿉니다.
PUBCHEM_STRUCTURE_FETCH_ENV = "PRIME_ALLOW_PUBCHEM_STRUCTURE_FETCH"


def pubchem_structure_fetch_enabled() -> bool:
    """구조 이미지의 PubChem 원격 조회 허용 여부 (기본 False)."""
    return os.environ.get(PUBCHEM_STRUCTURE_FETCH_ENV, "").strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }


# 2026-09-03 수정: 좁은 PDF 폰트 / latin-1 코어 폰트에서 깨지거나 예외를 내는
#   문자를 실제 콘텐츠 기준으로 전부 치환합니다. 기존 버전은 U+2600-27BF 와
#   U+10000 이상만 지워서 U+2500 U+2265 U+2026 U+2460 U+03B1 그리고 한글이
#   전부 그대로 통과했습니다.
_SANITIZE_MAP = {
    "\u2014": "-",     # em dash
    "\u2013": "-",     # en dash
    "\u2012": "-",
    "\u2010": "-",
    "\u2011": "-",
    "\u2212": "-",     # minus sign
    "\u00b5": "u",     # micro sign
    "\u03bc": "u",     # greek small mu
    "\u2192": "->",
    "\u2190": "<-",
    "\u2194": "<->",
    "\u21d2": "=>",
    "\u2605": "*",     # black star
    "\u2606": "*",
    "\u25b6": ">",     # black right-pointing triangle
    "\u25c0": "<",
    "\u2022": "-",     # bullet
    "\u2026": "...",   # horizontal ellipsis
    "\u2264": "<=",
    "\u2265": ">=",
    "\u2260": "!=",
    "\u2248": "~",
    "\u00d7": "x",
    "\u00f7": "/",
    "\u2018": "'",
    "\u2019": "'",
    "\u201c": '"',
    "\u201d": '"',
    "\u00a0": " ",
    "\u200b": "",
    "\ufeff": "",
    "\u2070": "^0", "\u00b9": "^1", "\u00b2": "^2", "\u00b3": "^3",
    "\u2074": "^4", "\u2075": "^5", "\u2076": "^6", "\u2077": "^7",
    "\u2078": "^8", "\u2079": "^9", "\u207a": "+", "\u207b": "-",
    "\u03b1": "alpha", "\u03b2": "beta", "\u03b3": "gamma", "\u03b4": "delta",
    "\u0394": "Delta", "\u03bb": "lambda", "\u03c0": "pi", "\u03c3": "sigma",
    "\u03a9": "Ohm", "\u00c5": "A", "\u212b": "A",
}

# 박스 그리기 문자(U+2500-257F) — 세로선 계열은 '|', 나머지는 '-'.
_BOX_VERTICAL = frozenset("\u2502\u2503\u2551\u250a\u250b\u254e\u254f")


def _translate_char(ch: str) -> str:
    mapped = _SANITIZE_MAP.get(ch)
    if mapped is not None:
        return mapped
    code = ord(ch)
    if 0x2500 <= code <= 0x257F:          # box drawing
        return "|" if ch in _BOX_VERTICAL else "-"
    if 0x2580 <= code <= 0x259F:          # block elements
        return "#"
    if 0x2460 <= code <= 0x2473:          # circled 1..20
        return f"({code - 0x2460 + 1})"
    if 0x2474 <= code <= 0x2487:          # parenthesized 1..20
        return f"({code - 0x2474 + 1})"
    if 0x2776 <= code <= 0x277F:          # dingbat negative circled 1..10
        return f"({code - 0x2776 + 1})"
    return ch


def pdf_latin1_fallback(text: str) -> str:
    """latin-1 코어 폰트(Helvetica) 최종 방어선 — 어떤 문자도 예외를 못 내게."""
    return text.encode("latin-1", "replace").decode("latin-1")


def pdf_sanitize(text: str, family: str | None = None) -> str:
    """PDF 안전 문자열.

    2026-09-03 수정:
      (1) 치환표를 실제 콘텐츠(U+2500 U+2265 U+2026 U+2460 U+03B1 등)까지 확대.
      (2) family == "Helvetica" (폰트 미발견 fallback) 이면 latin-1 로 강제
          인코딩해 fpdf2 의 FPDFUnicodeEncodingException 으로 죽지 않게 함.
    """
    cleaned = "".join(_translate_char(ch) for ch in str(text))
    # 이모지 / 보조평면 문자는 TTF 서브셋에 임베드되지 않는 경우가 많아 제거.
    cleaned = re.sub(r"[\U00010000-\U0010ffff]", "", cleaned)
    cleaned = re.sub(r"[\u2600-\u27BF\ufe0f]", "", cleaned)
    if family == "Helvetica":
        cleaned = pdf_latin1_fallback(cleaned)
    return cleaned


def pdf_set_font(pdf: FPDF, family: str, bold: bool, size: float) -> None:
    if family == "Helvetica":
        pdf.set_font("Helvetica", "B" if bold else "", size)
    else:
        pdf.set_font(family, "B" if bold else "", size)


def pdf_output_bytes(pdf: FPDF) -> bytes:
    """FPDF.output() 결과를 항상 bytes 로.

    2026-09-03 수정: p81·p86·p87 에 축자 동일본이 3벌 있었습니다 → 여기로 통합.
    """
    out = pdf.output()
    if isinstance(out, bytearray):
        return bytes(out)
    if isinstance(out, bytes):
        return out
    return out.encode("latin-1")


def structure_png_bytes(
    smiles: str,
    *,
    size: str | None = None,
    allow_pubchem: bool | None = None,
) -> bytes | None:
    """2D 구조 PNG — RDKit 우선, PubChem 조회는 명시적 opt-in 일 때만.

    2026-09-03 수정: p86/p87 판본과 p88/p89 판본(빈 SMILES 가드 유무)이 달라
    보고서 종류별로 구조 렌더링이 갈렸습니다 → 여기로 통합.
    또한 RDKit 실패 시 SMILES 를 PubChem REST 로 전송하던 동작을
    PRIME_ALLOW_PUBCHEM_STRUCTURE_FETCH=1 일 때만으로 제한합니다 (기본 생략).
    """
    from p84_structure_images import (
        DEFAULT_REPORT_IMAGE_SIZE,
        rdkit_clean_structure_png,
        structure_png,
    )

    smi = (smiles or "").strip()
    if not smi:
        return None
    img_size = size or DEFAULT_REPORT_IMAGE_SIZE
    png = rdkit_clean_structure_png(smi, size=img_size)
    if png:
        return png
    if allow_pubchem is None:
        allow_pubchem = pubchem_structure_fetch_enabled()
    if not allow_pubchem:
        logger.info(
            "구조 이미지 생략 — RDKit 렌더링 실패. PubChem 원격 조회는 기본 "
            "비활성입니다 (%s=1 로 명시적 허용).",
            PUBCHEM_STRUCTURE_FETCH_ENV,
        )
        return None
    return structure_png(None, smi, size=img_size)


MOLECULAR_WEIGHT_LABEL = "g/mol (RDKit average MW, Descriptors.MolWt)"


def average_mw_text(smiles: str) -> str:
    """평균 분자량(Descriptors.MolWt) 문자열.

    2026-09-03 수정: p88·p89 는 ExactMolWt(단동위원소), p86 은 MolWt(평균)를
    써서 라니티딘 분자량이 보고서에 따라 314.14 / 314.41 로 갈렸습니다.
    규제 문서 관례인 평균 분자량으로 통일합니다 (라벨 MOLECULAR_WEIGHT_LABEL).
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return "N/A"
    try:
        return f"{Descriptors.MolWt(mol):.2f}"
    except Exception:
        return "N/A"


def report_control_header_rows(
    *,
    report_id: str,
    generated: str,
    platform: str,
    document_version: str,
    assessment_type: str,
) -> list[tuple[str, str]]:
    """p88 · p89 문서관리표 공통 앞부분 (2026-09-03 수정: 중복 정의 제거)."""
    return [
        ("Report ID", report_id),
        ("Generated", generated),
        ("Software platform", platform),
        ("Document version", document_version),
        ("Assessment type", assessment_type),
    ]


def pdf_line_height(font_size_pt: float) -> float:
    """Line height in mm \u2014 scales with font size to avoid vertical clipping.

    2026-09-03 수정: p81 에 `max(4.2, size * 0.45)` 라는 두 번째 정의가 있어
    측정(p82)과 그리기(p81)가 최대 2.2% 어긋났습니다. p81 정의를 삭제하고
    이 함수 하나만 사용합니다.
    """
    return max(4.4, font_size_pt * 0.46)


@dataclass(frozen=True)
class PdfTextBlock:
    font_size: float
    line_height: float
    line_count: int
    block_height: float


def pdf_measure_wrapped_text(
    pdf: FPDF,
    family: str,
    *,
    bold: bool,
    width_mm: float,
    text: str,
    base_size_pt: float,
    min_size_pt: float = MIN_FONT_PT,
    family_has_bold: bool = True,
    max_lines: int = MAX_WRAP_LINES,
) -> PdfTextBlock:
    """Measure wrapped lines; shrink font until content fits or min size reached."""
    if width_mm <= 0:
        return PdfTextBlock(base_size_pt, pdf_line_height(base_size_pt), 1, pdf_line_height(base_size_pt))

    # 2026-09-03 수정: family 를 넘겨 Helvetica fallback 시 latin-1 로 강제 인코딩.
    content = pdf_sanitize(str(text).strip() or "N/A", family)
    size = float(base_size_pt)
    min_size = float(min_size_pt)

    while size >= min_size:
        lh = pdf_line_height(size)
        pdf_set_font(pdf, family, bold and family_has_bold, size)
        # 2026-09-03 수정: fpdf2 deprecated split_only=True -> dry_run/output.
        lines = pdf.multi_cell(width_mm, lh, content, dry_run=True, output="LINES")
        count = max(1, len(lines))
        if count <= max_lines or size <= min_size:
            return PdfTextBlock(size, lh, count, count * lh)
        size -= SHRINK_STEP_PT

    lh = pdf_line_height(min_size)
    pdf_set_font(pdf, family, bold and family_has_bold, min_size)
    # 2026-09-03 수정: fpdf2 deprecated split_only=True -> dry_run/output.
    lines = pdf.multi_cell(width_mm, lh, content, dry_run=True, output="LINES")
    count = max(1, len(lines))
    return PdfTextBlock(min_size, lh, count, count * lh)


def pdf_draw_wrapped_text(
    pdf: FPDF,
    family: str,
    *,
    bold: bool,
    x: float,
    y: float,
    width_mm: float,
    text: str,
    block: PdfTextBlock,
    family_has_bold: bool = True,
    align: str = "L",
    text_color: tuple[int, int, int] = (0, 0, 0),
) -> None:
    pdf.set_text_color(*text_color)
    pdf_set_font(pdf, family, bold and family_has_bold, block.font_size)
    pdf.set_xy(x, y)
    pdf.multi_cell(
        width_mm,
        block.line_height,
        # 2026-09-03 수정: family 를 넘겨 Helvetica fallback 시 latin-1 로 강제 인코딩.
        pdf_sanitize(str(text).strip() or "N/A", family),
        align=align,
    )


@dataclass(frozen=True)
class PdfTwoColumnRowLayout:
    label_block: PdfTextBlock
    value_block: PdfTextBlock
    row_height: float
    emphasize: bool


def pdf_measure_two_column_row(
    pdf: FPDF,
    family: str,
    *,
    label: str,
    value: str,
    label_width: float,
    value_width: float,
    label_size_pt: float = 10.0,
    value_size_pt: float = 10.0,
    value_bold: bool = False,
    emphasize: bool = False,
    family_has_bold: bool = True,
    row_pad_mm: float = CELL_PAD_MM,
    value_min_size_pt: float = MIN_FONT_PT,
) -> PdfTwoColumnRowLayout:
    label_block = pdf_measure_wrapped_text(
        pdf,
        family,
        bold=True,
        width_mm=label_width,
        text=label,
        base_size_pt=label_size_pt,
        family_has_bold=family_has_bold,
    )
    value_block = pdf_measure_wrapped_text(
        pdf,
        family,
        bold=value_bold,
        width_mm=value_width,
        text=value,
        base_size_pt=value_size_pt,
        min_size_pt=value_min_size_pt,
        family_has_bold=family_has_bold,
    )
    content_h = max(label_block.block_height, value_block.block_height)
    row_h = content_h + 2 * row_pad_mm
    return PdfTwoColumnRowLayout(
        label_block=label_block,
        value_block=value_block,
        row_height=row_h,
        emphasize=emphasize,
    )


def pdf_draw_two_column_row(
    pdf: FPDF,
    family: str,
    *,
    x: float,
    y: float,
    label: str,
    value: str,
    label_width: float,
    value_width: float,
    layout: PdfTwoColumnRowLayout,
    family_has_bold: bool = True,
    label_color: tuple[int, int, int] = (40, 40, 40),
    value_color: tuple[int, int, int] = (40, 40, 40),
    value_fill: tuple[int, int, int] | None = None,
    row_pad_mm: float = CELL_PAD_MM,
) -> None:
    inner_h = layout.row_height - 2 * row_pad_mm
    label_y = y + row_pad_mm + max(0.0, (inner_h - layout.label_block.block_height) / 2)
    value_y = y + row_pad_mm + max(0.0, (inner_h - layout.value_block.block_height) / 2)

    if value_fill is not None:
        pdf.set_fill_color(*value_fill)
        pdf.rect(x + label_width, y, value_width, layout.row_height, style="F")

    pdf_draw_wrapped_text(
        pdf,
        family,
        bold=True,
        x=x,
        y=label_y,
        width_mm=label_width,
        text=label,
        block=layout.label_block,
        family_has_bold=family_has_bold,
        text_color=label_color,
    )
    pdf_draw_wrapped_text(
        pdf,
        family,
        bold=layout.emphasize,
        x=x + label_width,
        y=value_y,
        width_mm=value_width,
        text=value,
        block=layout.value_block,
        family_has_bold=family_has_bold,
        text_color=value_color,
    )


def pdf_blend_rgb(
    a: tuple[int, int, int],
    b: tuple[int, int, int],
    *,
    weight_a: float = 0.18,
) -> tuple[int, int, int]:
    w = max(0.0, min(1.0, weight_a))
    return tuple(int(a[i] * w + b[i] * (1.0 - w)) for i in range(3))
