"""
SMIRKS multi-step degradation tree — PDF report with 2D structures per node.
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import Any

from fpdf import FPDF

from p41_degradation_smirks_rules import flatten_degradation_tree
from p80_formal_assessment_report import (
    DOCUMENT_VERSION,
    ICH_M7_REGULATORY,
    PLATFORM_FULL,
    PLATFORM_NAME,
)
from p83_m7_executive_visuals import build_m7_executive_rows, pdf_add_m7_executive_callout
from p81_report_export import (
    _PDF_SECTION_GAP_LINES,
    _register_pdf_fonts,
    _pdf_add_structure_figure,
    _pdf_blank_lines,
    _pdf_kv_rows,
    _pdf_section,
    _pdf_write_wrapped,
)
from p82_pdf_typesetting import pdf_output_bytes, structure_png_bytes
from p84_structure_images import DEFAULT_REPORT_IMAGE_SIZE

_TREE_SAFETY = (
    "DISCLAIMER: This multi-step SMIRKS degradation tree is a computational pathway "
    "simulation for ICH Q1A/Q1B planning and ICH M7(R2) hazard screening support. "
    "Predicted nodes do not replace analytical identification (LC-MS), structure "
    "elucidation, or confirmatory mutagenicity testing."
)


# 2026-09-03 수정 (F-4): `_structure_png_bytes`(p86/p87 판 vs p88/p89 판) 와
#   `_pdf_output_bytes`(p81·p86·p87 축자 3벌) 를 p82 로 통합. 여기는 재export.
#   RDKit 실패 시의 PubChem 원격 조회는 p82 에서 명시적 opt-in 이 됐습니다
#   (PRIME_ALLOW_PUBCHEM_STRUCTURE_FETCH).
def _structure_png_bytes(smiles: str, *, size: str = DEFAULT_REPORT_IMAGE_SIZE) -> bytes | None:
    return structure_png_bytes(smiles, size=size)


def _pdf_output_bytes(pdf: FPDF) -> bytes:
    return pdf_output_bytes(pdf)


def _node_detail_rows(node: dict[str, Any], *, parent_smiles: str) -> list[tuple[str, str]]:
    smiles = str(node.get("name") or "")
    gen = node.get("generation", 0)
    triggered = str(node.get("triggered_by") or "")
    rows: list[tuple[str, str]] = [
        ("SMILES", smiles),
        ("Generation", str(gen)),
        ("Parent SMILES", parent_smiles or "(root)"),
        ("Transformation rule", triggered or "Parent API (root node)"),
    ]
    tox = node.get("toxicity_analysis") or {}
    if tox:
        rows.extend(
            [
                ("Ames AI probability", str(tox.get("Ames_AI_Probability", "N/A"))),
                ("AI prediction class", str(tox.get("Prediction_Class", "N/A"))),
                ("GNN model", str(tox.get("GNN_Model", "N/A"))),
                ("ICH M7 status", str(tox.get("ICH_M7_Status", "N/A"))),
            ]
        )
        rationale = str(tox.get("ICH_M7_Rationale") or "").strip()
        if rationale:
            rows.append(("ICH M7 rationale", rationale[:200]))
    return rows


# 2026-09-03 수정 (E-5): 노드 ICH M7 등급 복구 + 실제 규제 문구 조회.
_M7_CLASS_RE = re.compile(r"Class\s*([1-5])")


def _node_m7_class(tox: dict[str, Any]) -> int | None:
    """노드 toxicity_analysis 에서 ICH M7 등급(1-5)을 복구.

    p91._slim_toxicity_analysis 를 거친 트리에는 `ICH_M7_Class` 키가 없고
    `ICH_M7_Status` 문자열만 남습니다. 예전에는 int 검사만 해서 callout 이
    통째로 사라졌으므로 상태 문자열에서도 등급을 읽습니다.
    """
    raw = tox.get("ICH_M7_Class")
    if isinstance(raw, bool):
        raw = None
    if isinstance(raw, int) and 1 <= raw <= 5:
        return raw
    if isinstance(raw, str) and raw.strip().isdigit() and 1 <= int(raw.strip()) <= 5:
        return int(raw.strip())
    m = _M7_CLASS_RE.search(str(tox.get("ICH_M7_Status") or ""))
    if m:
        return int(m.group(1))
    return None


def _node_executive_rows(
    tox: dict[str, Any],
    m7_cls: int,
    *,
    node_index: int,
) -> list[tuple[str, str]]:
    """p83.build_degradant_executive_rows 와 동일한 ICH_M7_REGULATORY 조회."""
    m7 = ICH_M7_REGULATORY.get(m7_cls, ICH_M7_REGULATORY[5])
    status = str(tox.get("ICH_M7_Status") or f"Class {m7_cls}")
    rationale = str(tox.get("ICH_M7_Rationale") or "").strip()
    basis = (
        f"ICH M7 hazard class assignment logic: Class {m7_cls} - {m7['definition']}"
    )
    if rationale:
        basis = f"{basis} Node rationale: {rationale[:220]}"
    else:
        basis = f"{basis} Reported status: {status}"
    return build_m7_executive_rows(
        m7_cls,
        classification=m7["name"],
        definition=m7["definition"],
        primary_route="Hybrid Ensemble - Layer A/B/C + GNN",
        classification_basis=basis,
        recommended_control=m7["control"],
        hazard_label=f"NODE {node_index} ICH M7 OUTCOME",
    )


def _collect_tree_nodes(tree: dict[str, Any] | None) -> list[tuple[dict[str, Any], str]]:
    """DFS — (node, parent_smiles)."""
    out: list[tuple[dict[str, Any], str]] = []

    def _walk(node: dict[str, Any] | None, parent: str) -> None:
        if not node:
            return
        smi = str(node.get("name") or "")
        out.append((node, parent))
        for child in node.get("children") or []:
            _walk(child, smi)

    _walk(tree, "")
    return out


def export_degradation_tree_pdf(
    tree: dict[str, Any] | None,
    parent_smiles: str,
    *,
    max_depth: int = 3,
    include_ai: bool = True,
    report_id: str | None = None,
) -> bytes:
    """PDF — one section per tree node with 2D structure depiction."""
    if not tree:
        raise ValueError("Empty degradation tree")

    nodes = _collect_tree_nodes(tree)
    flat = flatten_degradation_tree(tree)
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    rid = report_id or f"PRIME-DS-TREE-{datetime.now().strftime('%Y%m%dT%H%M%S')}"
    mode = "AI toxicity + ICH M7 hybrid ensemble" if include_ai else "Structure-only SMIRKS"

    pdf = FPDF()
    pdf.set_margins(18, 18, 18)
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.add_page()
    w = pdf.epw
    family, family_has_bold = _register_pdf_fonts(pdf)

    def kv(rows: list[tuple[str, str]]) -> None:
        _pdf_kv_rows(
            pdf,
            family,
            rows,
            w,
            family_has_bold=family_has_bold,
            gap_between_lines=_PDF_SECTION_GAP_LINES,
        )

    _pdf_write_wrapped(
        pdf,
        family,
        "MULTI-STEP SMIRKS DEGRADATION TREE REPORT",
        w,
        bold=True,
        size=14,
    )
    _pdf_write_wrapped(
        pdf,
        family,
        "Recursive pathway simulation — 2D structure per node",
        w,
        size=10,
    )
    pdf.ln(3)

    _pdf_section(pdf, family, "DOCUMENT CONTROL", w)
    kv(
        [
            ("Report Number", rid),
            ("Document Version", DOCUMENT_VERSION),
            ("Generated", generated),
            ("Platform", f"{PLATFORM_NAME} — {PLATFORM_FULL}"),
            ("Parent SMILES", str(tree.get("name") or parent_smiles)),
            ("Maximum depth", str(max_depth)),
            ("Analysis mode", mode),
            ("Total nodes", str(len(nodes))),
        ]
    )
    _pdf_blank_lines(pdf, 2)

    _pdf_section(pdf, family, "TREE OVERVIEW", w)
    _pdf_write_wrapped(
        pdf,
        family,
        "Flat node index (Generation 0 = Parent API):",
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 1)
    for idx, row in enumerate(flat, start=1):
        gen = row.get("Generation", "")
        trig = row.get("Triggered_By") or "Parent"
        smi_short = str(row.get("SMILES", ""))[:48]
        # 2026-09-03 수정: U+2026(...)을 직접 생성하면 latin-1 코어 폰트로
        #   떨어졌을 때 PDF 생성이 깨집니다. ASCII "..." 로 씁니다.
        suffix = "..." if len(str(row.get("SMILES", ""))) > 48 else ""
        line = f"  {idx:02d}. Gen {gen} | {trig} | {smi_short}{suffix}"
        if include_ai and row.get("ICH_M7_Status"):
            line += f" | {row.get('ICH_M7_Status')}"
        _pdf_write_wrapped(pdf, family, line, w, size=8)
    _pdf_blank_lines(pdf, 2)

    _pdf_section(pdf, family, "NODE-BY-NODE STRUCTURE APPENDIX", w)
    _pdf_write_wrapped(
        pdf,
        family,
        "Each node includes a 2D chemical structure depiction and pathway metadata.",
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 2)

    for idx, (node, parent_smi) in enumerate(nodes, start=1):
        gen = node.get("generation", 0)
        triggered = str(node.get("triggered_by") or "Parent API")
        smiles = str(node.get("name") or "")
        figure_ref = f"Figure {idx}"

        _pdf_write_wrapped(
            pdf,
            family,
            f"Node {idx} — Generation {gen} — {triggered}",
            w,
            bold=True,
            size=11,
        )
        _pdf_blank_lines(pdf, 1)

        png = _structure_png_bytes(smiles, size="420x320")
        if png:
            _pdf_add_structure_figure(
                pdf,
                family,
                png,
                w,
                # 2026-09-03 수정: U+2026 -> ASCII "..." (위와 동일한 이유).
                caption=f"{figure_ref}: {smiles[:60]}{'...' if len(smiles) > 60 else ''}",
            )
        else:
            _pdf_write_wrapped(
                pdf,
                family,
                f"{figure_ref}: Structure depiction unavailable.",
                w,
                size=9,
            )

        kv(_node_detail_rows(node, parent_smiles=parent_smi))

        tox = node.get("toxicity_analysis") or {}
        # 2026-09-03 수정 (E-5): 세 가지를 함께 고칩니다.
        #   (1) definition 칸에 classification 과 같은 문자열을, recommended_control
        #       칸에 "Per ICH M7(R2) control strategy for assigned class" 라는
        #       플레이스홀더를 넣고 있었습니다. p80/p86 은 같은 칸에 실제 관리
        #       전략을 인쇄하므로 보고서마다 정보 수준이 달랐습니다.
        #       p83.build_degradant_executive_rows 와 동일한 ICH_M7_REGULATORY
        #       조회로 교체합니다.
        #   (2) `isinstance(m7_cls, int)` 만 보던 탓에 트리가
        #       p91._slim_toxicity_analysis 를 거치면 callout 전체가 예외 없이
        #       사라졌습니다. ICH_M7_Class → ICH_M7_Status 순으로 복구합니다.
        #   (3) 그래도 등급을 못 구하면 조용히 넘어가지 않고 사유를 인쇄합니다.
        if include_ai and tox:
            m7_cls = _node_m7_class(tox)
            if m7_cls is None:
                _pdf_write_wrapped(
                    pdf,
                    family,
                    f"Node {idx}: ICH M7 executive outcome omitted - no ICH M7 class "
                    "present in this node's toxicity analysis (the tree may have been "
                    "slimmed for API export).",
                    w,
                    size=9,
                )
            else:
                pdf_add_m7_executive_callout(
                    pdf,
                    family,
                    w,
                    m7_cls,
                    _node_executive_rows(tox, m7_cls, node_index=idx),
                    family_has_bold=family_has_bold,
                    title=f"NODE {idx} - ICH M7 Class {m7_cls}",
                )

        _pdf_blank_lines(pdf, 2)

    _pdf_section(pdf, family, "SAFETY DISCLAIMER", w)
    _pdf_write_wrapped(pdf, family, _TREE_SAFETY, w, size=9)

    return _pdf_output_bytes(pdf)
