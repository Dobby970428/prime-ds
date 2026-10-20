"""
가혹 분해산물 종합 리포트 출력 — regulatory-style TXT · PDF · Excel.

Layout hierarchy (Parent-First, Degradants-Second):
  1. Document header & control (+ stress conditions)
  2. Primary — Parent API identification & executive hazard assessment
  3. Secondary appendix — Predicted stress degradants screening
  4. Safety disclaimer
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO
from typing import Any, Callable

from fpdf import FPDF
from openpyxl import Workbook
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from rdkit import Chem

from p23_chem_search_engine import safe_canonical_smiles
from p24_db_smiles_lookup import classify_smiles_against_db
from p43_degradation_predictor import (
    DegradantProduct,
    abundance_en,
    ich_m7_signal_label_en,
    mechanism_en,
    stress_environment_en,
    translate_degradant_text,
)
from p80_formal_assessment_report import DOCUMENT_VERSION, ICH_M7_REGULATORY, PLATFORM_FULL, PLATFORM_NAME
from p66_ich_m7 import classify_ich_m7
from p83_m7_executive_visuals import (
    build_degradant_executive_rows,
    build_m7_executive_rows,
    format_m7_executive_text_block,
    pdf_add_m7_executive_callout,
)
from p50_physicochemical import calculate_my_advanced_properties, calculate_opensource_properties
from p81_report_export import (
    _PDF_SECTION_GAP_LINES,
    _register_pdf_fonts,
    _pdf_add_structure_figure,
    _pdf_blank_lines,
    _pdf_kv_rows,
    _pdf_section,
    _pdf_write_wrapped,
)
from p82_pdf_typesetting import (
    MOLECULAR_WEIGHT_LABEL,
    pdf_output_bytes,
    structure_png_bytes,
)
from p84_structure_images import DEFAULT_REPORT_IMAGE_SIZE
from p62_toxicity_alerts import scan_ashby_tennant_alerts

_REPORT_LINE_WIDTH = 88
_SAFETY_DISCLAIMER = (
    "REGULATORY SUBMISSION DISCLAIMER: This in silico report supports ICH Q1A/Q1B forced "
    "degradation planning and ICH M7(R2) mutagenic impurity screening. Computational "
    "predictions — including RDKit physicochemical descriptors, Ashby-Tennant structural "
    "alerts, Tanimoto DB RAG confidence scores, and degradant pathway simulations — do not "
    "replace qualified analytical chemistry (HPLC/LC-MS), structure elucidation (NMR), or "
    "confirmatory bacterial mutagenicity testing. All results require expert toxicological "
    "review before use in regulatory filings or GMP decision-making."
)


@dataclass
class ParentSubstanceAssessment:
    """Parent API — primary report subject (executive hazard outcome)."""

    smiles: str
    canonical_smiles: str
    molecular_weight: float
    logp: str
    logd: str
    intrinsic_solubility: str
    tpsa: str
    db_status: str
    db_ref: str
    alerts: list[str] = field(default_factory=list)
    ashby_tennant_status: str = "NEGATIVE"
    ich_m7_class: int = 5
    ich_m7_name: str = ""
    ich_m7_definition: str = ""
    ich_m7_control: str = ""
    mutagenicity_signal: str = ""
    rdkit_note: str = ""


@dataclass
class DegradationReportContext:
    """PDF/Excel/TXT 생성기에 전달할 컨텍스트."""

    parent_smiles: str
    stress_environment: str
    temperature_c: float
    generated_at: str
    products: list[DegradantProduct]
    parent_assessment: ParentSubstanceAssessment
    platform: str = PLATFORM_NAME

    @property
    def report_id(self) -> str:
        stamp = self.generated_at.replace(" ", "T").replace(":", "").replace("-", "")
        return f"PRIME-DS-DEG-{stamp}"

    def to_dict(self) -> dict[str, Any]:
        pa = self.parent_assessment
        return {
            "platform": self.platform,
            "report_id": self.report_id,
            "parent_smiles": self.parent_smiles,
            "parent_canonical_smiles": pa.canonical_smiles,
            "parent_assessment": {
                "molecular_weight": pa.molecular_weight,
                "logp": pa.logp,
                "logd": pa.logd,
                "intrinsic_solubility": pa.intrinsic_solubility,
                "tpsa": pa.tpsa,
                "db_status": pa.db_status,
                "db_ref": pa.db_ref,
                "ich_m7_class": pa.ich_m7_class,
                "ich_m7_name": pa.ich_m7_name,
                "alerts": pa.alerts,
            },
            "stress_environment": self.stress_environment,
            "stress_environment_en": stress_environment_en(self.stress_environment),
            "temperature_c": self.temperature_c,
            "generated_at": self.generated_at,
            "product_count": len(self.products),
            "products": [p.to_report_row_en() for p in self.products],
        }


def _estimate_intrinsic_solubility(
    props: dict[str, float | int] | None,
    mol: Chem.Mol | None = None,
) -> str:
    """Delaney ESOL estimate of intrinsic aqueous solubility (non-experimental).

    2026-09-03 수정 (E-1) — 세 가지 오류를 함께 고쳤습니다.

    1) 단위가 1,000배 틀렸습니다.
       `mol/L × g/mol = g/L = mg/mL` 인데 `(10**log_s) * mw / 1000.0` 으로
       1000을 더 나눠 실제로는 g/mL 값을 만들어 놓고 문자열에는 "mg/mL"
       라고 적었습니다. (logP 2 · MW 300 · TPSA 90 → 0.0063 mg/mL 로 출력,
       같은 식의 올바른 값은 6.27 mg/mL). `/1000.0` 을 제거했습니다.

    2) `min(1.0, ...)` 클램프는 logS = 1.0, 즉 10 mol/L 라는 물리적으로
       불가능한 상한을 만들었습니다. 상한을 ESOL 학습영역 상단인 +0.5
       (약 3 mol/L)로 낮췄습니다.

    3) 계수 (0.5, −1, −0.003, +0.008) 은 ESOL 을 포함한 어떤 공표 모델과도
       일치하지 않는 출처 불명 값이었습니다. Delaney ESOL 로 교체했습니다.

       출처: J. S. Delaney, "ESOL: Estimating Aqueous Solubility Directly
       from Molecular Structure", J. Chem. Inf. Comput. Sci. 2004, 44(3),
       1000-1005.  logS = 0.16 - 0.63·clogP - 0.0062·MW + 0.066·RB - 0.74·AP
       (RB = rotatable bonds, AP = aromatic proportion = 방향족 원자수 / 중원자수)
    """
    if not props:
        return "N/A — structure parse failure"
    try:
        logp = float(props.get("LogP", 0))
        mw = float(props.get("MW", 0))
        if mol is None:
            return (
                "N/A — ESOL requires rotatable-bond count and aromatic proportion "
                "(RDKit mol unavailable)"
            )
        from rdkit.Chem import Descriptors, Lipinski

        rb = int(Lipinski.NumRotatableBonds(mol))
        heavy = int(mol.GetNumHeavyAtoms())
        aromatic = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        ap = (aromatic / heavy) if heavy else 0.0
        if mw <= 0:
            mw = float(Descriptors.MolWt(mol))

        log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap
        log_s = max(-8.0, min(0.5, log_s))
        # mol/L × g/mol = g/L = mg/mL (추가 나눗셈 없음)
        mg_ml = (10**log_s) * mw
        return (
            f"logS ~ {log_s:.2f} (log mol/L, Delaney ESOL 2004); "
            f"intrinsic solubility ~ {mg_ml:.4f} mg/mL at 25 degC "
            f"(in silico estimate — not experimental)"
        )
    except (TypeError, ValueError):
        return "N/A — solubility estimation error"


def build_parent_substance_assessment(parent_smiles: str) -> ParentSubstanceAssessment:
    """Parent API hybrid characterization + expert rule / DB cross-reference."""
    raw = (parent_smiles or "").strip()
    canon, parse_err = safe_canonical_smiles(raw)
    canon = canon or raw
    mol = Chem.MolFromSmiles(canon) if canon else None
    rdkit_note = parse_err or ("RDKit parse failure" if mol is None and raw else "")

    props = calculate_opensource_properties(canon) if canon else None
    adv = calculate_my_advanced_properties(canon, skip_deepchem=True) if canon else None

    alerts: list[str] = []
    if mol is not None:
        alerts, _ = scan_ashby_tennant_alerts(mol)
    m7 = classify_ich_m7(mol, alerts) if mol is not None else classify_ich_m7(None, [])
    m7_info = ICH_M7_REGULATORY.get(m7.class_number, ICH_M7_REGULATORY[5])

    cls, db_rec = classify_smiles_against_db(canon) if canon else ("theoretical", None)
    db_status = (
        "Validated — match to in-house database (experimental reference)"
        if cls == "validated"
        else "Theoretical — computational prediction only"
    )
    db_ref = "N/A"
    if db_rec:
        db_ref = (
            f"{db_rec.get('drug_name') or db_rec.get('name', '—')} · "
            f"CID {db_rec.get('cid', '—')}"
        )

    logp_val = props.get("LogP", "N/A") if props else "N/A"
    logp = f"{logp_val} (RDKit MolLogP)" if props else "N/A"
    logd_raw = (adv or {}).get("LogD (My Algo)", "N/A")
    if logd_raw in ("추후 업데이트 예정", "pending", "", "N/A"):
        logd = (
            f"~{logp_val} (approx. from LogP at pH 7.4; pKa unavailable)"
            if props
            else "N/A"
        )
    else:
        logd = str(logd_raw)

    mw = float(props["MW"]) if props else exact_mw_fallback(canon, mol)
    tpsa = f"{props['TPSA']} A^2" if props else "N/A"

    return ParentSubstanceAssessment(
        smiles=raw,
        canonical_smiles=canon,
        molecular_weight=mw,
        logp=logp,
        logd=logd,
        # 2026-09-03 수정: ESOL 은 회전가능결합 수·방향족 비율이 필요하므로 mol 전달.
        intrinsic_solubility=_estimate_intrinsic_solubility(props, mol),
        tpsa=tpsa,
        db_status=db_status,
        db_ref=db_ref,
        alerts=list(alerts),
        ashby_tennant_status="POSITIVE" if alerts else "NEGATIVE",
        ich_m7_class=m7.class_number,
        ich_m7_name=m7_info["name"],
        ich_m7_definition=m7_info["definition"],
        ich_m7_control=m7_info["control"],
        mutagenicity_signal=ich_m7_signal_label_en(m7.class_number),
        rdkit_note=rdkit_note,
    )


def exact_mw_fallback(smiles: str, mol: Chem.Mol | None) -> float:
    # 2026-09-03 수정: 이름과 달리 평균 분자량(Descriptors.MolWt)을 돌려줍니다.
    #   규제 문서 관례가 평균 분자량이므로 p88·p89 도 여기에 맞췄습니다.
    #   (공개 API 유지를 위해 함수명은 그대로 둡니다.)
    if mol is not None:
        from rdkit.Chem import Descriptors

        return round(float(Descriptors.MolWt(mol)), 2)
    try:
        from p42_degradation_rag_engine import exact_mw_from_smiles

        return exact_mw_from_smiles(smiles)
    except Exception:
        return 0.0


def build_degradation_report_context(
    parent_smiles: str,
    stress_environment: str,
    temperature_c: float,
    products: list[DegradantProduct],
) -> DegradationReportContext:
    parent = (parent_smiles or "").strip()
    return DegradationReportContext(
        parent_smiles=parent,
        stress_environment=stress_environment,
        temperature_c=float(temperature_c),
        generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        products=products,
        parent_assessment=build_parent_substance_assessment(parent),
    )


def degradation_report_ready(ctx: DegradationReportContext) -> tuple[bool, str]:
    if ctx.stress_environment == "선택 안함":
        return False, "가혹 조건을 선택한 뒤 분해산물 예측을 실행해 주세요."
    if not ctx.products:
        return False, "예측된 분해산물이 없습니다. SMILES를 확인해 주세요."
    return True, "TXT·PDF·Excel 다운로드가 준비되었습니다."


def safe_degradation_filename_base(parent_smiles: str) -> str:
    short = re.sub(r"[^A-Za-z0-9]+", "_", parent_smiles[:24]).strip("_") or "compound"
    return f"PRIME_DS_Stress_Degradant_{short}"


# 2026-09-03 수정 (F-4): `_structure_png_bytes` 는 p86/p87 판본과 p88/p89 판본이
#   서로 달라(빈 SMILES 가드 유무) 보고서 종류별로 구조 렌더링이 갈렸고,
#   `_pdf_output_bytes` 는 p81·p86·p87 에 축자 3벌이었습니다.
#   둘 다 p82 로 통합하고 여기서는 얇은 재export 만 남깁니다.
#   아울러 RDKit 실패 시의 PubChem 원격 조회는 p82 에서 명시적 opt-in
#   (PRIME_ALLOW_PUBCHEM_STRUCTURE_FETCH) 으로 바뀌었습니다.
def _structure_png_bytes(smiles: str, *, size: str = DEFAULT_REPORT_IMAGE_SIZE) -> bytes | None:
    return structure_png_bytes(smiles, size=size)


def _pdf_output_bytes(pdf: FPDF) -> bytes:
    return pdf_output_bytes(pdf)


def _parent_executive_rows(pa: ParentSubstanceAssessment) -> list[tuple[str, str]]:
    alerts_text = ", ".join(pa.alerts) if pa.alerts else "None detected"
    basis = (
        f"ICH M7 hazard class assignment logic: Class {pa.ich_m7_class} — "
        f"{pa.ich_m7_definition} Ashby-Tennant alerts: {alerts_text}."
    )
    return build_m7_executive_rows(
        pa.ich_m7_class,
        classification=pa.ich_m7_name,
        definition=pa.ich_m7_definition,
        primary_route="Option 4 — Expert rule-based assessment (Ashby-Tennant + DB cross-reference)",
        classification_basis=basis,
        recommended_control=pa.ich_m7_control,
        hazard_label="FINAL HAZARD CLASS (Parent API)",
    )


def _parent_assessment_rows(pa: ParentSubstanceAssessment) -> list[tuple[str, str]]:
    alerts_text = ", ".join(pa.alerts) if pa.alerts else "None detected"
    rows = [
        ("Parent SMILES", pa.canonical_smiles or pa.smiles),
        ("Structure representation", "See Figure 1 (2D structure depiction)"),
        # 2026-09-03 수정: 어느 분자량인지 라벨에 명기(평균 분자량 Descriptors.MolWt).
        ("Molecular weight", f"{pa.molecular_weight} {MOLECULAR_WEIGHT_LABEL}"),
        ("LogP", pa.logp),
        ("LogD", pa.logd),
        ("Intrinsic solubility", pa.intrinsic_solubility),
        ("TPSA", pa.tpsa),
        ("In-house DB status", pa.db_status),
        ("In-house DB reference", pa.db_ref),
        ("Ashby-Tennant screening", pa.ashby_tennant_status),
        ("Ashby-Tennant alerts", alerts_text),
        ("Mutagenicity signal", pa.mutagenicity_signal),
    ]
    if pa.rdkit_note:
        rows.append(("RDKit parse note", pa.rdkit_note))
    return rows


def _document_control_rows(ctx: DegradationReportContext) -> list[tuple[str, str]]:
    return [
        ("Report Number", ctx.report_id),
        ("Document Version", DOCUMENT_VERSION),
        ("Issue Date", ctx.generated_at),
        ("Software Platform", f"{PLATFORM_NAME} — {PLATFORM_FULL}"),
        ("Regulatory Basis", "ICH M7(R2); ICH Q1A(R2); ICH Q1B(R2)"),
        ("Stress condition", stress_environment_en(ctx.stress_environment)),
        ("Stress temperature", f"{ctx.temperature_c:.1f} degC"),
        ("Study design reference", "ICH Q1A/Q1B forced degradation (in silico complement)"),
        ("Prediction engine", "PRIME-DS Smart Conditional Degradation Engine + DB RAG"),
    ]


def _degradant_summary_table_rows(products: list[DegradantProduct]) -> list[tuple[str, str, str, str, str]]:
    rows: list[tuple[str, str, str, str, str]] = []
    for product in products:
        rows.append(
            (
                product.degradant_id,
                translate_degradant_text(product.pathway_name),
                f"{product.molecular_weight}",
                f"Class {product.ich_m7_class}",
                f"{product.confidence_score:.1f}%",
            )
        )
    return rows


def _degradant_detail_rows(product: DegradantProduct, figure_ref: str) -> list[tuple[str, str]]:
    return [
        ("Degradant ID", product.degradant_id),
        ("Degradant name", translate_degradant_text(product.name)),
        ("SMILES", product.smiles),
        ("Structure representation", figure_ref),
        ("Degradation pathway", translate_degradant_text(product.pathway_name)),
        # 2026-09-03 수정: 분해산물 MW 는 p43/p42 의 exact_mw_from_smiles
        #   (= Descriptors.ExactMolWt, 단동위원소)에서 옵니다. Parent 행의
        #   평균 분자량과 정의가 다르므로 어느 쪽인지 라벨에 명기합니다.
        #   (p42·p43 은 이번 작업의 편집 대상이 아니라 값 자체는 그대로 둡니다.)
        ("Molecular weight", f"{product.molecular_weight} g/mol (RDKit monoisotopic, ExactMolWt)"),
        ("Predicted relative abundance", abundance_en(product.relative_abundance)),
        ("Confidence score (%)", f"{product.confidence_score:.1f}%"),
        ("Confidence tier", translate_degradant_text(product.confidence_label)),
        ("DB RAG cross-reference", product.db_match_detail),
        (
            "Ashby-Tennant alerts",
            ", ".join(product.alerts) if product.alerts else "None detected",
        ),
        ("Reaction mechanism", mechanism_en(product.mechanism)),
        ("BioTransformer route", translate_degradant_text(product.biotransformer_pathway)),
    ]


def generate_degradation_report_text(ctx: DegradationReportContext) -> str:
    """Plain-text regulatory report — Parent-First, Degradants-Second."""
    pa = ctx.parent_assessment
    alerts_text = ", ".join(pa.alerts) if pa.alerts else "None detected"
    sep = "=" * _REPORT_LINE_WIDTH
    dash = "-" * _REPORT_LINE_WIDTH

    lines = [
        sep,
        "     IN-SILICO STRESS DEGRADATION ASSESSMENT REPORT",
        "  ICH Q1A/Q1B Stress Study — Parent API Primary Assessment",
        sep,
        "",
        "DOCUMENT CONTROL",
        dash,
    ]
    for label, value in _document_control_rows(ctx):
        lines.append(f"  {label:<28} {value}")
    lines.extend(["", sep, "1. PRIMARY ASSESSMENT — TARGET PARENT SUBSTANCE (API)", sep, ""])

    lines.extend(
        [
            f"  Parent SMILES              {pa.canonical_smiles or pa.smiles}",
            # 2026-09-03 수정: 어느 분자량인지 명기 (평균 분자량 Descriptors.MolWt).
            f"  Molecular weight           {pa.molecular_weight} {MOLECULAR_WEIGHT_LABEL}",
            f"  LogP                       {pa.logp}",
            f"  LogD                       {pa.logd}",
            f"  Intrinsic solubility       {pa.intrinsic_solubility}",
            f"  TPSA                       {pa.tpsa}",
            f"  In-house DB status         {pa.db_status}",
            f"  In-house DB reference      {pa.db_ref}",
            f"  Ashby-Tennant screening    {pa.ashby_tennant_status}",
            f"  Ashby-Tennant alerts       {alerts_text}",
        ]
    )
    if pa.rdkit_note:
        lines.append(f"  RDKit parse note             {pa.rdkit_note}")
    lines.extend(format_m7_executive_text_block(pa.ich_m7_class, _parent_executive_rows(pa)))

    lines.extend(
        [
            "",
            sep,
            "APPENDIX A — PREDICTED STRESS DEGRADATION PRODUCTS SCREENING",
            sep,
            f"  Predicted degradant count    {len(ctx.products)}",
            f"  Stress condition             {stress_environment_en(ctx.stress_environment)}",
            f"  Stress temperature           {ctx.temperature_c:.1f} degC",
            "",
            "  Summary table:",
            f"  {'ID':<12} {'Pathway':<28} {'MW':<10} {'ICH M7':<10} {'Conf.':<8}",
            f"  {'-' * 12} {'-' * 28} {'-' * 10} {'-' * 10} {'-' * 8}",
        ]
    )
    for did, pathway, mw, m7c, conf in _degradant_summary_table_rows(ctx.products):
        lines.append(f"  {did:<12} {pathway[:28]:<28} {mw:<10} {m7c:<10} {conf:<8}")

    for idx, product in enumerate(ctx.products, start=1):
        lines.extend(["", dash, f"  A.{idx}  {product.degradant_id} — {translate_degradant_text(product.name)}", dash])
        for label, value in _degradant_detail_rows(product, f"Appendix Figure A.{idx}"):
            lines.append(f"  {label:<28} {value}")
        lines.extend(
            format_m7_executive_text_block(
                product.ich_m7_class,
                build_degradant_executive_rows(product),
                title=f"DEGRADANT EXECUTIVE OUTCOME — {product.degradant_id}",
            )
        )

    lines.extend(
        [
            "",
            sep,
            "SAFETY DISCLAIMER",
            sep,
            f"  {_SAFETY_DISCLAIMER}",
            "",
        ]
    )
    return "\n".join(lines)


def export_degradation_txt(ctx: DegradationReportContext) -> bytes:
    return generate_degradation_report_text(ctx).encode("utf-8")


def _write_excel_summary_sheet(ws, ctx: DegradationReportContext, *, header_fill, bold) -> None:
    pa = ctx.parent_assessment
    summary_rows: list[tuple[str, object]] = [
        ("Report Number", ctx.report_id),
        ("Platform", f"{PLATFORM_NAME} — {PLATFORM_FULL}"),
        ("Document Version", DOCUMENT_VERSION),
        ("Generated", ctx.generated_at),
        ("Regulatory Basis", "ICH Q1A(R2); ICH Q1B(R2); ICH M7(R2)"),
        ("Stress Condition", stress_environment_en(ctx.stress_environment)),
        ("Temperature (degC)", ctx.temperature_c),
        ("--- PARENT SUBSTANCE (PRIMARY) ---", ""),
        ("Parent SMILES", pa.canonical_smiles or pa.smiles),
        ("Molecular weight (g/mol)", pa.molecular_weight),
        ("LogP", pa.logp),
        ("LogD", pa.logd),
        ("Intrinsic solubility", pa.intrinsic_solubility),
        ("In-house DB status", pa.db_status),
        ("FINAL HAZARD CLASS (Parent)", f"Class {pa.ich_m7_class}"),
        ("ICH M7 (Parent)", pa.ich_m7_name),
        ("Ashby-Tennant alerts (Parent)", ", ".join(pa.alerts) if pa.alerts else "None detected"),
        ("--- APPENDIX: DEGRADANTS ---", ""),
        ("Degradant Count", len(ctx.products)),
        ("Prediction Engine", "PRIME-DS Smart Conditional Degradation + DB RAG"),
    ]
    for row_idx, (label, value) in enumerate(summary_rows, start=1):
        ws.cell(row=row_idx, column=1, value=label).font = bold
        ws.cell(row=row_idx, column=1).fill = header_fill
        ws.cell(row=row_idx, column=2, value=value)


def export_degradation_excel(ctx: DegradationReportContext) -> bytes:
    """openpyxl — Summary (Parent-first) + Degradants appendix sheet."""
    wb = Workbook()
    ws_sum = wb.active
    ws_sum.title = "Summary"
    header_fill = PatternFill("solid", fgColor="E3F2FD")
    bold = Font(bold=True)

    _write_excel_summary_sheet(ws_sum, ctx, header_fill=header_fill, bold=bold)
    ws_sum.column_dimensions["A"].width = 32
    ws_sum.column_dimensions["B"].width = 72

    ws = wb.create_sheet("Appendix_Degradants")
    if ctx.products:
        columns = list(ctx.products[0].to_report_row_en().keys())
    else:
        columns = ["ID"]

    for col_idx, col_name in enumerate(columns, start=1):
        cell = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font = bold
        cell.fill = header_fill
        cell.alignment = Alignment(wrap_text=True, vertical="top")

    for row_idx, product in enumerate(ctx.products, start=2):
        row_data = product.to_report_row_en()
        for col_idx, col_name in enumerate(columns, start=1):
            ws.cell(row=row_idx, column=col_idx, value=row_data.get(col_name, ""))
            ws.cell(row=row_idx, column=col_idx).alignment = Alignment(
                wrap_text=True, vertical="top"
            )

    for col_idx in range(1, len(columns) + 1):
        ws.column_dimensions[get_column_letter(col_idx)].width = 20
    if "SMILES" in columns:
        ws.column_dimensions[get_column_letter(columns.index("SMILES") + 1)].width = 36

    buf = BytesIO()
    wb.save(buf)
    return buf.getvalue()


def export_degradation_excel_with_structures(ctx: DegradationReportContext) -> bytes:
    """openpyxl — Parent-first Summary + degradant structures appendix."""
    wb = Workbook()
    ws_sum = wb.active
    ws_sum.title = "Summary"
    header_fill = PatternFill("solid", fgColor="E3F2FD")
    bold = Font(bold=True)

    _write_excel_summary_sheet(ws_sum, ctx, header_fill=header_fill, bold=bold)
    ws_sum.column_dimensions["A"].width = 32
    ws_sum.column_dimensions["B"].width = 72

    ws = wb.create_sheet("Appendix_Degradants")
    columns = [
        "ID",
        "Degradation pathway",
        "Predicted degradant",
        "SMILES",
        "Molecular weight (g/mol)",
        "ICH M7 Class",
        "Confidence score (%)",
        "Confidence tier",
        "DB RAG match",
        "2D Structure",
    ]
    for col_idx, col_name in enumerate(columns, start=1):
        cell = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font = bold
        cell.fill = header_fill
        cell.alignment = Alignment(wrap_text=True, vertical="top")

    for row_idx, product in enumerate(ctx.products, start=2):
        row_data = [
            product.degradant_id,
            translate_degradant_text(product.pathway_name),
            translate_degradant_text(product.name),
            product.smiles,
            product.molecular_weight,
            product.ich_m7_class,
            product.confidence_score,
            product.confidence_label,
            product.db_match_detail,
            "",
        ]
        for col_idx, val in enumerate(row_data, start=1):
            ws.cell(row=row_idx, column=col_idx, value=val)
            ws.cell(row=row_idx, column=col_idx).alignment = Alignment(
                wrap_text=True, vertical="top"
            )
        ws.row_dimensions[row_idx].height = 95
        png = _structure_png_bytes(product.smiles, size="220x165")
        if png:
            img_buf = BytesIO(png)
            img_buf.seek(0)
            img = XLImage(img_buf)
            img.width = 160
            img.height = 120
            ws.add_image(img, f"J{row_idx}")

    widths = [10, 28, 22, 36, 14, 10, 12, 28, 32, 24]
    for col_idx, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(col_idx)].width = w

    buf = BytesIO()
    wb.save(buf)
    return buf.getvalue()


def _pdf_kv(
    pdf: FPDF,
    family: str,
    w: float,
    rows: list[tuple[str, str]],
    *,
    family_has_bold: bool,
) -> None:
    _pdf_kv_rows(
        pdf,
        family,
        rows,
        w,
        family_has_bold=family_has_bold,
        gap_between_lines=_PDF_SECTION_GAP_LINES,
    )


def _pdf_document_control_section(
    pdf: FPDF,
    family: str,
    w: float,
    ctx: DegradationReportContext,
    kv: Callable[[list[tuple[str, str]]], None],
) -> None:
    _pdf_section(pdf, family, "DOCUMENT CONTROL", w)
    kv(_document_control_rows(ctx))


def _pdf_parent_primary_section(
    pdf: FPDF,
    family: str,
    w: float,
    ctx: DegradationReportContext,
    kv: Callable[[list[tuple[str, str]]], None],
    *,
    family_has_bold: bool,
) -> None:
    pa = ctx.parent_assessment
    parent_png = _structure_png_bytes(ctx.parent_smiles)
    if parent_png:
        _pdf_add_structure_figure(
            pdf,
            family,
            parent_png,
            w,
            caption="Figure 1: Parent Drug Substance (API) — Primary Assessment Subject",
        )

    _pdf_section(
        pdf,
        family,
        "1. PRIMARY ASSESSMENT — TARGET PARENT SUBSTANCE (API)",
        w,
    )
    _pdf_write_wrapped(
        pdf,
        family,
        "Executive hazard classification and hybrid physicochemical characterization of the "
        "input Parent drug substance precedes degradant screening in this report.",
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 1)
    kv(_parent_assessment_rows(pa))
    pdf_add_m7_executive_callout(
        pdf,
        family,
        w,
        pa.ich_m7_class,
        _parent_executive_rows(pa),
        family_has_bold=family_has_bold,
        title=f"PARENT API EXECUTIVE OUTCOME — ICH M7 Class {pa.ich_m7_class}",
    )


def _pdf_degradants_appendix_section(
    pdf: FPDF,
    family: str,
    w: float,
    ctx: DegradationReportContext,
    kv: Callable[[list[tuple[str, str]]], None],
    *,
    family_has_bold: bool,
) -> None:
    _pdf_section(
        pdf,
        family,
        "APPENDIX A — PREDICTED STRESS DEGRADATION PRODUCTS SCREENING",
        w,
    )
    kv(
        [
            ("Predicted degradant count", str(len(ctx.products))),
            ("Stress condition", stress_environment_en(ctx.stress_environment)),
            ("Stress temperature", f"{ctx.temperature_c:.1f} degC"),
        ]
    )
    _pdf_write_wrapped(
        pdf,
        family,
        "Summary table — degradant IDs, pathways, and preliminary ICH M7 classes:",
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 1)
    for did, pathway, mw, m7c, conf in _degradant_summary_table_rows(ctx.products):
        _pdf_write_wrapped(
            pdf,
            family,
            f"  {did} | {pathway} | MW {mw} | {m7c} | Confidence {conf}",
            w,
            size=9,
        )
    _pdf_blank_lines(pdf, 2)
    _pdf_write_wrapped(
        pdf,
        family,
        "Detailed degradant entries (structure, transformation path, alerts, ICH M7):",
        w,
        size=9,
    )
    _pdf_blank_lines(pdf, 1)

    for idx, product in enumerate(ctx.products, start=1):
        figure_ref = f"Appendix Figure A.{idx}"
        _pdf_write_wrapped(
            pdf,
            family,
            f"A.{idx}  {product.degradant_id} — {translate_degradant_text(product.name)}",
            w,
            bold=True,
            size=10,
        )
        _pdf_blank_lines(pdf, 1)

        degradant_png = _structure_png_bytes(product.smiles, size="400x300")
        if degradant_png:
            _pdf_add_structure_figure(
                pdf,
                family,
                degradant_png,
                w,
                caption=(
                    f"{figure_ref}: Predicted degradant — {product.degradant_id} "
                    f"(MW {product.molecular_weight} g/mol)"
                ),
            )
        else:
            _pdf_write_wrapped(
                pdf,
                family,
                f"{figure_ref}: Structure depiction unavailable for {product.degradant_id}.",
                w,
                size=9,
            )
        kv(_degradant_detail_rows(product, figure_ref))
        pdf_add_m7_executive_callout(
            pdf,
            family,
            w,
            product.ich_m7_class,
            build_degradant_executive_rows(product),
            family_has_bold=family_has_bold,
            title=f"DEGRADANT EXECUTIVE OUTCOME — {product.degradant_id} (Class {product.ich_m7_class})",
        )
        _pdf_blank_lines(pdf, 2)


def _pdf_safety_disclaimer_section(pdf: FPDF, family: str, w: float) -> None:
    _pdf_section(pdf, family, "SAFETY DISCLAIMER", w)
    _pdf_write_wrapped(pdf, family, _SAFETY_DISCLAIMER, w, size=9)


def export_degradation_pdf(ctx: DegradationReportContext) -> bytes:
    """Regulatory PDF — Parent-First primary assessment, degradants in Appendix A."""
    pdf = FPDF()
    pdf.set_margins(18, 18, 18)
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.add_page()
    w = pdf.epw
    family, family_has_bold = _register_pdf_fonts(pdf)

    kv = lambda rows: _pdf_kv(pdf, family, w, rows, family_has_bold=family_has_bold)

    _pdf_write_wrapped(
        pdf,
        family,
        "IN-SILICO STRESS DEGRADATION ASSESSMENT REPORT",
        w,
        bold=True,
        size=14,
    )
    _pdf_write_wrapped(
        pdf,
        family,
        "Parent API Primary Assessment — ICH Q1A/Q1B Stress Degradant Screening Appendix",
        w,
        size=10,
    )
    pdf.ln(3)

    _pdf_document_control_section(pdf, family, w, ctx, kv)
    _pdf_parent_primary_section(
        pdf, family, w, ctx, kv, family_has_bold=family_has_bold
    )
    _pdf_degradants_appendix_section(
        pdf, family, w, ctx, kv, family_has_bold=family_has_bold
    )
    _pdf_safety_disclaimer_section(pdf, family, w)

    return _pdf_output_bytes(pdf)
