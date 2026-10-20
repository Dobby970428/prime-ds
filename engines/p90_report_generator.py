"""
PRIME-DS Process Report — Excel · PDF 리포트 생성 엔진 (Phase 8 통합 양식).

공정 유연물질 + ICH M7 Class 1~5 + AI 독성 스코어를 규제 제출용 Excel/PDF로보낸다.
"""

from __future__ import annotations

import base64
import html
import io
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from engines.p03_branding import UI_ENGINE_TITLE, log_process_warning
from engines.p13_chemical_registry import format_reagent_context_summary
from engines.p33_process_impurity_engine import (
    IMP_COL_CONFIDENCE,
    IMP_COL_ICH_M7,
    IMP_COL_ICH_M7_CLASS,
    IMP_COL_MW,
    IMP_COL_PATHWAY,
    IMP_COL_REGULATORY,
    IMP_COL_SMILES,
    IMP_COL_SOURCE,
    ProcessImpurityEngine,
)
from engines.p37_retro_synthesis_process import MW_DISPLAY_DECIMALS
from engines.p93_weasyprint_support import configure_weasyprint_dll_paths
from p84_structure_images import rdkit_png

configure_weasyprint_dll_paths()
try:
    import weasyprint
except Exception:
    weasyprint = None  # type: ignore[assignment,misc]

from fpdf import FPDF

from p80_formal_assessment_report import (
    PLATFORM_FULL,
    PLATFORM_NAME,
    REPORT_TITLE_LINE1,
    REPORT_TITLE_LINE2,
    AssessmentContext,
    build_process_product_context,
)
from p81_report_export import (
    _pdf_add_structure_figure,
    _pdf_blank_lines,
    _pdf_ensure_space,
    _pdf_kv_rows,
    _pdf_output_bytes,
    _pdf_section,
    _pdf_write_centered,
    _pdf_write_wrapped,
    _register_pdf_fonts,
    report_to_regulatory_pdf,
)
from p84_structure_images import rdkit_clean_structure_png

# ── Phase 8 통합 보고서 상수 ───────────────────────────────────────────────────
REPORT_BRAND_TITLE = "PRIME-DS Process Report"
INTEGRATED_REPORT_TITLE = (
    "PRIME-DS 인실리코 공정 유연물질 및 유전독성 가이드라인(ICH M7) 통합 평가 보고서"
)
EXCEL_SHEET_NAME = "PRIME-DS ICH M7 Report"
EXCEL_FILE_NAME = "PRIME_DS_ICH_M7_Integrated_Report.xlsx"
PDF_FILE_NAME = "PRIME_DS_ICH_M7_Integrated_Report.pdf"

RPT_COL_IMAGE = "2D Structure"
RPT_COL_SMILES = "1. 생성 유연물질 구조식(SMILES)"
RPT_COL_NAME_RXN = "2. 기원 명명 반응(Name Reaction)"
RPT_COL_TRIGGER = "3. 트리거 시약/조건"
RPT_COL_STRUCTURAL_ALERT = "4. Structural Alert 유무"
RPT_COL_AI_SCORE = "5. AI 독성 스코어"
RPT_COL_ICH_M7_CLASS = "6. 최종 ICH M7 Class (1~5)"
RPT_COL_REMARKS = "규제 대응 가이드 (비고)"

INTEGRATED_EXPORT_COLUMNS: list[str] = [
    RPT_COL_IMAGE,
    RPT_COL_SMILES,
    RPT_COL_NAME_RXN,
    RPT_COL_TRIGGER,
    RPT_COL_STRUCTURAL_ALERT,
    RPT_COL_AI_SCORE,
    RPT_COL_ICH_M7_CLASS,
    RPT_COL_REMARKS,
]

TTC_RELIEF_REMARK = (
    "임상 및 상용 생산 단계에서 일반 불순물 기준(TTC 해제)으로 관리 가능"
)

EXCEL_IMAGE_COL_IDX = 0
EXCEL_ROW_HEIGHT = 100.0
EXCEL_IMAGE_COL_WIDTH = 18.0
EXCEL_MOL_IMAGE_SIZE: tuple[int, int] = (120, 120)
EXCEL_IMAGE_X_OFFSET = 5
EXCEL_IMAGE_Y_OFFSET = 5
EXCEL_COL_WIDTHS: tuple[float, ...] = (16.0, 38.0, 22.0, 28.0, 20.0, 16.0, 14.0, 42.0)
EXCEL_HEADER_ROW = 2
EXCEL_DATA_START_ROW = 3


def _regulatory_remarks(ich_class: int | None, existing: str | None) -> str:
    """Class 4/5 → TTC 해제 규제 대응 가이드 자동 삽입."""
    base = str(existing or "").strip()
    if ich_class in (4, 5):
        if TTC_RELIEF_REMARK in base:
            return base
        return f"{base} | {TTC_RELIEF_REMARK}" if base else TTC_RELIEF_REMARK
    return base


def _format_structural_alert_cell(record: dict[str, Any]) -> str:
    if record.get("alert_present") is False:
        return "No"
    alerts = record.get("structural_alerts")
    if isinstance(alerts, list) and alerts:
        return "Yes — " + ", ".join(str(a) for a in alerts[:5])
    if record.get("alert_present"):
        return "Yes"
    ich_detail = str(record.get("ich_m7_detail") or "")
    if "Alerts:" in ich_detail:
        return ich_detail.split("·")[0].strip() if "·" in ich_detail else ich_detail
    return "No"


def _format_ai_score_cell(record: dict[str, Any]) -> str:
    prob = record.get("ames_mutagenicity_prob")
    if prob is None:
        label = record.get("ames_prediction_label")
        return str(label) if label else "N/A"
    try:
        p = float(prob)
        verdict = "Positive" if p >= 0.5 else "Negative"
        return f"P(Mutagenic)={p:.3f} ({verdict})"
    except (TypeError, ValueError):
        return "N/A"


def _format_ich_m7_class_cell(record: dict[str, Any]) -> str:
    cls = record.get("ich_m7_class")
    if cls is not None:
        risk = "High Risk" if int(cls) <= 3 else ("De-risked" if int(cls) == 4 else "Safe")
        return f"Class {int(cls)} ({risk})"
    legacy = str(record.get("ich_m7") or record.get(IMP_COL_ICH_M7_CLASS, ""))
    return legacy or "—"


def _name_reaction_from_record(record: dict[str, Any]) -> str:
    name_rxn = record.get("name_reaction")
    if name_rxn:
        return str(name_rxn)
    pathway = str(record.get("pathway") or "")
    if "Name Rxn" in pathway or "Catalog ·" in pathway:
        for token in pathway.replace("·", "|").split("|"):
            token = token.strip()
            if token and token not in ("Process", "Catalog", "Name Rxn"):
                return token[:60]
    return "—"


def _trigger_context(
    record: dict[str, Any],
    process_conditions: dict[str, Any] | None,
) -> str:
    parts: list[str] = []
    pathway = str(record.get("pathway") or "").strip()
    if pathway:
        parts.append(pathway[:100])
    source = str(record.get("source_reactant") or "").strip()
    if source and source not in pathway:
        parts.append(source[:60])
    if process_conditions:
        ctx = process_conditions.get("process_reagents") or {}
        summary = format_reagent_context_summary(ctx)
        if summary and summary != "시약 미선택":
            parts.append(summary)
        stress = process_conditions.get("stress_condition")
        if stress and "선택안함" not in str(stress):
            parts.append(f"Stress: {stress}")
        temp = process_conditions.get("temperature_c")
        ph = process_conditions.get("ph")
        if temp is not None:
            parts.append(f"{temp}°C")
        if ph is not None:
            parts.append(f"pH {ph}")
    return " · ".join(parts) if parts else "—"


def _record_display_name(rec: dict[str, Any]) -> str:
    for key in (
        "mechanism_ui_text",
        "reaction_name",
        "pathway",
        "source_reactant",
        "name",
    ):
        val = str(rec.get(key) or "").strip()
        if val:
            return val[:120]
    return "Process impurity candidate"


def impurity_records_to_appendix_species(
    records: list[dict[str, Any]],
    *,
    parent_api_smiles: str = "",
) -> list[dict[str, Any]]:
    """Process Mode 예측 레코드 → 다단계 트리형 Appendix 노드."""
    species: list[dict[str, Any]] = []
    for rec in records:
        smi = str(rec.get("smiles") or "").strip()
        if not smi:
            continue
        alerts_raw = rec.get("structural_alerts")
        if isinstance(alerts_raw, list):
            alerts = [str(a) for a in alerts_raw if str(a).strip()]
        elif alerts_raw:
            alerts = [str(alerts_raw)]
        else:
            alerts = []
        logp = rec.get("logp") or rec.get("LogP")
        if logp is None and smi:
            try:
                from p50_physicochemical import calculate_opensource_properties

                props = calculate_opensource_properties(smi)
                logp = props.get("LogP", "N/A") if props else "N/A"
            except Exception:
                logp = "N/A"
        triggered = str(
            rec.get("reaction_name")
            or rec.get("mechanism_ui_text")
            or rec.get("pathway")
            or "Side reaction"
        ).strip()[:100]
        parent = str(rec.get("source_reactant") or parent_api_smiles or "").strip()
        role = str(rec.get("impurity_role") or rec.get("impurity_classification") or "")
        is_main = role == "Main_Product" or "main" in role.casefold()
        species.append(
            {
                "name": _record_display_name(rec),
                "smiles": smi,
                "classification": str(
                    rec.get("impurity_classification")
                    or rec.get("classification")
                    or "Theoretical — computational prediction"
                ),
                "ich_m7_class": rec.get("ich_m7_class") or 5,
                "logp": str(logp) if logp is not None else "N/A",
                "stability": str(rec.get("stability") or "See process conditions"),
                "alerts": alerts,
                "triggered_by": triggered,
                "parent_smiles": parent,
                "generation": 0 if is_main else 1,
            }
        )
    return species


def _process_mechanism_summary(process_conditions: dict[str, Any] | None) -> str:
    if not process_conditions:
        return "PRIME-DS Universal Side-Reaction + ICH M7 hybrid engine"
    parts: list[str] = []
    name_rxn = str(process_conditions.get("name_reaction") or "").strip()
    if name_rxn and name_rxn not in ("선택안함", "—", "None"):
        parts.append(name_rxn)
    stress = str(process_conditions.get("stress_condition") or "").strip()
    if stress and "선택안함" not in stress:
        parts.append(stress)
    ctx = process_conditions.get("process_reagents") or {}
    summary = format_reagent_context_summary(ctx)
    if summary and summary != "시약 미선택":
        parts.append(summary)
    return (
        " · ".join(parts)
        if parts
        else "PRIME-DS Universal Side-Reaction + ICH M7 hybrid engine"
    )


def build_process_mode_regulatory_context(
    target_props: dict[str, Any],
    records: list[dict[str, Any]],
    process_conditions: dict[str, Any] | None = None,
    *,
    reactant1: str = "",
    reactant2: str = "",
    reagent_smiles: str = "",
    reference_api: str = "",
    indexed_count: int = 13040,
) -> AssessmentContext:
    """Process Mode — 단일화합물 검색과 동일한 AssessmentContext + Appendix A."""
    from p66_ich_m7 import classify_ich_m7
    from p50_physicochemical import (
        calculate_my_advanced_properties,
        calculate_opensource_properties,
        predict_chemical_stability,
    )
    from p62_toxicity_alerts import scan_ashby_tennant_alerts

    target_smi = str(target_props.get("smiles") or "").strip()
    mol = Chem.MolFromSmiles(target_smi) if target_smi else None
    if mol is None:
        raise ValueError("Invalid target SMILES for regulatory report context.")

    p_alerts, _ = scan_ashby_tennant_alerts(mol)
    ref_mol = Chem.MolFromSmiles(reference_api) if reference_api else None
    p_m7 = classify_ich_m7(mol, p_alerts, drug_mol=ref_mol)
    p_props = calculate_opensource_properties(target_smi)
    p_adv = calculate_my_advanced_properties(target_smi, skip_deepchem=True)
    p_stab = predict_chemical_stability(target_smi, p_props)
    p_mw = target_props.get("mol_wt") or (p_props.get("MW") if p_props else "—")

    species = impurity_records_to_appendix_species(
        records,
        parent_api_smiles=target_smi,
    )
    primary: dict[str, Any] = {}
    if species:
        primary = min(
            species,
            key=lambda row: int(row.get("ich_m7_class") or 5),
        )

    ctx = build_process_product_context(
        product_name="Target API (Drug Substance)",
        product_smiles=target_smi,
        mw=p_mw,
        reactant1=reactant1 or "",
        reactant2=reactant2 or "",
        reagent=reagent_smiles or "",
        nitrite=bool(process_conditions and process_conditions.get("nitrite")),
        reference_structure=reference_api or "",
        primary_impurity_name=str(primary.get("name") or ""),
        primary_impurity_smiles=str(primary.get("smiles") or ""),
        mechanism=_process_mechanism_summary(process_conditions),
        stability=p_stab,
        m7=p_m7,
        alerts=p_alerts,
        db_status="Theoretical — computational prediction only",
        db_ref="",
        pubchem_cid=None,
        properties=p_props,
        advanced=p_adv,
        indexed_count=indexed_count,
    )
    return replace(ctx, appendix_species=species)


def build_integrated_report_dataframe(
    records: list[dict[str, Any]],
    process_conditions: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Phase 5~7 레코드 → Phase 8 통합 보고서 6+1열 DataFrame."""
    if not records:
        return pd.DataFrame(columns=INTEGRATED_EXPORT_COLUMNS)

    rows: list[dict[str, str]] = []
    for rec in records:
        ich_cls = rec.get("ich_m7_class")
        try:
            ich_int = int(ich_cls) if ich_cls is not None else None
        except (TypeError, ValueError):
            ich_int = None
        rows.append(
            {
                RPT_COL_IMAGE: "",
                RPT_COL_SMILES: str(rec.get("smiles") or ""),
                RPT_COL_NAME_RXN: _name_reaction_from_record(rec),
                RPT_COL_TRIGGER: _trigger_context(rec, process_conditions),
                RPT_COL_STRUCTURAL_ALERT: _format_structural_alert_cell(rec),
                RPT_COL_AI_SCORE: _format_ai_score_cell(rec),
                RPT_COL_ICH_M7_CLASS: _format_ich_m7_class_cell(rec),
                RPT_COL_REMARKS: _regulatory_remarks(
                    ich_int,
                    str(rec.get("regulatory_action_text") or rec.get(IMP_COL_REGULATORY) or ""),
                ),
            }
        )
    return pd.DataFrame(rows)


def _resolve_export_dataframe(
    pred_df: pd.DataFrame | None,
    records: list[dict[str, Any]] | None,
    process_conditions: dict[str, Any] | None,
) -> pd.DataFrame:
    if records:
        return build_integrated_report_dataframe(records, process_conditions)
    if pred_df is not None and not pred_df.empty:
        raw_records: list[dict[str, Any]] = []
        norm = ProcessImpurityEngine.normalize_prediction_dataframe(pred_df)
        for _, row in norm.iterrows():
            raw_records.append(row.to_dict())
        return build_integrated_report_dataframe(raw_records, process_conditions)
    return pd.DataFrame(columns=INTEGRATED_EXPORT_COLUMNS)


def compose_pdf_base_conditions(
    target_props: dict[str, Any],
    process_conditions: dict[str, Any] | None,
) -> dict[str, Any]:
    """PDF executive overview용 공정·목적물 메타데이터."""
    base: dict[str, Any] = {
        "target_smiles": str(target_props.get("smiles") or ""),
        "mol_wt": target_props.get("mol_wt"),
        "exact_mol_wt": target_props.get("exact_mol_wt"),
        "solvent": "—",
        "reagent": "—",
        "temperature_c": "—",
        "ph": "—",
        "reaction_time_h": "—",
        "stress_condition": "—",
        "report_title": INTEGRATED_REPORT_TITLE,
    }
    if process_conditions:
        ctx = process_conditions.get("process_reagents") or {}
        base.update(
            {
                "solvent": process_conditions.get("solvent", "—"),
                "reagent": process_conditions.get("reagent", "—"),
                "solvents": ", ".join(ctx.get("solvents") or []) or "—",
                "acids": ", ".join(ctx.get("acids") or []) or "—",
                "bases": ", ".join(ctx.get("bases") or []) or "—",
                "catalysts": ", ".join(ctx.get("catalysts") or []) or "—",
                "temperature_c": process_conditions.get("temperature_c", "—"),
                "ph": process_conditions.get("ph", "—"),
                "reaction_time_h": process_conditions.get("reaction_time_h", "—"),
                "stress_condition": process_conditions.get("stress_condition", "—"),
                "name_reaction": process_conditions.get("name_reaction", "—"),
            }
        )
    return base


def _excel_row_smiles(row: pd.Series) -> str:
    raw = row.get(RPT_COL_SMILES, "")
    if raw is not None and str(raw).strip() and str(raw).lower() != "nan":
        return str(raw).strip()
    return ""


def _smiles_to_png_bytes(
    smiles: str,
    *,
    size: tuple[int, int] = EXCEL_MOL_IMAGE_SIZE,
) -> bytes | None:
    """RDKit Cairo 드로어 우선 · MolToImage 폴백."""
    cleaned = (smiles or "").strip()
    if not cleaned:
        return None
    png = rdkit_png(cleaned, size=f"{size[0]}x{size[1]}")
    if png:
        return png
    mol = Chem.MolFromSmiles(cleaned)
    if mol is None:
        return None
    try:
        pil_image = Draw.MolToImage(mol, size=size)
        image_bytes_io = io.BytesIO()
        pil_image.save(image_bytes_io, format="PNG")
        return image_bytes_io.getvalue()
    except Exception:
        return None


def _mol_to_image_bytes_io(smiles: str) -> io.BytesIO | None:
    raw = _smiles_to_png_bytes(smiles, size=EXCEL_MOL_IMAGE_SIZE)
    if not raw:
        return None
    image_bytes_io = io.BytesIO(raw)
    image_bytes_io.seek(0)
    return image_bytes_io


def generate_primeds_process_excel(
    pred_df: pd.DataFrame | None = None,
    *,
    records: list[dict[str, Any]] | None = None,
    process_conditions: dict[str, Any] | None = None,
) -> bytes:
    """Phase 8 통합 평가 보고서 Excel 바이너리."""
    export_df = _resolve_export_dataframe(pred_df, records, process_conditions)
    buffer = io.BytesIO()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
        workbook = writer.book
        worksheet = workbook.add_worksheet(EXCEL_SHEET_NAME)

        title_fmt = workbook.add_format(
            {
                "bold": True,
                "font_size": 14,
                "font_color": "#1a365d",
                "align": "left",
                "valign": "vcenter",
            }
        )
        meta_fmt = workbook.add_format({"font_size": 9, "font_color": "#4a5568"})
        header_fmt = workbook.add_format(
            {
                "bold": True,
                "bg_color": "#E3F2FD",
                "border": 1,
                "text_wrap": True,
                "valign": "vcenter",
                "align": "center",
            }
        )
        cell_fmt = workbook.add_format(
            {"border": 1, "text_wrap": True, "valign": "top"}
        )
        remarks_fmt = workbook.add_format(
            {"border": 1, "text_wrap": True, "valign": "top", "font_color": "#1565c0"}
        )

        n_cols = len(INTEGRATED_EXPORT_COLUMNS)
        worksheet.merge_range(0, 0, 0, n_cols - 1, INTEGRATED_REPORT_TITLE, title_fmt)
        worksheet.merge_range(
            1,
            0,
            1,
            n_cols - 1,
            f"Generated: {timestamp} · PRIME-DS Process Mode · ICH M7 Class 1~5 Hybrid",
            meta_fmt,
        )

        for col_i, width in enumerate(EXCEL_COL_WIDTHS):
            fmt = remarks_fmt if col_i == n_cols - 1 else cell_fmt
            worksheet.set_column(col_i, col_i, width, fmt)
        worksheet.set_column(EXCEL_IMAGE_COL_IDX, EXCEL_IMAGE_COL_IDX, EXCEL_IMAGE_COL_WIDTH)

        for col_i, title in enumerate(INTEGRATED_EXPORT_COLUMNS):
            worksheet.write(EXCEL_HEADER_ROW, col_i, title, header_fmt)
        worksheet.set_row(EXCEL_HEADER_ROW, 36, header_fmt)
        worksheet.freeze_panes(EXCEL_DATA_START_ROW, 0)

        if export_df.empty:
            buffer.seek(0)
            return buffer.getvalue()

        for row_i, (_, row) in enumerate(export_df.iterrows()):
            excel_row = EXCEL_DATA_START_ROW + row_i
            worksheet.set_row(excel_row, EXCEL_ROW_HEIGHT)
            for col_i, col_name in enumerate(INTEGRATED_EXPORT_COLUMNS):
                val = row.get(col_name, "")
                fmt = remarks_fmt if col_name == RPT_COL_REMARKS else cell_fmt
                worksheet.write(excel_row, col_i, str(val), fmt)
            smiles = _excel_row_smiles(row)
            if smiles:
                image_bytes_io = _mol_to_image_bytes_io(smiles)
                if image_bytes_io is not None:
                    worksheet.insert_image(
                        excel_row,
                        EXCEL_IMAGE_COL_IDX,
                        "molecule.png",
                        {
                            "image_data": image_bytes_io,
                            "x_offset": EXCEL_IMAGE_X_OFFSET,
                            "y_offset": EXCEL_IMAGE_Y_OFFSET,
                        },
                    )

    buffer.seek(0)
    return buffer.getvalue()


def _impurity_mol_to_base64_data_uri(
    smiles: str,
    *,
    size: tuple[int, int] = (150, 150),
) -> str:
    png = _smiles_to_png_bytes(smiles, size=size)
    if not png:
        return ""
    b64_str = base64.b64encode(png).decode("ascii")
    return f"data:image/png;base64,{b64_str}"


def _structure_png_for_report(smiles: str) -> bytes | None:
    """단일화합물 보고서와 동일한 RDKit 구조식 렌더링."""
    cleaned = (smiles or "").strip()
    if not cleaned:
        return None
    png = rdkit_clean_structure_png(cleaned)
    if png:
        return png
    return _smiles_to_png_bytes(cleaned, size=(480, 360))


def _impurity_export_detail_rows(row: pd.Series) -> list[tuple[str, str]]:
    """유연물질 1건 — 단일화합물 보고서 스타일 KV 테이블 행."""
    return [
        (RPT_COL_SMILES, str(row.get(RPT_COL_SMILES, "") or "")),
        (RPT_COL_NAME_RXN, str(row.get(RPT_COL_NAME_RXN, "") or "—")),
        (RPT_COL_TRIGGER, str(row.get(RPT_COL_TRIGGER, "") or "—")),
        (RPT_COL_STRUCTURAL_ALERT, str(row.get(RPT_COL_STRUCTURAL_ALERT, "") or "—")),
        (RPT_COL_AI_SCORE, str(row.get(RPT_COL_AI_SCORE, "") or "—")),
        (RPT_COL_ICH_M7_CLASS, str(row.get(RPT_COL_ICH_M7_CLASS, "") or "—")),
        (RPT_COL_REMARKS, str(row.get(RPT_COL_REMARKS, "") or "—")),
    ]


def _process_conditions_kv_rows(
    process_conditions: dict[str, Any] | None,
) -> list[tuple[str, str]]:
    if not process_conditions:
        return []
    ctx = process_conditions.get("process_reagents") or {}
    return [
        ("Name Reaction", str(process_conditions.get("name_reaction", "—"))),
        ("Stress Pathway", str(process_conditions.get("stress_condition", "—"))),
        ("Solvent(s)", ", ".join(ctx.get("solvents") or []) or "—"),
        ("Acid(s)", ", ".join(ctx.get("acids") or []) or "—"),
        ("Base(s)", ", ".join(ctx.get("bases") or []) or "—"),
        ("Catalyst(s)", ", ".join(ctx.get("catalysts") or []) or "—"),
        ("Temperature (°C)", str(process_conditions.get("temperature_c", "—"))),
        ("pH", str(process_conditions.get("ph", "—"))),
        ("Runtime (h)", str(process_conditions.get("reaction_time_h", "—"))),
    ]


def _build_regulatory_grade_html(
    export_df: pd.DataFrame,
    base_conditions_dict: dict[str, Any],
) -> str:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report_id = datetime.now().strftime("PRIME-DS-%Y%m%d-%H%M%S")
    doc_title = str(
        base_conditions_dict.get("report_title") or INTEGRATED_REPORT_TITLE
    )

    target_smiles = html.escape(str(base_conditions_dict.get("target_smiles") or "N/A"))
    avg_mw = base_conditions_dict.get("mol_wt")
    exact_mw = base_conditions_dict.get("exact_mol_wt")
    if avg_mw is not None and exact_mw is not None:
        mw_line = (
            f"Average MW: {float(avg_mw):.{MW_DISPLAY_DECIMALS}f} g/mol · "
            f"Exact mass: {float(exact_mw):.{MW_DISPLAY_DECIMALS}f} g/mol"
        )
    elif exact_mw is not None:
        mw_line = f"Exact mass: {float(exact_mw):.{MW_DISPLAY_DECIMALS}f} g/mol"
    else:
        mw_line = "—"

    target_img_uri = _impurity_mol_to_base64_data_uri(
        str(base_conditions_dict.get("target_smiles") or ""),
        size=(480, 360),
    )
    target_img_tag = (
        '<div class="structure-block">'
        f'<img class="mol-target" src="{target_img_uri}" alt="Target API 2D structure"/>'
        '<p class="structure-caption">Figure 1: Target API structure</p>'
        "</div>"
        if target_img_uri
        else "<em class='no-img'>Structure image unavailable</em>"
    )

    stress_display = html.escape(
        ProcessImpurityEngine.normalize_stress_condition(
            str(base_conditions_dict.get("stress_condition", "—"))
        )
    )
    name_rxn = html.escape(str(base_conditions_dict.get("name_reaction", "—")))
    solvents = html.escape(str(base_conditions_dict.get("solvents", "—")))
    acids = html.escape(str(base_conditions_dict.get("acids", "—")))
    bases = html.escape(str(base_conditions_dict.get("bases", "—")))
    catalysts = html.escape(str(base_conditions_dict.get("catalysts", "—")))
    temperature = html.escape(str(base_conditions_dict.get("temperature_c", "—")))
    ph_value = html.escape(str(base_conditions_dict.get("ph", "—")))
    runtime_h = html.escape(str(base_conditions_dict.get("reaction_time_h", "—")))

    impurity_cards_html = ""
    if export_df is not None and not export_df.empty:
        for idx, (_, row) in enumerate(export_df.iterrows(), start=1):
            smi = str(row.get(RPT_COL_SMILES, "") or "")
            base64_str = _impurity_mol_to_base64_data_uri(smi, size=(400, 300))
            img_block = (
                '<div class="structure-block">'
                f'<img class="mol-impurity" src="{base64_str}" alt="Impurity 2D"/>'
                '<p class="structure-caption">2D structure depiction</p>'
                "</div>"
                if base64_str
                else '<p class="no-img">Structure image unavailable</p>'
            )
            detail_rows = "".join(
                f"<tr><th>{html.escape(label)}</th>"
                f"<td>{html.escape(str(row.get(col, '') or '—'))}</td></tr>"
                for label, col in (
                    (RPT_COL_SMILES, RPT_COL_SMILES),
                    (RPT_COL_NAME_RXN, RPT_COL_NAME_RXN),
                    (RPT_COL_TRIGGER, RPT_COL_TRIGGER),
                    (RPT_COL_STRUCTURAL_ALERT, RPT_COL_STRUCTURAL_ALERT),
                    (RPT_COL_AI_SCORE, RPT_COL_AI_SCORE),
                    (RPT_COL_ICH_M7_CLASS, RPT_COL_ICH_M7_CLASS),
                    (RPT_COL_REMARKS, RPT_COL_REMARKS),
                )
            )
            impurity_cards_html += (
                f'<div class="impurity-card">'
                f'<h3 class="impurity-heading">I.{idx:03d} — Process impurity candidate</h3>'
                f"{img_block}"
                f'<table class="overview-grid impurity-detail">{detail_rows}</table>'
                f"</div>"
            )
        impurity_section = (
            f'<p class="summary-count">{len(export_df)} process impurity candidate(s) · '
            f"ICH M7 Class 1~5 hybrid assessment.</p>"
            f"{impurity_cards_html}"
        )
    else:
        impurity_section = (
            "<p class='empty-notice'>No predicted impurity data available. "
            "Run PRIME-DS Process Mode prediction first.</p>"
        )

    engine_note = (
        f"{UI_ENGINE_TITLE} · Structural Alert + DeepChem Ames QSAR · "
        "ICH M7 Class 1~5 hybrid classifier"
    )

    malgun_path = Path(r"C:\Windows\Fonts\malgun.ttf")
    font_face_css = ""
    if malgun_path.is_file():
        font_uri = malgun_path.as_uri()
        font_face_css = f"""
  @font-face {{
    font-family: "PrimeKorean";
    src: url("{font_uri}");
  }}"""

    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="utf-8"/>
<title>{html.escape(doc_title)}</title>
<style>
  {font_face_css}
  @page {{ size: A4; margin: 16mm; }}
  body {{ font-family: "PrimeKorean", "Malgun Gothic", "Segoe UI", Arial, sans-serif; font-size: 9pt; color: #1a202c; }}
  .doc-header {{ border-bottom: 3px solid #1a365d; padding-bottom: 10px; margin-bottom: 14px; }}
  .doc-header h1 {{ font-size: 14pt; color: #1a365d; margin: 0 0 4px 0; }}
  .doc-meta {{ font-size: 8pt; color: #4a5568; }}
  .executive-card {{ background: #f7fafc; border: 1.5px solid #1a365d; padding: 12px; margin-bottom: 14px; }}
  .executive-card h2 {{ font-size: 11pt; color: #1a365d; margin: 0 0 8px 0; }}
  table.overview-grid {{ width: 100%; border-collapse: collapse; font-size: 8.5pt; }}
  table.overview-grid th, table.overview-grid td {{ border: 0.5px solid #cbd5e0; padding: 5px 7px; }}
  table.overview-grid th {{ background: #ebf4ff; width: 26%; text-align: left; }}
  h2.section-title {{ font-size: 11pt; color: #1a365d; border-bottom: 1.5px solid #1a365d; margin: 16px 0 8px; }}
  .structure-block {{ text-align: center; margin: 10px 0 14px; page-break-inside: avoid; }}
  .structure-caption {{ font-size: 8pt; color: #4a5568; font-style: italic; margin: 6px 0 0; }}
  .impurity-card {{ border: 1px solid #cbd5e0; border-radius: 4px; padding: 10px 12px; margin-bottom: 14px; page-break-inside: avoid; }}
  .impurity-heading {{ font-size: 10pt; color: #1a365d; margin: 0 0 8px 0; }}
  table.impurity-detail {{ margin-top: 6px; }}
  td.smiles {{ font-family: Consolas, monospace; font-size: 7pt; word-break: break-all; }}
  td.remarks {{ color: #1565c0; font-size: 7pt; }}
  img.mol-target {{ width: 200px; height: auto; max-height: 150px; display: inline-block; }}
  img.mol-impurity {{ width: 180px; height: auto; max-height: 135px; display: inline-block; }}
  .no-img {{ color: #a0aec0; font-size: 7pt; text-align: center; }}
  .doc-footer {{ margin-top: 16px; font-size: 7.5pt; color: #718096; border-top: 1px solid #cbd5e0; padding-top: 8px; }}
</style>
</head>
<body>
  <div class="doc-header">
    <h1>{html.escape(doc_title)}</h1>
    <p class="doc-meta"><strong>ID:</strong> {html.escape(report_id)} · <strong>Generated:</strong> {html.escape(timestamp)}</p>
  </div>
  <div class="executive-card">
    <h2>공정 조건 Executive Summary</h2>
    <table class="overview-grid">
      <tr><th>Name Reaction</th><td>{name_rxn}</td></tr>
      <tr><th>Stress Pathway</th><td>{stress_display}</td></tr>
      <tr><th>Solvent(s)</th><td>{solvents}</td></tr>
      <tr><th>Acid(s)</th><td>{acids}</td></tr>
      <tr><th>Base(s)</th><td>{bases}</td></tr>
      <tr><th>Catalyst(s)</th><td>{catalysts}</td></tr>
      <tr><th>Temperature (°C)</th><td>{temperature}</td></tr>
      <tr><th>pH</th><td>{ph_value}</td></tr>
      <tr><th>Runtime (h)</th><td>{runtime_h}</td></tr>
    </table>
  </div>
  <h2 class="section-title">Target API (Drug Substance)</h2>
  <table class="overview-grid">
    <tr><th>SMILES</th><td class="smiles">{target_smiles}</td></tr>
    <tr><th>Molecular weight</th><td>{html.escape(mw_line)}</td></tr>
  </table>
  {target_img_tag}
  <h2 class="section-title">공정 유연물질 · ICH M7 통합 평가</h2>
  {impurity_section}
  <div class="doc-footer"><p>{html.escape(engine_note)}</p></div>
</body>
</html>"""


def generate_primeds_process_pdf(
    pred_df: pd.DataFrame | None = None,
    base_conditions_dict: dict[str, Any] | None = None,
    *,
    records: list[dict[str, Any]] | None = None,
    process_conditions: dict[str, Any] | None = None,
) -> bytes:
    """WeasyPrint 규제 등급 PDF (Phase 8 통합 양식)."""
    if weasyprint is None:
        raise RuntimeError(
            "WeasyPrint is not available. Use generate_primeds_process_comparative_pdf() "
            "for FPDF fallback."
        )
    export_df = _resolve_export_dataframe(pred_df, records, process_conditions)
    base = base_conditions_dict or {}
    html_content = _build_regulatory_grade_html(export_df, base)
    return weasyprint.HTML(string=html_content).write_pdf()


def _pdf_write_process_title_block(pdf: FPDF, family: str, w: float) -> None:
    """단일화합물 보고서와 동일한 표지 타이틀."""
    _pdf_ensure_space(pdf, 28)
    pdf.set_text_color(26, 54, 93)
    _pdf_write_centered(pdf, family, REPORT_TITLE_LINE1, w, bold=True, size=18)
    pdf.ln(1)
    _pdf_write_centered(pdf, family, REPORT_TITLE_LINE2, w, bold=True, size=16)
    pdf.ln(2)
    pdf.set_text_color(0, 0, 0)
    _pdf_write_centered(
        pdf,
        family,
        "ICH M7(R2) — Process Impurity Integrated Assessment",
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


def _target_mw_line(target_props: dict[str, Any]) -> str:
    avg_mw = target_props.get("mol_wt")
    exact_mw = target_props.get("exact_mol_wt")
    if avg_mw is not None and exact_mw is not None:
        return (
            f"{float(avg_mw):.{MW_DISPLAY_DECIMALS}f} g/mol (average) · "
            f"{float(exact_mw):.{MW_DISPLAY_DECIMALS}f} g/mol (exact)"
        )
    if exact_mw is not None:
        return f"{float(exact_mw):.{MW_DISPLAY_DECIMALS}f} g/mol (exact)"
    if avg_mw is not None:
        return f"{float(avg_mw):.{MW_DISPLAY_DECIMALS}f} g/mol (average)"
    return "N/A"


def _generate_primeds_process_pdf_fpdf(
    target_props: dict[str, Any],
    export_df: pd.DataFrame,
    process_conditions: dict[str, Any] | None = None,
) -> bytes:
    """FPDF PDF — 단일화합물 보고서와 동일한 중앙 구조식 + KV 테이블 레이아웃."""
    pdf = FPDF()
    pdf.set_margins(18, 18, 18)
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.add_page()
    w = pdf.epw
    family, family_has_bold = _register_pdf_fonts(pdf)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    _pdf_write_process_title_block(pdf, family, w)
    pdf.ln(2)
    _pdf_write_centered(pdf, family, INTEGRATED_REPORT_TITLE, w, bold=True, size=11)
    _pdf_write_centered(pdf, family, f"Generated: {timestamp}", w, size=8)
    pdf.ln(3)

    cond_rows = _process_conditions_kv_rows(process_conditions)
    if cond_rows:
        _pdf_section(pdf, family, "PROCESS CONDITIONS — EXECUTIVE SUMMARY", w)
        _pdf_kv_rows(pdf, family, cond_rows, w, family_has_bold=family_has_bold)
        _pdf_blank_lines(pdf, 1)

    target_smi = str(target_props.get("smiles") or "N/A")
    _pdf_section(pdf, family, "TARGET API (DRUG SUBSTANCE)", w)
    _pdf_kv_rows(
        pdf,
        family,
        [
            ("SMILES", target_smi),
            ("Molecular weight", _target_mw_line(target_props)),
        ],
        w,
        family_has_bold=family_has_bold,
    )
    target_png = _structure_png_for_report(target_smi)
    if target_png:
        _pdf_add_structure_figure(
            pdf,
            family,
            target_png,
            w,
            caption="Figure 1: Target API structure",
        )

    _pdf_section(pdf, family, "PROCESS IMPURITIES — ICH M7 INTEGRATED ASSESSMENT", w)
    if export_df.empty:
        _pdf_write_wrapped(pdf, family, "No predicted impurity data available.", w, size=9)
    else:
        _pdf_write_wrapped(
            pdf,
            family,
            (
                f"{len(export_df)} process impurity candidate(s) · "
                "ICH M7 Class 1~5 hybrid assessment."
            ),
            w,
            size=8,
        )
        _pdf_blank_lines(pdf, 1)
        for index, (_, row) in enumerate(export_df.iterrows(), start=1):
            smi = str(row.get(RPT_COL_SMILES, "") or "")
            _pdf_ensure_space(pdf, 90)
            _pdf_write_wrapped(
                pdf,
                family,
                f"I.{index:03d} — Process impurity candidate",
                w,
                bold=True,
                size=10,
            )
            imp_png = _structure_png_for_report(smi)
            if imp_png:
                _pdf_add_structure_figure(
                    pdf,
                    family,
                    imp_png,
                    w,
                    caption="2D structure depiction",
                )
            _pdf_kv_rows(
                pdf,
                family,
                _impurity_export_detail_rows(row),
                w,
                family_has_bold=family_has_bold,
            )
            _pdf_blank_lines(pdf, 2)

    _pdf_blank_lines(pdf, 1)
    _pdf_write_wrapped(
        pdf,
        family,
        f"{UI_ENGINE_TITLE} · Structural Alert + DeepChem Ames QSAR · "
        "ICH M7 Class 1~5 hybrid classifier",
        w,
        size=7.5,
    )
    return _pdf_output_bytes(pdf)


def generate_primeds_process_comparative_pdf(
    target_props: dict[str, Any],
    pred_df: pd.DataFrame | None = None,
    process_conditions: dict[str, Any] | None = None,
    *,
    records: list[dict[str, Any]] | None = None,
    reactant1: str = "",
    reactant2: str = "",
    reagent_smiles: str = "",
    reference_api: str = "",
) -> bytes:
    """단일화합물 검색 보고서와 동일한 규제형 PDF (구조식 중앙 + KV 테이블 + Appendix A)."""
    resolved_records = records or []
    if not resolved_records and pred_df is not None and not pred_df.empty:
        norm = ProcessImpurityEngine.normalize_prediction_dataframe(pred_df)
        resolved_records = [row.to_dict() for _, row in norm.iterrows()]

    try:
        ctx = build_process_mode_regulatory_context(
            target_props,
            resolved_records,
            process_conditions,
            reactant1=reactant1,
            reactant2=reactant2,
            reagent_smiles=reagent_smiles,
            reference_api=reference_api,
        )
        return report_to_regulatory_pdf(ctx)
    except Exception as reg_exc:
        log_process_warning(
            f"Regulatory process PDF failed → legacy FPDF layout: {reg_exc}"
        )

    export_df = _resolve_export_dataframe(pred_df, records, process_conditions)
    try:
        return _generate_primeds_process_pdf_fpdf(
            target_props, export_df, process_conditions
        )
    except Exception as fpdf_exc:
        log_process_warning(f"FPDF process report failed → WeasyPrint: {fpdf_exc}")

    if weasyprint is None:
        raise RuntimeError(f"PDF generation failed: {fpdf_exc}") from fpdf_exc

    base_conditions = compose_pdf_base_conditions(target_props, process_conditions)
    html_content = _build_regulatory_grade_html(export_df, base_conditions)
    return weasyprint.HTML(string=html_content).write_pdf()
