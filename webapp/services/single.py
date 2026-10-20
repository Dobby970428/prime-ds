"""모드 1 — 단일 화합물 독성 · 물성 스캔.

Streamlit 판의 _run_single_toxicity_scan + _render_single_scan_core_results에
해당합니다. 계산만 하고 결과 dict를 돌려줍니다.
"""

from __future__ import annotations

from typing import Any

from . import chem


def run_scan(
    smiles: str,
    *,
    reference_smiles: str = "",
    allow_slow_db_scan: bool = False,
    with_read_across: bool = True,
) -> dict[str, Any]:
    """단일 화합물 스캔 전체를 수행합니다.

    반환 dict:
        identity   구조 신원 (분자식·분자량·canonical)
        alerts     Layer B 구조 경보
        m7         ICH M7 판정 (dict)
        db         Validated / Theoretical 판정
        properties 물성 (baseline · advanced · stability)
        hybrid     3중 대조 결과 (자체 DB / 오픈소스 표)
        analogues  Read-Across 유사 구조
        errors     계산 중 발생한 오류 목록
    """
    result: dict[str, Any] = {
        "smiles": (smiles or "").strip(),
        "reference_smiles": (reference_smiles or "").strip(),
        "identity": None,
        "alerts": {"labels": [], "atoms": [], "error": None},
        "m7": None,
        "db": {"status": "", "record": None},
        "properties": {},
        "hybrid": None,
        "analogues": {"analogues": [], "threshold": None, "error": None},
        "pubchem": {},
        "errors": [],
    }

    identity = chem.describe(smiles)
    result["identity"] = identity
    if not identity["valid"]:
        result["errors"].append(identity["error"] or "구조를 해석하지 못했습니다")
        return result

    canonical = identity["canonical"]

    # --- Layer B 구조 경보 -------------------------------------------------
    result["alerts"] = chem.alerts(canonical)

    # --- ICH M7 판정 (Layer A + B + C 합산) --------------------------------
    m7_raw, m7_error = chem.classify(canonical, reference_smiles, result["alerts"]["labels"])
    result["m7"] = chem.m7_view(m7_raw)
    result["_m7_raw"] = m7_raw
    result["m7_error"] = m7_error
    if m7_raw is None:
        result["errors"].append(f"ICH M7 classification failed — {m7_error or 'unknown error'}")

    # --- 자체 DB 대조 (Validated / Theoretical) ----------------------------
    try:
        from p24_db_smiles_lookup import classify_smiles_against_db

        status, record = classify_smiles_against_db(
            canonical, None, allow_slow_scan=allow_slow_db_scan
        )
        result["db"] = {"status": status, "record": record}
    except Exception as exc:
        result["errors"].append(f"DB 대조: {exc}")

    # --- 물성 -------------------------------------------------------------
    props = chem.properties(canonical)
    result["properties"] = props
    result["errors"].extend(props.get("errors", []))

    # --- PubChem ----------------------------------------------------------
    try:
        from p25_pubchem_client import fetch_compound_by_smiles

        result["pubchem"] = fetch_compound_by_smiles(canonical) or {}
    except Exception:
        result["pubchem"] = {}

    # --- 3중 대조 (자체 DB · PubChem · 구조 경보) --------------------------
    try:
        from p68_hybrid_toxicity import run_hybrid_toxicity_assessment

        hybrid = run_hybrid_toxicity_assessment(
            canonical, local_df=None, drug_smiles=reference_smiles or None
        )
        result["hybrid"] = _hybrid_view(hybrid)
        result["_hybrid_raw"] = hybrid
    except Exception as exc:
        result["errors"].append(f"3중 대조: {exc}")

    # --- Read-Across ------------------------------------------------------
    if with_read_across:
        result["analogues"] = chem.read_across(canonical)

    return result


def _hybrid_view(assessment: Any) -> dict[str, Any] | None:
    if assessment is None:
        return None
    get = lambda n, d=None: getattr(assessment, n, d)  # noqa: E731
    classification = get("classification")
    return {
        "canonical": get("canon_smiles", ""),
        "classification": str(classification) if classification is not None else "",
        "verdict_label": get("verdict_label", ""),
        "is_validated": "valid" in str(get("verdict_label", "")).lower()
        or "실측" in str(get("verdict_label", "")),
        "db_record": get("db_record"),
        "alerts": list(get("alerts", []) or []),
        "highlight_atoms": list(get("highlight_atoms", []) or []),
        "alert_status": get("alert_status", ""),
        "alert_emoji": get("alert_emoji", ""),
        "pubchem": dict(get("pubchem", {}) or {}),
        "inhouse_table": list(get("inhouse_table", []) or []),
        "opensource_table": list(get("opensource_table", []) or []),
        "errors": list(get("errors", []) or []),
    }


def build_report_context(scan: dict[str, Any], *, name: str = "") -> Any:
    """규제 보고서 컨텍스트를 만듭니다 (p80)."""
    identity = scan.get("identity") or {}
    props = scan.get("properties") or {}
    m7_raw = scan.get("_m7_raw")
    if m7_raw is None:
        return None

    db = scan.get("db") or {}
    record = db.get("record") or {}
    pubchem = scan.get("pubchem") or {}

    try:
        from p80_formal_assessment_report import build_single_compound_context

        return build_single_compound_context(
            name=name or record.get("name") or pubchem.get("name") or "Unnamed compound",
            smiles=identity.get("canonical") or scan.get("smiles", ""),
            mw=identity.get("mw") or "",
            properties=props.get("baseline"),
            advanced=props.get("advanced"),
            stability=props.get("stability") or ("", ""),
            m7=m7_raw,
            alerts=list((scan.get("alerts") or {}).get("labels") or []),
            db_status=db.get("status", ""),
            db_ref=str(record.get("drug_file") or record.get("source") or ""),
            pubchem_cid=pubchem.get("cid"),
            reference_structure=scan.get("reference_smiles", ""),
        )
    except Exception:
        return None


def render_report_text(scan: dict[str, Any], *, name: str = "") -> str:
    ctx = build_report_context(scan, name=name)
    if ctx is None:
        return ""
    try:
        from p80_formal_assessment_report import generate_formal_assessment_report

        return generate_formal_assessment_report(ctx)
    except Exception as exc:
        return f"보고서 생성 실패: {exc}"
