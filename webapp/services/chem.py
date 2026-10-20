"""공용 화학 헬퍼 — SMILES 검증, 물성, 구조 정보."""

from __future__ import annotations

from typing import Any


def parse(smiles: str | None):
    """SMILES → RDKit Mol. 실패하면 None."""
    text = (smiles or "").strip()
    if not text:
        return None
    try:
        from rdkit import Chem

        return Chem.MolFromSmiles(text)
    except Exception:
        return None


def canonical(smiles: str | None) -> str:
    mol = parse(smiles)
    if mol is None:
        return ""
    try:
        from rdkit import Chem

        return Chem.MolToSmiles(mol)
    except Exception:
        return ""


def describe(smiles: str | None) -> dict[str, Any]:
    """구조의 기본 신원 정보. 판정에 쓰이지 않고 화면 표시용입니다."""
    out: dict[str, Any] = {
        "smiles": (smiles or "").strip(),
        "valid": False,
        "canonical": "",
        "formula": "",
        "mw": None,
        "heavy_atoms": None,
        "error": None,
    }
    text = out["smiles"]
    if not text:
        out["error"] = "구조가 비어 있습니다"
        return out

    mol = parse(text)
    if mol is None:
        out["error"] = "RDKit이 해석하지 못한 SMILES입니다"
        return out

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        out["valid"] = True
        out["canonical"] = Chem.MolToSmiles(mol)
        out["formula"] = rdMolDescriptors.CalcMolFormula(mol)
        out["mw"] = round(float(Descriptors.ExactMolWt(mol)), 4)
        out["heavy_atoms"] = mol.GetNumHeavyAtoms()
    except Exception as exc:
        out["error"] = str(exc)
    return out


def properties(smiles: str) -> dict[str, Any]:
    """물성 계산.

    baseline(p50 RDKit 표준)은 신뢰할 수 있습니다.
    advanced(p51 HybridPropertyEngine)의 pKa·LogD·끓는점은 고정 상수에서
    고른 추정값이라 표시용입니다. 화면에서 그렇게 안내합니다.
    """
    out: dict[str, Any] = {"baseline": None, "advanced": None, "stability": None, "errors": []}

    try:
        from p50_physicochemical import calculate_opensource_properties

        out["baseline"] = calculate_opensource_properties(smiles)
    except Exception as exc:
        out["errors"].append(f"기본 물성: {exc}")

    try:
        from p50_physicochemical import calculate_my_advanced_properties

        out["advanced"] = calculate_my_advanced_properties(smiles)
    except Exception as exc:
        out["errors"].append(f"고도화 물성: {exc}")

    try:
        from p50_physicochemical import predict_chemical_stability

        out["stability"] = predict_chemical_stability(smiles, out["baseline"])
    except Exception as exc:
        out["errors"].append(f"안정성: {exc}")

    return out


def alerts(smiles: str) -> dict[str, Any]:
    """Ashby-Tennant 구조 경보 (Layer B)."""
    mol = parse(smiles)
    if mol is None:
        return {"labels": [], "atoms": [], "error": "구조 해석 실패"}
    try:
        from p62_toxicity_alerts import scan_ashby_tennant_alerts

        labels, atoms = scan_ashby_tennant_alerts(mol)
        return {"labels": list(labels or []), "atoms": list(atoms or []), "error": None}
    except Exception as exc:
        return {"labels": [], "atoms": [], "error": str(exc)}


def classify(smiles: str, reference_smiles: str = "", alert_labels=None) -> tuple[Any, str | None]:
    """ICH M7 판정 (p66 — 상향 전용 합산).

    실패하면 (None, 사유)를 돌려줍니다. 예전에는 예외를 조용히 삼켜서
    화면에 «판정 실패»만 뜨고 원인을 알 수 없었습니다.
    """
    mol = parse(smiles)
    if mol is None:
        return None, "RDKit could not parse the structure"

    ref_mol = parse(reference_smiles) if reference_smiles else None
    try:
        from p66_ich_m7 import classify_ich_m7

        return classify_ich_m7(mol, alert_labels, ref_mol), None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def read_across(smiles: str, top_n: int = 3) -> dict[str, Any]:
    """유사 구조 Read-Across (p51 ReadAcrossEngine)."""
    try:
        from p51_prime_backend_core import READ_ACROSS_MIN_SIMILARITY, ReadAcrossEngine

        found = ReadAcrossEngine.find_analogues(
            smiles, top_n=top_n, min_similarity=READ_ACROSS_MIN_SIMILARITY
        )
        return {"analogues": list(found or []), "threshold": READ_ACROSS_MIN_SIMILARITY, "error": None}
    except Exception as exc:
        return {"analogues": [], "threshold": None, "error": str(exc)}


def oecd_assessment(smiles: str) -> dict[str, Any]:
    """OECD QSAR (AMBIT2) 조회.

    ⚠️ AMBIT 서버가 떠 있든 아니든 실제 계산은 RDKit입니다.
    run_ambit_bridge.bat이 띄우는 것도 진짜 OECD Toolbox가 아니라
    p55가 만든 흉내내기 서버입니다. 화면에서 그렇게 안내합니다.
    """
    out: dict[str, Any] = {"rows": [], "status": "", "endpoint": "", "is_proxy": True, "error": None}
    text = (smiles or "").strip()
    if not text:
        out["error"] = "구조가 비어 있습니다"
        return out

    try:
        from p01_prime_config import get_ambit_api_url
        from p52_ambit_client import build_oecd_assessment

        data = build_oecd_assessment(text) or {}
        out["endpoint"] = get_ambit_api_url()
        out["status"] = data.get("oecd_connection_status", "N/A")
        out["is_proxy"] = not str(out["status"]).startswith("Connected")
        out["rows"] = [
            ("연결 상태", data.get("oecd_connection_status", "N/A")),
            ("가수분해 (Hydrolysis)", data.get("oecd_hydrolysis", "N/A")),
            ("생분해 (Biodegradability)", data.get("oecd_biodegradability", "N/A")),
            ("변이원성 알러트", data.get("oecd_mutagenicity_alert", "N/A")),
            ("(Q)SAR 엔진", data.get("qsar_engine", "N/A")),
        ]
        out["raw"] = data
    except Exception as exc:
        out["error"] = str(exc)
    return out


def m7_view(result: Any) -> dict[str, Any] | None:
    """ICHM7Result를 템플릿이 쓰기 쉬운 dict로."""
    if result is None:
        return None
    get = lambda name, default=None: getattr(result, name, default)  # noqa: E731
    return {
        "class_number": get("class_number"),
        "class_name_ko": get("class_name_ko", ""),
        "definition": get("definition", ""),
        "control_strategy": get("control_strategy", ""),
        "rationale": get("rationale", ""),
        "alerts": list(get("alerts", []) or []),
        "shared_alerts": list(get("shared_alerts_with_drug", []) or []),
        "alert_present": bool(get("alert_present", False)),
        "empirical_db_match": get("empirical_db_match", "None"),
        "empirical_compound_name": get("empirical_compound_name", ""),
        "ames_evidence": get("ames_evidence", ""),
        "carcinogenicity_evidence": get("carcinogenicity_evidence", ""),
        "recommended_ttc": get("recommended_ttc", ""),
        "consensus_sources": list(get("consensus_sources", []) or []),
        "classifier_engine": get("classifier_engine", ""),
        "layer_c_probability": get("gnn_mutagenic_probability"),
        "layer_c_model": get("gnn_model_name", ""),
        "layer_c_log": get("gnn_fallback_log", ""),
        "ai_inferred": bool(get("ai_inferred_liability", False)),
        "assigned_label": get("assigned_class_label", ""),
        "rule_based_class": get("rule_based_class"),
        "ensemble_log": list(get("ensemble_log", []) or []),
        "is_preliminary": bool(get("is_preliminary", True)),
    }
