"""PubChem REST API 연동 — 화합물 정보 · 독성/Ames BioAssay 스캔."""

# 외부 조회

from __future__ import annotations

import urllib.parse
from typing import Any

import requests

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "PRIME-DS/1.0 (pharmaceutical-informatics)"})

_GENOTOX_KEYWORDS = (
    "ames",
    "mutagen",
    "mutagenic",
    "genotox",
    "genotoxic",
    "salmonella",
    "bacterial reverse mutation",
    "micronucleus",
    "chromosome aberration",
)


def _safe_get(url: str, *, timeout: int = 14) -> tuple[dict | list | None, str | None]:
    """(parsed_json, error_message)"""
    try:
        r = SESSION.get(url, timeout=timeout)
        if r.status_code == 404:
            return None, "PubChem에 해당 구조가 등록되어 있지 않습니다."
        if r.status_code != 200:
            return None, f"PubChem API 오류 (HTTP {r.status_code})"
        return r.json(), None
    except requests.Timeout:
        return None, "PubChem API 응답 시간 초과"
    except requests.RequestException as exc:
        return None, f"PubChem 네트워크 오류: {exc}"
    except ValueError:
        return None, "PubChem JSON 파싱 실패"


def fetch_compound_by_smiles(smiles: str) -> dict | None:
    """SMILES → PubChem CID, Name, Formula, MW."""
    if not smiles or not str(smiles).strip():
        return None
    encoded = urllib.parse.quote(smiles.strip(), safe="")
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}"
        f"/property/MolecularFormula,MolecularWeight,IUPACName,Title/JSON"
    )
    data, err = _safe_get(url)
    if err or not data:
        return None
    try:
        props = data["PropertyTable"]["Properties"][0]
        cid = props.get("CID")
        return {
            "cid": cid,
            "name": props.get("Title") or props.get("IUPACName") or "",
            "formula": props.get("MolecularFormula", ""),
            "mw": props.get("MolecularWeight"),
            "link": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else "",
            "api_status": "success",
            "api_message": f"🟢 PubChem CID {cid} 조회 성공",
        }
    except (KeyError, IndexError, TypeError):
        return None


def fetch_cid_by_smiles(smiles: str) -> int | None:
    info = fetch_compound_by_smiles(smiles)
    if info and info.get("cid"):
        try:
            return int(info["cid"])
        except (TypeError, ValueError):
            pass
    return None


def fetch_smiles_by_compound_name(name: str) -> str | None:
    """PubChem 화합물명 → Canonical SMILES (부형제 오픈 DB 연동)."""
    q = (name or "").strip()
    if not q:
        return None
    encoded = urllib.parse.quote(q, safe="")
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}"
        f"/property/CanonicalSMILES,IsomericSMILES/JSON"
    )
    data, err = _safe_get(url, timeout=16)
    if err or not data:
        return None
    try:
        props = data["PropertyTable"]["Properties"][0]
        smi = props.get("IsomericSMILES") or props.get("CanonicalSMILES")
        return str(smi).strip() if smi else None
    except (KeyError, IndexError, TypeError):
        return None


def fetch_compound_by_name(name: str) -> dict | None:
    """PubChem 화합물명 → CID, SMILES, MW."""
    smi = fetch_smiles_by_compound_name(name)
    if not smi:
        return None
    base = fetch_compound_by_smiles(smi) or {}
    return {
        **base,
        "smiles": smi,
        "query_name": name.strip(),
        "source": "PubChem",
    }


def _assay_text(assay: dict) -> str:
    parts = [
        str(assay.get("Name") or ""),
        str(assay.get("Description") or ""),
        str(assay.get("ProtocolName") or ""),
    ]
    return " ".join(parts).lower()


def fetch_genotoxicity_assays(cid: int | str, *, max_hits: int = 8) -> tuple[list[dict], str | None]:
    """
    PubChem BioAssay 요약에서 Ames / 변이원성 관련 시험 검색.
    Returns: (hits, error_message)
    """
    try:
        cid_int = int(cid)
    except (TypeError, ValueError):
        return [], "유효하지 않은 PubChem CID"

    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_int}"
        f"/assaysummary/JSON"
    )
    data, err = _safe_get(url, timeout=18)
    if err:
        return [], err
    if not data:
        return [], None

    try:
        summaries = data.get("AssaySummaries", {}).get("AssaySummary", [])
        if isinstance(summaries, dict):
            summaries = [summaries]
    except Exception:
        return [], "BioAssay 요약 파싱 실패"

    hits: list[dict] = []
    for assay in summaries:
        if not isinstance(assay, dict):
            continue
        text = _assay_text(assay)
        if not any(kw in text for kw in _GENOTOX_KEYWORDS):
            continue
        hits.append(
            {
                "aid": assay.get("AID"),
                "name": assay.get("Name") or assay.get("Description") or "BioAssay",
                "activity_outcome": assay.get("ActivityOutcome") or assay.get("Outcome") or "—",
                "source": assay.get("SourceName") or "PubChem",
            }
        )
        if len(hits) >= max_hits:
            break
    return hits, None


def fetch_full_toxicity_profile(smiles: str) -> dict[str, Any]:
    """
    SMILES 기반 PubChem 통합 프로필.
    compound 정보 + Ames/변이원성 BioAssay 매핑 여부.
    """
    result: dict[str, Any] = {
        "compound": None,
        "genotox_assays": [],
        "genotox_found": False,
        "api_status": "not_found",
        "api_message": "⚪ PubChem 검색 기록 없음",
        "error": None,
    }
    if not smiles or not str(smiles).strip():
        result["error"] = "SMILES가 비어 있습니다."
        result["api_message"] = "⚪ SMILES 없음 — PubChem 조회 불가"
        return result

    compound = fetch_compound_by_smiles(smiles)
    if not compound:
        result["api_status"] = "not_found"
        result["api_message"] = "⚪ PubChem에 등록된 동일 구조 없음"
        return result

    result["compound"] = compound
    result["api_status"] = "success"
    result["api_message"] = compound.get("api_message", "🟢 PubChem 조회 성공")

    cid = compound.get("cid")
    if cid:
        assays, assay_err = fetch_genotoxicity_assays(cid)
        result["genotox_assays"] = assays
        result["genotox_found"] = bool(assays)
        if assay_err:
            result["error"] = assay_err
        if assays:
            result["api_message"] = (
                f"🟢 PubChem CID {cid} — Ames/변이원성 관련 BioAssay {len(assays)}건 매핑"
            )

    return result
