"""OECD QSAR Toolbox (AMBIT2) REST API 클라이언트 + RDKit (Q)SAR offline proxy."""

# AMBIT2 서버에 요청

from __future__ import annotations

import json
from typing import Any

import requests

from p53_oecd_rdkit_fallback import rdkit_oecd_qsar_proxy
from p01_prime_config import get_ambit_api_url, get_ambit_query_timeout, use_rdkit_oecd_fallback

DEFAULT_AMBIT_API_URL = "http://localhost:8080/ambit2"
_QSAR_CACHE_VERSION = "v2-rdkit-proxy"
_QSAR_CACHE: dict[str, dict[str, Any]] = {}


def _refresh_ambit_config() -> tuple[str, float]:
    return get_ambit_api_url(), get_ambit_query_timeout()


AMBIT_API_URL = get_ambit_api_url()
AMBIT_QUERY_TIMEOUT = get_ambit_query_timeout()


def _cache_key(smiles: str) -> str:
    return f"{_QSAR_CACHE_VERSION}:{smiles}"


def _offline_fallback(smiles: str, *, reason: str = "connection") -> dict[str, Any]:
    """AMBIT 미연결 시 RDKit (Q)SAR proxy로 실측 가능한 예측값 반환."""
    if not use_rdkit_oecd_fallback():
        return {
            "Hydrolysis_Stability": "Data Missing (Server Offline)",
            "Biodegradability": "Data Missing (Server Offline)",
            "OECD_Alert_Found": "AMBIT unavailable — enable use_rdkit_oecd_fallback",
            "source": "offline",
            "fallback_reason": reason,
        }
    result = dict(rdkit_oecd_qsar_proxy(smiles))
    result["fallback_reason"] = reason
    return result


def query_ambit_qsar_properties(smiles: str) -> dict[str, Any]:
    """
    SMILES → OECD QSAR (AMBIT2) 조회.
    AMBIT 미가동·오류 시 RDKit (Q)SAR proxy로 hydrolysis/biodegradability/alert 산출.
    """
    smi = (smiles or "").strip()
    if not smi:
        return {"error": "Empty SMILES"}

    key = _cache_key(smi)
    if key in _QSAR_CACHE:
        return dict(_QSAR_CACHE[key])

    payload = {
        "smiles": smi,
        "features": ["hydrolysis_rate", "biodegradability", "mutagenicity_alert"],
    }
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }

    try:
        api_url, timeout = _refresh_ambit_config()
        response = requests.post(
            f"{api_url}/query/property",
            data=json.dumps(payload),
            headers=headers,
            timeout=timeout,
        )
        if response.status_code == 200:
            api_result = response.json()
            # 2026-09-03 수정: HTTP 200이면 무조건 "AMBIT2"로 덮어쓰던 것을,
            # 응답 body의 source/engine을 신뢰하도록 변경했습니다.
            # p55_ambit_bridge_server는 자신을 "PRIME-DS RDKit (Q)SAR proxy"로
            # 정직하게 표시하는데도 보고서에는 "OECD QSAR Toolbox (AMBIT2) REST
            # API / Connected"로 찍혀 정식 평가보고서로 직행했습니다.
            reported_source = str(api_result.get("source") or "")
            reported_engine = str(api_result.get("engine") or "")
            is_local_bridge = "PRIME-DS" in f"{reported_source} {reported_engine}"
            refined_result: dict[str, Any] = {
                "Hydrolysis_Stability": api_result.get("hydrolysis_rate", "N/A"),
                "Biodegradability": api_result.get("biodegradability", "N/A"),
                "OECD_Alert_Found": api_result.get("mutagenicity_alert", "None detected"),
                "source": "rdkit_proxy_bridge" if is_local_bridge else "AMBIT2",
                "server_reported_source": reported_source,
                "server_reported_engine": reported_engine,
            }
            _QSAR_CACHE[key] = refined_result
            return refined_result

        result = _offline_fallback(
            smi,
            reason=f"HTTP {response.status_code}",
        )
        _QSAR_CACHE[key] = result
        return result
    except requests.exceptions.ConnectionError:
        result = _offline_fallback(smi, reason="connection")
        _QSAR_CACHE[key] = result
        return result
    except requests.RequestException as exc:
        result = _offline_fallback(smi, reason=str(exc)[:120])
        _QSAR_CACHE[key] = result
        return result


BRIDGE_ENGINE_LABEL = (
    "PRIME-DS local bridge (RDKit computation, not OECD Toolbox)"
)


def _connection_status_label(data: dict[str, Any]) -> str:
    source = data.get("source", "")
    if source == "AMBIT2":
        return "Connected (AMBIT2 REST API)"
    # 2026-09-03 추가: 흉내 서버(p55) 응답을 «연결됨(OECD)»으로 표기하지 않습니다.
    if source == "rdkit_proxy_bridge":
        engine = str(data.get("server_reported_engine") or "").strip()
        return f"{BRIDGE_ENGINE_LABEL}" + (f" — {engine}" if engine else "")
    if source == "rdkit_proxy":
        reason = data.get("fallback_reason", "AMBIT unavailable")
        return (
            "PRIME-DS RDKit (Q)SAR proxy — AMBIT offline "
            f"({reason}; values computed locally)"
        )
    return "Unavailable"


def build_oecd_assessment(smiles: str) -> dict[str, str]:
    """보고서·UI 공통 OECD QSAR 결과 번들."""
    data = query_ambit_qsar_properties(smiles)
    if data.get("error"):
        proxy = rdkit_oecd_qsar_proxy((smiles or "").strip())
        data = {**proxy, "fallback_reason": data["error"]}

    hydro = str(data.get("Hydrolysis_Stability", "N/A"))
    bio = str(data.get("Biodegradability", "N/A"))
    alert = str(data.get("OECD_Alert_Found", "N/A"))
    # 2026-09-03 수정: 로컬 브리지(p55) 응답은 OECD Toolbox가 아니므로
    # qsar_engine·oecd_connection_status에 «AMBIT2 / Connected»를 쓰지 않습니다.
    source = data.get("source")
    connected = source == "AMBIT2"
    if connected:
        engine = "OECD QSAR Toolbox (AMBIT2) REST API"
    elif source == "rdkit_proxy_bridge":
        engine = (
            "PRIME-DS RDKit (Q)SAR proxy via local AMBIT2-compatible bridge "
            "(not OECD QSAR Toolbox)"
        )
    else:
        engine = "PRIME-DS RDKit (Q)SAR proxy (OECD 301 / ICH M7 read-across)"
    ames = f"OECD mutagenicity alert: {alert}; hydrolysis={hydro}; biodegradability={bio}"
    return {
        "qsar_engine": engine,
        "qsar_ames": ames,
        "oecd_connection_status": _connection_status_label(data),
        "oecd_hydrolysis": hydro,
        "oecd_biodegradability": bio,
        "oecd_mutagenicity_alert": alert,
    }


def ambit_qsar_report_fields(smiles: str) -> tuple[str, str]:
    """규제 보고서용 (Q)SAR engine / Ames prediction 문자열."""
    oecd = build_oecd_assessment(smiles)
    return oecd["qsar_engine"], oecd["qsar_ames"]


def ambit_opensource_rows(smiles: str) -> list[dict[str, str]]:
    """하이브리드 독성 평가 UI 표용 AMBIT / RDKit proxy 결과 행."""
    data = query_ambit_qsar_properties(smiles)
    if data.get("error") and data.get("source") != "rdkit_proxy":
        proxy = rdkit_oecd_qsar_proxy((smiles or "").strip())
        data = proxy

    if data.get("source") == "AMBIT2":
        status = "Connected (AMBIT2)"
    # 2026-09-03 추가: 로컬 브리지 케이스를 별도 표기.
    elif data.get("source") == "rdkit_proxy_bridge":
        status = BRIDGE_ENGINE_LABEL
    elif data.get("source") == "rdkit_proxy":
        status = "RDKit (Q)SAR proxy (AMBIT offline — local computation)"
    else:
        status = "Unavailable"

    return [
        {"항목": "OECD QSAR (AMBIT2)", "결과": status},
        {"항목": "  └ Hydrolysis stability", "결과": str(data.get("Hydrolysis_Stability", "N/A"))},
        {"항목": "  └ Biodegradability", "결과": str(data.get("Biodegradability", "N/A"))},
        {"항목": "  └ OECD mutagenicity alert", "결과": str(data.get("OECD_Alert_Found", "N/A"))},
    ]


if __name__ == "__main__":
    test_smiles = "CN(C)N=O"
    qsar_data = query_ambit_qsar_properties(test_smiles)
    print(f"[{test_smiles}] OECD QSAR result:")
    print(json.dumps(qsar_data, indent=4, ensure_ascii=False))
    print("\nAssessment bundle:")
    print(json.dumps(build_oecd_assessment(test_smiles), indent=4, ensure_ascii=False))
