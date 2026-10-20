"""PRIME-DS 로컬 설정 — 환경변수 없이 프로젝트 폴더에서 자동 로드."""

# 설정값 읽기

# ① 환경변수를 본다        (PRIME_LAYER_C_MODEL 같은 것)
# ② 없으면 prime_settings.json 을 본다
# ③ 그것도 없으면 기본값

from __future__ import annotations

import json
from pathlib import Path

_SETTINGS_PATH = Path(__file__).resolve().parent / "prime_settings.json"

_DEFAULTS = {
    "ambit_api_url": "http://localhost:8080/ambit2",
    "ambit_query_timeout_sec": 3,
    "use_rdkit_oecd_fallback": True,
    "ambit_bridge_host": "127.0.0.1",
    "ambit_bridge_port": 8080,
}


def load_prime_settings() -> dict:
    settings = dict(_DEFAULTS)
    if _SETTINGS_PATH.is_file():
        try:
            with open(_SETTINGS_PATH, encoding="utf-8") as handle:
                loaded = json.load(handle)
            if isinstance(loaded, dict):
                settings.update(loaded)
        except Exception:
            pass
    return settings


def get_ambit_api_url() -> str:
    import os

    env = os.environ.get("AMBIT_API_URL", "").strip()
    if env:
        return env.rstrip("/")
    return str(load_prime_settings().get("ambit_api_url", _DEFAULTS["ambit_api_url"])).rstrip("/")


def get_ambit_query_timeout() -> float:
    import os

    env = os.environ.get("AMBIT_QUERY_TIMEOUT", "").strip()
    if env:
        try:
            return float(env)
        except ValueError:
            pass
    try:
        return float(load_prime_settings().get("ambit_query_timeout_sec", 3))
    except (TypeError, ValueError):
        return 3.0


def use_rdkit_oecd_fallback() -> bool:
    """AMBIT 미연결 시 RDKit (Q)SAR proxy 사용 여부 (기본 True)."""
    import os

    env = os.environ.get("PRIME_USE_RDKIT_OECD_FALLBACK", "").strip().lower()
    if env in ("0", "false", "no", "off"):
        return False
    if env in ("1", "true", "yes", "on"):
        return True
    val = load_prime_settings().get("use_rdkit_oecd_fallback", True)
    return bool(val) if val is not None else True


def get_ambit_bridge_host() -> str:
    return str(load_prime_settings().get("ambit_bridge_host", _DEFAULTS["ambit_bridge_host"]))


def get_ambit_bridge_port() -> int:
    try:
        return int(load_prime_settings().get("ambit_bridge_port", _DEFAULTS["ambit_bridge_port"]))
    except (TypeError, ValueError):
        return int(_DEFAULTS["ambit_bridge_port"])


def get_ambit_bridge_bind() -> str:
    import os

    return os.environ.get("AMBIT_BRIDGE_HOST", get_ambit_bridge_host())
