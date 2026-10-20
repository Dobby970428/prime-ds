"""PRIME-DS 공정 유연물질 엔진 — 로컬 캐시 경로."""

from __future__ import annotations

from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = _PROJECT_ROOT / ".prime_process_cache"
VERIFICATION_FP_CACHE_PATH = CACHE_DIR / "verification_fp_cache.pkl"


def ensure_cache_dir() -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR
