"""템플릿 전역 컨텍스트 — 사이드바 통계 등."""

from __future__ import annotations

import threading
import time
from typing import Any

from flask import Flask

KETCHER_URL = "https://ketcher.mireklzicar.com/"

_lock = threading.Lock()
_cache: dict[str, Any] = {}
_cached_at = 0.0
_TTL = 120.0


def _collect_stats() -> dict[str, int]:
    stats = {"indexed": 0, "empirical": 0, "excipients": 0}

    try:
        from p20_db_index import count_indexed

        stats["indexed"] = int(count_indexed() or 0)
    except Exception:
        pass

    try:
        from p60_m7_empirical_db import empirical_db_stats

        raw = empirical_db_stats() or {}
        for key in ("total", "records", "rows", "count"):
            if isinstance(raw.get(key), int):
                stats["empirical"] = int(raw[key])
                break
    except Exception:
        pass

    try:
        from p71_excipient_opensource_registry import opensource_catalog_stats

        stats["excipients"] = int((opensource_catalog_stats() or {}).get("total", 0))
    except Exception:
        pass

    return stats


def stats(force: bool = False) -> dict[str, int]:
    """사이드바 통계. 무거운 조회라 2분간 캐시합니다."""
    global _cached_at
    now = time.time()
    with _lock:
        if force or not _cache or now - _cached_at > _TTL:
            _cache.clear()
            _cache.update(_collect_stats())
            _cached_at = now
        return dict(_cache)


def invalidate() -> None:
    global _cached_at
    with _lock:
        _cached_at = 0.0


def register(app: Flask) -> None:  # 하위 호환
    register_context(app)


def register_context(app: Flask) -> None:
    from ..i18n import LANGUAGES, current_lang, option_label, t
    from p16_ich_m7_terms import (
        CLASS_ACTIONS,
        CLASS_ACTIONS_KO,
        CLASS_DEFINITIONS,
        CLASS_DEFINITIONS_KO,
        CLASS_SHORT,
        CLASS_SHORT_KO,
        GUIDELINE_ID,
        GUIDELINE_TITLE,
    )

    @app.context_processor
    def _inject() -> dict[str, Any]:
        lang = current_lang()
        return {
            "g_stats": stats(),
            "g_ketcher_url": KETCHER_URL,
            "t": t,
            "opt": option_label,
            "lang": lang,
            "languages": LANGUAGES,
            # ICH M7 공식 용어 — 임의 번역 금지, p16 상수를 그대로 씁니다
            "m7_guideline": GUIDELINE_ID,
            "m7_guideline_title": GUIDELINE_TITLE,
            "m7_short": CLASS_SHORT_KO if lang == "ko" else CLASS_SHORT,
            "m7_definition": CLASS_DEFINITIONS_KO if lang == "ko" else CLASS_DEFINITIONS,
            "m7_action": CLASS_ACTIONS_KO if lang == "ko" else CLASS_ACTIONS,
        }
