"""Jinja 템플릿 필터 및 전역 함수."""

from __future__ import annotations

import html
import json
from typing import Any

from flask import Flask


def _shorten(value: Any, limit: int = 48) -> str:
    text = "" if value is None else str(value)
    if len(text) <= limit:
        return text
    return text[: max(limit - 1, 1)] + "…"


def _fmt(value: Any, digits: int = 2, dash: str = "—") -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return dash


def _pct(value: Any, digits: int = 1, dash: str = "—") -> str:
    try:
        return f"{float(value):.{digits}f}%"
    except (TypeError, ValueError):
        return dash


def _nl2br(value: Any) -> str:
    text = "" if value is None else str(value)
    return html.escape(text).replace("\n", "<br>")


def _tojson_safe(value: Any) -> str:
    """</script> 로 템플릿이 깨지지 않도록 이스케이프한 JSON."""
    try:
        raw = json.dumps(value, ensure_ascii=False, default=str)
    except (TypeError, ValueError):
        raw = "null"
    return raw.replace("</", "<\\/")


_M7_TONE = {
    1: "danger",
    2: "danger",
    3: "warn",
    4: "info",
    5: "ok",
}


def _m7_tone(value: Any) -> str:
    """ICH M7 Class 번호를 색상 톤 이름으로."""
    try:
        return _M7_TONE.get(int(value), "muted")
    except (TypeError, ValueError):
        return "muted"


def _confidence_tone(value: Any) -> str:
    try:
        score = float(value)
    except (TypeError, ValueError):
        return "muted"
    if score >= 85.0:
        return "danger"
    if score >= 50.0:
        return "warn"
    return "info"


def register(app: Flask) -> None:
    app.jinja_env.filters["shorten"] = _shorten
    app.jinja_env.filters["fmt"] = _fmt
    app.jinja_env.filters["pct"] = _pct
    app.jinja_env.filters["nl2br"] = _nl2br
    app.jinja_env.filters["tojson_safe"] = _tojson_safe
    app.jinja_env.filters["m7_tone"] = _m7_tone
    app.jinja_env.filters["confidence_tone"] = _confidence_tone
