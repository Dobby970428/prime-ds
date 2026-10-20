"""PRIME-DS Flask — 서버측 세션 저장소.

Streamlit의 st.session_state를 대체합니다.

Streamlit은 브라우저 탭 하나에 파이썬 세션 하나가 붙어 있어서
st.session_state에 무엇이든 담을 수 있었습니다. Flask는 요청이 끝나면
메모리가 사라지고, 쿠키 세션은 4KB 제한이 있어 분해산물 리스트나
리포트 컨텍스트 객체를 담을 수 없습니다.

그래서 쿠키에는 짧은 식별자(sid)만 넣고, 실제 데이터는 이 모듈의
프로세스 메모리에 보관합니다. 단일 사용자 데스크톱 실행이 전제입니다.
"""

from __future__ import annotations

import threading
import time
import uuid
from typing import Any

from flask import session

# 세션 하나가 유지되는 시간 (초). 기본 8시간.
_TTL_SECONDS = 8 * 60 * 60

# 세션 정리를 시도하는 최소 간격 (초).
_SWEEP_INTERVAL = 300

_lock = threading.RLock()
_store: dict[str, dict[str, Any]] = {}
_touched: dict[str, float] = {}
_last_sweep = 0.0

SESSION_COOKIE_KEY = "prime_sid"


def _sweep_locked(now: float) -> None:
    """만료된 세션을 지웁니다. 반드시 _lock을 잡은 상태에서 호출하세요."""
    global _last_sweep
    if now - _last_sweep < _SWEEP_INTERVAL:
        return
    _last_sweep = now
    dead = [sid for sid, ts in _touched.items() if now - ts > _TTL_SECONDS]
    for sid in dead:
        _store.pop(sid, None)
        _touched.pop(sid, None)


def current_sid() -> str:
    """현재 브라우저 세션의 식별자. 없으면 새로 만듭니다."""
    sid = session.get(SESSION_COOKIE_KEY)
    if not sid:
        sid = uuid.uuid4().hex
        session[SESSION_COOKIE_KEY] = sid
        session.permanent = True
    return sid


def bucket() -> dict[str, Any]:
    """현재 세션의 저장 공간을 돌려줍니다."""
    sid = current_sid()
    now = time.time()
    with _lock:
        _sweep_locked(now)
        _touched[sid] = now
        return _store.setdefault(sid, {})


def get(key: str, default: Any = None) -> Any:
    return bucket().get(key, default)


def put(key: str, value: Any) -> Any:
    bucket()[key] = value
    return value


def pop(key: str, default: Any = None) -> Any:
    return bucket().pop(key, default)


def has(key: str) -> bool:
    return key in bucket()


def clear_keys(*keys: str) -> None:
    b = bucket()
    for key in keys:
        b.pop(key, None)


def clear_prefix(prefix: str) -> None:
    """접두어로 시작하는 키를 전부 지웁니다."""
    b = bucket()
    for key in [k for k in b if k.startswith(prefix)]:
        b.pop(key, None)


def reset() -> None:
    """현재 세션을 완전히 비웁니다."""
    sid = current_sid()
    with _lock:
        _store.pop(sid, None)
        _touched.pop(sid, None)


# ---------------------------------------------------------------------------
# 표시 설정 — Streamlit에서 전역 변수였던 MAX_SHOW / COLS / IMG_SIZE / SHOW_SMI
# ---------------------------------------------------------------------------

DISPLAY_DEFAULTS: dict[str, Any] = {
    "max_show": 12,
    "cols": 3,
    "img_size": "300x300",
    "show_smiles": False,
}


def display_prefs() -> dict[str, Any]:
    prefs = bucket().setdefault("_display", dict(DISPLAY_DEFAULTS))
    for key, value in DISPLAY_DEFAULTS.items():
        prefs.setdefault(key, value)
    return prefs


def update_display_prefs(**kwargs: Any) -> dict[str, Any]:
    prefs = display_prefs()
    for key, value in kwargs.items():
        if key in DISPLAY_DEFAULTS and value is not None:
            prefs[key] = value
    return prefs


# ---------------------------------------------------------------------------
# 캔버스 SMILES — Streamlit의 mol_smiles_{key_suffix}
# ---------------------------------------------------------------------------

CANVAS_KEYS = (
    "single_search",
    "single_reference_drug",
    "reactant_1",
    "reactant_2",
    "reagent",
    "product",
    "reference_drug",
)


def canvas() -> dict[str, str]:
    return bucket().setdefault("_canvas", {})


def get_canvas(name: str) -> str:
    return (canvas().get(name) or "").strip()


def set_canvas(name: str, smiles: str | None) -> str:
    value = (smiles or "").strip()
    canvas()[name] = value
    return value
