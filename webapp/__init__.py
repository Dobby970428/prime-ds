"""PRIME-DS Flask 애플리케이션 팩토리.

Streamlit 판(p00_app.py)을 대체하는 웹 프론트엔드입니다.
화학 계산 엔진(p20~p91)은 그대로 재사용하고 화면만 새로 만들었습니다.

실행:
    python p00_flask_app.py
"""

from __future__ import annotations

import datetime as _dt
import os
import sys
from pathlib import Path

from flask import Flask

# 프로젝트 루트를 import 경로에 올립니다. webapp/ 안에서 pXX_* 모듈을
# 그대로 import할 수 있게 하기 위함입니다.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

PLATFORM_NAME = "PRIME-DS"
PLATFORM_TAGLINE = "Predictive Risk & Impurity Modeling Engine for Drug Substance"


def _secret_key() -> bytes:
    """세션 서명 키. 재시작 때마다 바뀌어도 무방합니다(단일 사용자)."""
    env = os.environ.get("PRIME_FLASK_SECRET", "").strip()
    if env:
        return env.encode("utf-8")
    return os.urandom(32)


def create_app() -> Flask:
    app = Flask(
        __name__,
        template_folder="templates",
        static_folder="static",
    )
    app.config.update(
        SECRET_KEY=_secret_key(),
        PERMANENT_SESSION_LIFETIME=_dt.timedelta(hours=8),
        MAX_CONTENT_LENGTH=32 * 1024 * 1024,
        JSON_AS_ASCII=False,
        TEMPLATES_AUTO_RELOAD=True,
        PRIME_ROOT=str(_ROOT),
    )
    app.jinja_env.trim_blocks = True
    app.jinja_env.lstrip_blocks = True

    from .views import register_blueprints

    register_blueprints(app)

    from . import filters

    filters.register(app)

    return app
