"""PRIME-DS — Flask 실행 진입점.

실행:
    python p00_flask_app.py

옵션 (환경변수):
    PRIME_FLASK_HOST    기본 127.0.0.1
    PRIME_FLASK_PORT    기본 8501
    PRIME_FLASK_DEBUG   1 이면 디버그 모드 (코드 수정 시 자동 재시작)
    PRIME_FLASK_SECRET  세션 서명 키 (미설정 시 실행할 때마다 새로 생성)
"""

from __future__ import annotations

import os
import sys
import threading
import webbrowser
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from webapp import create_app  # noqa: E402

app = create_app()


def _boot_checks() -> None:
    """실행 직전에 필수 구성요소를 확인하고 안내합니다."""
    problems: list[str] = []

    try:
        import rdkit  # noqa: F401
    except Exception:
        problems.append("RDKit이 설치되어 있지 않습니다 — pip install rdkit")

    if not (_ROOT / "models" / "ames_ecfp_rf.pkl").exists():
        problems.append(
            "Layer C 모델(models/ames_ecfp_rf.pkl)이 없습니다 — "
            "python s1_setup_layer_c_ames.py 로 학습하세요"
        )

    if not (_ROOT / "cache" / "m7_empirical.db").exists():
        problems.append(
            "실측 DB(cache/m7_empirical.db)가 없습니다 — "
            "python p61_build_m7_empirical_db.py 로 생성하세요"
        )

    forced = os.environ.get("PRIME_LAYER_C_MODEL", "").strip()
    if forced:
        problems.append(
            f"환경변수 PRIME_LAYER_C_MODEL='{forced}' 이(가) 설정되어 있습니다. "
            "ablation 실험용이 아니라면 해제하세요."
        )

    if problems:
        print("\n  ⚠️  확인이 필요한 항목")
        for item in problems:
            print(f"     · {item}")
        print()


def main() -> None:
    host = os.environ.get("PRIME_FLASK_HOST", "127.0.0.1")
    port = int(os.environ.get("PRIME_FLASK_PORT", "8501"))
    debug = os.environ.get("PRIME_FLASK_DEBUG", "").strip() in ("1", "true", "True")

    _boot_checks()

    url = f"http://{host}:{port}/"
    print("=" * 62)
    print("  PRIME-DS  ·  Predictive Risk & Impurity Modeling Engine")
    print("=" * 62)
    print(f"  브라우저에서 열기 :  {url}")
    print("  종료             :  Ctrl+C")
    print("=" * 62)
    print()

    # 리로더가 켜져 있으면 자식 프로세스에서만 브라우저를 엽니다.
    if not debug or os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        threading.Timer(1.2, lambda: webbrowser.open(url)).start()

    app.run(host=host, port=port, debug=debug, use_reloader=debug, threaded=True)


if __name__ == "__main__":
    main()
