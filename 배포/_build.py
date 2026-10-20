"""PRIME-DS 포터블 배포본 — 앱 파일 복사 · 런처 작성 · 검증 · 압축.

배포_빌드.bat 가 내려받은 embeddable 런타임으로 실행합니다. 직접 실행하지 마십시오.

    runtime\\python.exe _build.py <프로젝트폴더> <패키지폴더> [키트폴더]

키트폴더는 _verify.py 가 있는 곳(= 이 파일이 있는 곳)입니다. 생략하면 이
파일의 위치를 씁니다. 빌드 키트를 프로젝트 밖에 두어도 동작하게 하려고
경로를 인자로 받습니다.

하는 일
  1. 배포본에 들어갈 파일만 골라 <패키지폴더>\\app\\ 로 복사
  2. «PRIME-DS 실행.bat» 과 README.txt 작성
  3. 복사된 배포본을 **실제로 기동**해 검증 (_verify.py)
  4. ZIP 으로 묶기

들어가는 것 / 빠지는 것은 아래 INCLUDE·EXCLUDE 표가 전부입니다.
새 모듈을 추가하면 이 표도 같이 고쳐야 합니다.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

# ── 들어갈 것 ────────────────────────────────────────────────────────
ROOT_PY_GLOB = "p*.py"          # 엔진 모듈 전체 (p00 ~ p96)
ROOT_PY_SKIP = {
    # deepchem 이 있어야 돌아가는 체크포인트 생성 스크립트. 배포본에 둘 이유가 없다.
    "p97_build_attentivefp_model.py",
}
ROOT_EXTRA = [
    "prime_settings.json",
    "README.md",
    "PRIME-DS_쉬운설명_제제연구자용.md",
]
DIR_ALL = ["engines", "webapp", "templates"]     # __pycache__ 는 제외
MODELS_KEEP = ["ames_ecfp_rf.pkl", "tox21_ecfp_rf_sklearn.pkl", "tox21.csv.gz"]
CACHE_KEEP = [
    "m7_empirical.db",          # Layer A 실측 기록
    "ad_train_fps.pkl",         # 적용영역 지문
    "substruct_primary.pkl",    # 부분구조 검색 색인
    "substruct_manifest.json",
]

# ── 반드시 있어야 하는 것 (없으면 빌드 중단) ─────────────────────────
REQUIRED = [
    "p00_flask_app.py",
    "webapp/__init__.py",
    "models/ames_ecfp_rf.pkl",
    "cache/m7_empirical.db",
]

# ── 빠지는 것 ────────────────────────────────────────────────────────
#   data\ · drug_index.db   → 검색 데이터 추가팩 (검색데이터_추가팩_빌드.bat)
#   s1~s26, 논문자료\        → 논문 재현용. 배포본에 넣지 않기로 함
#   RxnMetrics\             → 옛 Streamlit 앱. 자체 venv 를 들고 있음
#   toxtree_sources\        → Toxtree 원본 추출물. 참고자료이지 실행에 쓰이지 않음
#   _archive\ _backup*      → 백업
EXCLUDE_DIRS = {"__pycache__", ".git", "venv", ".venv", ".streamlit"}


def log(msg: str) -> None:
    print(f"       {msg}", flush=True)


def copy_tree(src: Path, dst: Path) -> int:
    n = 0
    for root, dirs, files in os.walk(src):
        dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS]
        rel = Path(root).relative_to(src)
        (dst / rel).mkdir(parents=True, exist_ok=True)
        for f in files:
            if f.endswith((".pyc", ".pyo")):
                continue
            shutil.copy2(Path(root) / f, dst / rel / f)
            n += 1
    return n


LAUNCHER = """@echo off
chcp 65001 >nul
title PRIME-DS

echo.
echo  ==============================================================
echo   PRIME-DS  ·  불순물 예측 및 ICH M7 유전독성 평가 플랫폼
echo  ==============================================================
echo.
echo   잠시 뒤 브라우저가 자동으로 열립니다.
echo   이 창을 닫거나 Ctrl+C 를 누르면 종료됩니다.
echo.

cd /d "%~dp0app"
"%~dp0runtime\\python.exe" p00_flask_app.py

if errorlevel 1 (
    echo.
    echo  ------------------------------------------------------------
    echo   실행에 실패했습니다. 위 메시지를 그대로 전달해 주세요.
    echo  ------------------------------------------------------------
    echo.
    pause
)
"""

SELFCHECK_BAT = """@echo off
chcp 65001 >nul
title PRIME-DS 자가점검

rem  배포본이 제대로 동작하는지 확인합니다.
rem  다른 PC 로 옮긴 뒤, 또는 무언가 이상할 때 돌려 보세요.

cd /d "%~dp0app"
"%~dp0runtime\\python.exe" _self_check.py
echo.
pause
"""

README = """PRIME-DS 포터블 배포본
=====================================================================

실행
---------------------------------------------------------------------
  «PRIME-DS 실행.bat» 을 더블클릭하세요.
  브라우저가 http://127.0.0.1:8501/ 로 자동으로 열립니다.
  끌 때는 검은 창을 닫거나 Ctrl+C 를 누르세요.

  Python 을 설치할 필요가 없습니다. 필요한 것은 전부 runtime\\ 폴더
  안에 들어 있습니다. 레지스트리나 시스템 설정은 건드리지 않습니다.
  지울 때는 이 폴더를 통째로 삭제하면 끝입니다. USB 에 담아 옮겨도
  그대로 동작합니다.

요구 사항
---------------------------------------------------------------------
  · Windows 10 (64비트) 이상
  · 디스크 여유 공간 {size}
  · 인터넷은 필요하지 않습니다.
    다만 PubChem 조회 기능만은 연결되어 있을 때 동작합니다.

들어 있는 것
---------------------------------------------------------------------
  · 평가 4개 모듈 — 단일 화합물 · 공정 · 분해 · 첨가제 상호작용
  · 4계층 ICH M7 판정 엔진 (실측 기록 · 구조 경보 · 통계 · 읽기교차)
  · 보고서 출력 (Word · PDF · Excel)

들어 있지 않은 것
---------------------------------------------------------------------
  · 실측 기록 검색(Search) 모듈의 데이터
    왼쪽 메뉴의 Search 는 열리지만 결과가 비어 있습니다.
    쓰시려면 «PRIME-DS_검색데이터팩.zip» 을 이 폴더에 풀어 주세요.
  · 논문 재현용 스크립트 (s1 ~ s26)

문제가 생기면
---------------------------------------------------------------------
  · 8501 포트를 다른 프로그램이 쓰고 있으면 실행되지 않습니다.
    그때는 명령 프롬프트에서
        set PRIME_FLASK_PORT=8600
    을 먼저 실행한 뒤 «PRIME-DS 실행.bat» 을 실행하세요.
  · 백신이 runtime\\python.exe 를 막는 경우가 있습니다. 예외로 등록해
    주세요.

빌드 정보
---------------------------------------------------------------------
{buildinfo}
"""


def main() -> int:
    if len(sys.argv) < 3:
        print("사용법: _build.py <프로젝트폴더> <패키지폴더>")
        return 2
    proj = Path(sys.argv[1]).resolve()
    pkg = Path(sys.argv[2]).resolve()
    kit = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else Path(__file__).resolve().parent
    app = pkg / "app"
    rt = pkg / "runtime"

    if not (kit / "_verify.py").exists():
        print(f"\n  [실패] 자가점검 스크립트를 찾지 못했습니다: {kit / '_verify.py'}")
        print("         빌드 키트 7개 파일이 같은 폴더에 있어야 합니다.")
        return 1

    missing = [r for r in REQUIRED if not (proj / r).exists()]
    if missing:
        print("\n  [실패] 배포본에 반드시 필요한 파일이 없습니다:")
        for m in missing:
            print(f"     · {m}")
        if "models/ames_ecfp_rf.pkl" in missing:
            print("\n     Layer C 모델이 없습니다. s1_setup_layer_c_ames.py 로 만들 수")
            print("     있지만, 이 스크립트는 학습·검증 분할을 다시 생성하므로")
            print("     논문 수치와 어긋나게 됩니다. 기존 모델 파일을 복구하세요.")
        return 1

    app.mkdir(parents=True, exist_ok=True)
    n = 0

    # 1) 루트 엔진 모듈
    for f in sorted(proj.glob(ROOT_PY_GLOB)):
        if f.name in ROOT_PY_SKIP:
            continue
        shutil.copy2(f, app / f.name)
        n += 1
    log(f"엔진 모듈 {n}개")

    # 2) 루트 부속 파일
    for name in ROOT_EXTRA:
        src = proj / name
        if src.exists():
            shutil.copy2(src, app / name)
            n += 1

    # 3) 폴더 통째로
    for d in DIR_ALL:
        src = proj / d
        if not src.is_dir():
            continue
        c = copy_tree(src, app / d)
        n += c
        log(f"{d}\\ {c}개")

    # 4) 모델 · 캐시는 고른 것만
    for folder, keep in (("models", MODELS_KEEP), ("cache", CACHE_KEEP)):
        (app / folder).mkdir(parents=True, exist_ok=True)
        got = 0
        for name in keep:
            src = proj / folder / name
            if src.exists():
                shutil.copy2(src, app / folder / name)
                got += 1
                n += 1
        log(f"{folder}\\ {got}/{len(keep)}개")

    # 5) 자가점검 스크립트 · 런처 · README
    #    _verify.py 는 app\ 안에 놓아야 한다. 자기 위치를 기준으로 models\ · cache\
    #    를 찾고, 스크립트가 있는 폴더가 sys.path 에 들어가야 p62·webapp 을 import
    #    할 수 있기 때문이다. 배포본에 남겨 두면 사용자가 나중에 다시 돌릴 수 있다.
    shutil.copy2(kit / "_verify.py", app / "_self_check.py")
    n += 1
    (pkg / "PRIME-DS 실행.bat").write_text(LAUNCHER, encoding="utf-8", newline="\r\n")
    (pkg / "자가점검.bat").write_text(SELFCHECK_BAT, encoding="utf-8", newline="\r\n")

    total = sum(f.stat().st_size for f in pkg.rglob("*") if f.is_file())
    buildinfo = (
        f"  Python   : {sys.version.split()[0]}\n"
        f"  빌드 파일 : {n}개\n"
        f"  크기     : {total / 1e9:.2f} GB"
    )
    (pkg / "README.txt").write_text(
        README.format(size=f"약 {total / 1e9:.1f} GB", buildinfo=buildinfo),
        encoding="utf-8", newline="\r\n",
    )
    log(f"앱 파일 {n}개 · 전체 {total / 1e9:.2f} GB")

    # 6) 검증 — 복사된 배포본을 실제로 기동한다
    print()
    print("       ── 검증 ──────────────────────────────────────────")
    r = subprocess.run([str(rt / "python.exe"), "_self_check.py"], cwd=str(app))
    if r.returncode != 0:
        print("\n  [실패] 검증을 통과하지 못했습니다. 배포본을 배포하지 마십시오.")
        return 1

    # 7) 압축
    print()
    log("압축 중 (몇 분 걸립니다)")
    zpath = pkg.parent / "PRIME-DS_포터블.zip"
    if zpath.exists():
        zpath.unlink()
    files = [f for f in pkg.rglob("*") if f.is_file()]
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for i, f in enumerate(files, 1):
            z.write(f, Path(pkg.name) / f.relative_to(pkg))
            if i % 2000 == 0:
                log(f"  {i}/{len(files)}")
    log(f"ZIP {zpath.stat().st_size / 1e9:.2f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
