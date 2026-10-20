"""검색 데이터 추가팩 — drug_index.db + data\\ 를 ZIP 으로 묶습니다.

검색데이터_추가팩_빌드.bat 가 호출합니다.

    python _build_searchpack.py <프로젝트폴더> <출력폴더>

ZIP 안의 경로를 포터블 배포본의 구조에 맞춰 두었으므로,
받는 쪽은 **PRIME-DS 폴더에 그대로 압축을 풀면** 끝납니다.

    PRIME-DS\\
      app\\drug_index.db      ← 이 팩이 넣는 것
      app\\data\\CID*.csv      ← 이 팩이 넣는 것
      검색데이터_재색인.bat     ← 이 팩이 넣는 것

왜 재색인이 필요한가
    drug_index.db 는 각 CSV 의 **절대 경로**를 저장합니다(p20:58 resolve()).
    다른 PC 의 다른 폴더에 풀면 그 경로가 가리키는 파일이 없습니다.
    파일명으로 다시 찾는 보조 경로가 p22:118 에 있어 이름 검색은 살아
    있지만, 구조(부분구조) 검색은 어긋납니다. 한 번 재색인하면 전부
    맞아떨어집니다.
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path

REINDEX_BAT = """@echo off
chcp 65001 >nul
title PRIME-DS 검색 색인 다시 만들기

echo.
echo   data\\ 폴더의 CSV 를 다시 훑어 drug_index.db 를 새로 만듭니다.
echo   파일이 많아 5 ~ 15 분쯤 걸립니다. 한 번만 하면 됩니다.
echo.

cd /d "%~dp0app"
"%~dp0runtime\\python.exe" p26_import_index.py
echo.
pause
"""

NOTE = """검색 데이터 추가팩
=====================================================================

이 팩에 들어 있는 것
  · drug_index.db      의약품 파일 색인
  · data\\CID*.csv       PubChem 실측 기록

설치
  1. PRIME-DS 포터블 폴더(«PRIME-DS 실행.bat» 이 있는 그 폴더)에
     이 ZIP 을 그대로 풉니다. 폴더 구조가 맞춰져 있습니다.
  2. «검색데이터_재색인.bat» 을 한 번 실행합니다. (5 ~ 15 분)
  3. «PRIME-DS 실행.bat» 으로 띄우면 왼쪽 메뉴의 Search 가 동작합니다.

2번을 건너뛰면
  이름 검색은 되지만 구조(부분구조) 검색이 비어 보입니다.
  색인에 들어 있는 경로가 만든 PC 기준의 절대 경로이기 때문입니다.
"""


def main() -> int:
    if len(sys.argv) < 3:
        print("사용법: _build_searchpack.py <프로젝트폴더> <출력폴더>")
        return 2
    proj = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)

    db = proj / "drug_index.db"
    data = proj / "data"
    if not db.exists() and not data.is_dir():
        print("  [실패] drug_index.db 도 data\\ 폴더도 없습니다.")
        return 1

    files: list[tuple[Path, str]] = []
    if db.exists():
        files.append((db, "app/drug_index.db"))
    if data.is_dir():
        for f in sorted(data.rglob("*")):
            if f.is_file() and f.suffix.lower() in (".csv", ".xlsx", ".xlsm"):
                files.append((f, "app/data/" + f.relative_to(data).as_posix()))

    raw = sum(f.stat().st_size for f, _ in files)
    print(f"       대상 {len(files)}개 · 원본 {raw / 1e9:.2f} GB")

    zpath = out / "PRIME-DS_검색데이터팩.zip"
    if zpath.exists():
        zpath.unlink()
    # CSV 는 압축이 잘 듣습니다. 대신 파일 수가 많아 시간이 걸립니다.
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        z.writestr("검색데이터팩_읽어주세요.txt", NOTE.replace("\n", "\r\n"))
        z.writestr("검색데이터_재색인.bat", REINDEX_BAT.replace("\n", "\r\n"))
        for i, (src, arc) in enumerate(files, 1):
            z.write(src, arc)
            if i % 2000 == 0:
                print(f"         {i}/{len(files)}")
    print(f"       ZIP {zpath.stat().st_size / 1e9:.2f} GB  →  {zpath}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
