"""
data/ 폴더의 CSV·엑셀을 스캔해 drug_index.db 에 등록합니다.

사용:
  python p26_import_index.py
"""
# 2026-09-03 수정: 위 사용법의 구식 파일명(import_index.py)을 실제 파일명으로 정정.

# 색인 만들기(별도 실행 파일) - ⚙️인덱스 관리 버튼
# data/CID*.csv 16,823개를 읽는다
#    ↓
# drug_index.db 에 등록

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

from p20_db_index import DB_PATH, connect, count_indexed, init_db, upsert_record
from p22_drug_data import file_drug_hint, normalize_name

APP_DIR = Path(__file__).resolve().parent
DEFAULT_DATA = APP_DIR / "data"
GLOB_PATTERNS = ("CID*.csv", "CID*.xlsx", "CID*.xlsm")
BATCH_COMMIT = 300


def parent_cid_from_name(path: Path) -> str | None:
    m = re.match(r"CID(\d+)", path.stem, re.I)
    return m.group(1) if m else None


def fast_csv_row_count(path: Path) -> int:
    try:
        with open(path, "rb") as f:
            n = sum(1 for _ in f)
        return max(0, n - 1)
    except OSError:
        return 0


def read_table_meta(path: Path) -> tuple[str | None, int, bool]:
    suffix = path.suffix.lower()
    try:
        if suffix == ".csv":
            cols = pd.read_csv(path, nrows=0).columns.tolist()
            has_smiles = "SMILES" in cols
            if not has_smiles:
                return None, 0, False
            df_head = pd.read_csv(path, nrows=1, usecols=lambda c: c in ("Name", "SMILES"))
            row_count = fast_csv_row_count(path)
        else:
            df_head = pd.read_excel(path, nrows=1)
            has_smiles = "SMILES" in df_head.columns
            row_count = 0 if not has_smiles else max(1, len(pd.read_excel(path, usecols=["SMILES"])))
    except Exception as e:
        raise ValueError(f"읽기 실패: {e}") from e

    drug_from_row = None
    if "Name" in df_head.columns and len(df_head):
        val = df_head["Name"].iloc[0]
        if pd.notna(val):
            drug_from_row = str(val)
    return drug_from_row, row_count, has_smiles


def collect_files(data_dir: Path) -> list[Path]:
    # 2026-09-03 수정: glob(1단계) → rglob(재귀)로 통일.
    # 탐색 쪽(p22:31-40)은 d.rglob 을 쓰는데 인덱싱만 data_dir.glob 이라,
    # data/A/CID123_X.csv 같은 하위 폴더 파일이 목록에는 보이지만
    # drug_index.db 에는 등록되지 않아 구조 검색·실측 판정에서만 빠졌습니다.
    found: dict[str, Path] = {}
    for pattern in GLOB_PATTERNS:
        for p in data_dir.rglob(pattern):
            if p.is_file():
                found[str(p.resolve())] = p
    return sorted(found.values(), key=lambda p: p.name.lower())


def import_all(data_dir: Path, verbose: bool = True) -> tuple[int, int, int]:
    init_db()
    files = collect_files(data_dir)
    ok, skip, fail = 0, 0, 0
    conn = connect()

    for i, path in enumerate(files, start=1):
        try:
            row_name, row_count, has_smiles = read_table_meta(path)
            if not has_smiles:
                skip += 1
                continue

            hint = file_drug_hint(path)
            drug_name = row_name or hint.replace("_", " ")
            upsert_record(
                path,
                drug_name=drug_name,
                name_normalized=normalize_name(drug_name) + "|" + normalize_name(hint),
                parent_cid=parent_cid_from_name(path),
                row_count=row_count,
                has_smiles=True,
                conn=conn,
            )
            ok += 1
            if i % BATCH_COMMIT == 0:
                conn.commit()
            if verbose and (i % 500 == 0 or i == len(files)):
                print(f"  ... {i}/{len(files)} (등록 {ok}, skip {skip}, fail {fail})")
        except Exception as e:
            fail += 1
            if verbose and fail <= 20:
                print(f"[FAIL] {path.name}: {e}")

    conn.commit()
    conn.close()
    return ok, skip, fail


def main() -> int:
    parser = argparse.ArgumentParser(description="의약품 데이터 파일 인덱스 구축")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA, help="데이터 폴더")
    args = parser.parse_args()

    data_dir = args.data.resolve()
    if not data_dir.is_dir():
        print(f"폴더 없음: {data_dir}")
        return 1

    print(f"스캔: {data_dir}")
    print(f"DB: {DB_PATH}")
    ok, skip, fail = import_all(data_dir)
    print(f"\n완료 — 등록: {ok}, skip: {skip}, fail: {fail}, DB 총: {count_indexed()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
