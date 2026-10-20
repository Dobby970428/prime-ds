"""의약품명으로 CSV 데이터 파일 찾기·로드."""

# 찾은 약의 CSV 열기

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

APP_DIR = Path(__file__).resolve().parent
DATA_DIR = APP_DIR / "data"
SEARCH_DIRS = [DATA_DIR]

GLOB_PATTERNS = ("CID*.csv", "CID*.xlsx", "CID*.xlsm")


def normalize_name(text: str) -> str:
    return re.sub(r"[^a-z0-9가-힣]+", "", text.lower())


def file_drug_hint(path: Path) -> str:
    # CID1645_3-AMINOBENZAMIDE.csv -> 3-AMINOBENZAMIDE
    stem = path.stem
    if "_" in stem:
        return stem.split("_", 1)[1]
    return stem


def discover_data_files() -> list[Path]:
    found: dict[str, Path] = {}
    for d in SEARCH_DIRS:
        if not d.is_dir():
            continue
        for pattern in GLOB_PATTERNS:
            for p in d.rglob(pattern):
                if p.is_file():
                    found[str(p.resolve())] = p
    return sorted(found.values(), key=lambda p: p.name.lower())


_FILE_LIST_CACHE: list[Path] | None = None


def get_data_files_cached() -> list[Path]:
    """16k+ 파일 목록 — 프로세스당 1회 스캔."""
    global _FILE_LIST_CACHE
    if _FILE_LIST_CACHE is None:
        _FILE_LIST_CACHE = discover_data_files()
    return _FILE_LIST_CACHE


def discover_csv_files() -> list[Path]:
    return discover_data_files()


def match_score_for_labels(query: str, labels: list[str]) -> float:
    qn = normalize_name(query)
    q_tokens = set(re.findall(r"[a-z0-9]+", query.lower()))
    best = 0.0
    for label in labels:
        if not label:
            continue
        ln = normalize_name(label)
        if not ln:
            continue
        if qn == ln:
            best = max(best, 100.0)
        elif qn in ln or ln in qn:
            best = max(best, 80.0)
        else:
            lt = set(re.findall(r"[a-z0-9]+", label.lower()))
            overlap = len(q_tokens & lt) / max(len(q_tokens), 1)
            best = max(best, overlap * 60.0)
    return best


def _row_drug_name(path: Path) -> str | None:
    try:
        df = pd.read_csv(path, usecols=["Name"], nrows=5)
        if "Name" in df.columns and len(df):
            return str(df["Name"].iloc[0])
    except Exception:
        pass
    return None


def match_drug_files(query: str, files: list[Path] | None = None) -> list[tuple[Path, float, str]]:
    """(경로, 점수, 표시명) — 점수 높을수록 유사."""
    if not query.strip():
        return []
    files = files or discover_data_files()

    scored: list[tuple[Path, float, str]] = []
    for path in files:
        hint = file_drug_hint(path)
        best = match_score_for_labels(query, [hint, path.stem])
        if best >= 30.0:
            scored.append((path, best, hint))

    scored.sort(key=lambda x: (-x[1], x[2]))
    return scored


def load_related_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path)
    if "SMILES" not in df.columns:
        raise ValueError(f"SMILES 열이 없습니다: {path.name}")
    return df


def resolve_data_path(stored_path: Path, file_name: str) -> Path:
    """DB에 저장된 경로가 없으면 data/ 폴더에서 파일명으로 재탐색."""
    if stored_path.exists():
        return stored_path
    alt = DATA_DIR / file_name
    if alt.exists():
        return alt
    return stored_path
