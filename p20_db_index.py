"""의약품 파일 메타데이터 SQLite 인덱스 (2만 파일 검색용)."""

# 이름으로 찾기

from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
DB_PATH = APP_DIR / "drug_index.db"


def connect() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def init_db() -> None:
    with connect() as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS drug_files (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                file_path TEXT NOT NULL UNIQUE,
                file_name TEXT NOT NULL,
                drug_name TEXT,
                name_normalized TEXT,
                parent_cid TEXT,
                row_count INTEGER DEFAULT 0,
                has_smiles INTEGER DEFAULT 0,
                indexed_at TEXT NOT NULL
            )
            """
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_drug_name_norm ON drug_files(name_normalized)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_drug_file_name ON drug_files(file_name)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_drug_parent_cid ON drug_files(parent_cid)"
        )


def upsert_record(
    file_path: Path,
    drug_name: str,
    name_normalized: str,
    parent_cid: str | None,
    row_count: int,
    has_smiles: bool,
    conn: sqlite3.Connection | None = None,
) -> None:
    now = datetime.now(timezone.utc).isoformat()
    params = (
        str(file_path.resolve()),
        file_path.name,
        drug_name,
        name_normalized,
        parent_cid,
        row_count,
        1 if has_smiles else 0,
        now,
    )
    sql = """
        INSERT INTO drug_files (
            file_path, file_name, drug_name, name_normalized,
            parent_cid, row_count, has_smiles, indexed_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(file_path) DO UPDATE SET
            file_name=excluded.file_name,
            drug_name=excluded.drug_name,
            name_normalized=excluded.name_normalized,
            parent_cid=excluded.parent_cid,
            row_count=excluded.row_count,
            has_smiles=excluded.has_smiles,
            indexed_at=excluded.indexed_at
    """
    if conn is not None:
        conn.execute(sql, params)
    else:
        with connect() as c:
            c.execute(sql, params)


def count_indexed() -> int:
    # 2026-09-03 수정: init_db() 호출 추가. 빈 폴더 최초 실행에서 테이블이 없어
    # sqlite3.OperationalError: no such table: drug_files 가 났고, 이 예외가
    # p42:875 에 삼켜지며 이미 계산한 Tanimoto 결과까지 버려졌습니다.
    init_db()
    with connect() as conn:
        row = conn.execute("SELECT COUNT(*) AS c FROM drug_files").fetchone()
        return int(row["c"]) if row else 0


def _search_index_once(query: str, conn: sqlite3.Connection) -> list[sqlite3.Row]:
    import re

    from p22_drug_data import normalize_name

    q = query.strip()
    qn = normalize_name(q)
    if not qn and not q:
        return []

    ql = q.lower()
    like_qn = f"%{qn}%" if qn else f"%{ql}%"
    like_ql = f"%{ql}%"

    # 2026-09-03 수정: 도달 불가 elif 삭제. 첫 정규식 ^cid\s*(\d+)$ 의 \s* 가
    # 공백 0개도 허용하므로 "CID14219" 를 이미 첫 분기가 잡습니다
    # (실측: 첫 정규식 True → elif 분기는 절대 실행되지 않음).
    cid_m = re.match(r"^cid\s*(\d+)$", q, re.I)
    if cid_m:
        parent_cid = cid_m.group(1)
    elif q.isdigit():
        parent_cid = q
    else:
        parent_cid = "-1"

    rows = conn.execute(
        """
        SELECT * FROM drug_files
        WHERE name_normalized LIKE ?
           OR lower(file_name) LIKE ?
           OR lower(drug_name) LIKE ?
           OR parent_cid = ?
        LIMIT 800
        """,
        (like_qn, like_ql, like_ql, parent_cid),
    ).fetchall()

    if len(rows) < 3:
        rows = conn.execute(
            "SELECT * FROM drug_files WHERE lower(file_name) LIKE ? LIMIT 800",
            (like_ql,),
        ).fetchall()
    return rows


def search_index(query: str, limit: int = 30) -> list[dict]:
    # 2026-09-03 수정: 미사용 import 제거(`re`, `normalize_name`).
    # 두 이름 모두 이 함수 본문에서 쓰이지 않으며, 실제 사용처인
    # _search_index_once() 는 자체적으로 import 합니다.
    from p22_drug_data import match_score_for_labels
    from p21_drug_search_aliases import expand_search_queries

    q = query.strip()
    if not q:
        return []

    # 2026-09-03 수정: init_db() 호출 추가 (count_indexed() 와 같은 이유).
    # 빈 폴더 최초 실행에서 no such table: drug_files 로 터지며,
    # search_index_status() 도 함께 실패했습니다.
    init_db()

    queries = expand_search_queries(q)
    seen_paths: set[str] = set()
    scored: list[tuple[float, dict]] = []

    with connect() as conn:
        for term in queries:
            rows = _search_index_once(term, conn)
            for row in rows:
                r = dict(row)
                path_key = r["file_path"]
                if path_key in seen_paths:
                    continue
                seen_paths.add(path_key)
                labels = [r.get("drug_name") or "", r["file_name"], Path(r["file_path"]).stem]
                score = max(match_score_for_labels(q, labels), match_score_for_labels(term, labels))
                if score >= 18.0:
                    scored.append((score, r))

        if not scored and count_indexed() < 5000:
            all_rows = conn.execute("SELECT * FROM drug_files").fetchall()
            for row in all_rows:
                r = dict(row)
                path_key = r["file_path"]
                if path_key in seen_paths:
                    continue
                score = max(
                    match_score_for_labels(q, [r.get("drug_name") or "", r["file_name"], Path(r["file_path"]).stem]),
                    *[
                        match_score_for_labels(
                            term, [r.get("drug_name") or "", r["file_name"], Path(r["file_path"]).stem]
                        )
                        for term in queries
                    ],
                )
                if score >= 22.0:
                    scored.append((score, r))

    scored.sort(key=lambda x: (-x[0], x[1].get("drug_name") or ""))
    out = []
    for score, rec in scored[:limit]:
        rec["match_score"] = score
        out.append(rec)
    return out


def get_by_path(path: Path) -> dict | None:
    with connect() as conn:
        row = conn.execute(
            "SELECT * FROM drug_files WHERE file_path = ?",
            (str(path.resolve()),),
        ).fetchone()
        return dict(row) if row else None
