"""
ICH M7 empirical evidence index — Hansen / ISSCAN / Tox21 seed records.

SQLite-backed canonical-SMILES lookup (InChIKey primary key) for sub-10 ms
consensus classification. Built automatically from data/m7_empirical_seed.csv.
"""

from __future__ import annotations

import csv
import sqlite3
import time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from rdkit import Chem

_APP_DIR = Path(__file__).resolve().parent
_SEED_CSV = _APP_DIR / "data" / "m7_empirical_seed.csv"
_DB_PATH = _APP_DIR / "cache" / "m7_empirical.db"

_VALID_OUTCOMES = frozenset({"positive", "negative", "unknown"})


@dataclass(frozen=True)
class EmpiricalM7Hit:
    """Aggregated empirical outcome for one canonical structure."""

    matched: bool = False
    inchikey: str = ""
    canonical_smiles: str = ""
    compound_name: str = ""
    sources: tuple[str, ...] = ()
    ames_result: str = "unknown"
    carcinogenicity_result: str = "unknown"
    match_type: str = "none"
    latency_ms: float = 0.0

    @property
    def source_label(self) -> str:
        if not self.sources:
            return "None"
        return "/".join(sorted(set(self.sources)))

    def to_report_dict(self) -> dict[str, str]:
        # 2026-09-15 수정: 라벨이 «(Hansen/ISSCAN)» 으로 고정돼 있어, 실제로는
        #   TDC 가 컴파일한 Ames 학습 분할에서 나온 적중까지 1차 실험 기록처럼
        #   보였습니다. DB 5,324건 중 5,274건이 TDC-AMES 이고 1차 기록은 50건
        #   (Hansen 29 · Tox21 12 · ISSCAN 9) 뿐입니다. s2 가 적재 시 같은
        #   취지로 경고합니다.
        #
        #   고정 라벨을 떼고 source 문자열을 그대로 노출합니다. 어디서 온
        #   적중인지 보고서에서 바로 읽힙니다.
        return {
            "Empirical DB Source Match": self.source_label,
            "Ames Evidence": self.ames_result,
            "Carcinogenicity Evidence": self.carcinogenicity_result,
            "Empirical Match Type": self.match_type,
            "Empirical Compound Name": self.compound_name or "N/A",
        }


def _canonicalize(smiles: str) -> tuple[str, str, str | None]:
    """Return (canonical_smiles, inchikey, error)."""
    smi = (smiles or "").strip()
    if not smi:
        return "", "", "Empty SMILES"
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return "", "", "Invalid SMILES"
        can = Chem.MolToSmiles(mol, canonical=True)
        key = Chem.MolToInchiKey(mol)
        return can, key, None
    except Exception as exc:
        return "", "", str(exc)


def _merge_outcome(current: str, new: str) -> str:
    """Merge evidence: positive (worst) > negative (confirmed safe) > unknown."""
    rank = {"positive": 3, "negative": 2, "unknown": 1}
    cur = current if current in _VALID_OUTCOMES else "unknown"
    nxt = new if new in _VALID_OUTCOMES else "unknown"
    return cur if rank[cur] >= rank[nxt] else nxt


def _ensure_db() -> Path:
    _DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    if _DB_PATH.is_file() and _DB_PATH.stat().st_size > 0:
        return _DB_PATH
    build_m7_empirical_db(force=True)
    return _DB_PATH


def build_m7_empirical_db(*, force: bool = False) -> int:
    """Build SQLite index from seed CSV. Returns row count."""
    if _DB_PATH.is_file() and not force:
        return _count_rows()

    if not _SEED_CSV.is_file():
        raise FileNotFoundError(f"Seed file missing: {_SEED_CSV}")

    _DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    if _DB_PATH.is_file():
        _DB_PATH.unlink()

    conn = sqlite3.connect(_DB_PATH)
    try:
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.executescript(
            """
            CREATE TABLE m7_empirical (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                inchikey TEXT NOT NULL,
                canonical_smiles TEXT NOT NULL,
                compound_name TEXT,
                source TEXT NOT NULL,
                ames_result TEXT NOT NULL,
                carcinogenicity_result TEXT NOT NULL
            );
            CREATE INDEX idx_m7_inchikey ON m7_empirical(inchikey);
            CREATE INDEX idx_m7_canonical ON m7_empirical(canonical_smiles);
            CREATE VIEW empirical_toxicology_db AS
            SELECT
                canonical_smiles AS smiles,
                ames_result AS ames,
                carcinogenicity_result AS carc,
                source,
                compound_name,
                inchikey
            FROM m7_empirical;
            """
        )

        inserted = 0
        with _SEED_CSV.open(encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                raw_smi = (row.get("canonical_smiles") or "").strip()
                can, key, err = _canonicalize(raw_smi)
                if err or not can:
                    continue
                ames = (row.get("ames_result") or "unknown").strip().lower()
                carc = (row.get("carcinogenicity_result") or "unknown").strip().lower()
                if ames not in _VALID_OUTCOMES:
                    ames = "unknown"
                if carc not in _VALID_OUTCOMES:
                    carc = "unknown"
                conn.execute(
                    """
                    INSERT INTO m7_empirical
                    (inchikey, canonical_smiles, compound_name, source, ames_result, carcinogenicity_result)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (
                        key,
                        can,
                        (row.get("compound_name") or "").strip(),
                        (row.get("source") or "Hansen").strip(),
                        ames,
                        carc,
                    ),
                )
                inserted += 1
        conn.commit()
        inserted = conn.execute("SELECT COUNT(*) FROM m7_empirical").fetchone()[0]
        return int(inserted)
    finally:
        conn.close()
        _db_connection.cache_clear()


def _count_rows() -> int:
    conn = sqlite3.connect(_DB_PATH)
    try:
        return int(conn.execute("SELECT COUNT(*) FROM m7_empirical").fetchone()[0])
    finally:
        conn.close()


@lru_cache(maxsize=1)
def _db_connection() -> sqlite3.Connection:
    _ensure_db()
    conn = sqlite3.connect(_DB_PATH, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA query_only=ON")
    return conn


def lookup_empirical_m7(smiles: str) -> EmpiricalM7Hit:
    """
    Fast empirical lookup by canonical SMILES / InChIKey.
    Aggregates Hansen, ISSCAN, and Tox21 rows for the same structure.
    """
    t0 = time.perf_counter()
    can, key, err = _canonicalize(smiles)
    if err or not can:
        return EmpiricalM7Hit(
            matched=False,
            match_type="invalid_smiles",
            latency_ms=(time.perf_counter() - t0) * 1000,
        )

    try:
        conn = _db_connection()
        rows = conn.execute(
            """
            SELECT compound_name, source, ames_result, carcinogenicity_result
            FROM m7_empirical
            WHERE inchikey = ? OR canonical_smiles = ?
            """,
            (key, can),
        ).fetchall()
    except Exception:
        build_m7_empirical_db(force=True)
        _db_connection.cache_clear()
        return lookup_empirical_m7(smiles)

    if not rows:
        return EmpiricalM7Hit(
            matched=False,
            inchikey=key,
            canonical_smiles=can,
            match_type="none",
            latency_ms=(time.perf_counter() - t0) * 1000,
        )

    return _rows_to_hit(can, key, rows, (time.perf_counter() - t0) * 1000)


def get_default_db_connection() -> sqlite3.Connection:
    """Shared read-only SQLite handle for evaluate_ich_m7_class()."""
    return _db_connection()


def query_empirical_toxicology_db(
    db_connection: sqlite3.Connection,
    canonical_smiles: str,
) -> tuple[bool, str | None, str | None, dict[str, str]]:
    """
    Low-level lookup matching evaluate_ich_m7_class pseudocode:
      SELECT ames, carc FROM empirical_toxicology_db WHERE smiles=?

    Returns (is_in_db, ames_result, carc_result, metadata).
    """
    can = (canonical_smiles or "").strip()
    if not can:
        return False, None, None, {}

    db_connection.row_factory = sqlite3.Row
    rows = db_connection.execute(
        """
        SELECT ames, carc, source, compound_name
        FROM empirical_toxicology_db
        WHERE smiles = ?
        """,
        (can,),
    ).fetchall()

    if not rows:
        return False, None, None, {"canonical_smiles": can}

    ames = "unknown"
    carc = "unknown"
    sources: list[str] = []
    names: list[str] = []
    for row in rows:
        sources.append(str(row["source"]))
        if row["compound_name"]:
            names.append(str(row["compound_name"]))
        ames = _merge_outcome(ames, str(row["ames"]))
        carc = _merge_outcome(carc, str(row["carc"]))

    meta = {
        "canonical_smiles": can,
        "source": "/".join(sorted(set(sources))),
        "compound_name": names[0] if names else "",
    }
    return True, ames, carc, meta


def _rows_to_hit(can: str, key: str, rows: list, latency_ms: float) -> EmpiricalM7Hit:
    ames = "unknown"
    carc = "unknown"
    sources: list[str] = []
    names: list[str] = []
    for row in rows:
        sources.append(str(row["source"]))
        if row["compound_name"]:
            names.append(str(row["compound_name"]))
        ames = _merge_outcome(ames, str(row["ames_result"]))
        carc = _merge_outcome(carc, str(row["carcinogenicity_result"]))

    return EmpiricalM7Hit(
        matched=True,
        inchikey=key,
        canonical_smiles=can,
        compound_name=names[0] if names else "",
        sources=tuple(sources),
        ames_result=ames,
        carcinogenicity_result=carc,
        match_type="exact",
        latency_ms=latency_ms,
    )


def empirical_db_stats() -> dict[str, int | str]:
    _ensure_db()
    conn = _db_connection()
    total = conn.execute("SELECT COUNT(*) FROM m7_empirical").fetchone()[0]
    compounds = conn.execute("SELECT COUNT(DISTINCT inchikey) FROM m7_empirical").fetchone()[0]
    return {
        "records": int(total),
        "unique_structures": int(compounds),
        "db_path": str(_DB_PATH),
    }
