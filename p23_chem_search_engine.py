"""
RDKit 기반 화합물 검색 엔진 — Exact / Substructure / Tanimoto Similarity.

대용량 CSV(2천만+ 행)를 고려해:
  - drug_primary: 의약품 파일당 대표 화합물 + Morgan FP (유사도 Top-N, 빠름)
  - compounds: Canonical SMILES 역인덱스 (Exact Match, 선택적 전체 구축)
  - Substructure: SubstructLibrary (대표 화합물 고속) · 선택적 CSV 심층 스캔
"""

# 구조로 찾기
# ① 정확일치: canonical SMILES 대조-"이 구조가 DB에 있나?" (즉시)
# ② 부분구조: RDKit SubstractLibrary-"이 조각을 포함한 약은?" (밀리초)
# ③ 유사도: Morgan 지문 Tanimoto-"가장 비슷한 약 3개는?" (초 단위)

# 유사도 검색 엔진
# p42가 분해산물을 만든다
#    ↓
# "이 구조, DB에 비슷한 게 있나?" ← p23에게 물어봄
#    ↓
# Tanimoto 87% → 신뢰도 87점

from __future__ import annotations

import json
import os
import pickle
import sqlite3
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdSubstructLibrary

from p20_db_index import DB_PATH, connect, count_indexed, init_db
from p22_drug_data import load_related_table, resolve_data_path
from p30_reaction_engine import canonical_smiles, _mol_from_smiles

APP_DIR = Path(__file__).resolve().parent
CACHE_DIR = APP_DIR / "cache"
SUBSTRUCT_LIB_PATH = CACHE_DIR / "substruct_primary.pkl"
SUBSTRUCT_MANIFEST_PATH = CACHE_DIR / "substruct_manifest.json"
MORGAN_RADIUS = 2
MORGAN_BITS = 2048
ProgressFn = Callable[[str, float], None] | None

_LIB_CACHE: tuple[rdSubstructLibrary.SubstructLibrary, list[dict]] | None = None


def safe_canonical_smiles(smiles: str) -> tuple[str | None, str | None]:
    """
    SMILES → Canonical SMILES (RDKit MolFromSmiles / MolToSmiles).
    Returns (canonical, error_message).
    """
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return None, "SMILES가 비어 있습니다."
    raw = smiles.strip()
    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        return None, f"RDKit이 SMILES를 해석하지 못했습니다: {raw[:80]}"
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        return None, f"구조 정제(Sanitize) 실패: {exc}"
    try:
        canon = Chem.MolToSmiles(mol, canonical=True)
    except Exception as exc:
        return None, f"Canonical SMILES 변환 실패: {exc}"
    if not canon:
        return None, "Canonical SMILES 결과가 비어 있습니다."
    return canon, None


def mol_from_smiles_safe(smiles: str) -> tuple[Chem.Mol | None, str | None]:
    canon, err = safe_canonical_smiles(smiles)
    if err:
        return None, err
    mol = _mol_from_smiles(canon)
    if mol is None:
        return None, "분자 객체 생성 실패"
    return mol, None


def fp_to_bytes(fp: DataStructs.ExplicitBitVect) -> bytes:
    return DataStructs.BitVectToFPSText(fp).encode("utf-8")


def bytes_to_fp(data: bytes) -> DataStructs.ExplicitBitVect:
    return DataStructs.CreateFromFPSText(data.decode("utf-8"))


def morgan_fp(mol: Chem.Mol) -> DataStructs.ExplicitBitVect:
    return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_BITS)


def init_chem_search_tables() -> None:
    init_db()
    with connect() as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS drug_primary (
                file_path TEXT PRIMARY KEY,
                file_name TEXT NOT NULL,
                drug_name TEXT,
                canon_smiles TEXT NOT NULL,
                compound_name TEXT,
                compound_cid TEXT,
                fp_blob BLOB NOT NULL
            )
            """
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_drug_primary_canon ON drug_primary(canon_smiles)"
        )
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS compounds (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                canon_smiles TEXT NOT NULL,
                name TEXT,
                compound_cid TEXT,
                drug_file TEXT NOT NULL,
                drug_name TEXT,
                UNIQUE(canon_smiles, drug_file)
            )
            """
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_compounds_canon ON compounds(canon_smiles)"
        )


def count_drug_primary() -> int:
    init_chem_search_tables()
    with connect() as conn:
        row = conn.execute("SELECT COUNT(*) AS c FROM drug_primary").fetchone()
        return int(row["c"]) if row else 0


def count_compounds() -> int:
    init_chem_search_tables()
    with connect() as conn:
        row = conn.execute("SELECT COUNT(*) AS c FROM compounds").fetchone()
        return int(row["c"]) if row else 0


def _first_valid_row(path: Path) -> dict | None:
    try:
        df = pd.read_csv(path, usecols=lambda c: c in ("SMILES", "Name", "Compound_CID"), nrows=50)
    except Exception:
        return None
    if "SMILES" not in df.columns:
        return None
    for _, row in df.iterrows():
        smi = row.get("SMILES")
        if not isinstance(smi, str) or not smi.strip():
            continue
        canon, _ = safe_canonical_smiles(smi)
        if not canon:
            continue
        return {
            "canon_smiles": canon,
            "compound_name": str(row.get("Name", "") or ""),
            "compound_cid": str(row.get("Compound_CID", "") or ""),
        }
    return None


def build_drug_primary_index(progress: ProgressFn = None) -> tuple[int, int]:
    """의약품 파일당 대표 화합물 1건 + Morgan FP 구축 (유사도 검색용)."""
    init_chem_search_tables()
    with connect() as conn:
        files = conn.execute(
            "SELECT file_path, file_name, drug_name FROM drug_files ORDER BY file_name"
        ).fetchall()
    total = len(files)
    ok, skip = 0, 0
    with connect() as conn:
        for i, rec in enumerate(files, start=1):
            path = resolve_data_path(Path(rec["file_path"]), rec["file_name"])
            if progress and (i == 1 or i % 100 == 0 or i == total):
                progress(f"대표 화합물 인덱싱 {i}/{total}", i / max(total, 1))
            if not path.exists():
                skip += 1
                continue
            meta = _first_valid_row(path)
            if not meta:
                skip += 1
                continue
            mol = _mol_from_smiles(meta["canon_smiles"])
            if mol is None:
                skip += 1
                continue
            conn.execute(
                """
                INSERT INTO drug_primary (
                    file_path, file_name, drug_name, canon_smiles,
                    compound_name, compound_cid, fp_blob
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(file_path) DO UPDATE SET
                    file_name=excluded.file_name,
                    drug_name=excluded.drug_name,
                    canon_smiles=excluded.canon_smiles,
                    compound_name=excluded.compound_name,
                    compound_cid=excluded.compound_cid,
                    fp_blob=excluded.fp_blob
                """,
                (
                    str(path.resolve()),
                    rec["file_name"],
                    rec["drug_name"],
                    meta["canon_smiles"],
                    meta["compound_name"],
                    meta["compound_cid"],
                    fp_to_bytes(morgan_fp(mol)),
                ),
            )
            ok += 1
            if i % 200 == 0:
                conn.commit()
        conn.commit()
    return ok, skip


def build_compound_index(progress: ProgressFn = None, max_files: int | None = None) -> tuple[int, int, int]:
    """
    전체 CSV Canonical SMILES 역인덱스 구축 (Exact Match 고속화).
    대용량 — UI에서 진행률 표시 권장.
    """
    init_chem_search_tables()
    with connect() as conn:
        files = conn.execute(
            "SELECT file_path, file_name, drug_name FROM drug_files ORDER BY file_name"
        ).fetchall()
    if max_files:
        files = files[:max_files]
    total = len(files)
    inserted, skipped_files, rows_seen = 0, 0, 0
    with connect() as conn:
        for i, rec in enumerate(files, start=1):
            path = resolve_data_path(Path(rec["file_path"]), rec["file_name"])
            if progress and (i == 1 or i % 20 == 0 or i == total):
                progress(f"화합물 역인덱스 {i}/{total} (누적 {inserted:,}건)", i / max(total, 1))
            if not path.exists():
                skipped_files += 1
                continue
            try:
                usecols = ["SMILES", "Name", "Compound_CID"]
                for chunk in pd.read_csv(path, usecols=lambda c: c in usecols, chunksize=8000):
                    if "SMILES" not in chunk.columns:
                        break
                    batch: list[tuple] = []
                    for _, row in chunk.iterrows():
                        rows_seen += 1
                        smi = row.get("SMILES")
                        if not isinstance(smi, str) or not smi.strip():
                            continue
                        canon, _ = safe_canonical_smiles(smi)
                        if not canon:
                            continue
                        batch.append(
                            (
                                canon,
                                str(row.get("Name", "") or ""),
                                str(row.get("Compound_CID", "") or ""),
                                rec["file_name"],
                                rec["drug_name"],
                            )
                        )
                    if batch:
                        before = conn.execute("SELECT COUNT(*) FROM compounds").fetchone()[0]
                        conn.executemany(
                            """
                            INSERT OR IGNORE INTO compounds
                            (canon_smiles, name, compound_cid, drug_file, drug_name)
                            VALUES (?, ?, ?, ?, ?)
                            """,
                            batch,
                        )
                        after = conn.execute("SELECT COUNT(*) FROM compounds").fetchone()[0]
                        inserted += max(0, after - before)
                    if i % 10 == 0:
                        conn.commit()
            except Exception:
                skipped_files += 1
                continue
        conn.commit()
    return inserted, skipped_files, rows_seen


def exact_match_search(smiles: str) -> dict | None:
    """Canonical SMILES Exact Match — compounds → drug_primary → None."""
    canon, err = safe_canonical_smiles(smiles)
    if not canon:
        return None
    init_chem_search_tables()
    with connect() as conn:
        row = conn.execute(
            """
            SELECT canon_smiles, name, compound_cid, drug_file, drug_name
            FROM compounds WHERE canon_smiles = ? LIMIT 1
            """,
            (canon,),
        ).fetchone()
        if row:
            return {
                "match_type": "exact_compound_index",
                "canon_smiles": row["canon_smiles"],
                "name": row["name"],
                "cid": row["compound_cid"],
                "drug_file": row["drug_file"],
                "drug_name": row["drug_name"],
            }
        row = conn.execute(
            """
            SELECT canon_smiles, compound_name, compound_cid, file_name, drug_name
            FROM drug_primary WHERE canon_smiles = ? LIMIT 1
            """,
            (canon,),
        ).fetchone()
        if row:
            return {
                "match_type": "exact_drug_primary",
                "canon_smiles": row["canon_smiles"],
                "name": row["compound_name"],
                "cid": row["compound_cid"],
                "drug_file": row["file_name"],
                "drug_name": row["drug_name"],
            }
    return None


def clear_substruct_library_cache() -> None:
    """하부구조 라이브러리 메모리·디스크 캐시 무효화."""
    global _LIB_CACHE
    _LIB_CACHE = None
    for path in (SUBSTRUCT_LIB_PATH, SUBSTRUCT_MANIFEST_PATH):
        try:
            path.unlink(missing_ok=True)
        except OSError:
            pass


def substruct_library_ready() -> bool:
    """디스크 캐시가 현재 drug_primary 건수와 일치하는지."""
    if not SUBSTRUCT_LIB_PATH.is_file() or not SUBSTRUCT_MANIFEST_PATH.is_file():
        return False
    try:
        manifest = json.loads(SUBSTRUCT_MANIFEST_PATH.read_text(encoding="utf-8"))
        return int(manifest.get("drug_primary_count", -1)) == count_drug_primary()
    except Exception:
        return False


def build_substruct_library(progress: ProgressFn = None) -> tuple[int, str | None]:
    """
    drug_primary 기준 RDKit SubstructLibrary 구축 (멀티스레드 검색용).
    최초 1회 약 5~15초 · 이후 디스크 캐시에서 즉시 로드.
    """
    global _LIB_CACHE
    init_chem_search_tables()
    primary_n = count_drug_primary()
    if primary_n < 100:
        return 0, (
            "대표 화합물 인덱스가 비어 있습니다. "
            "⚙️ 설정 → 「대표 화합물 FP 구축」을 먼저 실행하세요."
        )

    with connect() as conn:
        rows = conn.execute(
            """
            SELECT canon_smiles, compound_name, compound_cid, file_name, drug_name
            FROM drug_primary ORDER BY file_name
            """
        ).fetchall()

    mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    meta: list[dict] = []
    total = len(rows)
    for i, rec in enumerate(rows, start=1):
        smi = rec["canon_smiles"]
        if not smi:
            continue
        mols.AddSmiles(smi)
        meta.append(
            {
                "drug_name": rec["drug_name"] or rec["file_name"],
                "drug_file": rec["file_name"],
                "compound_name": str(rec["compound_name"] or ""),
                "compound_cid": str(rec["compound_cid"] or ""),
                "canon_smiles": smi,
            }
        )
        if progress and (i == 1 or i % 2000 == 0 or i == total):
            progress(f"하부구조 라이브러리 구축 {i}/{total}", i / max(total, 1))

    if not meta:
        return 0, "대표 화합물 SMILES를 읽지 못했습니다."

    library = rdSubstructLibrary.SubstructLibrary(mols)
    rdSubstructLibrary.AddPatterns(library)

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(SUBSTRUCT_LIB_PATH, "wb") as f:
        pickle.dump((library, meta), f, protocol=pickle.HIGHEST_PROTOCOL)
    SUBSTRUCT_MANIFEST_PATH.write_text(
        json.dumps({"drug_primary_count": primary_n, "entries": len(meta)}),
        encoding="utf-8",
    )
    _LIB_CACHE = (library, meta)
    return len(meta), None


def _load_substruct_library(
    progress: ProgressFn = None,
) -> tuple[rdSubstructLibrary.SubstructLibrary | None, list[dict], str | None]:
    global _LIB_CACHE
    if _LIB_CACHE is not None:
        return _LIB_CACHE[0], _LIB_CACHE[1], None

    if not substruct_library_ready():
        if progress:
            progress("하부구조 검색 인덱스 준비 중 (최초 1회)…", 0.05)
        n, err = build_substruct_library(progress=progress)
        if err:
            return None, [], err
        if progress:
            progress(f"인덱스 준비 완료 — {n:,}건", 0.2)

    try:
        with open(SUBSTRUCT_LIB_PATH, "rb") as f:
            library, meta = pickle.load(f)
    except Exception as exc:
        clear_substruct_library_cache()
        return None, [], f"하부구조 인덱스 로드 실패: {exc}"

    _LIB_CACHE = (library, meta)
    return library, meta, None


def _scan_file_substructure(args: tuple) -> list[dict]:
    path_str, file_name, drug_name, query_smiles, max_hits = args
    path = Path(path_str)
    if not path.is_file():
        return []

    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        return []

    hits: list[dict] = []
    try:
        usecols = ["SMILES", "Name", "Compound_CID"]
        for chunk in pd.read_csv(path, usecols=lambda c: c in usecols, chunksize=8000):
            if "SMILES" not in chunk.columns:
                break
            for _, row in chunk.iterrows():
                if len(hits) >= max_hits:
                    return hits
                smi = row.get("SMILES")
                if not isinstance(smi, str) or not smi.strip():
                    continue
                db_mol = Chem.MolFromSmiles(smi.strip())
                if db_mol is None:
                    continue
                try:
                    if not db_mol.HasSubstructMatch(query_mol):
                        continue
                except Exception:
                    continue
                canon, _ = safe_canonical_smiles(smi)
                hits.append(
                    {
                        "drug_name": drug_name or file_name,
                        "drug_file": file_name,
                        "compound_name": str(row.get("Name", "") or ""),
                        "compound_cid": str(row.get("Compound_CID", "") or ""),
                        "canon_smiles": canon or smi.strip(),
                        "query_smiles": query_smiles,
                    }
                )
    except Exception:
        pass
    return hits


def substructure_search_deep(
    query_smiles: str,
    *,
    max_results: int = 500,
    progress: ProgressFn = None,
    max_workers: int | None = None,
) -> tuple[list[dict], str | None]:
    """전체 CSV 행 심층 스캔 (병렬). 수 분 소요 가능."""
    canon, err = safe_canonical_smiles(query_smiles)
    if err or not canon:
        return [], err or "쿼리 SMILES 오류"

    init_chem_search_tables()
    with connect() as conn:
        files = conn.execute(
            "SELECT file_path, file_name, drug_name FROM drug_files ORDER BY row_count ASC"
        ).fetchall()

    workers = max_workers or min(8, max(2, (os.cpu_count() or 4) - 1))
    tasks: list[tuple] = []
    for rec in files:
        path = resolve_data_path(Path(rec["file_path"]), rec["file_name"])
        tasks.append(
            (
                str(path),
                rec["file_name"],
                rec["drug_name"],
                canon,
                max(1, max_results // max(len(files) // workers, 1)),
            )
        )

    results: list[dict] = []
    total = len(tasks)
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_scan_file_substructure, t) for t in tasks]
        for fut in as_completed(futures):
            if len(results) >= max_results:
                break
            done += 1
            if progress and (done == 1 or done % 200 == 0 or done == total):
                progress(
                    f"심층 스캔 {done}/{total} · 매칭 {len(results)}건",
                    done / max(total, 1),
                )
            try:
                for hit in fut.result():
                    if len(results) >= max_results:
                        break
                    results.append(hit)
            except Exception:
                continue
    return results, None


def substructure_search(
    query_smiles: str,
    *,
    max_results: int = 500,
    progress: ProgressFn = None,
    mode: str = "fast",
) -> tuple[list[dict], str | None]:
    """
    하부구조 검색.
    - fast (기본): SubstructLibrary + drug_primary (~1만 종, 수 ms~수 초)
    - deep: 전체 CSV 병렬 심층 스캔 (느림)
    """
    query_mol, err = mol_from_smiles_safe(query_smiles)
    if err or query_mol is None:
        return [], err or "쿼리 분자 생성 실패"

    if mode == "deep":
        return substructure_search_deep(
            Chem.MolToSmiles(query_mol, canonical=True),
            max_results=max_results,
            progress=progress,
        )

    library, meta, load_err = _load_substruct_library(progress=progress)
    if load_err or library is None:
        return [], load_err

    if progress:
        progress("멀티스레드 SubstructLibrary 검색…", 0.85)

    try:
        indices = library.GetMatches(
            query_mol,
            maxResults=max_results,
            numThreads=-1,
        )
    except Exception as exc:
        return [], f"하부구조 검색 실패: {exc}"

    query_canon = Chem.MolToSmiles(query_mol, canonical=True)
    results: list[dict] = []
    for idx in indices:
        if idx < 0 or idx >= len(meta):
            continue
        rec = dict(meta[idx])
        rec["query_smiles"] = query_canon
        results.append(rec)

    if progress:
        progress(f"완료 — {len(results)}건", 1.0)
    return results, None


def similarity_search_top_n(
    query_smiles: str,
    *,
    top_n: int = 3,
) -> tuple[list[dict], str | None]:
    """Morgan FP + Tanimoto — drug_primary(13k) 기준 Top-N."""
    query_mol, err = mol_from_smiles_safe(query_smiles)
    if err or query_mol is None:
        return [], err or "쿼리 분자 생성 실패"

    init_chem_search_tables()
    primary_n = count_drug_primary()
    if primary_n < max(100, count_indexed() // 20):
        return (
            [],
            "대표 화합물 FP 인덱스가 비어 있습니다. "
            "상단 ⚙️ 표시 설정 → 「대표 화합물 FP 구축」을 먼저 실행하세요.",
        )

    query_fp = morgan_fp(query_mol)
    with connect() as conn:
        rows = conn.execute(
            "SELECT file_name, drug_name, canon_smiles, compound_name, compound_cid, fp_blob FROM drug_primary"
        ).fetchall()

    if not rows:
        return [], "대표 화합물 인덱스가 비어 있습니다. 인덱스 구축을 실행하세요."

    fps: list[DataStructs.ExplicitBitVect] = []
    meta: list[dict] = []
    for row in rows:
        try:
            fps.append(bytes_to_fp(row["fp_blob"]))
            meta.append(dict(row))
        except Exception:
            continue

    if not fps:
        return [], "지문 데이터를 읽을 수 없습니다."

    sims = DataStructs.BulkTanimotoSimilarity(query_fp, fps)
    ranked = sorted(enumerate(sims), key=lambda x: -x[1])[:top_n]

    out: list[dict] = []
    for rank, (idx, score) in enumerate(ranked, start=1):
        rec = meta[idx]
        out.append(
            {
                "rank": rank,
                "tanimoto": float(score),
                "tanimoto_pct": round(float(score) * 100, 1),
                "drug_name": rec.get("drug_name") or rec.get("file_name"),
                "drug_file": rec.get("file_name"),
                "compound_name": rec.get("compound_name") or "",
                "compound_cid": rec.get("compound_cid") or "",
                "canon_smiles": rec.get("canon_smiles"),
            }
        )
    return out, None


def search_index_status() -> dict:
    return {
        "drug_files": count_indexed(),
        "drug_primary": count_drug_primary(),
        "compounds_indexed": count_compounds(),
        "substruct_ready": substruct_library_ready(),
        "db_path": str(DB_PATH),
    }
