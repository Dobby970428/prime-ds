"""drug_index.db + CSV 기반 SMILES 실측 데이터 조회."""

# "이 구조가 실제로 존재가 확인된 화합물인가" 확인 (실측/이론 판정)
# 📌 실측 검출 이력 있음 (Validated)
# 🧪 이론적 예측 (Theoretical)

# 판정 순서
# ① 지금 열려 있는 CSV에 있나?
# ② drug_index.db 에 있나?
# ③ (느린 검색 허용 시) PubChem에 물어보고 그 약의 CSV를 뒤진다
#    ↓
# 하나라도 걸리면 Validated, 아니면 Theoretical

from __future__ import annotations

from pathlib import Path

from p23_chem_search_engine import exact_match_search, safe_canonical_smiles  # noqa: F401
from p20_db_index import search_index
from p22_drug_data import load_related_table, resolve_data_path
from p25_pubchem_client import fetch_cid_by_smiles
from p30_reaction_engine import build_db_smiles_index, canonical_smiles


def lookup_smiles_in_csv(smiles: str, df) -> dict | None:
    """단일 CSV DataFrame에서 Canonical SMILES Exact Match."""
    canon, err = safe_canonical_smiles(smiles)
    if not canon:
        return None
    idx = build_db_smiles_index(df)
    rec = idx.get(canon)
    if not rec:
        return None
    return {
        "name": rec.get("name", ""),
        "cid": rec.get("cid"),
        "drug_file": None,
        "drug_name": None,
        "canon_smiles": canon,
    }


def lookup_smiles_in_index_db(smiles: str, max_files: int = 4) -> dict | None:
    """
    drug_index.db를 통해 SMILES가 실측 DB에 있는지 조회.
    PubChem CID → 인덱스 drug 파일 → CSV SMILES 스캔.
    """
    canon, _ = safe_canonical_smiles(smiles)
    if not canon:
        return None

    cid = fetch_cid_by_smiles(canon)
    queries: list[str] = []
    if cid:
        queries.extend([f"CID{cid}", str(cid)])
    queries.append(canon[:20])

    seen_paths: set[str] = set()
    for q in queries:
        if not q:
            continue
        for hit in search_index(q, limit=max_files * 2):
            path_key = hit.get("file_path", "")
            if path_key in seen_paths:
                continue
            seen_paths.add(path_key)
            try:
                p = resolve_data_path(Path(path_key), hit["file_name"])
                df = load_related_table(p)
                rec = lookup_smiles_in_csv(canon, df)
                if rec:
                    rec["drug_file"] = hit["file_name"]
                    rec["drug_name"] = hit.get("drug_name") or hit["file_name"]
                    return rec
            except Exception:
                continue
            if len(seen_paths) >= max_files:
                break
    return None


def _exact_to_record(exact: dict) -> dict:
    return {
        "name": exact.get("name") or "",
        "cid": exact.get("cid"),
        "drug_file": exact.get("drug_file"),
        "drug_name": exact.get("drug_name"),
        "canon_smiles": exact.get("canon_smiles"),
        "match_type": exact.get("match_type"),
    }


def classify_smiles_against_db(
    smiles: str,
    local_df=None,
    *,
    allow_slow_scan: bool = False,
) -> tuple[str, dict | None]:
    """
    Returns: (classification, db_record)
    classification: 'validated' | 'theoretical'

    기본은 Canonical Exact 인덱스·로컬 DF만 조회 (UI 응답 속도).
    allow_slow_scan=True 시 PubChem CID 기반 파일 순회(느림) 추가.
    """
    canon, err = safe_canonical_smiles(smiles)
    if err:
        return "theoretical", None

    if local_df is not None:
        rec = lookup_smiles_in_csv(canon, local_df)
        if rec:
            return "validated", rec

    exact = exact_match_search(canon)
    if exact:
        return "validated", _exact_to_record(exact)

    if allow_slow_scan:
        rec = lookup_smiles_in_index_db(canon, max_files=3)
        if rec:
            rec["canon_smiles"] = canon
            return "validated", rec
    return "theoretical", None


def resolve_local_df_from_structure(smiles: str):
    """
    그린 구조(SMILES) → drug_index.db Exact Match → 해당 의약품 CSV(DataFrame).
    Returns (local_df | None, db_record | None).
    """
    if not smiles or not str(smiles).strip():
        return None, None
    try:
        exact = exact_match_search(smiles)
        db_rec = _exact_to_record(exact) if exact else None

        drug_file = (exact or {}).get("drug_file") if exact else None
        if not drug_file:
            slow = lookup_smiles_in_index_db(smiles, max_files=2)
            if slow:
                db_rec = slow
                drug_file = slow.get("drug_file")

        if drug_file:
            for hit in search_index(str(drug_file), limit=5):
                try:
                    if hit.get("file_name") != drug_file and drug_file not in str(
                        hit.get("file_name", "")
                    ):
                        continue
                    p = resolve_data_path(Path(hit["file_path"]), hit["file_name"])
                    return load_related_table(p), db_rec
                except Exception:
                    continue

        if db_rec:
            return None, db_rec
        cls, rec = classify_smiles_against_db(smiles)
        return None, rec if cls == "validated" else None
    except Exception:
        return None, None
