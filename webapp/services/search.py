"""검색 — 의약품명 · 부분구조 · 유사도.

주의할 점 두 가지가 있었습니다.

1. Streamlit 판 p00_app.py는 expand_search_queries와 has_hangul을
   import하지 않은 채 사용합니다 (검색 폴백 경로와 결과 0건 분기에서
   NameError). 여기서는 제대로 import합니다.

2. p11_structure_search_history는 st.session_state를 직접 조작하므로
   Streamlit 밖에서는 쓸 수 없습니다. 이력은 webapp.state로 대체합니다.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Any

from .. import state

HISTORY_KEY = "_search_history"
HISTORY_MAX = 20


# ---------------------------------------------------------------------------
# 검색 이력 — p11을 대체합니다
# ---------------------------------------------------------------------------


def record_history(smiles: str, *, label: str = "독성 스캔") -> None:
    from . import chem

    canon = chem.canonical(smiles) or (smiles or "").strip()
    if not canon:
        return
    items = [i for i in state.get(HISTORY_KEY, []) if i.get("smiles") != canon]
    items.insert(
        0,
        {
            "smiles": canon,
            "label": label,
            "at": datetime.now().strftime("%Y-%m-%d %H:%M"),
        },
    )
    state.put(HISTORY_KEY, items[:HISTORY_MAX])


def history(limit: int = 12) -> list[dict[str, Any]]:
    return list(state.get(HISTORY_KEY, []))[:limit]


def clear_history() -> None:
    state.pop(HISTORY_KEY)


# ---------------------------------------------------------------------------
# 의약품명 검색
# ---------------------------------------------------------------------------


def _has_hangul(text: str) -> bool:
    try:
        from p21_drug_search_aliases import has_hangul

        return bool(has_hangul(text))
    except Exception:
        return any("가" <= ch <= "힣" for ch in text or "")


def drug_search(query: str, *, limit: int = 20) -> dict[str, Any]:
    out: dict[str, Any] = {
        "query": (query or "").strip(),
        "matches": [],
        "has_hangul": False,
        "error": None,
    }
    if not out["query"]:
        return out
    out["has_hangul"] = _has_hangul(out["query"])

    matches: list[dict[str, Any]] = []

    try:
        from p20_db_index import search_index

        for row in search_index(out["query"], limit=limit) or []:
            path = row.get("file_path") or ""
            matches.append(
                {
                    "path": str(path),
                    "file_name": row.get("file_name") or (Path(str(path)).name if path else ""),
                    "drug_name": row.get("drug_name") or "",
                    "score": row.get("match_score"),
                    "cid": row.get("compound_cid") or row.get("cid"),
                }
            )
    except Exception as exc:
        out["error"] = f"색인 검색: {exc}"

    # 색인에서 못 찾으면 별칭 확장 + 파일명 매칭으로 폴백
    if not matches:
        try:
            from p21_drug_search_aliases import expand_search_queries
            from p22_drug_data import get_data_files_cached, match_drug_files

            files = get_data_files_cached()
            for term in expand_search_queries(out["query"]) or [out["query"]]:
                for item in (match_drug_files(term, files) or [])[:limit]:
                    path = item[0]
                    score = item[1] if len(item) > 1 else None
                    hint = item[2] if len(item) > 2 else ""
                    matches.append(
                        {
                            "path": str(path),
                            "file_name": Path(str(path)).name,
                            "drug_name": hint or Path(str(path)).stem,
                            "score": score,
                            "cid": None,
                        }
                    )
                if matches:
                    break
        except Exception as exc:
            out["error"] = (out["error"] or "") + f" · 파일명 검색: {exc}"

    seen: set[str] = set()
    unique = []
    for m in matches:
        if not m["path"] or m["path"] in seen:
            continue
        seen.add(m["path"])
        unique.append(m)
    out["matches"] = unique[:limit]
    return out


def drug_table(file_path: str, file_name: str = "", *, max_show: int = 12) -> dict[str, Any]:
    """선택한 의약품 파일의 관련 화합물 표."""
    out: dict[str, Any] = {
        "rows": [],
        "path": file_path,
        "file_name": file_name,
        "total": 0,
        "error": None,
    }
    try:
        from p22_drug_data import load_related_table, resolve_data_path

        stored = Path(file_path)
        resolved = resolve_data_path(stored, file_name or stored.name)
        frame = load_related_table(resolved)
        out["total"] = int(len(frame))
    except Exception as exc:
        out["error"] = str(exc)
        return out

    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen

        from p62_toxicity_alerts import scan_ashby_tennant_alerts
        from p66_ich_m7 import classify_ich_m7

        columns = {c.lower(): c for c in frame.columns}
        smi_col = columns.get("smiles", "SMILES")
        name_col = columns.get("name") or columns.get("cmpdname")
        cid_col = columns.get("cid")

        for _, record in frame.head(max(1, int(max_show))).iterrows():
            smiles = str(record.get(smi_col) or "").strip()
            mol = Chem.MolFromSmiles(smiles) if smiles else None

            labels: list[str] = []
            atoms: list[int] = []
            m7_class = None
            logp = None
            if mol is not None:
                try:
                    labels, atoms = scan_ashby_tennant_alerts(mol)
                except Exception:
                    labels, atoms = [], []
                try:
                    m7_class = getattr(classify_ich_m7(mol, labels, None), "class_number", None)
                except Exception:
                    m7_class = None
                try:
                    logp = round(float(Crippen.MolLogP(mol)), 2)
                except Exception:
                    logp = None

            out["rows"].append(
                {
                    "name": str(record.get(name_col) or "") if name_col else "",
                    "cid": record.get(cid_col) if cid_col else None,
                    "smiles": smiles,
                    "m7_class": m7_class,
                    "logp": logp,
                    "alerts": list(labels or []),
                    "highlight_atoms": list(atoms or []),
                }
            )
    except Exception as exc:
        out["error"] = f"화합물 평가: {exc}"

    return out


# ---------------------------------------------------------------------------
# DB 검색 리포트 (p89)
# ---------------------------------------------------------------------------


def _drug_label(table: dict[str, Any], query: str = "") -> str:
    """파일명에서 의약품 이름을 뽑아낸다. 실패하면 검색어를 쓴다."""
    name = str(table.get("file_name") or "")
    stem = Path(name).stem if name else ""
    stem = re.sub(r"[_-]?(related|impurit\w*|compounds?|list|data)$", "", stem, flags=re.I)
    stem = stem.replace("_", " ").strip()
    return stem or (query or "").strip() or "Unknown drug"


def _strip_emoji(text: str) -> str:
    """PDF 폰트에 없는 그림문자를 걷어낸다 (엔진 라벨에 ✅ ⚠ 등이 섞여 있다)."""
    out = "".join(ch for ch in str(text or "") if ord(ch) < 0x2190 or 0xAC00 <= ord(ch) <= 0xD7A3)
    return " ".join(out.split()) or "—"


def drug_db_report_context(table: dict[str, Any], *, query: str = ""):
    """
    화면 표를 p89 리포트 컨텍스트로 옮긴다.

    화면 쪽은 소문자 키(name · smiles · m7_class …)를 쓰고
    p89는 표시용 열 이름(Name · SMILES · ICH M7 …)을 기대하므로 여기서 맞춘다.
    """
    from p89_drug_db_report_export import build_drug_db_report_context

    try:
        from p50_physicochemical import stability_summary_for_table
    except Exception:
        stability_summary_for_table = None

    rows: list[dict[str, Any]] = []
    for record in table.get("rows") or []:
        smiles = str(record.get("smiles") or "").strip()

        stability = "—"
        if stability_summary_for_table is not None and smiles:
            try:
                stability = _strip_emoji(stability_summary_for_table(smiles)[0])
            except Exception:
                stability = "—"

        alerts = [str(a) for a in (record.get("alerts") or []) if str(a).strip()]
        m7 = record.get("m7_class")
        logp = record.get("logp")

        rows.append(
            {
                "Name": str(record.get("name") or ""),
                "SMILES": smiles,
                "CID": record.get("cid"),
                "LogP": "—" if logp is None else logp,
                "ICH M7": f"Class {m7}" if m7 else "—",
                "Ashby-Tennant": ", ".join(alerts) if alerts else "—",
                "공정 안정성": stability,
            }
        )

    return build_drug_db_report_context(
        drug_label=_drug_label(table, query),
        search_query=query,
        source_file=str(table.get("path") or ""),
        total_records=int(table.get("total") or 0),
        rows=rows,
    )


def drug_db_report_text(table: dict[str, Any], *, query: str = "") -> str:
    from p89_drug_db_report_export import generate_drug_db_report_text

    return generate_drug_db_report_text(drug_db_report_context(table, query=query))


# ---------------------------------------------------------------------------
# 구조 검색
# ---------------------------------------------------------------------------


def substructure(query_smiles: str, *, mode: str = "fast", max_results: int = 500) -> dict[str, Any]:
    out: dict[str, Any] = {
        "query": query_smiles,
        "canonical": "",
        "mode": mode,
        "hits": [],
        "error": None,
    }
    try:
        from p23_chem_search_engine import safe_canonical_smiles, substructure_search

        canon, err = safe_canonical_smiles(query_smiles)
        if err or not canon:
            out["error"] = err or "구조를 해석하지 못했습니다"
            return out
        out["canonical"] = canon

        hits, search_err = substructure_search(canon, max_results=max_results, mode=mode)
        out["hits"] = list(hits or [])
        if search_err:
            out["error"] = search_err
    except Exception as exc:
        out["error"] = str(exc)
    return out


def similarity(query_smiles: str, *, top_n: int = 3) -> dict[str, Any]:
    out: dict[str, Any] = {"query": query_smiles, "canonical": "", "hits": [], "error": None}
    try:
        from p23_chem_search_engine import safe_canonical_smiles, similarity_search_top_n

        canon, err = safe_canonical_smiles(query_smiles)
        if err or not canon:
            out["error"] = err or "구조를 해석하지 못했습니다"
            return out
        out["canonical"] = canon

        hits, search_err = similarity_search_top_n(canon, top_n=top_n)
        out["hits"] = list(hits or [])
        if search_err:
            out["error"] = search_err
    except Exception as exc:
        out["error"] = str(exc)
    return out
