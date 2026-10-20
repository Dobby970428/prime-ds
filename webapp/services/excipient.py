"""부형제-API 상호작용.

적용되는 반응은 세 갈래뿐입니다 — Maillard 이민, 포르밀화, 니트로소화.
아질산염 프로파일이 있는 부형제 7종만 니트로소화 경로를 탑니다.
"""

from __future__ import annotations

from typing import Any


def _msg(key: str) -> str:
    from ..i18n import t

    return t(key)

SKIP_LABEL = "선택 안함"

_OTHER = "Other"


def catalog(lang: str = "en") -> dict[str, Any]:
    """부형제 목록과 카테고리.

    영어 화면에서는 엔진의 한글 label_ko 대신 webapp/excipient_terms.py의
    영문 표시 이름을 씁니다. 엔진 데이터는 건드리지 않습니다.
    """
    out: dict[str, Any] = {
        "choices": [],
        "categories": [],
        "labels": {},
        "rows": [],
        "stats": {},
        "error": None,
    }
    try:
        from p72_excipient_interaction_engine import (
            EXCIPIENT_CATEGORIES,
            EXCIPIENT_CHOICES,
            EXCIPIENT_LIBRARY,
            excipient_catalog_rows,
            excipient_labels_for_ui,
        )

        out["choices"] = list(EXCIPIENT_CHOICES or [])
        out["categories"] = list(EXCIPIENT_CATEGORIES or [])
        rows = list(excipient_catalog_rows() or [])

        if lang == "en":
            from ..excipient_terms import alert_en, column_en, name_en

            out["labels"] = {
                key: f"{name_en(key)} [{(EXCIPIENT_LIBRARY.get(key) or {}).get('category', '')}]"
                for key in out["choices"]
            }
            out["rows"] = [
                {
                    column_en(col): (
                        name_en(row.get("key", ""), value)
                        if col in ("부형제", "부형제키")
                        else alert_en(row.get("key", ""), value)
                        if col == "주의보"
                        else value
                    )
                    for col, value in row.items()
                }
                for row in rows
            ]
        else:
            out["labels"] = dict(excipient_labels_for_ui() or {})
            out["rows"] = rows
    except Exception as exc:
        out["error"] = str(exc)

    try:
        from p71_excipient_opensource_registry import opensource_catalog_stats

        out["stats"] = dict(opensource_catalog_stats() or {})
    except Exception:
        out["stats"] = {}

    # 2026-09-30 수정: total 을 선택 가능한 항목 수로 덮어쓰면서 master·opensource 는
    #   레지스트리(76종) 기준 그대로 남아, 화면에 "catalogue 143 (master 47 · open
    #   source 29)" 처럼 47+29 가 143 의 내역인 것처럼 보였습니다. 두 집계가 서로 다른
    #   집합이므로 레지스트리 총계를 records 로 따로 보존해 화면에서 구분합니다.
    #     choices          = p72.EXCIPIENT_LIBRARY        — 선택 가능한 143종
    #     records          = p71.get_full_excipient_library() — 상세 레코드 76종
    out["stats"]["records"] = out["stats"].get("total", 0)
    out["stats"]["total"] = len(out["choices"])

    # 카테고리별 개수 (필터 칩에 표시)
    counts: dict[str, int] = {}
    try:
        from p72_excipient_interaction_engine import EXCIPIENT_LIBRARY

        for key in out["choices"]:
            cat_name = str((EXCIPIENT_LIBRARY.get(key) or {}).get("category") or "Other")
            counts[cat_name] = counts.get(cat_name, 0) + 1
    except Exception:
        pass
    out["category_counts"] = counts

    return out


def filter_keys(keys, *, categories=None, query: str = "", lang: str = "en") -> list[str]:
    """카테고리·검색어로 부형제 키를 걸러냅니다.

    ⚠️ `p72.filter_excipient_keys`는 **키워드 전용** 인자입니다
    (`def filter_excipient_keys(*, categories=None, search_query="")`).
    위치 인자로 부르면 TypeError가 나고, 예외를 삼키면 필터가 통째로
    무시됩니다. 실제로 카테고리 선택이 전혀 먹지 않고 있었습니다.

    엔진 쪽 검색은 `label_ko`만 보므로, 영어 화면에서 «macrogol» 같은
    영문명으로 찾을 수 있도록 여기서 한 번 더 걸러 줍니다.
    """
    cats = [c for c in (categories or []) if c and c not in ("전체", "All")]
    needle = (query or "").strip().lower()
    pool = list(keys or [])

    # ── 카테고리 ─────────────────────────────────────────────────────
    if cats:
        try:
            from p72_excipient_interaction_engine import EXCIPIENT_LIBRARY

            allowed = set(cats)
            pool = [
                k
                for k in pool
                if str((EXCIPIENT_LIBRARY.get(k) or {}).get("category") or "") in allowed
            ]
        except Exception:
            pass

    # ── 검색어 — 키 · 한글명 · 영문명 모두 대상 ──────────────────────
    if needle:
        try:
            from p72_excipient_interaction_engine import EXCIPIENT_LIBRARY
        except Exception:
            EXCIPIENT_LIBRARY = {}
        try:
            from ..excipient_terms import name_en
        except Exception:
            name_en = None  # type: ignore[assignment]

        matched = []
        for key in pool:
            entry = EXCIPIENT_LIBRARY.get(key) or {}
            parts = [
                key.replace("_", " "),
                str(entry.get("label_ko") or ""),
                str(entry.get("category") or ""),
            ]
            if name_en is not None:
                parts.append(name_en(key, ""))
            if needle in " ".join(parts).lower():
                matched.append(key)
        pool = matched

    return pool


def group_by_category(keys: list[str], catalog_data: dict[str, Any]) -> list[dict[str, Any]]:
    """카테고리별로 묶어 화면에 접어 보여주기 위한 구조.

    139종을 한 목록에 쏟으면 스크롤에 묻힙니다.
    """
    try:
        from p72_excipient_interaction_engine import EXCIPIENT_LIBRARY
    except Exception:
        EXCIPIENT_LIBRARY = {}

    order = list(catalog_data.get("categories") or [])
    rows_by_key = {r.get("key"): r for r in (catalog_data.get("rows") or [])}
    labels = catalog_data.get("labels") or {}

    buckets: dict[str, list[dict[str, Any]]] = {}
    for key in keys:
        category = str((EXCIPIENT_LIBRARY.get(key) or {}).get("category") or "Other")
        row = rows_by_key.get(key) or {}
        buckets.setdefault(category, []).append(
            {
                "key": key,
                "label": labels.get(key, key),
                "caution": row.get("Caution") or row.get("주의보") or "",
            }
        )

    grouped: list[dict[str, Any]] = []
    for category in order + [c for c in buckets if c not in order]:
        items = buckets.get(category)
        if not items:
            continue
        items.sort(key=lambda x: x["label"])
        grouped.append({"category": category, "items": items, "count": len(items)})
    return grouped


def run(api_smiles: str, excipient_keys: list[str], *, lang: str = "en") -> dict[str, Any]:
    """상호작용 예측을 실행합니다."""
    out: dict[str, Any] = {
        "mode": "matrix",
        "api_smiles": (api_smiles or "").strip(),
        "excipient_keys": list(excipient_keys or []),
        "hits": [],
        "amines": None,
        "counts": {"Matrix": 0, "NDSRI": 0, "NDMA": 0, _OTHER: 0},
        "errors": [],
        "_ctx": None,
    }

    if not out["api_smiles"]:
        out["errors"].append(_msg("flash.parse_failed"))
        return out
    if not out["excipient_keys"]:
        out["errors"].append(_msg("flash.select_excipient"))
        return out

    try:
        from p72_excipient_interaction_engine import AdvancedExcipientMatrixEngine

        engine = AdvancedExcipientMatrixEngine()
        hits = engine.evaluate_formulation_risk(out["api_smiles"], out["excipient_keys"])
        payload = [h.to_dict() if hasattr(h, "to_dict") else dict(h) for h in (hits or [])]
        out["hits"] = payload
    except Exception as exc:
        out["errors"].append(f"{_msg('excipient.results')}: {exc}")
        return out

    for hit in out["hits"]:
        _normalize_hit(hit)
        _localize_hit(hit, lang)
        risk = str(hit.get("risk_class") or _OTHER)
        out["counts"][risk] = out["counts"].get(risk, 0) + 1

    # 아민 자리 — 니트로소화 가능성의 전제
    try:
        from p30_reaction_engine import detect_amine_sites

        sites = detect_amine_sites(out["api_smiles"])
        out["amines"] = {
            "secondary": getattr(sites, "secondary_count", 0),
            "tertiary": getattr(sites, "tertiary_count", 0),
            "nitrosatable": bool(getattr(sites, "has_nitrosatable_amine", False)),
        }
    except Exception:
        out["amines"] = None

    # 보고서 컨텍스트 — 다운로드 라우트에서 씁니다
    try:
        from p72_excipient_interaction_engine import excipient_labels_for_ui
        from p88_excipient_report_export import build_excipient_report_context

        out["_ctx"] = build_excipient_report_context(
            out["api_smiles"],
            out["excipient_keys"],
            out["hits"],
            excipient_label_map=excipient_labels_for_ui() or {},
            amine_secondary=(out["amines"] or {}).get("secondary", 0),
            amine_tertiary=(out["amines"] or {}).get("tertiary", 0),
        )
    except Exception as exc:
        out["errors"].append(f"{_msg('reports.regulatory')}: {exc}")

    return out


def run_binary(
    api_smiles: str,
    excipient_keys: list[str],
    *,
    ratio: tuple[float, float] = (1.0, 1.0),
    lang: str = "en",
) -> dict[str, Any]:
    """1:1 (w/w) 이원 배합적합성 — 부형제를 한 쌍씩 평가합니다.

    매트릭스 모드와 달리 부형제마다 엔진을 새로 만들기 때문에, 같은
    불순물을 만드는 부형제가 여럿이어도 전부 잡힙니다.
    """
    out: dict[str, Any] = {
        "mode": "binary",
        "api_smiles": (api_smiles or "").strip(),
        "excipient_keys": list(excipient_keys or []),
        "ratio": list(ratio),
        "pairs": [],
        "heatmap": None,
        "summary": {},
        "amines": None,
        "errors": [],
        "_ctx": None,
        "_all_hits": [],
    }

    if not out["api_smiles"]:
        out["errors"].append(_msg("flash.parse_failed"))
        return out
    if not out["excipient_keys"]:
        out["errors"].append(_msg("flash.select_excipient"))
        return out

    try:
        from p17_binary_compatibility import evaluate_binary_compatibility, heatmap
    except Exception as exc:
        out["errors"].append(f"Binary compatibility engine unavailable: {exc}")
        return out

    result = evaluate_binary_compatibility(
        out["api_smiles"], out["excipient_keys"], ratio=ratio
    )
    out["errors"].extend(result.errors)
    out["api_mw"] = result.api_mw
    out["api_formula"] = result.api_formula

    payload = result.to_dict()
    out["pairs"] = payload["pairs"]
    hm = heatmap(result)
    if lang == "en":
        from ..excipient_terms import name_en

        for row in hm["rows"]:
            row["excipient_name"] = name_en(row["excipient_key"], row["excipient_name"])
    out["heatmap"] = hm

    all_hits: list[dict[str, Any]] = []
    for pair in out["pairs"]:
        for hit in pair["hits"]:
            _normalize_hit(hit)
            _localize_hit(hit, lang)
            all_hits.append(hit)
        if lang == "en":
            from ..excipient_terms import name_en

            pair["excipient_name"] = name_en(pair["excipient_key"], pair["excipient_name"])
    out["_all_hits"] = all_hits

    counts = {"Matrix": 0, "NDSRI": 0, "NDMA": 0, _OTHER: 0}
    for hit in all_hits:
        counts[str(hit.get("risk_class") or _OTHER)] = (
            counts.get(str(hit.get("risk_class") or _OTHER), 0) + 1
        )
    bands = {"high": 0, "moderate": 0, "low": 0, "none": 0}
    for pair in out["pairs"]:
        bands[pair["risk_band"]] = bands.get(pair["risk_band"], 0) + 1

    out["counts"] = counts
    out["summary"] = {
        "total_pairs": payload["total_pairs"],
        "interacting": payload["interacting_pairs"],
        "clean": payload["clean_pairs"],
        "total_hits": len(all_hits),
        "bands": bands,
        "worst_band": next(
            (b for b in ("high", "moderate", "low") if bands.get(b)), "none"
        ),
    }

    # 아민 자리 — 니트로소화 가능성의 전제
    try:
        from p30_reaction_engine import detect_amine_sites

        sites = detect_amine_sites(out["api_smiles"])
        out["amines"] = {
            "secondary": getattr(sites, "secondary_count", 0),
            "tertiary": getattr(sites, "tertiary_count", 0),
            "nitrosatable": bool(getattr(sites, "has_nitrosatable_amine", False)),
        }
    except Exception:
        out["amines"] = None

    # 보고서 컨텍스트 — 전체 hit을 모아 넘깁니다
    try:
        from p72_excipient_interaction_engine import excipient_labels_for_ui
        from p88_excipient_report_export import build_excipient_report_context

        out["_ctx"] = build_excipient_report_context(
            out["api_smiles"],
            out["excipient_keys"],
            all_hits,
            excipient_label_map=excipient_labels_for_ui() or {},
            amine_secondary=(out["amines"] or {}).get("secondary", 0),
            amine_tertiary=(out["amines"] or {}).get("tertiary", 0),
        )
    except Exception as exc:
        out["errors"].append(f"{_msg('reports.regulatory')}: {exc}")

    return out


def _as_bool(value: Any) -> bool | None:
    """엔진이 bool을 문자열 'True'/'False'로 돌려주는 경우가 있습니다."""
    if isinstance(value, bool):
        return value
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in ("true", "1", "yes", "y"):
        return True
    if text in ("false", "0", "no", "n"):
        return False
    return None


def _localize_hit(hit: dict[str, Any], lang: str) -> None:
    """영어 화면에서 부형제 이름을 영문으로 바꿉니다."""
    if lang != "en":
        return
    from ..excipient_terms import name_en

    raw = str(hit.get("excipient_name") or "").strip()
    if raw:
        hit["excipient_name"] = name_en(raw, raw)


def _normalize_hit(hit: dict[str, Any]) -> None:
    """화면에서 쓰기 쉽게 자료형을 다듬습니다.

    ich_m7_class_hint는 'Class 1-3 (alpha-H present)' 같은 문자열이라
    숫자 배지로 바로 못 씁니다. 표시용 텍스트와 대표 등급을 분리합니다.
    """
    hit["alpha_hydrogen_present"] = _as_bool(hit.get("alpha_hydrogen_present"))
    hit["steric_hindrance"] = _as_bool(hit.get("steric_hindrance"))

    raw = str(hit.get("ich_m7_class_hint") or "").strip()
    hit["m7_hint_text"] = raw
    numbers = [int(ch) for ch in raw if ch.isdigit() and ch != "0"]
    hit["m7_hint_class"] = min(numbers) if numbers else None
    hit["m7_hint_tone"] = (
        "danger" if hit["m7_hint_class"] in (1, 2)
        else "warn" if hit["m7_hint_class"] == 3
        else "info" if hit["m7_hint_class"] == 4
        else "ok" if hit["m7_hint_class"] == 5
        else "muted"
    )


def nitrite_profiles(lang: str = "en") -> list[dict[str, Any]]:
    """아질산염 잔류 프로파일이 등록된 부형제 (7종)."""
    try:
        # 2026-08-26 — p73이 p72로 흡수되면서 이리로 옮겨 왔습니다.
        from p72_excipient_interaction_engine import nitrite_matrix_catalog_rows

        rows = list(nitrite_matrix_catalog_rows() or [])
    except Exception:
        return []

    if lang != "en":
        return rows

    from ..excipient_terms import column_en, name_en, nitrite_note_en

    out = []
    for row in rows:
        key = str(row.get("부형제키") or row.get("key") or "")
        out.append(
            {
                column_en(col): (
                    name_en(key, value)
                    if col in ("부형제키", "부형제")
                    else nitrite_note_en(key, value)
                    if col == "가이드"
                    else value
                )
                for col, value in row.items()
            }
        )
    return out
