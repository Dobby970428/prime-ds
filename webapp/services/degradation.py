"""분해산물 예측 · 분해 트리.

Streamlit 판(p95)은 캐시 검증·재예측·정제·에러 표시가 한 함수에 섞여
있었습니다. 여기서는 계산만 하고 경고·오류를 목록으로 돌려줍니다.
"""

from __future__ import annotations

import json
from typing import Any

STRESS_NONE = "선택 안함"


def stress_choices() -> tuple[str, ...]:
    try:
        from p43_degradation_predictor import STRESS_ENV_CHOICES

        return tuple(STRESS_ENV_CHOICES)
    except Exception:
        return (STRESS_NONE,)


def predict(
    parent_smiles: str,
    stress_environment: str,
    temperature_c: float = 25.0,
) -> dict[str, Any]:
    """분해산물을 예측합니다.

    반환:
        products        화면용 dict 목록
        _products_raw   보고서 생성에 쓸 원본 객체 목록
        warnings/errors 사용자에게 보여줄 메시지
    """
    out: dict[str, Any] = {
        "parent_smiles": (parent_smiles or "").strip(),
        "parent_canonical": "",
        "parent_formula": "",
        "parent_mw": None,
        "stress_environment": stress_environment,
        "temperature_c": float(temperature_c),
        "products": [],
        "_products_raw": [],
        "warnings": [],
        "errors": [],
        "engine_version": "",
    }

    if not out["parent_smiles"]:
        out["errors"].append("구조가 비어 있습니다")
        return out

    if not stress_environment or stress_environment == STRESS_NONE:
        out["warnings"].append(
            "가혹 조건을 선택하지 않았습니다. 조건을 골라야 분해 경로가 적용됩니다."
        )
        return out

    try:
        from p02_degradation_constants import ENGINE_CACHE_VERSION

        out["engine_version"] = ENGINE_CACHE_VERSION
    except Exception:
        pass

    try:
        from p43_degradation_predictor import predict_advanced_degradants

        outcome = predict_advanced_degradants(
            out["parent_smiles"], stress_environment, float(temperature_c)
        )
    except Exception as exc:
        out["errors"].append(f"분해 예측 실패: {exc}")
        return out

    out["parent_canonical"] = getattr(outcome, "parent_canonical_smiles", "") or ""
    out["parent_formula"] = getattr(outcome, "parent_formula", "") or ""
    out["parent_mw"] = getattr(outcome, "parent_exact_mw", None)

    if getattr(outcome, "error_message", None):
        out["errors"].append(str(outcome.error_message))
        return out

    products = list(getattr(outcome, "products", []) or [])

    # 보관 중에는 일어나지 않을 대사성 반응을 걸러냅니다 (p44 → p43 경유).
    try:
        from p43_degradation_predictor import sanitize_degradant_products

        before = len(products)
        products = list(sanitize_degradant_products(products, out["parent_smiles"]) or [])
        removed = before - len(products)
        if removed > 0:
            out["warnings"].append(
                f"보관 조건에서 일어나기 어려운 후보 {removed}건을 제외했습니다 "
                "(간 대사에서만 일어나는 방향족 수산화 등)."
            )
        if before and not products:
            out["warnings"].append(
                "예측된 후보가 모두 제외되었습니다. 이 구조·조건에서는 유의미한 분해산물이 없습니다."
            )
    except Exception as exc:
        out["warnings"].append(f"후보 정제 단계 생략: {exc}")

    out["_products_raw"] = products
    out["products"] = [_product_view(p) for p in products]
    return out


def _product_view(product: Any) -> dict[str, Any]:
    get = lambda n, d=None: getattr(product, n, d)  # noqa: E731
    return {
        "id": get("degradant_id", ""),
        "pathway": get("pathway_name", ""),
        "name": get("name", ""),
        "smiles": get("smiles", ""),
        "mw": get("molecular_weight", 0.0),
        "m7_class": get("ich_m7_class"),
        "m7_label": get("ich_m7_label", ""),
        "signal_tier": get("signal_tier", ""),
        "signal_label": get("signal_label_ko", ""),
        "mechanism": get("mechanism", ""),
        "biotransformer": get("biotransformer_pathway", ""),
        "abundance": get("relative_abundance", ""),
        "confidence": get("confidence_score", 0.0),
        "confidence_label": get("confidence_label", ""),
        "db_match": get("db_match_detail", ""),
        "alerts": list(get("alerts", []) or []),
    }


def summary(run: dict[str, Any]) -> dict[str, Any]:
    """상단 지표용 집계."""
    products = run.get("products") or []
    scores = [p["confidence"] for p in products if p.get("confidence") is not None]
    tiers = {"danger": 0, "caution": 0, "safe": 0}
    for p in products:
        tier = p.get("signal_tier")
        if tier in tiers:
            tiers[tier] += 1
    high_risk = sum(1 for p in products if (p.get("m7_class") or 9) <= 3)
    return {
        "count": len(products),
        "avg_confidence": (sum(scores) / len(scores)) if scores else 0.0,
        "tiers": tiers,
        "high_risk": high_risk,
    }


# ---------------------------------------------------------------------------
# 분해 트리
# ---------------------------------------------------------------------------


def build_tree(
    parent_smiles: str,
    *,
    max_depth: int = 3,
    include_ai: bool = True,
) -> dict[str, Any]:
    """다세대 분해 트리를 만들고 직렬화 형식을 한 번에 준비합니다.

    Streamlit 판은 JSON·CSV·Cytoscape·PDF를 화면 그릴 때마다 전부
    만들었습니다. 여기서는 직렬화 텍스트만 만들어 두고 PDF는 실제로
    다운로드를 누를 때 생성합니다.
    """
    out: dict[str, Any] = {
        "parent_smiles": (parent_smiles or "").strip(),
        "max_depth": int(max_depth),
        "include_ai": bool(include_ai),
        "nodes": [],
        "json_text": "{}",
        "cytoscape_text": "{}",
        "cytoscape": None,
        "csv_text": "",
        "error": None,
    }
    if not out["parent_smiles"]:
        out["error"] = "구조가 비어 있습니다"
        return out

    try:
        from p41_degradation_smirks_rules import (
            degradation_tree_cytoscape_json,
            degradation_tree_to_cytoscape,
            degradation_tree_to_json,
            flatten_degradation_tree,
            generate_degradation_tree,
            generate_toxic_degradation_tree,
        )

        if include_ai:
            tree = generate_toxic_degradation_tree(out["parent_smiles"], max_depth=int(max_depth))
        else:
            tree = generate_degradation_tree(out["parent_smiles"], max_depth=int(max_depth))

        out["nodes"] = list(flatten_degradation_tree(tree) or [])
        out["json_text"] = degradation_tree_to_json(tree)
        out["cytoscape_text"] = degradation_tree_cytoscape_json(tree)
        try:
            out["cytoscape"] = degradation_tree_to_cytoscape(tree)
        except Exception:
            out["cytoscape"] = json.loads(out["cytoscape_text"] or "{}")
    except Exception as exc:
        out["error"] = f"분해 트리 생성 실패: {exc}"
        return out

    # CSV — pandas 없이도 되도록 직접 씁니다.
    try:
        import csv
        import io

        rows = out["nodes"]
        if rows:
            buf = io.StringIO()
            keys: list[str] = []
            for row in rows:
                for k in row:
                    if k not in keys:
                        keys.append(k)
            writer = csv.DictWriter(buf, fieldnames=keys, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({k: row.get(k, "") for k in keys})
            out["csv_text"] = "﻿" + buf.getvalue()
    except Exception:
        out["csv_text"] = ""

    return out


def backend_core_api(
    parent_smiles: str,
    *,
    formulation_ph: float = 7.4,
    max_depth: int = 3,
    include_structure_b64: bool = True,
) -> dict[str, Any]:
    """물성 + Read-Across가 붙은 Backend Core API 응답."""
    try:
        from p51_prime_backend_core import build_degradation_api_response

        payload = build_degradation_api_response(
            parent_smiles,
            formulation_ph=float(formulation_ph),
            max_depth=int(max_depth),
            read_across_top_n=1,
            include_structure_b64=bool(include_structure_b64),
        )
        return {"ok": True, "payload": payload, "error": None}
    except Exception as exc:
        return {"ok": False, "payload": None, "error": str(exc)}
