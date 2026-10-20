"""
PRIME-DS unified JSON — Process_Synthesis_Analysis + Degradation_MultiStep_Tree.

Frontend / regulatory API binding format for hybrid process + toxic degradation tree.
"""

from __future__ import annotations

import json
from typing import Any

from p31_process_reaction_rules import predict_process_reactions


def _format_ich_m7_status_export(tox: dict[str, Any]) -> str:
    status = str(tox.get("ICH_M7_Status") or "").strip()
    cls = tox.get("ICH_M7_Class")

    if "AI Inferred" in status or "3-Alert" in status:
        return "Class 3-Alert (AI Inferred Liability)"

    if cls == 1 or status.startswith("Class 1"):
        return "Class 1 (Known mutagenic carcinogen)"
    if cls == 2 or status.startswith("Class 2"):
        return "Class 2 (Known mutagen with unknown carcinogenicity)"
    if cls == 3 or (status.startswith("Class 3") and "3-Alert" not in status):
        return "Class 3 (Mutagenic impurity, spec control required)"
    if cls == 4 or status.startswith("Class 4"):
        return "Class 4 (Alerting structure — confirmatory Ames recommended)"
    if cls == 5 or status.startswith("Class 5"):
        return "Class 5 (Non-mutagenic substance)"

    if status:
        return status
    return "Class 5 (Non-mutagenic substance)"


def _slim_toxicity_analysis(tox: dict[str, Any] | None) -> dict[str, Any]:
    """API/차트 바인딩용 축약 — 단, 보고서 생성에 필요한 키는 버리지 않습니다.

    2026-09-03 수정 (E-5): 기존 축약은 `ICH_M7_Class` 와 `GNN_Model` 을
    버렸습니다. 그 결과 slim 트리를 p87 로 넘기면
      - `isinstance(m7_cls, int)` 가 항상 False 라 노드별 ICH M7 executive
        callout 이 예외 없이 통째로 사라지고,
      - `_node_detail_rows` 의 "GNN model" 이 항상 N/A 로 인쇄되었습니다.
    두 키를 보존합니다(기존 키는 그대로이므로 공개 스키마는 하위 호환).
    """
    tox = tox or {}
    slim: dict[str, Any] = {
        "Ames_AI_Probability": str(tox.get("Ames_AI_Probability") or "N/A"),
        "Prediction_Class": str(tox.get("Prediction_Class") or "Unknown"),
        "ICH_M7_Status": _format_ich_m7_status_export(tox),
        "GNN_Model": str(tox.get("GNN_Model") or "N/A"),
    }
    raw_cls = tox.get("ICH_M7_Class")
    if isinstance(raw_cls, bool):
        raw_cls = None
    if isinstance(raw_cls, int) and 1 <= raw_cls <= 5:
        slim["ICH_M7_Class"] = raw_cls
    elif isinstance(raw_cls, str) and raw_cls.strip().isdigit() and 1 <= int(raw_cls.strip()) <= 5:
        slim["ICH_M7_Class"] = int(raw_cls.strip())
    rationale = str(tox.get("ICH_M7_Rationale") or "").strip()
    if rationale:
        slim["ICH_M7_Rationale"] = rationale
    return slim


def format_process_synthesis_analysis(proc_report: dict[str, Any]) -> dict[str, Any]:
    """Map predict_process_reactions() output to Process_Synthesis_Analysis schema."""
    mechanisms = proc_report.get("Triggered_Mechanisms") or []
    main_products = proc_report.get("Main_Products") or []
    side_rows = proc_report.get("Side_Products_Impurities") or []

    if len(mechanisms) == 1:
        triggered = mechanisms[0]
    elif mechanisms:
        triggered = ", ".join(mechanisms)
    else:
        triggered = ""

    main_smi = ""
    if main_products:
        main_smi = str(main_products[0].get("SMILES") or "")

    impurities = [
        {
            "Side_Product_SMILES": str(row.get("SMILES") or ""),
            "Alert_Description": str(row.get("Alert_Description") or ""),
        }
        for row in side_rows
        if row.get("SMILES")
    ]

    return {
        "Triggered_Mechanism": triggered,
        "Main_Product_SMILES": main_smi,
        "Process_Impurities": impurities,
    }


def format_degradation_multistep_tree(tree: dict[str, Any] | None) -> dict[str, Any]:
    """Slim degradation tree for API / chart binding (toxicity fields only)."""
    if not tree:
        return {}

    node: dict[str, Any] = {
        "name": tree.get("name", ""),
        "generation": tree.get("generation", 0),
        "children": [
            format_degradation_multistep_tree(child)
            for child in (tree.get("children") or [])
        ],
    }
    if tree.get("triggered_by"):
        node["triggered_by"] = tree["triggered_by"]
    if tree.get("toxicity_analysis"):
        node["toxicity_analysis"] = _slim_toxicity_analysis(tree["toxicity_analysis"])
    return node


def build_prime_integrated_report(
    *,
    reactant_smiles_list: list[str],
    degradation_parent_smiles: str,
    max_depth: int = 3,
    include_toxicity: bool = True,
) -> dict[str, Any]:
    """
    Build unified PRIME-DS report:
    Process_Synthesis_Analysis + Degradation_MultiStep_Tree.
    """
    from p41_degradation_smirks_rules import (
        generate_degradation_tree,
        generate_toxic_degradation_tree,
    )

    proc_report = predict_process_reactions(reactant_smiles_list)
    parent = (degradation_parent_smiles or "").strip()
    if not parent:
        mains = proc_report.get("Main_Products") or []
        if mains:
            parent = str(mains[0].get("SMILES") or "")

    if parent:
        if include_toxicity:
            raw_tree = generate_toxic_degradation_tree(parent, max_depth=max_depth)
        else:
            raw_tree = generate_degradation_tree(
                parent, max_depth=max_depth, include_toxicity=False
            )
    else:
        raw_tree = None

    return {
        "Process_Synthesis_Analysis": format_process_synthesis_analysis(proc_report),
        "Degradation_MultiStep_Tree": format_degradation_multistep_tree(raw_tree),
    }


def integrated_report_to_json(report: dict[str, Any], *, indent: int = 2) -> str:
    return json.dumps(report, indent=indent, ensure_ascii=False)


if __name__ == "__main__":
    sample = build_prime_integrated_report(
        reactant_smiles_list=[
            "O=C(O)Cc1ccccc1",
            "CCN(CC)",
            "CC(=O)O",
        ],
        degradation_parent_smiles="CCOC(=O)CCSCCc1ccccc1",
        max_depth=2,
    )
    print(integrated_report_to_json(sample))
