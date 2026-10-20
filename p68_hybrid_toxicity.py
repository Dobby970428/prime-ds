"""
[상태 메모] 2026-09-03 확인: 이 저장소의 어떤 `p*.py`/`s*.py`도 이 모듈을
import 하지 않습니다. 공개 심볼 `HybridToxicityAssessment` ·
`run_hybrid_toxicity_assessment` 는 문자열 검색을 포함해도 다른 .py 파일에서
이름조차 등장하지 않습니다(참조는 문서 2건뿐: `코드스키마.md` 범위 표기,
`60번대_설명.md:351`). 이 스냅샷에는 `webapp/` 폴더가 없어 화면 쪽 참조 여부를
확인할 수 없으므로 삭제하지 않고 보존합니다.
→ `webapp/` 쪽 참조 여부를 확인한 뒤 유지·삭제를 결정하세요.

하이브리드 독성 평가 — In-House DB + PubChem + Ashby-Tennant (Derek/Sarah style).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from rdkit import Chem

from p23_chem_search_engine import safe_canonical_smiles
from p24_db_smiles_lookup import classify_smiles_against_db
from p66_ich_m7 import ICHM7Result, classify_ich_m7
from p25_pubchem_client import fetch_full_toxicity_profile
from p62_toxicity_alerts import alert_status_label, scan_ashby_tennant_alerts

try:
    from p52_ambit_client import ambit_opensource_rows
except ImportError:
    ambit_opensource_rows = None

Classification = Literal["validated", "theoretical"]


@dataclass
class HybridToxicityAssessment:
    canon_smiles: str
    classification: Classification
    verdict_label: str
    db_record: dict | None
    alerts: list[str] = field(default_factory=list)
    highlight_atoms: list[int] = field(default_factory=list)
    alert_status: str = ""
    alert_emoji: str = ""
    pubchem: dict[str, Any] = field(default_factory=dict)
    m7: ICHM7Result | None = None
    inhouse_table: list[dict] = field(default_factory=list)
    opensource_table: list[dict] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)


def _verdict_label(cls: Classification) -> str:
    if cls == "validated":
        return "📌 판정: [실측 데이터 기반 유연물질 (Validated)]"
    return "🧪 판정: [오픈소스 규칙 기반 예측 불순물 (Theoretical)]"


def run_hybrid_toxicity_assessment(
    smiles: str,
    *,
    local_df=None,
    drug_smiles: str | None = None,
) -> HybridToxicityAssessment | None:
    """
    단일 SMILES에 대한 하이브리드 독성 평가.
    실패 시 None (SMILES 오류).
    """
    canon, err = safe_canonical_smiles(smiles)
    if err or not canon:
        return None

    mol = Chem.MolFromSmiles(canon)
    if mol is None:
        return None

    errors: list[str] = []
    cls_str, db_rec = classify_smiles_against_db(canon, local_df=local_df)
    classification: Classification = "validated" if cls_str == "validated" else "theoretical"

    alerts, hl = scan_ashby_tennant_alerts(mol)
    emoji, alert_status = alert_status_label(alerts)

    pub_profile: dict[str, Any] = {}
    try:
        pub_profile = fetch_full_toxicity_profile(canon)
        if pub_profile.get("error"):
            errors.append(str(pub_profile["error"]))
    except Exception as exc:
        errors.append(f"PubChem 연동 예외: {exc}")
        pub_profile = {"api_status": "error", "api_message": "⚪ PubChem 조회 실패"}

    drug_mol = None
    if drug_smiles:
        try:
            drug_mol = Chem.MolFromSmiles(drug_smiles.strip())
        except Exception:
            errors.append("원료의약품 SMILES 해석 실패")

    m7 = classify_ich_m7(mol, alerts, drug_mol=drug_mol)

    # --- In-house table ---
    inhouse_rows: list[dict] = [
        {"항목": "자체 DB (13,040종)", "결과": f"{cls_str.upper()} — {classification}"},
    ]
    if db_rec:
        inhouse_rows.extend(
            [
                {"항목": "매칭 화합물명", "결과": db_rec.get("name") or db_rec.get("drug_name") or "—"},
                {"항목": "의약품 파일", "결과": db_rec.get("drug_file") or db_rec.get("drug_name") or "—"},
                {"항목": "Compound CID", "결과": str(db_rec.get("cid") or "—")},
                {"항목": "Canonical SMILES", "결과": db_rec.get("canon_smiles") or canon},
            ]
        )
    else:
        inhouse_rows.append({"항목": "매칭 결과", "결과": "자체 DB에 동일 Canonical SMILES 없음"})

    # --- Open-source table ---
    compound = pub_profile.get("compound") or {}
    opensource_rows: list[dict] = [
        {"항목": "PubChem API", "결과": pub_profile.get("api_message", "—")},
        {"항목": "PubChem CID", "결과": str(compound.get("cid") or "—")},
        {"항목": "화합물명", "결과": compound.get("name") or "—"},
        {"항목": "분자식 / MW", "결과": f"{compound.get('formula', '—')} / {compound.get('mw', '—')}"},
        {
            "항목": "Ashby-Tennant 스캔",
            "결과": f"{emoji} {alert_status}",
        },
        {
            "항목": "검출 알러트",
            "결과": ", ".join(alerts) if alerts else "없음",
        },
        {
            "항목": "ICH M7 Class (Consensus)",
            "결과": (
                f"Class {m7.class_number} — {m7.empirical_db_match} · Ames {m7.ames_evidence}"
                if m7
                else "—"
            ),
        },
        {
            "항목": "Empirical DB (Hansen/ISSCAN/Tox21)",
            "결과": (
                (
                    f"{m7.empirical_db_match} · Carc {m7.carcinogenicity_evidence} · "
                    f"{m7.recommended_ttc[:60]}{'…' if len(m7.recommended_ttc) > 60 else ''}"
                )
                if m7 and m7.empirical_db_match != "None"
                else "No exact match in empirical M7 index"
            ),
        },
        {
            "항목": "Layer C GNN P(Mutagenic)",
            "결과": (
                f"{m7.gnn_mutagenic_probability:.3f} · {m7.gnn_model_name}"
                if m7 and m7.gnn_mutagenic_probability is not None
                else "GNN unavailable (RDKit/DB standalone mode)"
            ),
        },
        {
            "항목": "Ensemble Assigned Class",
            "결과": (
                f"{m7.assigned_class_label or f'Class {m7.class_number}'} "
                f"(rule Class {m7.rule_based_class})"
                if m7 and m7.rule_based_class is not None
                else (f"Class {m7.class_number}" if m7 else "—")
            ),
        },
    ]

    if ambit_opensource_rows is not None:
        try:
            opensource_rows.extend(ambit_opensource_rows(canon))
        except Exception as exc:
            opensource_rows.append({"항목": "OECD QSAR (AMBIT2)", "결과": f"조회 실패: {exc}"})

    genotox = pub_profile.get("genotox_assays") or []
    if genotox:
        opensource_rows.append(
            {
                "항목": "Ames/변이원성 BioAssay",
                "결과": f"{len(genotox)}건 PubChem 매핑",
            }
        )
        for i, assay in enumerate(genotox[:3], start=1):
            opensource_rows.append(
                {
                    "항목": f"  └ BioAssay #{i}",
                    "결과": f"AID {assay.get('aid')} · {assay.get('name', '')[:60]} · {assay.get('activity_outcome', '')}",
                }
            )
    else:
        opensource_rows.append(
            {"항목": "Ames/변이원성 BioAssay", "결과": "PubChem 요약에 관련 시험 미검출"}
        )

    return HybridToxicityAssessment(
        canon_smiles=canon,
        classification=classification,
        verdict_label=_verdict_label(classification),
        db_record=db_rec,
        alerts=alerts,
        highlight_atoms=hl,
        alert_status=alert_status,
        alert_emoji=emoji,
        pubchem=pub_profile,
        m7=m7,
        inhouse_table=inhouse_rows,
        opensource_table=opensource_rows,
        errors=errors,
    )
