"""
하이브리드 공정 예측 엔진 — Reaction + Degradation + 3-way Cross-Validation.

① drug_index.db 실측 대조
② PubChem API (Ames/변이원성 BioAssay)
③ Ashby-Tennant Rule-based 알러트
"""

# 경로 A

# 지휘자 + 3중 대조

# 생성된 각 물질마다:
#   ① p24_db_smiles_lookup   자체 DB에 있나? (실측/이론)
#   ② p25_pubchem_client     PubChem에 Ames 데이터 있나?
#   ③ p62_toxicity_alerts    구조 경보 걸리나?

# 8개 범주 표
# cross_reaction · ndsri · reagent_side · reagent_ndsri
# reagent_cross · degradant · unreacted · degradation_risk

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal

from rdkit import Chem

from p24_db_smiles_lookup import classify_smiles_against_db
from p40_degradation_engine import DegradationRiskProfile, scan_degradation_vulnerabilities
from p66_ich_m7 import ICHM7Result, classify_ich_m7
from p50_physicochemical import stability_summary_for_table
from p32_process_reaction_engine import ProcessAnalysisResult, analyze_process_reaction
from p25_pubchem_client import fetch_full_toxicity_profile
from p30_reaction_engine import ImpurityClassification, ImpurityResult
from p62_toxicity_alerts import alert_status_label, scan_ashby_tennant_alerts

Category = Literal[
    "cross_reaction",
    "ndsri",
    "reagent_side",
    "reagent_ndsri",
    "reagent_cross",
    "degradant",
    "unreacted",
    "degradation_risk",
]


@dataclass
class HybridProcessRow:
    """종합 공정 보고서 1행."""
    category: Category
    category_label: str
    smiles: str
    rule_label: str
    verdict: str
    is_validated: bool
    db_name: str
    db_cid: str
    pubchem_cid: str
    pubchem_link: str
    pubchem_status: str
    genotox_summary: str
    alert_status: str
    alert_emoji: str
    ashby_tennant: str
    ich_m7_class: int
    ich_m7_label: str
    logp: str
    stability: str
    highlight_atoms: list[int] = field(default_factory=list)
    impurity: ImpurityResult | None = None
    pubchem_profile: dict[str, Any] = field(default_factory=dict)
    degradation_flags: list[str] = field(default_factory=list)


@dataclass
class HybridProcessResult:
    reactant1: str
    reactant2: str
    reagent: str
    product: str
    reference_drug: str
    nitrite_in_process: bool
    analysis: ProcessAnalysisResult
    product_risk: DegradationRiskProfile
    rows: list[HybridProcessRow] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    validated_count: int = 0
    theoretical_count: int = 0
    reagent_roles: list[str] = field(default_factory=list)


def _safe_mol(smiles: str) -> Chem.Mol | None:
    try:
        if not smiles or not str(smiles).strip():
            return None
        mol = Chem.MolFromSmiles(str(smiles).strip())
        if mol is None:
            return None
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def _enrich_species(
    imp: ImpurityResult,
    *,
    category: Category,
    category_label: str,
    local_df,
    drug_mol: Chem.Mol | None,
    degradation_flags: list[str] | None = None,
) -> tuple[HybridProcessRow | None, str | None]:
    """단일 불순물/분해산물에 3-way 하이브리드 검증 적용."""
    try:
        mol = _safe_mol(imp.smiles)
        if mol is None:
            return None, f"SMILES 해석 실패: {imp.smiles[:40]}"

        alerts = imp.ashby_tennant_alerts
        hl = imp.highlight_atoms
        if not alerts:
            try:
                alerts, hl = scan_ashby_tennant_alerts(mol)
            except Exception:
                alerts, hl = [], []

        emoji, alert_status = alert_status_label(alerts)

        cls_str, db_rec = imp.classification.value, None
        try:
            cls_str, db_rec = classify_smiles_against_db(imp.smiles, local_df=local_df)
        except Exception:
            pass

        is_validated = cls_str == "validated" or imp.classification == ImpurityClassification.VALIDATED
        if is_validated:
            verdict = "📌 [실측 데이터 기반 유연물질 (Validated)]"
        else:
            verdict = "🧪 [오픈소스 규칙 기반 예측 불순물 (Theoretical)]"

        pub_profile: dict[str, Any] = {}
        pub_err: str | None = None
        try:
            pub_profile = fetch_full_toxicity_profile(imp.smiles)
            if pub_profile.get("error"):
                pub_err = str(pub_profile["error"])
        except Exception as exc:
            pub_err = str(exc)
            pub_profile = {"api_status": "error", "api_message": "⚪ PubChem 조회 실패"}

        compound = pub_profile.get("compound") or {}
        cid = compound.get("cid") or imp.db_cid or (db_rec or {}).get("cid")
        link = compound.get("link") or (f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else "")
        genotox = pub_profile.get("genotox_assays") or []
        if genotox:
            genotox_summary = f"PubChem Ames/변이원성 BioAssay {len(genotox)}건"
        elif pub_profile.get("genotox_found"):
            genotox_summary = "변이원성 시험 매핑 있음"
        else:
            genotox_summary = "변이원성 BioAssay 미검출"

        try:
            m7 = classify_ich_m7(mol, alerts, drug_mol=drug_mol)
        except Exception:
            m7 = ICHM7Result(
                class_number=5,
                class_name_ko="Class 5",
                definition="",
                control_strategy="",
                rationale="분류 실패",
            )

        try:
            stab, logp = stability_summary_for_table(imp.smiles)
        except Exception:
            stab, logp = "—", "—"

        db_name = (
            (db_rec or {}).get("name")
            or (db_rec or {}).get("drug_name")
            or imp.db_name
            or compound.get("name")
            or "—"
        )

        row = HybridProcessRow(
            category=category,
            category_label=category_label,
            smiles=imp.smiles,
            rule_label=imp.rule_label,
            verdict=verdict,
            is_validated=is_validated,
            db_name=str(db_name),
            db_cid=str((db_rec or {}).get("cid") or imp.db_cid or "—"),
            pubchem_cid=str(cid or "—"),
            pubchem_link=str(link or ""),
            pubchem_status=str(pub_profile.get("api_message", "—")),
            genotox_summary=genotox_summary,
            alert_status=alert_status,
            alert_emoji=emoji,
            ashby_tennant=", ".join(alerts) if alerts else "없음",
            ich_m7_class=m7.class_number,
            ich_m7_label=f"Class {m7.class_number}",
            logp=str(logp),
            stability=str(stab),
            highlight_atoms=hl,
            impurity=imp,
            pubchem_profile=pub_profile,
            degradation_flags=degradation_flags or [],
        )
        return row, pub_err
    except Exception as exc:
        return None, f"종합 검증 오류 ({imp.smiles[:30]}…): {exc}"


def run_hybrid_process_analysis(
    reactant1: str,
    reactant2: str,
    product: str,
    *,
    reagent: str = "",
    reference_drug_smiles: str = "",
    nitrite_in_process: bool = False,
    local_df=None,
) -> HybridProcessResult:
    """
    하이브리드 공정 예측 — Reaction + Degradation + 3-way 매칭.
    모든 RDKit/API 단계는 try-except로 보호됩니다.
    """
    errors: list[str] = []
    product_risk = DegradationRiskProfile(canon_smiles="")
    analysis = ProcessAnalysisResult()
    rows: list[HybridProcessRow] = []
    ref_smi = (reference_drug_smiles or "").strip()
    effective_product = (product or ref_smi or "").strip()

    try:
        try:
            product_risk = (
                scan_degradation_vulnerabilities(effective_product)
                if effective_product
                else product_risk
            )
        except Exception as exc:
            errors.append(f"분해 취약성 스캔 오류: {exc}")

        def _safe_cls(s: str) -> tuple[str, dict | None]:
            try:
                return classify_smiles_against_db(s, local_df=local_df)
            except Exception:
                return "theoretical", None

        def _safe_alerts(mol: Chem.Mol) -> tuple[list[str], list[int]]:
            try:
                return scan_ashby_tennant_alerts(mol)
            except Exception:
                return [], []

        try:
            analysis = analyze_process_reaction(
                reactant1,
                reactant2,
                effective_product,
                reagent=reagent or "",
                classification_fn=_safe_cls,
                alert_fn=_safe_alerts,
                nitrite_in_process=nitrite_in_process,
            )
        except Exception as exc:
            errors.append(f"반응 엔진 오류: {exc}")

        drug_mol = _safe_mol(ref_smi) if ref_smi else (_safe_mol(effective_product) if effective_product else None)

        category_map: list[tuple[list[ImpurityResult], Category, str]] = [
            (getattr(analysis, "reagent_reaction_products", []), "reagent_side", "시약 부반응 / NDSRI"),
            (analysis.cross_reaction_products, "cross_reaction", "교차반응 / NDSRI"),
            (analysis.degradants, "degradant", "목적물 분해산물"),
            (analysis.unreacted_starting_materials, "unreacted", "미반응 출발물질·시약"),
        ]

        for items, cat, cat_label in category_map:
            for imp in items or []:
                try:
                    flags = product_risk.vulnerability_labels if cat == "degradant" else []
                    row, err = _enrich_species(
                        imp,
                        category=cat,
                        category_label=cat_label,
                        local_df=local_df,
                        drug_mol=drug_mol,
                        degradation_flags=flags,
                    )
                    if row:
                        rows.append(row)
                    if err:
                        errors.append(err)
                except Exception as exc:
                    errors.append(f"종목 검증 스킵: {exc}")
    except Exception as exc:
        errors.append(f"공정 예측 엔진 치명 오류: {exc}")

    validated = sum(1 for r in rows if r.is_validated)
    theoretical = len(rows) - validated

    reagent_roles: list[str] = []
    profile = getattr(analysis, "reagent_profile", None)
    if profile:
        reagent_roles = list(getattr(profile, "role_labels_ko", ()))

    return HybridProcessResult(
        reactant1=reactant1 or "",
        reactant2=reactant2 or "",
        reagent=reagent or "",
        product=effective_product,
        reference_drug=ref_smi,
        nitrite_in_process=nitrite_in_process,
        analysis=analysis,
        product_risk=product_risk,
        rows=rows,
        errors=errors,
        validated_count=validated,
        theoretical_count=theoretical,
        reagent_roles=reagent_roles,
    )


def hybrid_rows_to_report_dicts(rows: list[HybridProcessRow]) -> list[dict]:
    """ICH M7 부록용 dict 목록."""
    out: list[dict] = []
    for r in rows:
        out.append(
            {
                "name": r.rule_label or r.db_name or r.smiles[:40],
                "smiles": r.smiles,
                "classification": r.verdict,
                "ich_m7_class": r.ich_m7_class,
                "logp": r.logp,
                "stability": r.stability,
                "alerts": r.ashby_tennant.split(", ") if r.ashby_tennant != "없음" else [],
                "result": r.impurity,
                "m7": None,
            }
        )
    return out
