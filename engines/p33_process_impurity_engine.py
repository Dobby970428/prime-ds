"""
PRIME-DS — Process Impurity Engine (독립 구동형 공정 유연물질 예측)
====================================================================
RDKit AllChem 반응 SMARTS 기반 가혹 분해(degradant) 생성 및
Morgan Fingerprint Tanimoto 유사도 매칭(내부 10K 검증 DB 시뮬레이션).
"""

# 경로 B
# 가혹조건·촉매
# 공정 조건 패널에서 고른 가혹 조건(산성/염기성/산화)과 촉매를 반영

from __future__ import annotations

import pickle
from functools import lru_cache
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors

# ---------------------------------------------------------------------------
# Stress degradation pathways — UI selectbox mapping
# ---------------------------------------------------------------------------
STRESS_NONE = "선택안함 (No Accelerated Stress)"
STRESS_ACIDIC_HYDROLYSIS = "Acidic Hydrolysis (Harsh)"
STRESS_ALKALINE_HYDROLYSIS = "Alkaline Hydrolysis (Harsh)"
STRESS_OXIDATIVE = "Oxidative Stress Conditions"

STRESS_CONDITION_OPTIONS: list[str] = [
    STRESS_NONE,
    STRESS_ACIDIC_HYDROLYSIS,
    STRESS_ALKALINE_HYDROLYSIS,
    STRESS_OXIDATIVE,
]

# Legacy Korean labels → English (backward compatibility)
_STRESS_CONDITION_ALIASES: dict[str, str] = {
    "선택안함": STRESS_NONE,
    "선택 안함": STRESS_NONE,
    "No Stress": STRESS_NONE,
    "산성 가혹 가수분해 (Acidic Hydrolysis)": STRESS_ACIDIC_HYDROLYSIS,
    "염기성 가혹 가수분해 (Alkaline Hydrolysis)": STRESS_ALKALINE_HYDROLYSIS,
    "과산화 산화 조건 (Oxidative Stress)": STRESS_OXIDATIVE,
}


def normalize_stress_condition(condition: str) -> str:
    """Map legacy Korean stress labels to current English constants."""
    return _STRESS_CONDITION_ALIASES.get(condition, condition)

# Executive Dashboard / Export column names (English)
IMP_COL_SOURCE = "Source Structure"
IMP_COL_PATHWAY = "Stress Degradation Pathway"
IMP_COL_SMILES = "Predicted Impurity / Scaffold (SMILES)"
IMP_COL_MW = "Exact MW (g/mol)"
IMP_COL_ICH_M7 = "ICH M7 Genotoxicity Assessment"
IMP_COL_ICH_M7_CLASS = "ICH M7 Class"
IMP_COL_REGULATORY = "Regulatory Guidance"
IMP_COL_CONFIDENCE = "🎯 Confidence & DB Verification"
IMP_COL_IMPURITY_CLASS = "유연물질 분류"
IMP_COL_MECHANISM_UI = "메커니즘 설명"
IMP_COL_MECHANISM_RULE = "매칭 메커니즘"
IMP_COL_IMAGE = "🧪 2D Structure Image"

ICH_M7_CLASS_EMOJI: dict[int, str] = {
    1: "🔴",
    2: "🟠",
    3: "🟧",
    4: "🔵",
    5: "🟢",
}

IMPURITY_RESULT_COLUMNS: list[str] = [
    IMP_COL_SOURCE,
    IMP_COL_IMPURITY_CLASS,
    IMP_COL_MECHANISM_RULE,
    IMP_COL_MECHANISM_UI,
    IMP_COL_PATHWAY,
    IMP_COL_SMILES,
    IMP_COL_MW,
    IMP_COL_ICH_M7_CLASS,
    IMP_COL_REGULATORY,
    IMP_COL_ICH_M7,
    IMP_COL_CONFIDENCE,
]

# ---------------------------------------------------------------------------
# 반응 SMARTS — 기전론적(Mechanistic) 변환 정의
# ---------------------------------------------------------------------------
ESTER_CLEAVAGE_SMARTS = "[#6:1](=[O:2])-[O:3]-[#6:4]>>[#6:1](=[O:2])O.[O:3]-[#6:4]"
AMIDE_CLEAVAGE_SMARTS = "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])O.[N:3]"
SULFUR_OXIDATION_SMARTS = "[S:1]>>[S:1]=O"
NITROGEN_OXIDATION_SMARTS = "[#7:1]>>[#7:1]=O"

HYDROLYSIS_REACTIONS: tuple[str, ...] = (ESTER_CLEAVAGE_SMARTS, AMIDE_CLEAVAGE_SMARTS)
OXIDATION_REACTIONS: tuple[str, ...] = (SULFUR_OXIDATION_SMARTS, NITROGEN_OXIDATION_SMARTS)

# 촉매 매개 homo-coupling 부반응 SMARTS (2× 동일 기질 → 대칭 이량체)
HOMO_COUPLING_REACTIONS: tuple[tuple[str, str], ...] = (
    (
        "aryl_halide",
        "[c:1]-[Cl,Br,I:2].[c:3]-[Cl,Br,I:4]>>[c:1]-[c:3]",
    ),
    (
        "vinyl_halide",
        "[C:1]=[C:2]-[Cl,Br,I:3].[C:4]=[C:5]-[Cl,Br,I:6]>>[C:1]=[C:2]-[C:4]=[C:5]",
    ),
    (
        "terminal_alkyne",
        "[C:1]#[CH:2].[C:3]#[CH:4]>>[C:1]#[C:2]-[C:3]#[C:4]",
    ),
    (
        "aryl_boronic",
        "[c:1]-B(O)O.[c:2]-B(O)O>>[c:1]-[c:2]",
    ),
)

# 공정 시약 C — 촉매 분류 (homo-coupling 예측 활성화 대상)
Pd_CATALYSTS: frozenset[str] = frozenset(
    {
        "Pd(PPh3)4",
        "Pd(PPh3)4 (Tetrakis)",
        "Pd/C",
        "Palladium on Carbon (Pd/C)",
        "Pd(OAc)2",
        "Palladium(II) Acetate (Pd(OAc)2)",
        "Pd(dppf)Cl2",
        "Pd2(dba)3",
        "Pd",
    }
)
Ni_CATALYSTS: frozenset[str] = frozenset(
    {"NiCl2", "Ni(dppp)Cl2", "Nickel(II) Chloride (NiCl2)", "Ni(PPh3)4", "Ni"}
)
Cu_CATALYSTS: frozenset[str] = frozenset(
    {"CuI", "Copper(I) Iodide (CuI)", "CuBr", "Cu"}
)
CATALYST_REAGENTS: frozenset[str] = Pd_CATALYSTS | Ni_CATALYSTS | Cu_CATALYSTS

# ICH M7 유전독성 구조 경보 (Structural Alerts) — SMARTS 기반 휴리스틱
M7_ALERT_SMARTS: list[tuple[str, str, str]] = [
    ("[N;$(N-c1ccccc1)]", "ICH M7 Class 1", "🚨 High risk — aromatic amine (known mutagenicity concern)"),
    ("[N+](=O)[O-]", "ICH M7 Class 2", "🚨 High risk — nitro group"),
    ("[nH]1cccc1", "ICH M7 Class 3", "⚠️ Alert — heteroaromatic N (structural alert)"),
    ("C=CC=O", "ICH M7 Class 3", "⚠️ Alert — α,β-unsaturated carbonyl (Michael acceptor)"),
    ("[Cl,Br,I]", "ICH M7 Class 4", "ℹ️ Alkyl halide — potential alkylating agent (verify)"),
]

M7_DEFAULT_CLASS = ("ICH M7 Class 5", "✅ No structural alert detected — low risk estimate")

DB_TAG_VERIFIED = "🔥 ≥90% scaffold match with verified detection records"
DB_TAG_MECHANISTIC = "⚠️ Mechanistic pathway & analogous functional group validated"
DB_TAG_MODEL = "💡 Chemoinformatics model-based prediction"

MORGAN_RADIUS = 2
MORGAN_NBITS = 2048
VERIFICATION_DB_TARGET_SIZE = 10_000
TEMPERATURE_BOOST_THRESHOLD_C = 60.0


# ---------------------------------------------------------------------------
# 분자 · 반응 유틸리티
# ---------------------------------------------------------------------------
def canonical_smiles(mol: Chem.Mol) -> str:
    """RDKit canonical SMILES (중복·원본 비교용)."""
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        log_invalid_smiles_bypass(context="canonical_smiles")
        return ""


def safe_sanitize(mol: Chem.Mol, *, log_bypass: bool = True) -> Chem.Mol | None:
    """SanitizeMol 성공 시 Mol 반환, 실패 시 None."""
    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        if log_bypass:
            log_invalid_smiles_bypass(context="sanitize")
        return None


def mol_from_smiles(smiles: str, *, log_bypass: bool = True) -> Chem.Mol | None:
    """SMILES → Mol (실패·원자가 오류 시 None, bypass 로그)."""
    cleaned = (smiles or "").strip()
    if not cleaned:
        return None
    try:
        mol = Chem.MolFromSmiles(cleaned)
        if mol is None:
            if log_bypass:
                log_invalid_smiles_bypass(context=cleaned[:80])
            return None
        sanitized = safe_sanitize(mol)
        return sanitized
    except Exception:
        if log_bypass:
            log_invalid_smiles_bypass(context=cleaned[:80])
        return None


def unpack_reaction_product_matrix(product_matrix: tuple) -> list[Chem.Mol]:
    """
    RunReactants 반환 Product Matrix Tuple 전체 순회.

  CRITICAL: 각 outcome 튜플 내 모든 생성물(소분자 Fragment + 잔존 Core Scaffold)을
    독립 Mol로 언패킹. 단일 첫 번째 생성물만 취하는 유기화학적 버그 방지.
    """
    unpacked: list[Chem.Mol] = []
    if not product_matrix:
        return unpacked

    for outcome in product_matrix:
        if outcome is None:
            continue
        # outcome: 생성물 Mol 튜플 (다분자 cleavage 시 2개 이상)
        if isinstance(outcome, Chem.Mol):
            mol = safe_sanitize(outcome)
            if mol is not None:
                unpacked.append(mol)
            continue
        for prod in outcome:
            if prod is None:
                continue
            mol = safe_sanitize(prod)
            if mol is not None:
                unpacked.append(mol)
    return unpacked


def run_reaction_smarts(mol: Chem.Mol, smarts: str) -> list[Chem.Mol]:
    """단일 SMARTS 반응 실행 후 모든 생성물 반환."""
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
    except Exception:
        return []
    if rxn is None:
        return []
    try:
        outcomes = rxn.RunReactants((mol,))
    except Exception:
        return []
    return unpack_reaction_product_matrix(outcomes)


def run_bimolecular_reaction_smarts(mol: Chem.Mol, smarts: str) -> list[Chem.Mol]:
    """동일 기질 2분자 반응 (homo-coupling 등) — 모든 생성물 반환."""
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
    except Exception:
        return []
    if rxn is None:
        return []
    try:
        outcomes = rxn.RunReactants((mol, mol))
    except Exception:
        return []
    return unpack_reaction_product_matrix(outcomes)


def is_catalyst_reagent(reagent: str) -> bool:
    """시약 C가 Pd/Ni/Cu 촉매인지 여부 (homo-coupling 예측 트리거)."""
    from engines.p13_chemical_registry import is_catalyst_reagent_name, normalize_reagent_name

    cleaned = normalize_reagent_name(reagent)
    if is_catalyst_reagent_name(cleaned):
        return True
    return cleaned in CATALYST_REAGENTS


def catalyst_family_label(reagent: str) -> str:
    """촉매 금속 계열 표시명."""
    from engines.p13_chemical_registry import catalyst_metal_label, normalize_reagent_name

    cleaned = normalize_reagent_name(reagent)
    label = catalyst_metal_label(cleaned)
    if label != "Catalyst":
        return label
    if cleaned in Pd_CATALYSTS:
        return "Pd"
    if cleaned in Ni_CATALYSTS:
        return "Ni"
    if cleaned in Cu_CATALYSTS:
        return "Cu"
    return "Catalyst"


def deduplicate_molecules(
    mols: list[Chem.Mol],
    *,
    exclude_canonical: set[str] | None = None,
    ctx: ProcessRunContext | None = None,
) -> list[Chem.Mol]:
    """Canonical SMILES 기준 중복 제거 및 exclude 집합 필터 (라디칼 중간체 제외)."""
    try:
        from p44_degradation_shelf_life_filter import is_spin_stabilized_neutral_mol as _neutral_mol
    except ImportError:
        def _neutral_mol(mol: Chem.Mol | None) -> bool:
            if mol is None:
                return False
            try:
                if sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()) >= 1:
                    return False
                Chem.SanitizeMol(mol)
                return True
            except Exception:
                return False

    seen: set[str] = set(exclude_canonical or ())
    unique: list[Chem.Mol] = []
    for mol in mols:
        if ctx is not None and not ctx.check_continue():
            break
        try:
            if not _neutral_mol(mol):
                continue
            smi = canonical_smiles(mol)
            if not smi or smi in seen:
                continue
            seen.add(smi)
            unique.append(mol)
        except Exception:
            log_invalid_smiles_bypass(context="deduplicate_molecules")
            continue
    return unique


def predict_process_impurities(
    input_smiles: str,
    condition: str,
    temperature: float,
    *,
    db_fps: list | None = None,
    db_smiles: list[str] | None = None,
    ctx: ProcessRunContext | None = None,
) -> list[dict[str, Any]]:
    """
    핵심 가혹 분해·산화 예측 엔진.

    Parameters
    ----------
    input_smiles : str
        입력 기질 SMILES (주원료 A/B 또는 Reagent).
    condition : str
        가혹 분해 경로 (STRESS_* 상수).
    temperature : float
        공정 온도 (°C) — 신뢰도 가중치 보정에 사용.

    Returns
    -------
    list[dict]
        degradant 레코드 (pathway, smiles, exact_mw, ich_m7, db_match, ...).
    """
    try:
        mol = mol_from_smiles(input_smiles)
        if mol is None:
            return []

        parent_canonical = canonical_smiles(mol)
        if not parent_canonical:
            return []

        exclude = {parent_canonical}
        candidate_mols: list[Chem.Mol] = []
        pathway_label = normalize_stress_condition(condition)
        condition = pathway_label

        if condition == STRESS_NONE:
            return []
        if condition in (STRESS_ACIDIC_HYDROLYSIS, STRESS_ALKALINE_HYDROLYSIS):
            for smarts in HYDROLYSIS_REACTIONS:
                if ctx and not ctx.check_continue():
                    break
                candidate_mols.extend(run_reaction_smarts(mol, smarts))
        elif condition == STRESS_OXIDATIVE:
            for smarts in OXIDATION_REACTIONS:
                if ctx and not ctx.check_continue():
                    break
                candidate_mols.extend(run_reaction_smarts(mol, smarts))
        else:
            return []

        products = deduplicate_molecules(
            candidate_mols, exclude_canonical=exclude, ctx=ctx
        )
        if not products:
            return []

        if db_fps is None or db_smiles is None:
            db_fps, db_smiles = get_verification_fingerprint_db()
        records: list[dict[str, Any]] = []

        for prod in products:
            if ctx and not ctx.check_continue():
                break
            try:
                smi = canonical_smiles(prod)
                if not smi or smi in exclude:
                    continue

                exact_mw = round(Descriptors.ExactMolWt(prod), 2)
                ich_class, ich_detail = classify_ich_m7_genotoxicity(prod)
                match_pct, status_tag, reliability_pct = compute_db_verification_score(
                    prod,
                    db_fps,
                    db_smiles,
                    temperature_c=temperature,
                )

                records.append(
                    {
                        "pathway": pathway_label,
                        "smiles": smi,
                        "exact_mw": exact_mw,
                        "ich_m7": ich_class,
                        "ich_m7_detail": ich_detail,
                        "db_match_pct": match_pct,
                        "db_status_tag": status_tag,
                        "reliability_pct": reliability_pct,
                        "confidence_display": (
                            f"{reliability_pct:.1f}% — {status_tag}"
                        ),
                    }
                )
                if ctx:
                    ctx.register_product()
            except Exception:
                log_invalid_smiles_bypass(context="predict_process_impurities")
                continue

        records.sort(key=lambda r: r["db_match_pct"], reverse=True)
        return ctx.truncate_records(records) if ctx else records
    except Exception as exc:
        log_process_warning(f"predict_process_impurities aborted: {exc}")
        return []


def predict_homo_coupling_impurities(
    input_smiles: str,
    reagent: str,
    temperature: float,
    *,
    db_fps: list | None = None,
    db_smiles: list[str] | None = None,
    ctx: ProcessRunContext | None = None,
) -> list[dict[str, Any]]:
    """
    촉매(Pd/Ni/Cu) 매개 homo-coupling 부반응 유연물질 예측.

    할로겐화 방향족·비닐 할라이드·말단 알킨·보론산 기질에서
    2× 동일 분자 결합 이량체(대칭 dimer)를 RDKit SMARTS로 생성한다.
    """
    if not is_catalyst_reagent(reagent):
        return []

    try:
        mol = mol_from_smiles(input_smiles)
        if mol is None:
            return []

        parent_canonical = canonical_smiles(mol)
        if not parent_canonical:
            return []

        exclude = {parent_canonical}
        candidate_mols: list[Chem.Mol] = []
        metal = catalyst_family_label(reagent)
        pathway_label = f"{metal}-Catalyzed Homo-Coupling (Side Product)"

        for _tag, smarts in HOMO_COUPLING_REACTIONS:
            if ctx and not ctx.check_continue():
                break
            candidate_mols.extend(run_bimolecular_reaction_smarts(mol, smarts))

        products = deduplicate_molecules(
            candidate_mols, exclude_canonical=exclude, ctx=ctx
        )
        if not products:
            return []

        if db_fps is None or db_smiles is None:
            db_fps, db_smiles = get_verification_fingerprint_db()

        records: list[dict[str, Any]] = []
        for prod in products:
            if ctx and not ctx.check_continue():
                break
            try:
                smi = canonical_smiles(prod)
                if not smi or smi in exclude:
                    continue

                exact_mw = round(Descriptors.ExactMolWt(prod), 2)
                ich_class, ich_detail = classify_ich_m7_genotoxicity(prod)
                match_pct, status_tag, reliability_pct = compute_db_verification_score(
                    prod,
                    db_fps,
                    db_smiles,
                    temperature_c=temperature,
                )

                records.append(
                    {
                        "pathway": pathway_label,
                        "smiles": smi,
                        "exact_mw": exact_mw,
                        "ich_m7": ich_class,
                        "ich_m7_detail": ich_detail,
                        "db_match_pct": match_pct,
                        "db_status_tag": status_tag,
                        "reliability_pct": reliability_pct,
                        "confidence_display": f"{reliability_pct:.1f}% — {status_tag}",
                    }
                )
                if ctx:
                    ctx.register_product()
            except Exception:
                log_invalid_smiles_bypass(context="predict_homo_coupling_impurities")
                continue

        records.sort(key=lambda r: r["db_match_pct"], reverse=True)
        return ctx.truncate_records(records) if ctx else records
    except Exception as exc:
        log_process_warning(f"predict_homo_coupling_impurities aborted: {exc}")
        return []


def classify_ich_m7_genotoxicity(mol: Chem.Mol) -> tuple[str, str]:
    """ICH M7 유전독성 구조 경보 기반 휴리스틱 판정."""
    try:
        for smarts, m7_class, detail in M7_ALERT_SMARTS:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None and mol.HasSubstructMatch(pattern):
                return m7_class, detail
    except Exception:
        log_invalid_smiles_bypass(context="classify_ich_m7")
    return M7_DEFAULT_CLASS


def morgan_fp(mol: Chem.Mol):
    """Morgan Fingerprint (radius=2, 2048 bits)."""
    try:
        return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_NBITS)
    except Exception:
        log_invalid_smiles_bypass(context="morgan_fp")
        return None


def compute_db_verification_score(
    mol: Chem.Mol,
    db_fps: list,
    db_smiles: list[str],
    *,
    temperature_c: float,
) -> tuple[float, str, float]:
    """
    10K 내부 검증 DB 대비 Tanimoto 유사도 매칭 + 온도 가중 신뢰도.

    Returns
    -------
    (match_pct, status_tag, reliability_pct)
    """
    query_fp = morgan_fp(mol)
    if query_fp is None or not db_fps:
        base_score = 42.0
    else:
        try:
            sims = DataStructs.BulkTanimotoSimilarity(query_fp, db_fps)
            best_sim = max(sims) if sims else 0.0
            base_score = round(best_sim * 100.0, 2)
        except Exception:
            log_invalid_smiles_bypass(context="compute_db_verification_score")
            base_score = 42.0

    # Linear Weight Adjustment: T > 60°C 시 신뢰도 비례 상승
    reliability = base_score
    if temperature_c > TEMPERATURE_BOOST_THRESHOLD_C:
        boost = min(15.0, (temperature_c - TEMPERATURE_BOOST_THRESHOLD_C) * 0.25)
        reliability = min(100.0, base_score + boost)

    if base_score >= 85.0:
        tag = DB_TAG_VERIFIED
    elif base_score >= 50.0:
        tag = DB_TAG_MECHANISTIC
    else:
        tag = DB_TAG_MODEL

    return base_score, tag, round(reliability, 2)


from engines.p04_cache_paths import CACHE_DIR, VERIFICATION_FP_CACHE_PATH, ensure_cache_dir
from engines.p03_branding import (
    UI_ENGINE_TITLE,
    UI_PROCESS_MODE_TITLE,
    log_process_info,
    log_process_warning,
)
from engines.p05_safeguards import (
    MAX_PROCESS_PRODUCTS,
    PROCESS_REQUEST_TIMEOUT_SEC,
    ProcessRunContext,
    log_invalid_smiles_bypass,
)


def _seed_smiles_corpus() -> list[str]:
    """의약 시드 + 인덱스 기반 체인 SMILES로 10,000건 고유 코퍼스를 빠르게 구축."""
    pharma_seeds = [
        "CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "COc1ccc(CC2SC(=O)NC2=O)cc1",
        "CC(=O)NC1=CC=C(O)C=C1", "c1ccc2[nH]cnc2c1",
        "O=C(O)c1ccccc1O", "CCN(CC)CC", "c1ccncc1", "C1CCOC1",
        "CC(=O)O", "NC(=O)c1ccccc1", "COC(=O)c1ccccc1", "c1ccc(S)cc1",
        "CN(C)C=O", "ClCCl", "CC#N", "O", "c1ccccc1", "Cc1ccccc1",
        "Oc1ccccc1", "Nc1ccccc1", "Fc1ccccc1", "Clc1ccccc1", "Brc1ccccc1",
        "c1ccc(C(=O)O)cc1", "c1ccc(OC)cc1", "C1CCNCC1", "C1COCCO1",
        "CC(C)O", "CCO", "CC(=O)OC", "COC", "CS(C)=O", "c1ccoc1",
        "c1ccsc1", "NCC(=O)O", "O=C(O)C(O)C(O)C(=O)O",
    ]

    seen: set[str] = set()
    corpus: list[str] = []

    def _add(smi: str) -> None:
        if len(corpus) >= VERIFICATION_DB_TARGET_SIZE:
            return
        mol = mol_from_smiles(smi)
        if mol is None:
            return
        canon = canonical_smiles(mol)
        if canon not in seen:
            seen.add(canon)
            corpus.append(canon)

    for smi in pharma_seeds:
        _add(smi)

    # 100×100 비대칭 아미드 체인 격자 = 10,000건 고유 SMILES
    for a in range(1, 101):
        for b in range(1, 101):
            if len(corpus) >= VERIFICATION_DB_TARGET_SIZE:
                break
            _add("C" * a + "C(=O)N" + "C" * b)
        if len(corpus) >= VERIFICATION_DB_TARGET_SIZE:
            break

    return corpus[:VERIFICATION_DB_TARGET_SIZE]


@lru_cache(maxsize=1)
def get_verification_fingerprint_db() -> tuple[list, list[str]]:
    """
    10K 검증 DB Morgan Fingerprint 매트릭스 (프로세스·디스크 캐시).

    Returns
    -------
    (fp_list, smiles_list)
    """
    if VERIFICATION_FP_CACHE_PATH.is_file():
        try:
            with VERIFICATION_FP_CACHE_PATH.open("rb") as fh:
                cached = pickle.load(fh)
            if (
                isinstance(cached, tuple)
                and len(cached) == 2
                and len(cached[1]) >= int(VERIFICATION_DB_TARGET_SIZE * 0.9)
            ):
                return cached
        except Exception:
            pass

    log_process_info("10K verification fingerprint DB building — cache miss or stale")
    smiles_list = _seed_smiles_corpus()
    fps: list = []
    valid_smiles: list[str] = []
    for smi in smiles_list:
        mol = mol_from_smiles(smi)
        if mol is None:
            continue
        fps.append(morgan_fp(mol))
        valid_smiles.append(canonical_smiles(mol))

    payload = (fps, valid_smiles)
    try:
        ensure_cache_dir()
        with VERIFICATION_FP_CACHE_PATH.open("wb") as fh:
            pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
        log_process_info(f"Verification FP cache saved ({len(valid_smiles)} structures)")
    except Exception as exc:
        log_process_warning(f"Verification FP cache write skipped: {exc}")
    return payload


def run_batch_impurity_prediction(
    reactant_entries: list[tuple[str, str, dict[str, Any]]],
    condition: str,
    temperature_c: float,
    *,
    reagent: str = "None",
    process_reagents: dict[str, Any] | None = None,
    ctx: ProcessRunContext | None = None,
) -> list[dict[str, Any]]:
    """
    다중 반응물에 대해 predict_process_impurities 일괄 실행.

    Parameters
    ----------
    reactant_entries : list of (label, key, props_dict)
    condition : 가혹 분해 경로 (STRESS_NONE 이면 stress 시뮬레이션 생략)
    temperature_c : 공정 온도
    reagent : 레거시 단일 시약 C (촉매 homo-coupling)
    process_reagents : 4-카테고리 Registry 컨텍스트 (build_process_reagent_context 출력)
    """
    all_records: list[dict[str, Any]] = []
    global_seen: set[str] = set()
    db_fps, db_smiles = get_verification_fingerprint_db()
    run_ctx = ctx or ProcessRunContext()
    run_ctx.restart_timer()

    reagent_ctx = process_reagents or {}
    homo_catalysts: list[str] = list(reagent_ctx.get("homo_coupling_catalysts") or [])
    if not homo_catalysts and is_catalyst_reagent(reagent):
        homo_catalysts = [reagent]

    for label, _key, props in reactant_entries:
        if not run_ctx.check_continue():
            break
        try:
            if not props.get("valid") or not props.get("smiles"):
                continue
            parent_smi = str(props.get("smiles") or "").strip()
            if not parent_smi:
                continue
            parent_mol = mol_from_smiles(parent_smi)
            if parent_mol is None:
                continue
            parent_canon = canonical_smiles(parent_mol)
            if not parent_canon:
                continue
            global_seen.add(parent_canon)

            preds = predict_process_impurities(
                parent_smi,
                condition,
                temperature_c,
                db_fps=db_fps,
                db_smiles=db_smiles,
                ctx=run_ctx,
            )
            if homo_catalysts:
                for catalyst_name in homo_catalysts:
                    if not run_ctx.check_continue():
                        break
                    preds = preds + predict_homo_coupling_impurities(
                        parent_smi,
                        catalyst_name,
                        temperature_c,
                        db_fps=db_fps,
                        db_smiles=db_smiles,
                        ctx=run_ctx,
                    )
            for rec in preds:
                if not run_ctx.check_continue():
                    break
                smi = rec.get("smiles")
                if not smi or smi in global_seen:
                    continue
                global_seen.add(smi)
                all_records.append({"source_reactant": label, **rec})

            selected_materials = reagent_ctx.get("selected_materials") or {}
            if any(selected_materials.get(cat) for cat in ("Solvents", "Acids", "Bases", "Catalysts")):
                if not run_ctx.check_continue():
                    break
                try:
                    from engines.p13_chemical_registry import build_selected_materials_payload
                    from engines.p34_name_reaction_process import get_advanced_process_engine

                    materials = build_selected_materials_payload(reagent_ctx)
                    raw_preds = get_advanced_process_engine().run_raw_material_simulation(
                        parent_smi,
                        materials,
                        reagent_context=reagent_ctx,
                    )
                    for pred in raw_preds:
                        if not run_ctx.check_continue():
                            break
                        smi = str(pred.get("smiles") or "").strip()
                        if not smi or smi in global_seen:
                            continue
                        mol = mol_from_smiles(smi)
                        if mol is None:
                            continue
                        global_seen.add(smi)
                        origin = str(pred.get("origin_reaction", "RawMaterial"))
                        mechanism = str(pred.get("mechanism", ""))
                        pathway = f"Process: Raw Material · {origin}"
                        if mechanism:
                            pathway = f"{pathway} · {mechanism}"
                        ich_class, ich_detail = classify_ich_m7_genotoxicity(mol)
                        match_pct, status_tag, reliability_pct = compute_db_verification_score(
                            mol, db_fps, db_smiles, temperature_c=temperature_c
                        )
                        all_records.append(
                            {
                                "source_reactant": label,
                                "pathway": pathway,
                                "smiles": smi,
                                "exact_mw": round(Descriptors.ExactMolWt(mol), 2),
                                "ich_m7": ich_class,
                                "ich_m7_detail": ich_detail,
                                "db_match_pct": match_pct,
                                "db_status_tag": status_tag,
                                "reliability_pct": reliability_pct,
                                "confidence_display": f"{reliability_pct:.1f}% — {status_tag}",
                                "name_reaction": origin,
                                "is_name_reaction_byproduct": True,
                                "catalog_mechanism": mechanism,
                            }
                        )
                except Exception as exc:
                    log_process_warning(f"raw material simulation skipped ({label}): {exc}")
        except Exception as exc:
            log_process_warning(f"batch entry skipped ({label}): {exc}")
            log_invalid_smiles_bypass(context=label)
            continue

    return run_ctx.truncate_records(all_records)


def impurity_records_to_dataframe(records: list[dict[str, Any]]) -> "pd.DataFrame":
    """Executive Dashboard 테이블용 DataFrame."""
    import pandas as pd

    if not records:
        return pd.DataFrame(columns=IMPURITY_RESULT_COLUMNS)

    rows = [
        {
            IMP_COL_SOURCE: r.get("source_reactant", ""),
            IMP_COL_IMPURITY_CLASS: r.get("impurity_classification", "Process Impurity"),
            IMP_COL_MECHANISM_RULE: r.get("reaction_name") or r.get("name_reaction") or "",
            IMP_COL_MECHANISM_UI: r.get("mechanism_ui_text", ""),
            IMP_COL_PATHWAY: r["pathway"],
            IMP_COL_SMILES: r["smiles"],
            IMP_COL_MW: r["exact_mw"],
            IMP_COL_ICH_M7_CLASS: _format_ich_m7_class_cell(r),
            IMP_COL_REGULATORY: r.get("regulatory_action_text", ""),
            IMP_COL_ICH_M7: f"{r.get('ich_m7', '')} — {r.get('ich_m7_detail', '')}".strip(" —"),
            IMP_COL_CONFIDENCE: r["confidence_display"],
            "source_reactant": r.get("source_reactant", ""),
            "db_match_pct": r.get("db_match_pct"),
            "reliability_pct": r.get("reliability_pct"),
            "ich_m7_class": r.get("ich_m7_class"),
            "impurity_warn": r.get("impurity_warn", False),
            "impurity_role": r.get("impurity_role"),
        }
        for r in records
    ]
    return pd.DataFrame(rows)


def _format_ich_m7_class_cell(record: dict[str, Any]) -> str:
    cls = record.get("ich_m7_class")
    if cls is None:
        return record.get("ich_m7", "—")
    emoji = ICH_M7_CLASS_EMOJI.get(int(cls), "⚪")
    risk = "High Risk" if int(cls) <= 3 else ("De-risked" if int(cls) == 4 else "Safe")
    return f"{emoji} Class {int(cls)} ({risk})"


# Legacy Korean / raw-record column names → current English export columns
_LEGACY_IMPURITY_COLUMN_MAP: dict[str, str] = {
    "예측 기준 (출발 구조)": IMP_COL_SOURCE,
    "가혹 분해 경로": IMP_COL_PATHWAY,
    "예측 유연물질 / 모핵 잔존체 (SMILES)": IMP_COL_SMILES,
    "예측 유연물질 / 분해 산물 (SMILES)": IMP_COL_SMILES,
    "실시간 분자량 (g/mol)": IMP_COL_MW,
    "정확 분자량 (g/mol)": IMP_COL_MW,
    "ICH M7 유전독성": IMP_COL_ICH_M7,
    "ICH M7 유전독성 평가": IMP_COL_ICH_M7,
    "🎯 신뢰도 & DB 검증": IMP_COL_CONFIDENCE,
    "신뢰도 & DB 검증": IMP_COL_CONFIDENCE,
}


def normalize_prediction_dataframe(pred_df: "pd.DataFrame") -> "pd.DataFrame":
    """
    Normalize prediction DataFrame for dashboard / Excel export.

    Handles legacy Korean column headers and raw engine record field names.
    """
    import pandas as pd

    if pred_df is None or pred_df.empty:
        return pd.DataFrame(columns=IMPURITY_RESULT_COLUMNS)

    df = pred_df.copy()
    rename = {k: v for k, v in _LEGACY_IMPURITY_COLUMN_MAP.items() if k in df.columns}
    if rename:
        df = df.rename(columns=rename)

    fallbacks: tuple[tuple[str, str], ...] = (
        (IMP_COL_SOURCE, "source_reactant"),
        (IMP_COL_MECHANISM_RULE, "reaction_name"),
        (IMP_COL_MECHANISM_UI, "mechanism_ui_text"),
        (IMP_COL_PATHWAY, "pathway"),
        (IMP_COL_SMILES, "smiles"),
        (IMP_COL_MW, "exact_mw"),
        (IMP_COL_ICH_M7, "ich_m7"),
        (IMP_COL_ICH_M7_CLASS, "ich_m7_class"),
        (IMP_COL_REGULATORY, "regulatory_action_text"),
        (IMP_COL_CONFIDENCE, "confidence_display"),
    )
    for target, raw in fallbacks:
        if target not in df.columns and raw in df.columns:
            df[target] = df[raw]

    if IMP_COL_ICH_M7 in df.columns:
        ich = df[IMP_COL_ICH_M7].astype(str)
        if "ich_m7_detail" in df.columns:
            detail = df["ich_m7_detail"].fillna("").astype(str)
            needs_detail = ~ich.str.contains("—", regex=False) & detail.ne("")
            df.loc[needs_detail, IMP_COL_ICH_M7] = ich[needs_detail] + " — " + detail[needs_detail]

    for col in IMPURITY_RESULT_COLUMNS:
        if col not in df.columns:
            df[col] = ""

    return df[IMPURITY_RESULT_COLUMNS + [
        c for c in ("source_reactant", "smiles", "pathway", "exact_mw", "db_match_pct")
        if c in df.columns
    ]]


class ProcessImpurityEngine:
    """PRIME-DS 독립 구동형 공정 유연물질 예측 엔진 파사드."""

    STRESS_NONE = STRESS_NONE
    STRESS_ACIDIC_HYDROLYSIS = STRESS_ACIDIC_HYDROLYSIS
    STRESS_ALKALINE_HYDROLYSIS = STRESS_ALKALINE_HYDROLYSIS
    STRESS_OXIDATIVE = STRESS_OXIDATIVE
    STRESS_CONDITION_OPTIONS = STRESS_CONDITION_OPTIONS
    IMPURITY_RESULT_COLUMNS = IMPURITY_RESULT_COLUMNS
    IMP_COL_SOURCE = IMP_COL_SOURCE
    IMP_COL_PATHWAY = IMP_COL_PATHWAY
    IMP_COL_SMILES = IMP_COL_SMILES
    IMP_COL_MW = IMP_COL_MW
    IMP_COL_ICH_M7 = IMP_COL_ICH_M7
    IMP_COL_CONFIDENCE = IMP_COL_CONFIDENCE
    UI_ENGINE_TITLE = UI_ENGINE_TITLE
    UI_PROCESS_MODE_TITLE = UI_PROCESS_MODE_TITLE
    MAX_PROCESS_PRODUCTS = MAX_PROCESS_PRODUCTS
    PROCESS_REQUEST_TIMEOUT_SEC = PROCESS_REQUEST_TIMEOUT_SEC

    normalize_stress_condition = staticmethod(normalize_stress_condition)
    predict_impurities = staticmethod(predict_process_impurities)
    predict_homo_coupling_impurities = staticmethod(predict_homo_coupling_impurities)
    run_batch_impurity_prediction = staticmethod(run_batch_impurity_prediction)
    is_catalyst_reagent = staticmethod(is_catalyst_reagent)
    impurity_records_to_dataframe = staticmethod(impurity_records_to_dataframe)
    get_verification_fingerprint_db = staticmethod(get_verification_fingerprint_db)
    normalize_prediction_dataframe = staticmethod(normalize_prediction_dataframe)
