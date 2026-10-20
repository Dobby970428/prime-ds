"""
분해산물(Degradants) 추적 · 공정/보관 취약 구조 패턴 스캔 (RDKit).

목적물(Product) 기반 가수분해·산화·환원 취약성 평가 및
rdChemReactions 분해 경로 예측을 보조합니다.
"""

# 진단용- 화면의 "공정 안정성" 항목
# 패턴 12종을 SMARTS로 훑어서 "이 분자에 약한 고리가 있나" 확인

# 아스피린 예시
# 분자식     C9H8O4
# 분자량     180.16
# LogP       1.31
# 요약       취약 패턴 1건 — 에스터 결합 — 가수분해 취약

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from p30_reaction_engine import _mol_from_smiles, canonical_smiles

# ── 취약 구조 패턴 (SMARTS) — 분해 위험 징후 ──
DEGRADATION_VULNERABILITY: tuple[tuple[str, str, str], ...] = (
    ("ester_labile", "에스터 결합 — 가수분해 취약", "[#6:1](=O)O[#6:2]"),
    ("amide_labile", "아미드 결합 — 산/염기 가수분해 취약", "[#6:1](=O)[#7:2]"),
    # 2026-09-03 수정: 기존 "[#6;R]1[#6;R][#6;R][#6;R][#7;R]1=O"는 고리 질소에 =O가 붙은
    #   (즉 N=O를 가진) 5원 고리를 요구해 2-피롤리디논·δ-락탐·NMP 어느 것도 잡지 못했습니다.
    #   고리 아미드(N-C=O)를 직접 기술하도록 교체 — 고리 크기와 무관하게 락탐을 잡습니다.
    ("lactam_ring", "락탐 고리 — 개환 가수분해", "[NX3;R][CX3;R](=O)"),
    ("carbamate", "카르바메이트 — 가수분해", "[#7][#6](=O)O[#6]"),
    ("imide", "이미드 — 가수분해", "[#6:1](=O)[#7:2][#6:3](=O)"),
    ("sulfonamide", "설폰아미드 — 가수분해", "[#16](=O)(=O)[#7]"),
    # 2026-09-03 수정: 기존 "[#6;R]1[#6;R][#6;R][#6;R]1=O"에는 질소가 없어 β-락탐이 아니라
    #   시클로부타논을 잡았습니다(페니실린 G는 False). 4원 고리 아미드로 교체.
    ("beta_lactam", "β-락탐 — 고리 개환", "[NX3;r4]1[CX3;r4](=O)[CX4;r4][CX4;r4]1"),
    ("aromatic_amine_ox", "방향족 아민 — N-산화 취약", "[c][N;H2,H1]"),
    ("phenol_ox", "페놀 — 산화 취약", "[c][OH]"),
    ("nitro_reducible", "니트로기 — 환원 취약", "[N+](=O)[O-]"),
    ("disulfide", "이황화 결합 — 환원 분해", "[#16][#16]"),
    # 2026-09-03 수정: 기존 "[#6][OH][O][#6]"는 하이드록실 산소에 결합을 3개 요구해
    #   매치 가능한 분자 자체가 존재하지 않았습니다(원자가 위반). 실제 헤미아세탈 골격
    #   C(-OH)(-O-C)로 교체 — 글루코스·락톨에 매치됩니다.
    ("hemiacetal", "헤미아세탈 — 산성 가수분해", "[CX4]([OX2H])[OX2][#6]"),
)

# ── 추가 분해 SMIRKS (rdChemReactions) — 가수분해만 ──
# 2026-09-03 수정: 가수분해 산물을 카복실레이트 음이온([O-])으로 내보내면 중성형으로 저장된
#   PubChem DB canonical SMILES와 절대 일치하지 않아 항상 "이론적"으로 분류됩니다.
#   p30·p32·p40·p41 전부 중성형([OX2H])으로 통일했습니다.
EXTRA_DEGRADATION_SMIRKS: tuple[tuple[str, str, str], ...] = (
    # 2026-09-03 수정: 기존 락탐 SMIRKS는 반응물에 원자 맵이 하나도 없는데 산물이 :1,:2를
    #   참조해, 매치되더라도 기질과 무관한 조각만 나왔습니다(게다가 반응물 SMARTS 자체가
    #   아무것도 잡지 못했습니다). 고리 N-C(=O) 결합만 끊고 나머지 고리는 그대로 이어지는
    #   개환(하나의 산물 성분 `( ... )`)으로 재작성 — 2-피롤리디논 → GABA,
    #   페니실린 G → 페니실로산이 나옵니다.
    (
        "lactam_hydrolysis",
        "락탐 가수분해",
        "[NX3;R:1][CX3;R:2](=[OX1:3])>>([NX3:1].[CX3:2](=[OX1:3])[OX2H])",
    ),
    (
        "carbamate_hydrolysis",
        "카르바메이트 가수분해",
        "[#7:1][#6:2](=O)O[#6:3]>>[#7:1].[#6:2](=O)[OX2H].[#6:3][OH]",
    ),
    (
        "sulfonamide_hydrolysis",
        "설폰아미드 가수분해",
        "[#16:1](=O)(=O)[#7:2]>>[#16:1](=O)(=O)[OX2H].[#7:2]",
    ),
    (
        "imide_hydrolysis",
        "이미드 가수분해",
        "[#6:1](=O)[#7:2][#6:3](=O)>>[#6:1](=O)[OX2H].[#7:2].[#6:3](=O)[OX2H]",
    ),
)


@dataclass
class DegradationRiskProfile:
    """목적물 분해 취약성 프로필."""
    canon_smiles: str
    vulnerabilities: list[str] = field(default_factory=list)
    vulnerability_labels: list[str] = field(default_factory=list)
    molecular_formula: str = ""
    molecular_weight: float = 0.0
    logp: float | None = None
    tpsa: float | None = None
    rotatable_bonds: int = 0
    summary: str = ""


def _safe_scan_pattern(mol: Chem.Mol, smarts: str) -> bool:
    try:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            return False
        return bool(mol.HasSubstructMatch(pat))
    except Exception:
        return False


def scan_degradation_vulnerabilities(smiles: str) -> DegradationRiskProfile:
    """
    목적물 구조에서 분해 취약 패턴을 스캔합니다.
    RDKit 연산 오류 시 빈 프로필을 반환합니다.
    """
    empty = DegradationRiskProfile(canon_smiles="")
    try:
        canon = canonical_smiles(smiles)
        if not canon:
            return empty
        mol = _mol_from_smiles(canon)
        if mol is None:
            return empty

        vuln_ids: list[str] = []
        vuln_labels: list[str] = []
        for vid, label, smarts in DEGRADATION_VULNERABILITY:
            try:
                if _safe_scan_pattern(mol, smarts):
                    vuln_ids.append(vid)
                    vuln_labels.append(label)
            except Exception:
                continue

        try:
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            rot = Descriptors.NumRotatableBonds(mol)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mw = round(Descriptors.MolWt(mol), 2)
        except Exception:
            logp, tpsa, rot, formula, mw = None, None, 0, "", 0.0

        summary = (
            f"취약 패턴 {len(vuln_labels)}건 검출 — {', '.join(vuln_labels[:4])}"
            + (" …" if len(vuln_labels) > 4 else "")
            if vuln_labels
            else "주요 가수분해/산화 취약 패턴 미검출"
        )

        return DegradationRiskProfile(
            canon_smiles=canon,
            vulnerabilities=vuln_ids,
            vulnerability_labels=vuln_labels,
            molecular_formula=formula,
            molecular_weight=mw,
            logp=logp,
            tpsa=tpsa,
            rotatable_bonds=rot,
            summary=summary,
        )
    except Exception:
        return empty
