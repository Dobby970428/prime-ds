"""
PRIME-DS RDKit (Q)SAR proxy — OECD QSAR Toolbox / AMBIT2 미연결 시 로컬 예측.

AMBIT REST API 대신 RDKit 물성·Ashby-Tennant 알러트·ICH Q1A 가수분해 규칙으로
hydrolysis / biodegradability / mutagenicity alert 값을 산출합니다.
"""

# 서버 없을 때 RDKit으로 대신 계산

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from p50_physicochemical import calculate_opensource_properties, predict_chemical_stability
from p62_toxicity_alerts import scan_ashby_tennant_alerts

_HALOGEN_SMARTS = "[F,Cl,Br,I]"
_NITRO_AROM_SMARTS = "[c][N+](=O)[O-]"
_ESTER_SMARTS = "[CX3](=O)[OX2H0]"
_AMIDE_SMARTS = "[CX3](=O)[NX3]"


def _count_halogen_atoms(mol: Chem.Mol) -> int:
    pat = Chem.MolFromSmarts(_HALOGEN_SMARTS)
    if pat is None:
        return 0
    return len(mol.GetSubstructMatches(pat))


def predict_hydrolysis_proxy(smiles: str) -> str:
    """OECD hydrolysis stability — RDKit / ICH Q1A-style labile group rules."""
    # 2026-09-03 수정: 파싱 불가 SMILES를 진입부에서 걸러냅니다.
    #   전에는 Biodegradability·OECD_Alert만 "N/A — invalid SMILES"를 내고
    #   가수분해만 "Stable to hydrolysis under typical aqueous conditions"
    #   라는 **긍정 판정**이 나갔습니다 (예: rdkit_oecd_qsar_proxy("Q@@bad")).
    #   빈 문자열은 RDKit이 원자 0개짜리 유효 분자로 읽으므로 함께 봅니다.
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None or mol.GetNumAtoms() == 0:
        return "N/A — invalid SMILES (RDKit parse failure)"

    props = calculate_opensource_properties(smiles)
    status, detail = predict_chemical_stability(smiles, props)
    low = (status + " " + detail).lower()
    if status == "Unknown":
        return f"N/A — {detail} (RDKit proxy)"
    if "unstable" in low or "불안정" in status:
        return (
            "High hydrolysis susceptibility — labile ester/carbamate-type functionality "
            "(RDKit proxy; confirm per ICH Q1A stress study)"
        )
    if "moderate" in low or "보통" in status:
        return (
            "Moderate hydrolysis risk — amide or acid-labile groups under pH/temperature stress "
            "(RDKit proxy)"
        )
    ester = Chem.MolFromSmarts(_ESTER_SMARTS)
    amide = Chem.MolFromSmarts(_AMIDE_SMARTS)
    if ester and mol.HasSubstructMatch(ester):
        return "Ester present — monitor acidic/basic hydrolysis (RDKit substructure proxy)"
    if amide and mol.HasSubstructMatch(amide):
        return "Amide present — limited hydrolysis under mild conditions (RDKit proxy)"
    return "Stable to hydrolysis under typical aqueous conditions (RDKit proxy)"


def predict_biodegradability_proxy(smiles: str) -> str:
    """OECD 301-style biodegradability read-across proxy from RDKit descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "N/A — invalid SMILES (RDKit parse failure)"

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    halogen = _count_halogen_atoms(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)

    score = 0
    if mw <= 350:
        score += 1
    if logp <= 3.0:
        score += 1
    if halogen == 0:
        score += 1
    if arom <= 1:
        score += 1
    if tpsa >= 40:
        score += 1
    if rot >= 2:
        score += 1

    nitro = Chem.MolFromSmarts(_NITRO_AROM_SMARTS)
    if nitro and mol.HasSubstructMatch(nitro):
        return "Not readily biodegradable — aromatic nitro / persistent motif (RDKit proxy)"

    if halogen >= 3 or (mw > 500 and logp > 4):
        return "Not readily biodegradable — high MW/lipophilicity or multi-halogenated (RDKit proxy)"
    if score >= 5:
        return "Readily biodegradable — favorable MW/LogP/TPSA profile (OECD 301-like RDKit proxy)"
    if score >= 3:
        return "Inherently biodegradable — moderate structural complexity (RDKit proxy)"
    return "Not readily biodegradable — persistent aromatic/lipophilic profile (RDKit proxy)"


def predict_mutagenicity_alert_proxy(smiles: str) -> str:
    """OECD M7 / Ashby-Tennant structural alert proxy."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "N/A — invalid SMILES (RDKit parse failure)"
    alerts, _ = scan_ashby_tennant_alerts(mol)
    if alerts:
        shown = ", ".join(alerts[:4])
        extra = f" (+{len(alerts) - 4} more)" if len(alerts) > 4 else ""
        return f"Alert detected: {shown}{extra} (Ashby-Tennant / ICH M7 proxy)"
    return "None detected (Ashby-Tennant screening negative — RDKit proxy)"


def rdkit_oecd_qsar_proxy(smiles: str) -> dict[str, str]:
    """Full OECD QSAR result bundle for offline / fallback mode."""
    smi = (smiles or "").strip()
    return {
        "Hydrolysis_Stability": predict_hydrolysis_proxy(smi),
        "Biodegradability": predict_biodegradability_proxy(smi),
        "OECD_Alert_Found": predict_mutagenicity_alert_proxy(smi),
        "source": "rdkit_proxy",
    }
