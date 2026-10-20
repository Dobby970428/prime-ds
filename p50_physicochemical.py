"""RDKit 기반 물성 계산 · 공정 안정성 예측 · 자체 알고리즘 플러그인 구역."""

# RDKit의 표준 계산 그대로 사용
# RDKit: 파이썬 라이브러리를 활용하는 분자식 표현 프로그램
# 뒤에서 p51을 호출해 pKa·LogD·끓는점을 덧붙임.
# LogP: 기름과 물에 나뉘어 녹는 정도(Crippen 방법)]
# TPSA: 분자 표면 중 극성 부분의 넓이. 세포막 투과성 지표
# H-Bind Donor/Acceptor: 수소결합을 주는/받는 자리 수

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors

_ESTER_SMARTS = "[CX3](=O)[OX2H0]"
_AMIDE_SMARTS = "[CX3](=O)[NX3]"


def calculate_opensource_properties(smiles: str) -> dict[str, float | int] | None:
    """RDKit 오픈소스 기반 물성 (baseline)."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return {
            "LogP": round(Descriptors.MolLogP(mol), 2),
            "MW": round(Descriptors.MolWt(mol), 2),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "H-Bond Donors": Descriptors.NumHDonors(mol),
            "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        }
    except Exception:
        return None


def calculate_my_advanced_properties(
    smiles: str,
    *,
    skip_deepchem: bool = False,
) -> dict[str, str | float] | None:
    """
    연구원 고유 물성 알고리즘 이식 구역 (logD, pKa, bp 등).
    DeepChem Tox21 기반 Ames 유전독성 AI 예측 포함.
    skip_deepchem=True: 규제 보고서 생성 등 빠른 경로 (OECD/AMBIT만 조회).
    """
    baseline = calculate_opensource_properties(smiles)
    if baseline is None:
        return None

    result: dict[str, str | float] = {
        "LogP (My Algo)": baseline["LogP"],
        "LogD (My Algo)": "추후 업데이트 예정",
        "pKa (My Algo)": "추후 업데이트 예정",
        "Boiling Point": "추후 업데이트 예정",
    }

    try:
        from p51_prime_backend_core import HybridPropertyEngine

        hybrid = HybridPropertyEngine.calculate_properties(smiles)
        if hybrid:
            result["LogP (My Algo)"] = hybrid["LogP_Baseline"]
            result["LogD (My Algo)"] = hybrid["LogD_pH7.4"]
            result["pKa (My Algo)"] = hybrid["Estimated_pKa"]
            result["Boiling Point"] = f"{hybrid['Boiling_Point_C']} °C"
            result["TPSA (Hybrid)"] = hybrid["TPSA"]
    except Exception:
        pass

    try:
        from p52_ambit_client import query_ambit_qsar_properties

        ambit = query_ambit_qsar_properties(smiles)
        result["Hydrolysis (OECD/AMBIT)"] = ambit.get("Hydrolysis_Stability", "N/A")
        result["Biodegradability (OECD/AMBIT)"] = ambit.get("Biodegradability", "N/A")
        result["OECD Mutagenicity Alert"] = ambit.get("OECD_Alert_Found", "N/A")
        if ambit.get("error"):
            result["OECD Status"] = f"Unavailable — {ambit['error']}"
        elif ambit.get("source") == "rdkit_proxy":
            result["OECD Status"] = "RDKit (Q)SAR proxy (AMBIT offline)"
        elif ambit.get("source") == "AMBIT2":
            result["OECD Status"] = "Connected (AMBIT2)"
        else:
            result["OECD Status"] = "Offline (AMBIT server not reachable)"
    except Exception as exc:
        result["Hydrolysis (OECD/AMBIT)"] = "Unavailable"
        result["Biodegradability (OECD/AMBIT)"] = "Unavailable"
        result["OECD Mutagenicity Alert"] = "Unavailable"
        result["OECD Status"] = f"Error — {exc}"

    try:
        from p63_deepchem_predictor import predict_single_smiles

        if not skip_deepchem:
            dc_pred = predict_single_smiles(smiles)
            if dc_pred:
                prob = dc_pred["Ames_Mutagenicity_Prob"]
                result["Ames (DeepChem AI)"] = round(float(prob) * 100, 1)
                result["Ames Class (DeepChem)"] = dc_pred["Prediction_Class"]
    except Exception:
        if not skip_deepchem:
            result["Ames (DeepChem AI)"] = "예측 불가"
            result["Ames Class (DeepChem)"] = "Unavailable"

    return result


def predict_chemical_stability(
    smiles: str,
    properties: dict | None = None,
) -> tuple[str, str]:
    """구조·LogP 기반 공정 안정성(가수분해 등) 규칙 예측."""
    # 2026-09-03 수정: 파싱 불가 SMILES가 «✅ 안정»으로 보고되던 문제.
    #   mol이 None이면 has_ester·has_amide가 둘 다 False가 되어 마지막
    #   분기로 떨어져 "주요 가수분해·산화 취약 구조가 발견되지 않았으며"
    #   라는 **긍정 판정**이 나갔습니다. 구조를 읽지도 못한 상태에서
    #   안정하다고 말하는 셈입니다. (p53의 Biodegradability·OECD_Alert는
    #   이미 "N/A — invalid SMILES"를 올바르게 돌려주고 있었습니다.)
    #   진입부에서 먼저 걸러냅니다.
    #   빈 문자열은 RDKit이 «원자 0개짜리 유효한 분자»로 읽으므로
    #   None 검사만으로는 걸러지지 않습니다. 원자 수도 함께 봅니다.
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None or mol.GetNumAtoms() == 0:
        return "Unknown", "SMILES 해석 실패"

    props = properties or calculate_opensource_properties(smiles)
    if not props:
        return "Unknown", "물성 데이터 부족"

    try:
        ester_pat = Chem.MolFromSmarts(_ESTER_SMARTS)
        amide_pat = Chem.MolFromSmarts(_AMIDE_SMARTS)
        has_ester = bool(mol and ester_pat and mol.HasSubstructMatch(ester_pat))
        has_amide = bool(mol and amide_pat and mol.HasSubstructMatch(amide_pat))
        logp = float(props.get("LogP", 0))

        if has_ester and logp < 1.5:
            return (
                "⚠️ 불안정 (Unstable)",
                "수용액 공정 조건에서 친수성 에스테르 결합의 급격한 가수분해(Hydrolysis) 리스크가 높음.",
            )
        if has_amide:
            return (
                "ℹ️ 보통 (Moderately Stable)",
                "강산/강염기 혹은 고온 조건에서 제한적인 아미드 가수분해 가능성 있음.",
            )
        return (
            "✅ 안정 (Stable)",
            "주요 가수분해·산화 취약 구조가 발견되지 않았으며, 화학적 거동이 안정할 것으로 예측됨.",
        )
    except Exception:
        return "Unknown", "안정성 연산 오류"


def stability_summary_for_table(smiles: str) -> tuple[str, str | float]:
    """테이블용 LogP · 안정성 등급."""
    props = calculate_opensource_properties(smiles)
    status, _ = predict_chemical_stability(smiles, props)
    logp = props["LogP"] if props else "N/A"
    return status, logp
