"""
니트로사민 예측 — 화학적 타당성 검증 레이어 (ICH M7 / QbD).

1. 알파-수소 유전독성 검증
2. 입체장애 기반 Nitrite 노출 등급 감쇄
"""

# 생성된 니트로사민의 타당성 검사
# ① 알파-수소: 니트로소 질소 옆 탄소에 수소가 있는지 확인. 니트로사민이 대사활성화되려면 그 자리에서 반응이 시작됨.
# 알파-수소 있음-> Class 1~3(위험)
# 알파-수소 없음-> Class 5로 완화
# ② 입체 장애: 주변이 복잡해 반응이 어려우면 아질산염 노출 등급을 한 칸 낮춤.

from __future__ import annotations

import re
from typing import Any

from rdkit import Chem

# 사용자 지정 + N-N=O 토폴로지에 맞는 RDKit 확장 (동일 화학적 의미)
# 2026-09-03 수정: `[CX4H]`는 «수소 정확히 1개»를 뜻해 NDMA·NDEA·NDPA·NMBA·
# N-니트로소모르폴린이 전부 False로 나왔습니다(감사 메타데이터가 «사용자 지정
# 패턴 검출 실패»로 읽힘). H1~H3을 모두 허용하도록 수정.
ALPHA_HYDROGEN_NITROSO_SMARTS_USER = "[NX3](-[CX4;H1,H2,H3])N=O"
ALPHA_HYDROGEN_NITROSO_SMARTS = "[NX3](-[CX4H1,CX4H2,CX4H3])-[NX2]=O"
NITROSO_MOTIF_SMARTS = "[#7]-[#7]=[#8]"

# 2026-09-03 수정: 앵커가 전부 그냥 `[NX3]`이라 니트로사민 질소가 아닌 «아무 질소»의
# 입체장애를 읽었습니다(6개 탄소 떨어진 N-tert-butyl 때문에 아질산염 노출 등급이
# 부당하게 하향). 앵커를 니트로소 질소로 고정하고, 패턴 0·1이 완전히 포섭되던
# `[NX3]-[CX4D4]` 중복과 사슬 4번 패턴을 삭제했습니다.
STERIC_HINDRANCE_SMARTS: tuple[str, ...] = (
    # 니트로소 질소에 직결된 3급/4급 α-탄소 (tert-alkyl)
    "[NX3](-[NX2]=O)-[CX4D4]",
    # 니트로소 질소에 직결된 분지형 2차 α-탄소 (isopropyl 등)
    "[NX3](-[NX2]=O)-[CX4H1](-[CX4])-[CX4]",
    # 아지리딘형 3원 고리 니트로사민
    "[NX3;r3](-[NX2]=O)",
)

_NITRITE_LEVEL_LADDER: tuple[str, ...] = (
    "Negligible",
    "< 0.01 (Low)",
    "0.05 - 1.0 (Medium)",
    "0.1 - 2.0 (Medium-High)",
    "0.1 - 5.0 (High)",
)

_NON_MUTAGENIC_TYPE = "비변이원성 니트로사민 (ICH M7 Class 5 위험도 낮음)"
# 2026-09-03 수정: 니트로사민이 아니거나 파싱되지 않은 구조에 니트로사민 전용
# 알파-수소 판정을 적용하지 않기 위한 표시값.
_NOT_APPLICABLE_HINT = "Not applicable — not a nitrosamine / unparsable"
_ALPHA_H_ABSENT_NOTE = "알파-수소 부재로 인한 DNA 반응성 상실 구조"
_STERIC_KINETICS_NOTE = (
    "Kinetics 제한 요인(Steric Hindrance)으로 인해 실제 생성률은 현저히 낮을 수 있음"
)


def _compile_smarts(pattern: str) -> Chem.Mol | None:
    try:
        return Chem.MolFromSmarts(pattern)
    except Exception:
        return None


_ALPHA_PAT_USER = _compile_smarts(ALPHA_HYDROGEN_NITROSO_SMARTS_USER)
_ALPHA_PAT = _compile_smarts(ALPHA_HYDROGEN_NITROSO_SMARTS)
_NITROSO_PAT = _compile_smarts(NITROSO_MOTIF_SMARTS)
_STERIC_PATS = tuple(p for s in STERIC_HINDRANCE_SMARTS if (p := _compile_smarts(s)))


def has_alpha_hydrogen_on_nitrosamine(smiles: str) -> bool:
    """니트로소기 인접 탄소 알파-수소 존재 여부."""
    return alpha_hydrogen_audit(smiles)["alpha_hydrogen_present"]


def alpha_hydrogen_audit(smiles: str) -> dict[str, bool]:
    """
    알파-수소 검증 감사 정보.
    - parsed: RDKit 파싱 성공 여부
    - is_nitrosamine: N-N=O 모티프 보유 여부
    - user_pattern: 요청 SMARTS `[NX3](-[CX4;H1,H2,H3])N=O`
    - validated_pattern: N-N=O 니트로사민 표준 토폴로지 확장

    2026-09-03 수정: 「파싱 실패」·「니트로소기 없음」·「알파-H 부재」가 모두 같은
    False로 뭉개져 잘못된 SMILES·벤젠·아플라톡신까지 Class 5(비변이원성
    니트로사민)로 강등되던 문제. 세 상태를 분리해 반환합니다.
    """
    mol = Chem.MolFromSmiles((smiles or "").strip())
    out = {
        "parsed": False,
        "is_nitrosamine": False,
        "alpha_hydrogen_user_pattern": False,
        "alpha_hydrogen_validated_pattern": False,
        "alpha_hydrogen_present": False,
    }
    if mol is None:
        return out
    out["parsed"] = True
    if _NITROSO_PAT is None or not mol.HasSubstructMatch(_NITROSO_PAT):
        return out
    out["is_nitrosamine"] = True
    if _ALPHA_PAT_USER and mol.HasSubstructMatch(_ALPHA_PAT_USER):
        out["alpha_hydrogen_user_pattern"] = True
    if _ALPHA_PAT and mol.HasSubstructMatch(_ALPHA_PAT):
        out["alpha_hydrogen_validated_pattern"] = True
    out["alpha_hydrogen_present"] = (
        out["alpha_hydrogen_user_pattern"] or out["alpha_hydrogen_validated_pattern"]
    )
    return out


def has_steric_hindrance_near_nitrosamine_n(smiles: str) -> bool:
    """니트로소(N-N=O) 질소에 직결된 α-탄소의 입체장애 여부."""
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return False
    for pat in _STERIC_PATS:
        try:
            if mol.HasSubstructMatch(pat):
                return True
        except Exception:
            continue
    return False


def _normalize_nitrite_token(level: str) -> str:
    text = (level or "").strip()
    for token in _NITRITE_LEVEL_LADDER:
        if token in text:
            return token
    return text


def demote_nitrite_exposure_level(level: str) -> str:
    """Nitrite ppm 등급을 한 단계 하향."""
    text = (level or "").strip()
    if not text or text == "Negligible":
        return text
    core = _normalize_nitrite_token(text)
    try:
        idx = _NITRITE_LEVEL_LADDER.index(core)
    except ValueError:
        return f"{text} (Steric-adjusted ↓)"
    if idx <= 0:
        return text
    lowered = _NITRITE_LEVEL_LADDER[idx - 1]
    suffix = ""
    m = re.search(r"\(([^)]+)\)\s*$", text)
    if m and m.group(1) not in _NITRITE_LEVEL_LADDER:
        suffix = f" ({m.group(1)} → Steric-adjusted)"
    elif "(Steric-adjusted" not in text:
        suffix = " (Steric-adjusted ↓)"
    return f"{lowered}{suffix}" if suffix else lowered


def apply_chemical_validity_layer(hit: dict[str, Any]) -> dict[str, Any]:
    """
    NDSRI/NDMA hit dict에 규제 신뢰성 검증 메타데이터를 병합.
    부형제 QbD nitrite_exposure와 regulatory_alert와 융합.
    """
    out = dict(hit)
    smi = str(out.get("smiles") or "")
    base_alert = str(out.get("regulatory_alert") or out.get("alert_context") or "")
    notes: list[str] = []

    alpha_audit = alpha_hydrogen_audit(smi)
    alpha_present = alpha_audit["alpha_hydrogen_present"]
    parsed = alpha_audit["parsed"]
    is_nitrosamine = alpha_audit["is_nitrosamine"]
    out["alpha_hydrogen_present"] = alpha_present
    out["alpha_hydrogen_user_pattern"] = alpha_audit["alpha_hydrogen_user_pattern"]
    out["alpha_hydrogen_validated_pattern"] = alpha_audit["alpha_hydrogen_validated_pattern"]
    out["structure_parsed"] = parsed
    out["is_nitrosamine"] = is_nitrosamine

    # 2026-09-03 수정: Class 5 강등은 «파싱 성공 + 니트로사민 + 알파-H 부재»일 때만.
    # 이전에는 파싱 실패·비니트로사민(벤젠·아플라톡신 등)까지 전부
    # Class 5 / NDSRI-Low / 비변이원성 니트로사민으로 찍혔습니다.
    if parsed and is_nitrosamine and not alpha_present:
        out["impurity_type"] = _NON_MUTAGENIC_TYPE
        out["ich_m7_class_hint"] = "Class 5"
        out["risk_class"] = "NDSRI-Low"
        out["nitrosamine_label"] = _NON_MUTAGENIC_TYPE
        notes.append(_ALPHA_H_ABSENT_NOTE)
    elif parsed and is_nitrosamine:
        out["ich_m7_class_hint"] = out.get("ich_m7_class_hint") or "Class 1-3 (alpha-H present)"
    else:
        out["ich_m7_class_hint"] = out.get("ich_m7_class_hint") or _NOT_APPLICABLE_HINT
        out["chemical_validity_scope"] = _NOT_APPLICABLE_HINT

    steric = has_steric_hindrance_near_nitrosamine_n(smi)
    out["steric_hindrance"] = steric

    nitrite_raw = str(out.get("nitrite_exposure_raw") or out.get("nitrite_exposure") or "")
    out["nitrite_exposure_raw"] = nitrite_raw
    nitrite_effective = nitrite_raw
    if steric and nitrite_raw:
        nitrite_effective = demote_nitrite_exposure_level(nitrite_raw)
        out["nitrite_exposure_adjusted"] = nitrite_effective
        out["kinetics_comment"] = _STERIC_KINETICS_NOTE
        notes.append(_STERIC_KINETICS_NOTE)
    else:
        out["nitrite_exposure_adjusted"] = nitrite_effective
        out["kinetics_comment"] = ""

    if nitrite_effective != nitrite_raw and nitrite_raw:
        out["nitrite_exposure"] = nitrite_effective

    if notes:
        merged_alert = " | ".join(filter(None, [base_alert, *notes]))
        out["regulatory_alert"] = merged_alert
        out["alert_context"] = merged_alert
        out["chemical_validity_notes"] = " · ".join(notes)
    else:
        out["chemical_validity_notes"] = ""

    return out
