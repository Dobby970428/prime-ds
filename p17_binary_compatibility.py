"""이원 배합적합성 (binary drug–excipient compatibility) — 1:1 w/w 기준.

전처방(preformulation) 배합적합성 시험은 원료의약품과 부형제를 **한 쌍씩**
섞어 평가합니다. 흔히 1:1 (w/w) 이원 혼합물을 쓰고, 필요에 따라 실제
처방비를 씁니다.

기존 `p72.AdvancedExcipientMatrixEngine.evaluate_formulation_risk`는 선택한
부형제를 **한꺼번에** 넣고 돌립니다. 그 안에서 `seen_smiles` 집합을
부형제 사이에 공유하기 때문에, 같은 불순물을 만드는 부형제가 여러 개면
**먼저 처리된 하나만** 결과에 남습니다.

    아테놀올 + [유당, 덱스트로스, 자당, 옥수수전분]
        한꺼번에 → 4건  (덱스트로스·옥수수전분이 «상호작용 없음»으로 보임)
        한 쌍씩  → 7건  (넷 다 Maillard 이민 생성)

셋 다 환원당이라 같은 이민을 만드는데도 하나만 잡히던 것입니다. 선택
순서만 바꿔도 결과가 달라지므로 위양성이 아니라 **위음성**입니다.

이 모듈은 부형제마다 엔진을 새로 만들어 독립적으로 평가합니다.

2026-08-25 작성.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

__all__ = [
    "BinaryPair",
    "BinaryCompatibilityResult",
    "evaluate_binary_compatibility",
    "MECHANISM_GROUPS",
    "DEFAULT_RATIO",
]

# 이원 혼합물 기본 배합비 (원료의약품 : 부형제, 중량 기준)
DEFAULT_RATIO: tuple[float, float] = (1.0, 1.0)

# 한 쌍에서 같은 반응 규칙이 만드는 산물의 표시 상한.
#
# 폴리올 부형제(소르비톨·수크로스 등)는 수산기가 여러 개라, 에스터를 가진
# 원료와 만나면 위치만 다른 같은 변환이 수십 개씩 나옵니다. 아믈로디핀에서
# 에스터 교환만 657건이 나왔습니다. 정보량은 같으므로 대표 몇 개만 남기고
# 나머지는 건수로 접습니다.
MAX_PRODUCTS_PER_RULE = 3

# 히트맵 세로축 — 반응 기전 묶음
MECHANISM_GROUPS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("maillard", "Maillard / imine", ("maillard", "imine", "schiff")),
    ("formylation", "Formylation", ("formyl",)),
    ("ndsri", "Nitrosamine (NDSRI)", ("ndsri", "nitroso", "nitrosamine")),
    ("ndma", "NDMA", ("ndma",)),
    ("other", "Other", ()),
)


@dataclass
class BinaryPair:
    """원료의약품 + 부형제 한 쌍의 평가 결과."""

    excipient_key: str
    excipient_name: str
    category: str = ""

    hits: list[dict[str, Any]] = field(default_factory=list)
    collapsed: dict[str, int] = field(default_factory=dict)
    error: str | None = None

    # 배합비 정보
    excipient_mw: float | None = None
    mole_ratio: float | None = None          # 부형제 mol / API mol (1:1 w/w 기준)
    ratio_note: str = ""

    # 위험도
    risk_score: float = 0.0
    risk_band: str = "none"                  # high | moderate | low | none
    worst_m7_class: int | None = None        # 범위 힌트는 상한(덜 심각한 쪽)
    # 2026-09-03 추가: `"1-3"`처럼 범위 그대로 보여 확정 등급과 추정 범위를 구분.
    worst_m7_class_display: str = ""
    mechanisms: set[str] = field(default_factory=set)

    @property
    def has_interaction(self) -> bool:
        return bool(self.hits)

    def to_dict(self) -> dict[str, Any]:
        return {
            "excipient_key": self.excipient_key,
            "excipient_name": self.excipient_name,
            "category": self.category,
            "hits": self.hits,
            "collapsed": self.collapsed,
            "collapsed_total": sum(self.collapsed.values()),
            "error": self.error,
            "excipient_mw": self.excipient_mw,
            "mole_ratio": self.mole_ratio,
            "ratio_note": self.ratio_note,
            "risk_score": round(self.risk_score, 1),
            "risk_band": self.risk_band,
            "worst_m7_class": self.worst_m7_class,
            "worst_m7_class_display": self.worst_m7_class_display or (
                str(self.worst_m7_class) if self.worst_m7_class is not None else ""
            ),
            "mechanisms": sorted(self.mechanisms),
            "hit_count": len(self.hits),
            "has_interaction": self.has_interaction,
        }


@dataclass
class BinaryCompatibilityResult:
    api_smiles: str = ""
    api_mw: float | None = None
    api_formula: str = ""
    ratio: tuple[float, float] = DEFAULT_RATIO
    pairs: list[BinaryPair] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @property
    def interacting(self) -> list[BinaryPair]:
        return [p for p in self.pairs if p.has_interaction]

    @property
    def clean(self) -> list[BinaryPair]:
        return [p for p in self.pairs if not p.has_interaction and not p.error]

    def to_dict(self) -> dict[str, Any]:
        return {
            "api_smiles": self.api_smiles,
            "api_mw": self.api_mw,
            "api_formula": self.api_formula,
            "ratio": list(self.ratio),
            "pairs": [p.to_dict() for p in self.pairs],
            "errors": self.errors,
            "total_pairs": len(self.pairs),
            "interacting_pairs": len(self.interacting),
            "clean_pairs": len(self.clean),
        }


# ---------------------------------------------------------------------------
# 배합비
# ---------------------------------------------------------------------------


def _collapse_positional_isomers(
    hits: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """같은 규칙이 만든 위치 이성질체를 대표 몇 개로 접습니다.

    반환: (표시할 산물, {규칙명: 접힌 건수})
    화학적 의미가 같은 변환이 반복되는 것이라 정보 손실은 없습니다.
    """
    grouped: dict[str, list[dict[str, Any]]] = {}
    for hit in hits:
        grouped.setdefault(str(hit.get("rule_name") or "?"), []).append(hit)

    kept: list[dict[str, Any]] = []
    collapsed: dict[str, int] = {}
    for rule, items in grouped.items():
        if len(items) <= MAX_PRODUCTS_PER_RULE:
            kept.extend(items)
            continue
        shown = items[:MAX_PRODUCTS_PER_RULE]
        hidden = len(items) - len(shown)
        collapsed[rule] = hidden
        for hit in shown:
            hit["isomers_collapsed"] = hidden
        kept.extend(shown)
    return kept, collapsed


def _molecular_weight(smiles: str | None) -> float | None:
    text = (smiles or "").strip()
    if not text:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(text)
        if mol is None:
            return None
        return float(Descriptors.MolWt(mol))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# 부형제 분자량의 «종류»
#
# ⚠️ p70/p19의 고분자 부형제 SMILES는 실제 고분자가 아니라 **단량체 대용
#    구조**입니다. 예: Polyethylene_Glycol_4000 = `OCCOCCOCCO`
#    (트리에틸렌글리콜, MW 150.17). 이 값으로 몰비를 내면 PEG 4000의
#    실제 몰비(≈0.079) 대비 26배 과대평가됩니다. 크로스포비돈·MCC처럼
#    애초에 몰수가 정의되지 않는 부형제도 있습니다.
#
#    `_molecular_weight()`가 유효한 값을 돌려주므로 기존 `if not exc_mw`
#    가드는 절대 발동하지 않았습니다. p70을 고칠 수 없는 상황이므로
#    이 모듈에서 방어합니다.
# ---------------------------------------------------------------------------

# 2026-09-03 추가: 이름으로 판별하는 고분자 부형제 표지어.
_POLYMER_NAME_MARKERS: tuple[str, ...] = (
    "polyethylene glycol", "peg", "macrogol", "poloxamer", "povidone",
    "polyvinyl", "pvp", "pva", "polysorbate", "cellulose", "mcc", "avicel",
    "hpmc", "hypromellose", "methylcellulose", "ethylcellulose",
    "carmellose", "carboxymethyl", "croscarmellose", "crospovidone",
    "starch", "dextrin", "maltodextrin", "dextran", "carbomer", "carbopol",
    "eudragit", "shellac", "gelatin", "gum", "alginate", "alginic",
    "chitosan", "carrageenan", "pectin", "xanthan", "guar", "polymer",
    "폴리", "셀룰로오스", "셀룰로스", "전분", "검", "고분자",
)

# 몰수가 정의되지 않는 부형제 종류 (mw_kind 필드 또는 자동 판별 결과)
_UNDEFINED_MW_KINDS = frozenset(
    {"polymer", "polymeric", "inorganic", "mineral", "monomer_surrogate", "surrogate", "undefined"}
)

_RATIO_NOTE_UNDEFINED = (
    "polymer/inorganic — molar ratio undefined; "
    "monomer surrogate MW used for reactivity screening only"
)


def _excipient_mw_kind(entry: dict[str, Any] | None, key: str = "") -> str:
    """부형제 분자량의 종류를 돌려줍니다.

    `mw_kind` 필드가 있으면 그대로 존중하고, 없으면 category·이름·SMILES로
    고분자/무기염 여부를 판별합니다. 반환값: "polymer" | "inorganic" | "molecular".
    """
    entry = entry or {}
    declared = str(entry.get("mw_kind") or "").strip().lower()
    if declared:
        return declared

    blob = " ".join(
        str(entry.get(k) or "")
        for k in ("label_ko", "label", "name", "category")
    ).lower()
    blob = f"{blob} {str(key).replace('_', ' ').lower()}"
    if any(marker in blob for marker in _POLYMER_NAME_MARKERS):
        return "polymer"

    smiles = str(entry.get("smiles") or "").strip()
    if smiles:
        try:
            from rdkit import Chem

            mol = Chem.MolFromSmiles(smiles)
        except Exception:
            mol = None
        # 탄소가 하나도 없으면 무기염·무기 광물 (탈크·이산화규소·인산염 등)
        if mol is not None and not any(a.GetSymbol() == "C" for a in mol.GetAtoms()):
            return "inorganic"

    return "molecular"


def _mole_ratio(
    api_mw: float | None,
    exc_mw: float | None,
    ratio: tuple[float, float],
    *,
    mw_kind: str = "molecular",
) -> tuple[float | None, str]:
    """중량비에서 몰비를 구합니다.

    1:1 (w/w)이면 같은 무게를 섞는다는 뜻이므로, 분자량이 작은 쪽이
    몰 수가 많습니다. 반환값은 «API 1 몰당 부형제 몇 몰».
    """
    # 2026-09-03 수정: 고분자·무기 부형제는 단량체 대용 MW가 유효하게 계산되어
    # 기존 `if not exc_mw` 가드를 그대로 통과했습니다 (PEG 4000 몰비 26배 오차).
    # 분자량의 «종류»를 먼저 확인하고, 정의되지 않으면 몰비를 내지 않습니다.
    if str(mw_kind or "").strip().lower() in _UNDEFINED_MW_KINDS:
        return None, _RATIO_NOTE_UNDEFINED

    if not api_mw or not exc_mw or api_mw <= 0 or exc_mw <= 0:
        # 고분자·무기물처럼 단일 분자량이 없는 부형제가 30종 있습니다.
        # 반응은 reactive_probes로 평가되지만 몰비는 정의되지 않습니다.
        return None, "no defined molecular weight (polymer or inorganic)"

    w_api, w_exc = ratio
    if w_api <= 0:
        return None, "invalid ratio"

    mol_api = w_api / api_mw
    mol_exc = w_exc / exc_mw
    if mol_api <= 0:
        return None, "invalid ratio"

    value = mol_exc / mol_api
    if value >= 1.0:
        note = f"excipient in {value:.2f}-fold molar excess"
    else:
        note = f"excipient at {value:.2f} mol per mol of drug substance"
    return value, note


# ---------------------------------------------------------------------------
# 위험도 점수
#
# 규제 기준이 아니라 «어느 조합을 먼저 검토할지» 정하는 우선순위 지표입니다.
# 계산 근거를 화면에 그대로 보여줄 수 있도록 단순하게 만들었습니다.
# ---------------------------------------------------------------------------

_RISK_CLASS_WEIGHT = {"NDMA": 45.0, "NDSRI": 40.0, "Matrix": 15.0}
# 확정 등급(단일 Class) 가중치
_M7_WEIGHT = {1: 30.0, 2: 25.0, 3: 15.0, 4: 5.0, 5: 0.0}
# 2026-09-03 추가: `"Class 1-3"`처럼 **범위**로만 추정된 힌트는 확정 등급과
# 같은 무게로 더하면 안 됩니다. 상한(덜 심각한 쪽) 기준으로 낮춰 반영합니다.
_M7_RANGE_WEIGHT = {1: 20.0, 2: 18.0, 3: 12.0, 4: 4.0, 5: 0.0}
_NITRITE_WEIGHT = {"high": 12.0, "medium-high": 9.0, "medium": 6.0, "low": 2.0}

BAND_THRESHOLDS = ((60.0, "high"), (30.0, "moderate"), (0.1, "low"))


def _as_bool(value: Any) -> bool | None:
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


# 2026-09-03 수정: 기존 구현은 문자열 안의 **가장 작은 숫자**를 등급으로 읽어
# p74가 넘기는 `"Class 1-3 (alpha-H present)"`(범위 표시)를 확정 Class 1로
# 읽었습니다. Class 1은 정의상 실측 Ames 양성 + 발암성 양성을 요구하므로
# in-silico로 부여될 수 없는데도 `_M7_WEIGHT[1]=30.0`이 붙어 위험구간이
# `high`로 승격되고, 보고서 M7 컬럼과 논문 표의 위험 순위가 전부 여기서 나왔습니다.
_M7_HINT_RE = re.compile(r"class\s*([1-5])(?:\s*[-–—~]\s*([1-5]))?", re.IGNORECASE)


def _m7_class_range_from_hint(raw: Any) -> tuple[int | None, int | None]:
    """`"Class 1-3 (...)"` → (1, 3), `"Class 2"` → (2, 2), 없으면 (None, None)."""
    match = _M7_HINT_RE.search(str(raw or ""))
    if not match:
        return None, None
    low = int(match.group(1))
    high = int(match.group(2)) if match.group(2) else low
    if high < low:
        low, high = high, low
    return low, high


def _m7_class_from_hint(raw: Any) -> int | None:
    """등급 하나로 환산 — 범위이면 **상한(덜 심각한 쪽)** 을 채택합니다."""
    _low, high = _m7_class_range_from_hint(raw)
    return high


def _m7_class_display(raw: Any) -> str:
    """보고서 표시용 — 범위는 `"1-3"` 그대로 남겨 확정과 추정을 구분합니다."""
    low, high = _m7_class_range_from_hint(raw)
    if high is None:
        return ""
    return str(high) if low == high else f"{low}-{high}"


def _nitrite_level(raw: Any) -> str:
    text = str(raw or "").lower()
    if "high" in text and "medium" in text:
        return "medium-high"
    for level in ("high", "medium", "low"):
        if level in text:
            return level
    return ""


def _mechanism_group(hit: dict[str, Any]) -> str:
    blob = " ".join(
        str(hit.get(k) or "")
        for k in ("risk_class", "rule_name", "impurity_type", "mechanism_origin", "nitrosamine_label")
    ).lower()
    for key, _label, markers in MECHANISM_GROUPS:
        if markers and any(m in blob for m in markers):
            return key
    return "other"


def _hit_severity(hit: dict[str, Any]) -> float:
    """불순물 한 건의 심각도 (0~100)."""
    score = _RISK_CLASS_WEIGHT.get(str(hit.get("risk_class") or "Matrix"), 10.0)

    low, high = _m7_class_range_from_hint(hit.get("ich_m7_class_hint"))
    if high is not None:
        # 2026-09-03 수정: 확정 등급과 범위 추정에 서로 다른 가중치를 씁니다.
        table = _M7_WEIGHT if low == high else _M7_RANGE_WEIGHT
        score += table.get(high, 0.0)

    level = _nitrite_level(hit.get("nitrite_exposure") or hit.get("nitrite_exposure_raw"))
    if level:
        score += _NITRITE_WEIGHT.get(level, 0.0)

    # α-수소가 없는 니트로사민은 대사활성화가 어려워 위험이 낮습니다 (p74)
    # 2026-09-03 수정: 니트로사민이 아닌 hit(Maillard 이민 등)은 alpha_hydrogen_present가
    # 의미 없이 False로 찍히므로 감점 대상에서 제외합니다 (p74 `is_nitrosamine` 참조).
    if (
        _as_bool(hit.get("alpha_hydrogen_present")) is False
        and _as_bool(hit.get("is_nitrosamine")) is not False
    ):
        score -= 20.0
    if _as_bool(hit.get("steric_hindrance")) is True:
        score -= 5.0

    return max(0.0, score)


def score_pair(
    hits: list[dict[str, Any]],
    collapsed: dict[str, int] | None = None,
) -> tuple[float, str, int | None, set[str], str]:
    """쌍의 위험 점수(0~100)·구간·최악 등급·기전 집합·등급 표시문자열.

    건수를 그대로 더하면 부형제 하나가 산물을 많이 만든다는 이유만으로
    점수가 수백까지 올라갑니다. 심사자가 보는 방식대로
    «가장 심각한 한 건» 을 기준으로 삼고, 건수는 가산점으로만 반영합니다.

    `collapsed`는 `_collapse_positional_isomers`가 돌려주는 {규칙명: 접힌 건수}
    입니다. 2026-09-03 수정: 이전에는 `hits[:1]`의 `isomers_collapsed`만 더해서,
    첫 hit이 접히지 않은 규칙이면 보정이 0이 되었습니다 — 같은 데이터에서도
    순회 순서에 따라 점수가 달라지는 순서 의존 버그였습니다.
    """
    if not hits:
        return 0.0, "none", None, set(), ""

    severities = [_hit_severity(h) for h in hits]
    worst: int | None = None
    worst_display = ""
    mechanisms: set[str] = set()

    for hit in hits:
        hint = hit.get("ich_m7_class_hint")
        m7 = _m7_class_from_hint(hint)
        if m7 is not None and (worst is None or m7 < worst):
            worst = m7
            worst_display = _m7_class_display(hint)
        mechanisms.add(_mechanism_group(hit))

    score = max(severities)
    n_effective = len(hits) + sum(int(v or 0) for v in (collapsed or {}).values())
    if n_effective > 1:
        score += min(15.0, 3.0 * (n_effective - 1))
    score = max(0.0, min(100.0, score))

    band = "none"
    for threshold, name in BAND_THRESHOLDS:
        if score >= threshold:
            band = name
            break

    return score, band, worst, mechanisms, worst_display


# ---------------------------------------------------------------------------
# 본체
# ---------------------------------------------------------------------------


def evaluate_binary_compatibility(
    api_smiles: str,
    excipient_keys: list[str],
    *,
    ratio: tuple[float, float] = DEFAULT_RATIO,
) -> BinaryCompatibilityResult:
    """원료의약품과 각 부형제를 **한 쌍씩** 평가합니다.

    부형제마다 엔진 인스턴스를 새로 만듭니다. 그래야 `seen_smiles`
    중복 제거가 쌍 안에서만 적용되어, 같은 불순물을 만드는 부형제가
    여럿이어도 전부 잡힙니다.
    """
    result = BinaryCompatibilityResult(
        api_smiles=(api_smiles or "").strip(), ratio=ratio
    )

    if not result.api_smiles:
        result.errors.append("Drug substance structure is empty")
        return result
    if not excipient_keys:
        result.errors.append("No excipient selected")
        return result

    # API 신원
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        api_mol = Chem.MolFromSmiles(result.api_smiles)
        if api_mol is None:
            result.errors.append("RDKit could not parse the drug substance structure")
            return result
        result.api_mw = float(Descriptors.MolWt(api_mol))
        result.api_formula = rdMolDescriptors.CalcMolFormula(api_mol)
    except Exception as exc:
        result.errors.append(f"Structure analysis failed: {exc}")
        return result

    try:
        from p72_excipient_interaction_engine import (
            EXCIPIENT_LIBRARY,
            AdvancedExcipientMatrixEngine,
        )
    except Exception as exc:
        result.errors.append(f"Excipient engine unavailable: {exc}")
        return result

    for key in excipient_keys:
        entry = EXCIPIENT_LIBRARY.get(key) or {}
        pair = BinaryPair(
            excipient_key=key,
            excipient_name=str(entry.get("label_ko") or key),
            category=str(entry.get("category") or ""),
        )

        pair.excipient_mw = _molecular_weight(entry.get("smiles"))
        # 2026-09-03 수정: 단량체 대용 SMILES의 MW로 고분자 몰비를 내지 않도록
        # 분자량의 종류(mw_kind / category·이름·SMILES 판별)를 먼저 확인합니다.
        pair.mole_ratio, pair.ratio_note = _mole_ratio(
            result.api_mw,
            pair.excipient_mw,
            ratio,
            mw_kind=_excipient_mw_kind(entry, key),
        )

        # ⚠️ 부형제마다 엔진을 새로 만듭니다 — 중복 제거를 쌍 안으로 가둡니다
        try:
            engine = AdvancedExcipientMatrixEngine()
            hits = engine.evaluate_formulation_risk(result.api_smiles, [key])
            raw = [h.to_dict() if hasattr(h, "to_dict") else dict(h) for h in (hits or [])]
            pair.hits, pair.collapsed = _collapse_positional_isomers(raw)
        except Exception as exc:
            pair.error = f"{type(exc).__name__}: {exc}"

        # 2026-09-03 수정: 접힘 보정을 pair.collapsed 전체 합계로 계산 (순서 의존 제거)
        (
            pair.risk_score,
            pair.risk_band,
            pair.worst_m7_class,
            pair.mechanisms,
            pair.worst_m7_class_display,
        ) = score_pair(pair.hits, pair.collapsed)

        result.pairs.append(pair)

    # 위험도 높은 순 → 같으면 이름 순
    result.pairs.sort(key=lambda p: (-p.risk_score, p.excipient_name))
    return result


def heatmap(result: BinaryCompatibilityResult) -> dict[str, Any]:
    """부형제 × 반응기전 격자.

    셀 값은 그 기전으로 생긴 불순물 건수입니다.
    """
    groups = [(key, label) for key, label, _ in MECHANISM_GROUPS]
    rows = []
    for pair in result.pairs:
        counts = {key: 0 for key, _ in groups}
        for hit in pair.hits:
            counts[_mechanism_group(hit)] += 1
        rows.append(
            {
                "excipient_key": pair.excipient_key,
                "excipient_name": pair.excipient_name,
                "risk_band": pair.risk_band,
                "counts": counts,
                "total": len(pair.hits),
            }
        )
    return {
        "columns": [{"key": k, "label": lbl} for k, lbl in groups],
        "rows": rows,
        "max_count": max((r["total"] for r in rows), default=0),
    }
