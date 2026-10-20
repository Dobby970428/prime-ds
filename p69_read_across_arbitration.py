"""
Layer D — 유사도 기반 read-across 중재층 (ICH M7 전문가 검토 트리거).

왜 이 층이 필요한가
--------------------
ICH M7(R2)는 두 가지 상보적 (Q)SAR 방법론(전문가 규칙 · 통계)을 요구하고,
**두 결과가 상충하면 전문가 검토로 넘기라**고 규정합니다. 그런데 이 저장소의
합의 규칙(`p65`)은 Layer B(알럿)와 Layer C(통계)가 엇갈려도 임계값만 보고
등급을 정했습니다. 상충 자체를 «정보»로 쓰지 않았던 것입니다.

이 모듈은 그 빈자리를 채웁니다. 두 층이 엇갈릴 때에 한하여 실측 풀
(Hansen/ISSCAN/Tox21, `p60`의 SQLite)에서 **구조 유사체**를 찾아,
그 유사체들의 실측 Ames 결과를 제3의 근거로 제시합니다.

Layer A 와 무엇이 다른가
------------------------
    Layer A   질의 구조와 **정확히 같은** 구조를 조회 (InChIKey / canonical SMILES)
    Layer D   질의 구조와 **유사한** 구조를 조회 (Morgan ECFP-2048 Tanimoto)

따라서 Layer A 가 적중하면 Layer D 는 동작하지 않습니다. 실측이 있는데
유사체를 볼 이유가 없기 때문입니다. 같은 이유로 유사도 0.999 이상(사실상
동일 구조)은 read-across 근거에서 제외합니다 — 그것은 Layer A 의 일입니다.

판정 원칙 — 자동 하향은 하지 않습니다
--------------------------------------
유사체가 **양성**으로 몰리면 승급(escalation)을 제안합니다. 기존 합의 규칙이
승급 전용(upgrade-only)이므로 결이 같습니다.

유사체가 **음성**으로 몰리는 경우, 이 모듈은 등급을 **내리지 않습니다.**
대신 «전문가 검토 시 하향 가능»을 근거와 함께 기록합니다. 실측 유사체의
음성 근거는 하향의 정당한 사유가 될 수 있으나, ICH M7 에서 그 판단은 사람이
하는 것이지 임계값이 하는 것이 아니기 때문입니다. 시스템이 조용히 내리면
그 판단이 어디서 이루어졌는지 추적할 수 없습니다.

가중 방식
---------
유사체마다 실측 Ames 가 양성/음성이므로, 단순 다수결이 아니라 **유사도 가중**
양성 분율을 씁니다. 0.9 로 닮은 양성 하나가 0.5 로 닮은 음성 하나보다
무겁습니다.

    weighted_positive_fraction = Σ(sim · [양성]) / Σ(sim)

0.60 이상이면 양성 지지, 0.40 이하면 음성 지지, 그 사이는 판단 보류입니다.
경계를 둘로 나눈 것은 «애매한 유사체 집합»이 한쪽으로 반올림되지 않게
하기 위함입니다.

2026-09-08 작성.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# 유사체 검색에 쓰는 최소 Tanimoto — p51 의 기본값과 같습니다.
READ_ACROSS_MIN_SIMILARITY = 0.50

# 이 이상은 사실상 동일 구조로 보고 read-across 근거에서 제외 (Layer A 의 영역)
NEAR_IDENTITY_CUTOFF = 0.999

# 가중 양성 분율의 판정 경계
POSITIVE_SUPPORT_THRESHOLD = 0.60
NEGATIVE_SUPPORT_THRESHOLD = 0.40

# 유사체를 몇 개까지 근거로 삼을지
DEFAULT_TOP_N = 3

VERDICT_SUPPORTS_POSITIVE = "supports_positive"
VERDICT_SUPPORTS_NEGATIVE = "supports_negative"
VERDICT_INCONCLUSIVE = "inconclusive"
VERDICT_NOT_TRIGGERED = "not_triggered"


@dataclass
class ReadAcrossArbitration:
    """Layer B ↔ Layer C 상충에 대한 read-across 중재 결과."""

    triggered: bool = False
    disagreement: str = ""            # "alert_without_qsar" | "qsar_without_alert" | ""
    verdict: str = VERDICT_NOT_TRIGGERED
    analogues: list[dict[str, Any]] = field(default_factory=list)
    weighted_positive_fraction: float | None = None
    n_analogues: int = 0
    suggested_class: int | None = None      # 승급 제안 (하향은 제안만 하고 적용 안 함)
    downgrade_candidate_class: int | None = None
    expert_review_required: bool = False
    rationale: str = ""
    pool_size: int = 0
    pool_degraded: bool = False
    pool_degraded_reason: str = ""
    pool_source: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "Layer D Triggered": self.triggered,
            "Layer D Disagreement": self.disagreement or "None",
            "Layer D Verdict": self.verdict,
            "Layer D Analogues": "; ".join(
                f"{a.get('name', '?')} (Tanimoto {a.get('similarity', 0):.2f}, "
                f"Ames {a.get('ames_raw', 'unknown')})"
                for a in self.analogues
            )
            or "None",
            "Layer D Weighted Positive Fraction": self.weighted_positive_fraction,
            "Layer D Rationale": self.rationale,
            "Layer D Pool Size": self.pool_size,
            "Layer D Pool Degraded": self.pool_degraded,
        }


def detect_disagreement(
    *, alert_present: bool, p_mutagenic: float | None, threshold: float = 0.50
) -> str:
    """Layer B 와 Layer C 가 엇갈리는지 판정합니다.

    반환값
        "alert_without_qsar"  알럿은 있는데 통계 모델은 음성
        "qsar_without_alert"  알럿은 없는데 통계 모델은 양성
        ""                    일치하거나 판정 불가
    """
    if p_mutagenic is None:
        return ""
    qsar_positive = p_mutagenic >= threshold
    if alert_present and not qsar_positive:
        return "alert_without_qsar"
    if (not alert_present) and qsar_positive:
        return "qsar_without_alert"
    return ""


def _weighted_positive_fraction(analogues: list[dict[str, Any]]) -> float | None:
    """유사도 가중 양성 분율. 판정 가능한 유사체가 없으면 None."""
    num = 0.0
    den = 0.0
    for a in analogues:
        raw = str(a.get("ames_raw", "unknown")).lower()
        if raw not in ("positive", "negative"):
            continue
        sim = float(a.get("similarity", 0.0))
        den += sim
        if raw == "positive":
            num += sim
    if den <= 0.0:
        return None
    return num / den


def arbitrate(
    smiles: str,
    *,
    alert_present: bool,
    p_mutagenic: float | None,
    layer_a_hit: bool = False,
    rule_class: int = 5,
    min_similarity: float = READ_ACROSS_MIN_SIMILARITY,
    top_n: int = DEFAULT_TOP_N,
    qsar_threshold: float = 0.50,
) -> ReadAcrossArbitration:
    """Layer B ↔ Layer C 상충 시 read-across 로 중재합니다.

    Layer A 가 적중했거나 두 층이 일치하면 아무것도 하지 않고
    `triggered=False` 로 돌려줍니다 — 불필요한 지문 계산을 피하기 위함입니다.
    """
    out = ReadAcrossArbitration()

    if layer_a_hit:
        out.rationale = "Layer A hit — read-across not required (experimental record available)."
        return out

    disagreement = detect_disagreement(
        alert_present=alert_present, p_mutagenic=p_mutagenic, threshold=qsar_threshold
    )
    if not disagreement:
        return out

    out.triggered = True
    out.disagreement = disagreement

    try:
        from p51_prime_backend_core import ReadAcrossEngine

        envelope = ReadAcrossEngine.find_analogues_with_pool_info(
            smiles, top_n + 2, min_similarity=min_similarity
        )
    except Exception as exc:  # noqa: BLE001
        out.verdict = VERDICT_INCONCLUSIVE
        out.rationale = f"read-across unavailable ({type(exc).__name__}: {exc})"
        return out

    out.pool_size = int(envelope.get("pool_size", 0))
    out.pool_degraded = bool(envelope.get("pool_degraded", False))
    out.pool_degraded_reason = str(envelope.get("pool_degraded_reason", ""))
    out.pool_source = str(envelope.get("pool_source", ""))

    # 사실상 동일 구조는 read-across 근거가 아니라 Layer A 의 일입니다.
    analogues = [
        a
        for a in envelope.get("analogues", [])
        if float(a.get("similarity", 0.0)) < NEAR_IDENTITY_CUTOFF
    ][:top_n]
    out.analogues = analogues
    out.n_analogues = len(analogues)

    if not analogues:
        out.verdict = VERDICT_INCONCLUSIVE
        out.rationale = (
            f"no analogue at Tanimoto ≥ {min_similarity:.2f} in the empirical pool "
            f"(pool size {out.pool_size}"
            + (", pool degraded" if out.pool_degraded else "")
            + ") — read-across inconclusive; the disagreement stands unresolved."
        )
        out.expert_review_required = True
        return out

    frac = _weighted_positive_fraction(analogues)
    out.weighted_positive_fraction = None if frac is None else round(frac, 3)

    if frac is None:
        out.verdict = VERDICT_INCONCLUSIVE
        out.rationale = (
            f"{len(analogues)} analogue(s) found but none carries a usable Ames outcome — "
            "read-across inconclusive."
        )
        out.expert_review_required = True
        return out

    top = analogues[0]
    evidence = (
        f"{len(analogues)} analogue(s), closest {top.get('name', '?')} "
        f"(Tanimoto {float(top.get('similarity', 0.0)):.2f}, Ames {top.get('ames_raw', '?')}); "
        f"similarity-weighted positive fraction {frac:.2f}"
    )

    if frac >= POSITIVE_SUPPORT_THRESHOLD:
        out.verdict = VERDICT_SUPPORTS_POSITIVE
        # 승급 제안 — 기존 합의 규칙이 승급 전용이므로 결이 같습니다.
        out.suggested_class = 3 if rule_class >= 5 else min(rule_class, 3)
        out.rationale = (
            f"read-across supports the mutagenic call — {evidence}. "
            f"Escalation to Class {out.suggested_class} proposed "
            f"({'structural alert without a positive (Q)SAR call' if disagreement == 'alert_without_qsar' else 'positive (Q)SAR call without a structural alert'})."
        )
        return out

    if frac <= NEGATIVE_SUPPORT_THRESHOLD:
        out.verdict = VERDICT_SUPPORTS_NEGATIVE
        # 자동 하향은 하지 않습니다 — 근거만 제시하고 사람에게 넘깁니다.
        out.downgrade_candidate_class = 5
        out.expert_review_required = True
        out.rationale = (
            f"read-across supports the non-mutagenic call — {evidence}. "
            "The class is NOT lowered automatically; expert review may reclassify to "
            "Class 5 on this evidence (ICH M7 places that judgement with the reviewer, "
            "not with a threshold)."
        )
        return out

    out.verdict = VERDICT_INCONCLUSIVE
    out.expert_review_required = True
    out.rationale = (
        f"read-across inconclusive — {evidence}, between the "
        f"{NEGATIVE_SUPPORT_THRESHOLD:.2f} and {POSITIVE_SUPPORT_THRESHOLD:.2f} decision bounds. "
        "The Layer B / Layer C disagreement stands and is referred to expert review."
    )
    return out


__all__ = [
    "READ_ACROSS_MIN_SIMILARITY",
    "NEAR_IDENTITY_CUTOFF",
    "POSITIVE_SUPPORT_THRESHOLD",
    "NEGATIVE_SUPPORT_THRESHOLD",
    "VERDICT_SUPPORTS_POSITIVE",
    "VERDICT_SUPPORTS_NEGATIVE",
    "VERDICT_INCONCLUSIVE",
    "VERDICT_NOT_TRIGGERED",
    "ReadAcrossArbitration",
    "detect_disagreement",
    "arbitrate",
]
