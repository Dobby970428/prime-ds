"""
학습 endpoint ablation — 쌍체 검정 (표 2 비교 ③).

왜 있는가
    s11 의 (A-6) 표 2 는 report_ames.txt · report_tox21.txt 의 **집계
    혼동행렬만** 읽는다. 그래서 비율 지표의 Wilson CI 는 내지만
    «Ames판 대 Tox21판» 의 쌍체 McNemar 검정은 산출하지 못하고,
    다음과 같이 자리표시자를 남겨 두었다.

        쌍체 차이의 CI 는 화합물별 예측이 없어 산출하지 않았다
        — 지어내지 않고 자리표시자로 둔다.

    2026-09-15 에 s5 를 PRIME_LAYER_C_MODEL=tox21 로 재실행해 화합물별
    예측 CSV 를 보존했다. 이 스크립트가 그 빈칸을 채운다.

무엇을 내는가
    ① McNemar 정확검정 (이항검정, Holm 보정)
       Ames판과 Tox21판이 «같은 화합물을 맞히느냐» 를 쌍으로 비교한다.
    ② 민감도·특이도 쌍체 차이의 부트스트랩 95% CI
       s11 과 같은 조건 (2,000회 · seed 20260903).
    ③ 정합성 점검 — Layer B 는 두 판이 동일해야 한다.
       규칙 계층은 Layer C 모델과 무관하므로 b=0 · c=0 이 정상이다.
       0 이 아니면 두 CSV 가 다른 검증셋에서 나온 것이다.

입력
    data/pipeline_validation_results.csv         (Ames 판, s5)
    data/pipeline_validation_results_tox21.csv   (Tox21 판, s5 --tag tox21)

출력
    논문자료/02_검증결과/ablation_paired_tests_report.txt
    논문자료/02_검증결과/ablation_paired_tests.json

사용법
    py s18_ablation_paired_tests.py
    py s18_ablation_paired_tests.py --ames 다른.csv --tox21 다른.csv
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

APP = Path(__file__).resolve().parent
DATA = APP / "data"
OUT_DIR = APP / "논문자료" / "02_검증결과"
OUT_TXT = OUT_DIR / "ablation_paired_tests_report.txt"
OUT_JSON = OUT_DIR / "ablation_paired_tests.json"

SEED = 20260903          # s11 과 동일
BOOTSTRAP_N = 2000       # s11 과 동일

# (열 이름, 표시 이름) — 표 2 의 세 방식
METHODS = (
    ("pred_layer_c", "Layer C 단독 (통계)"),
    ("pred_ensemble", "4층 앙상블"),
    ("pred_layer_b", "Layer B 단독 (전문가 규칙)"),
)

_log: list[str] = []


def w(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


# ---------------------------------------------------------------------------
# 통계
# ---------------------------------------------------------------------------
def binom_cdf(k: int, n: int, p: float = 0.5) -> float:
    """P(X <= k), X ~ Bin(n, p). scipy 없이 계산한다."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    return sum(math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
               for i in range(k + 1))


def mcnemar_exact(b: int, c: int) -> float:
    """
    McNemar 정확검정 (양측).

    불일치쌍이 적을 수 있으므로 χ² 근사가 아닌 이항검정을 쓴다.
    귀무가설에서 b ~ Bin(b+c, 0.5) 이다. s11 · s12 · s14 와 같은 정의.
    """
    n = b + c
    if n == 0:
        return 1.0
    return min(1.0, 2.0 * binom_cdf(min(b, c), n, 0.5))


def holm(pvals: list[float]) -> list[float]:
    """Holm-Bonferroni 보정. 원래 순서로 돌려준다."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    out = [0.0] * m
    running = 0.0
    for rank, idx in enumerate(order):
        adj = min(1.0, (m - rank) * pvals[idx])
        running = max(running, adj)      # 단조 증가 강제
        out[idx] = running
    return out


def rates(y: list[int], p: list[int]) -> dict:
    tp = sum(1 for a, b_ in zip(y, p) if a == 1 and b_ == 1)
    tn = sum(1 for a, b_ in zip(y, p) if a == 0 and b_ == 0)
    fp = sum(1 for a, b_ in zip(y, p) if a == 0 and b_ == 1)
    fn = sum(1 for a, b_ in zip(y, p) if a == 1 and b_ == 0)
    return {
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "sens": tp / (tp + fn) if tp + fn else float("nan"),
        "spec": tn / (tn + fp) if tn + fp else float("nan"),
    }


def bootstrap_diff(y: list[int], pa: list[int], pb: list[int],
                   metric: str) -> tuple[float, float, float]:
    """
    (A − B) 쌍체 차이의 부트스트랩 백분위 95% CI.

    화합물 단위로 재표집하므로 두 방식이 같은 표본을 보고, 쌍 구조가
    유지된다. s11 의 Matthews 계수 CI 와 같은 조건이다.
    """
    rng = random.Random(SEED)
    n = len(y)
    idx_all = range(n)
    point = rates(y, pa)[metric] - rates(y, pb)[metric]

    diffs: list[float] = []
    for _ in range(BOOTSTRAP_N):
        idx = [rng.randrange(n) for _ in idx_all]
        yy = [y[i] for i in idx]
        aa = [pa[i] for i in idx]
        bb = [pb[i] for i in idx]
        da = rates(yy, aa)[metric]
        db = rates(yy, bb)[metric]
        if da == da and db == db:        # NaN 제외
            diffs.append(da - db)
    if not diffs:
        return (point, float("nan"), float("nan"))
    diffs.sort()
    lo = diffs[int(0.025 * len(diffs))]
    hi = diffs[min(len(diffs) - 1, int(0.975 * len(diffs)))]
    return (point, lo, hi)


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ames", default=str(DATA / "pipeline_validation_results.csv"))
    ap.add_argument("--tox21", default=str(DATA / "pipeline_validation_results_tox21.csv"))
    args = ap.parse_args()

    import pandas as pd

    p_ames = Path(args.ames)
    p_tox = Path(args.tox21)
    for p in (p_ames, p_tox):
        if not p.exists():
            w(f"[오류] 없는 파일: {p}")
            if p is p_tox:
                w("       아래를 한 줄로 실행해 만드십시오 (PowerShell).")
                w('       $env:PRIME_LAYER_C_MODEL="tox21"; '
                  "py s5_run_pipeline_validation.py --tag tox21; "
                  "Remove-Item Env:\\PRIME_LAYER_C_MODEL")
            return 1

    a = pd.read_csv(p_ames)
    t = pd.read_csv(p_tox)

    w("=" * 78)
    w("  학습 endpoint ablation — 쌍체 검정 (표 2 비교 ③)")
    w("=" * 78)
    w()
    w(f"  Ames 판  : {p_ames.name}  {len(a):,}행")
    w(f"  Tox21 판 : {p_tox.name}  {len(t):,}행")
    w()

    cols = ["smiles", "actual_ames"] + [c for c, _ in METHODS]
    missing = [c for c in cols if c not in a.columns or c not in t.columns]
    if missing:
        w(f"[오류] 필요한 열이 없습니다: {missing}")
        return 1

    m = a[cols].merge(t[cols], on=["smiles", "actual_ames"],
                      how="inner", suffixes=("_ames", "_tox21"),
                      validate="one_to_one")
    w(f"  화합물 단위 결합 : {len(m):,}건")
    if len(m) != len(a) or len(m) != len(t):
        w(f"  ⚠ 결합에서 빠진 행이 있다 (Ames {len(a) - len(m):,} · "
          f"Tox21 {len(t) - len(m):,}). 두 판이 같은 검증셋인지 확인할 것.")
    w(f"  실제 양성 {int((m['actual_ames'] == 1).sum()):,} · "
      f"음성 {int((m['actual_ames'] == 0).sum()):,}")
    w()

    y = m["actual_ames"].astype(int).tolist()

    # ── 먼저 McNemar 를 모두 모아 Holm 을 한 번에 적용한다 ──────────────
    raw: list[tuple[str, int, int, float]] = []
    for col, label in METHODS:
        pa = m[f"{col}_ames"].astype(int).tolist()
        pb = m[f"{col}_tox21"].astype(int).tolist()
        # b = Ames 판만 맞힌 건수, c = Tox21 판만 맞힌 건수
        b = sum(1 for yy, x, z in zip(y, pa, pb) if x == yy and z != yy)
        c = sum(1 for yy, x, z in zip(y, pa, pb) if z == yy and x != yy)
        raw.append((label, b, c, mcnemar_exact(b, c)))

    adj = holm([r[3] for r in raw])

    w("-" * 78)
    w("  ① McNemar 정확검정 — Ames 판 대 Tox21 판 (Holm 보정, 3 비교)")
    w("-" * 78)
    w()
    w(f"  {'방식':<28} {'b':>6} {'c':>6} {'p':>12} {'p(Holm)':>12}")
    w("  " + "-" * 68)
    for (label, b, c, p), pa_ in zip(raw, adj):
        w(f"  {label:<28} {b:6,} {c:6,} {p:12.4g} {pa_:12.4g}")
    w()
    w("  b = Ames 판만 맞힌 건수, c = Tox21 판만 맞힌 건수.")
    w()

    # ── 정합성 점검 ─────────────────────────────────────────────────────
    lb = next(r for r in raw if r[0].startswith("Layer B"))
    w("-" * 78)
    w("  ③ 정합성 점검 — Layer B 는 두 판이 동일해야 한다")
    w("-" * 78)
    w()
    if lb[1] == 0 and lb[2] == 0:
        w("  ✅ b = 0 · c = 0 — 규칙 계층 판정이 두 판에서 완전히 같다.")
        w("     Layer B 는 Layer C 모델과 무관하므로 이것이 정상이다.")
        w("     두 CSV 가 같은 검증셋에서 나왔다는 근거이기도 하다.")
    else:
        w(f"  ⚠ b = {lb[1]} · c = {lb[2]} — 규칙 계층이 두 판에서 다르다.")
        w("     Layer B 는 Layer C 모델에 의존하지 않으므로 이럴 수 없다.")
        w("     두 CSV 가 서로 다른 검증셋·다른 경보 집합에서 나온 것이다.")
        w("     **아래 수치를 쓰지 말고 원인을 먼저 찾을 것.**")
    w()

    # ── 쌍체 차이의 부트스트랩 CI ───────────────────────────────────────
    w("-" * 78)
    w("  ② 민감도·특이도 쌍체 차이 (Ames − Tox21)")
    w("-" * 78)
    w()
    w(f"  부트스트랩 백분위 95% CI · {BOOTSTRAP_N:,}회 · seed {SEED}")
    w("  (화합물 단위 재표집 — 두 방식이 같은 표본을 보므로 쌍 구조가 유지된다)")
    w()

    j_methods = {}
    for col, label in METHODS:
        pa = m[f"{col}_ames"].astype(int).tolist()
        pb = m[f"{col}_tox21"].astype(int).tolist()
        ra, rb = rates(y, pa), rates(y, pb)

        ds, dslo, dshi = bootstrap_diff(y, pa, pb, "sens")
        dp, dplo, dphi = bootstrap_diff(y, pa, pb, "spec")

        w(f"  ── {label}")
        w(f"     민감도  Ames {ra['sens'] * 100:5.1f}%  Tox21 {rb['sens'] * 100:5.1f}%"
          f"   차이 {ds * 100:+6.1f}%p  [{dslo * 100:+.1f} – {dshi * 100:+.1f}]")
        w(f"     특이도  Ames {ra['spec'] * 100:5.1f}%  Tox21 {rb['spec'] * 100:5.1f}%"
          f"   차이 {dp * 100:+6.1f}%p  [{dplo * 100:+.1f} – {dphi * 100:+.1f}]")
        sig = "0 을 포함하지 않는다" if (dslo > 0 or dshi < 0) else "0 을 포함한다"
        w(f"     → 민감도 차이의 CI 는 {sig}.")
        w()

        j_methods[col] = {
            "label": label,
            "ames": ra, "tox21": rb,
            "sens_diff": {"point": ds, "lo": dslo, "hi": dshi},
            "spec_diff": {"point": dp, "lo": dplo, "hi": dphi},
        }

    w("-" * 78)
    w("  해석 메모")
    w("-" * 78)
    w()
    w("  · 이 표는 «머신러닝을 썼는가» 가 아니라 «무엇으로 학습했는가» 가")
    w("    결정한다는 것을 보인다. 알고리즘(RandomForest) · 기술자(ECFP4)")
    w("    · 검증셋이 모두 같고 학습 endpoint 만 다르다.")
    w("  · 앙상블 쪽 차이가 Layer C 쪽보다 작다면, 상향 전용 합의 규칙에서")
    w("    규칙 계층이 등급을 끌고 가기 때문이다. 통계 계층이 기여하지")
    w("    못하면 4계층이 규칙 계층 하나로 주저앉는다.")
    w()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_TXT.write_text("\n".join(_log), encoding="utf-8")
    OUT_JSON.write_text(json.dumps({
        "script": "s18_ablation_paired_tests.py",
        "seed": SEED,
        "bootstrap_n": BOOTSTRAP_N,
        "n_paired": len(m),
        "mcnemar": [
            {"method": lab, "b": b, "c": c, "p": p, "p_holm": pa_}
            for (lab, b, c, p), pa_ in zip(raw, adj)
        ],
        "layer_b_identical": bool(lb[1] == 0 and lb[2] == 0),
        "methods": j_methods,
    }, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n요약 기록 : {OUT_TXT}")
    print(f"기계 판독 : {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
