"""
PRIME-DS 전체 파이프라인(3층 앙상블) 외부검증.

data/ames_validation.csv 의 화합물을 실제 PRIME-DS 판정 경로에 통과시켜
ICH M7 등급을 산출하고, 실제 Ames 결과와 대조해 성능 지표를 계산합니다.

이진 판정 기준 (ICH M7 운용 논리)
    Class 1 · 2 · 3  →  양성 (TTC 관리 또는 확인시험 필요)
    Class 4 · 5      →  음성 (ICH Q3A/Q3B 일반 기준)

동시에 다음 세 가지를 나란히 비교합니다
    Layer B 단독   — 구조 경보만 (expert rule-based, Derek 계열과 대응)
    Layer C 단독   — 통계 모델만 (statistical QSAR, Sarah/CASE Ultra 계열과 대응)
    전체 앙상블    — A + B + C 통합 후 최종 ICH M7 등급

사용법
    python s5_run_pipeline_validation.py
    python s5_run_pipeline_validation.py data\\다른검증셋.csv
    python s5_run_pipeline_validation.py --tag ames     # ablation 조건별 산출물 분리
    python s5_run_pipeline_validation.py --tag tox21

결과물
    data/pipeline_validation_results.csv          화합물별 상세 결과 (태그 없을 때)
    data/pipeline_validation_report.txt           요약 지표      (태그 없을 때)
    data/pipeline_validation_results_{tag}.csv    --tag 지정 시
    data/pipeline_validation_report_{tag}.txt     --tag 지정 시

2026-09-03 수정
    C-3c  리포트 헤더에 사용한 검증셋 파일의 SHA-256·건수를 기록합니다.
    C-9   ROC-AUC / PR-AUC / 임계값 스윕 / Wilson 95% CI 를 추가하고,
          Layer C 예측 실패를 «음성»으로 처리하던 fillna(0) 을 제거했습니다.
    C-12  --tag 로 산출물 파일명을 분리하고, 헤더에 PRIME_LAYER_C_MODEL 값과
          실제 사용한 모델 파일의 경로·mtime·크기를 남깁니다.
          (report_ames.txt 와 pipeline_validation_report.txt 가 바이트 동일해
           어느 실행에서 나온 산출물인지 추적할 수 없던 문제)
"""

from __future__ import annotations

import hashlib
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
FROZEN_DIR = DATA_DIR / "frozen"
# 2026-09-03 수정(C-3): 작업본(excluded_reason 표기본)을 우선 쓰고,
#                       없으면 동결본을 그대로 읽는다.
WORKING_INPUT = DATA_DIR / "ames_validation.csv"
FROZEN_INPUT = FROZEN_DIR / "ames_validation_v1.csv"
DEFAULT_INPUT = WORKING_INPUT if WORKING_INPUT.exists() else FROZEN_INPUT
OUT_CSV = DATA_DIR / "pipeline_validation_results.csv"
OUT_TXT = DATA_DIR / "pipeline_validation_report.txt"

POSITIVE_CLASSES = {1, 2, 3}  # 양성으로 간주할 ICH M7 등급

BAR = "=" * 70
_log: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# 2026-09-03 추가(C-9): Wilson score 95% 신뢰구간
def wilson_ci(successes: int, total: int, z: float = 1.959963984540054):
    if total <= 0:
        return float("nan"), float("nan")
    p = successes / total
    denom = 1.0 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    half = (z * ((p * (1 - p) / total + z * z / (4 * total * total)) ** 0.5)) / denom
    return max(0.0, center - half), min(1.0, center + half)


def fmt_ci(successes: int, total: int) -> str:
    lo, hi = wilson_ci(successes, total)
    if lo != lo:
        return "(CI 불가)"
    return f"[{lo * 100:.1f}–{hi * 100:.1f}]"


def threshold_sweep(y_true, prob, step: float = 0.05):
    """임계값 스윕 — ICH M7 맥락에서는 NPV 최대화 지점 선택이 핵심."""
    rows = []
    t = 0.05
    while t <= 0.9501:
        pred = [1 if p >= t else 0 for p in prob]
        TP = sum(1 for a, b in zip(y_true, pred) if a == 1 and b == 1)
        TN = sum(1 for a, b in zip(y_true, pred) if a == 0 and b == 0)
        FP = sum(1 for a, b in zip(y_true, pred) if a == 0 and b == 1)
        FN = sum(1 for a, b in zip(y_true, pred) if a == 1 and b == 0)
        rows.append({
            "threshold": round(t, 2), "TP": TP, "TN": TN, "FP": FP, "FN": FN,
            "sensitivity": TP / (TP + FN) if (TP + FN) else float("nan"),
            "specificity": TN / (TN + FP) if (TN + FP) else float("nan"),
            "npv": TN / (TN + FN) if (TN + FN) else float("nan"),
        })
        t += step
    return rows


def metrics(y_true, y_pred) -> dict:
    TP = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    TN = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 0)
    FP = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    FN = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)

    def s(a, b):
        return a / b if b else float("nan")

    sens = s(TP, TP + FN)
    spec = s(TN, TN + FP)
    denom = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
    return {
        "TP": TP, "TN": TN, "FP": FP, "FN": FN,
        # 2026-09-03 수정(C-9): 분모를 실제 평가 건수(TP+TN+FP+FN)로 잡는다.
        #   이전에는 len(y_true) 를 썼는데, 판정 실패분을 걸러낸 뒤에도
        #   전체 길이를 쓰면 일치율 분모가 어긋난다.
        "n_evaluated": TP + TN + FP + FN,
        "sensitivity": sens,
        "specificity": spec,
        "concordance": s(TP + TN, TP + TN + FP + FN),
        "npv": s(TN, TN + FN),
        "ppv": s(TP, TP + FP),
        "balanced_accuracy": (sens + spec) / 2,
        "mcc": ((TP * TN - FP * FN) / denom) if denom else float("nan"),
        "fn_rate": 1 - sens if sens == sens else float("nan"),
    }


def show_metrics(title: str, m: dict, n: int, note: str = "") -> None:
    say(f"  {title}")
    say(f"      혼동행렬        실제양성  실제음성")
    say(f"        예측 양성      {m['TP']:7,}  {m['FP']:8,}")
    say(f"        예측 음성      {m['FN']:7,}  {m['TN']:8,}")
    # 2026-09-03 수정(C-9): Wilson 95% CI 병기
    say(f"      Sensitivity  {m['sensitivity'] * 100:6.1f} % {fmt_ci(m['TP'], m['TP'] + m['FN']):>14s}   "
        f"Specificity  {m['specificity'] * 100:6.1f} % {fmt_ci(m['TN'], m['TN'] + m['FP']):>14s}")
    say(f"      Concordance  {m['concordance'] * 100:6.1f} % {fmt_ci(m['TP'] + m['TN'], m['n_evaluated']):>14s}   "
        f"NPV          {m['npv'] * 100:6.1f} % {fmt_ci(m['TN'], m['TN'] + m['FN']):>14s}")
    say(f"      PPV          {m['ppv'] * 100:6.1f} % {fmt_ci(m['TP'], m['TP'] + m['FP']):>14s}   "
        f"Balanced acc {m['balanced_accuracy'] * 100:6.1f} %")
    say(f"      MCC          {m['mcc']:6.3f}       "
        f"False neg    {m['fn_rate'] * 100:6.1f} %")
    say(f"      (n = {m['n_evaluated']:,}{f' / 판정 {n:,}' if n != m['n_evaluated'] else ''})")
    if note:
        say(f"      {note}")
    say()


def main() -> int:
    # 2026-09-03 수정(C-12): --tag 로 ablation 조건별 산출물을 분리한다.
    argv = sys.argv[1:]
    tag = ""
    positional: list[str] = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--tag" and i + 1 < len(argv):
            tag = argv[i + 1].strip()
            i += 2
            continue
        if a.startswith("--tag="):
            tag = a.split("=", 1)[1].strip()
            i += 1
            continue
        if a.startswith("--"):
            say(f"[오류] 알 수 없는 옵션: {a}   (사용 가능: --tag <이름>)")
            return 1
        positional.append(a)
        i += 1

    safe_tag = "".join(c for c in tag if c.isalnum() or c in "-_")
    if tag and safe_tag != tag:
        say(f"[오류] --tag 에는 영문자·숫자·'-'·'_' 만 쓸 수 있습니다: {tag!r}")
        return 1
    out_csv = DATA_DIR / (f"pipeline_validation_results_{safe_tag}.csv" if safe_tag
                          else OUT_CSV.name)
    out_txt = DATA_DIR / (f"pipeline_validation_report_{safe_tag}.txt" if safe_tag
                          else OUT_TXT.name)

    say(BAR)
    say("  PRIME-DS 전체 파이프라인 외부검증 (Layer A + B + C)")
    say(BAR)
    say()

    src = Path(positional[0]) if positional else DEFAULT_INPUT
    if not src.is_absolute():
        src = (BASE_DIR / src).resolve()
    if not src.exists():
        say(f"[오류] 검증셋을 찾을 수 없습니다: {src}")
        # 2026-09-03 수정(F-2): 실제 파일명으로 정정
        say("       먼저 python s1_setup_layer_c_ames.py 를 실행하세요.")
        return 1

    import pandas as pd
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")

    df_all = pd.read_csv(src)
    if "smiles" not in df_all.columns or "label" not in df_all.columns:
        say(f"[오류] 'smiles', 'label' 열이 필요합니다. 현재 열: {list(df_all.columns)}")
        return 1

    # --- 2026-09-03 추가(C-3c): 실행 출처 기록 --------------------------------
    say("  [실행 출처 기록]")
    say(f"  실행 시각        : {datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds')}")
    say(f"  태그             : {safe_tag or '(없음 — 기존 파일명 유지)'}")
    say(f"  검증셋 파일      : {src}")
    say(f"  검증셋 SHA-256   : {sha256_of(src)}")
    say(f"  검증셋 행 수     : {len(df_all):,}")
    if FROZEN_INPUT.exists():
        frozen_digest = sha256_of(FROZEN_INPUT)
        say(f"  동결본 SHA-256   : {frozen_digest}   (data/frozen/{FROZEN_INPUT.name})")
        if FROZEN_INPUT.resolve() == src.resolve():
            say("                     ← 동결본을 직접 사용")
    else:
        say("  동결본           : 없음 (s1_setup_layer_c_ames.py 로 동결을 권장)")
    forced = os.environ.get("PRIME_LAYER_C_MODEL", "")
    say(f"  PRIME_LAYER_C_MODEL = {forced or '(미지정 — 자동 선택)'}")

    # 2026-09-03 수정(C-3b): 행 삭제 대신 excluded_reason 컬럼으로 필터링한다.
    n_excluded = 0
    if "excluded_reason" in df_all.columns:
        reasons = df_all["excluded_reason"].fillna("").astype(str)
        excluded = df_all[reasons != ""]
        n_excluded = len(excluded)
        df = df_all[reasons == ""].copy()
    else:
        df = df_all.copy()
        excluded = df_all.iloc[0:0]
    say(f"  제외 표기(excluded_reason) : {n_excluded:,}건 → 분석 대상 {len(df):,}건")
    if n_excluded:
        counts: dict[str, int] = {}
        for r in excluded["excluded_reason"].astype(str):
            head = r.split(":")[0]
            counts[head] = counts.get(head, 0) + 1
        for k, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            say(f"      · {k}: {n:,}건")
    say()

    if df.empty:
        say("[오류] 분석 대상이 0건입니다. excluded_reason 을 확인하세요.")
        return 1

    say(f"검증셋: {src.name}   {len(df):,}건")
    say(f"        양성 {int((df['label'] == 1).sum()):,} / "
        f"음성 {int((df['label'] == 0).sum()):,}")
    say()

    from p60_m7_empirical_db import get_default_db_connection
    from p65_m7_ensemble_matrix import run_m7_ensemble_matrix

    try:
        conn = get_default_db_connection()
    except Exception as exc:  # noqa: BLE001
        say(f"(실측 DB 연결 실패: {exc} — Layer A 없이 진행)")
        conn = None

    import p63_deepchem_predictor as dp

    status = dp.model_status()
    say(f"Layer C 모델: {status['model_label']}")
    say(f"              {'Ames 학습 모델' if status.get('trained_on_ames') else 'Tox21 대체 모델'}")
    # 2026-09-03 추가(C-12): 실제 사용한 모델 파일을 추적 가능하게 기록한다.
    model_path = status.get("model_path")
    if model_path:
        mp = Path(str(model_path))
        if mp.exists():
            st = mp.stat()
            say(f"  모델 파일      : {mp}")
            say(f"  모델 mtime     : "
                f"{datetime.fromtimestamp(st.st_mtime, timezone.utc).astimezone().isoformat(timespec='seconds')}")
            say(f"  모델 크기      : {st.st_size:,} bytes")
        else:
            say(f"  모델 파일      : {mp}  (파일 없음)")
    say(f"  학습 엔드포인트: {status.get('endpoint') or '(불명)'}")
    if not status.get("trained_on_ames"):
        say("  ⚠️ 이 모델의 확률값은 Ames 세균 변이원성이 아닙니다 "
            "(Tox21 SR-p53 대리 지표).")
    say()
    say("판정 진행 중... (수 분 걸릴 수 있습니다)")

    rows = []
    t0 = time.time()
    failed = 0
    first_error: str | None = None
    for i, rec in enumerate(df.itertuples(index=False), start=1):
        smi = getattr(rec, "smiles")
        truth = int(getattr(rec, "label"))
        try:
            res = run_m7_ensemble_matrix(smi, db_connection=conn)
        except Exception as exc:  # noqa: BLE001
            failed += 1
            if first_error is None:
                import traceback
                first_error = f"{type(exc).__name__}: {exc}"
                say(f"\n[첫 실패 원인] {first_error}")
                say(traceback.format_exc().strip().splitlines()[-3])
            continue

        p = res.layer_c.p_mutagenic
        rows.append({
            "smiles": smi,
            "actual_ames": truth,
            "layer_a_db_hit": bool(res.layer_a.is_in_db),
            "layer_b_alert": bool(res.layer_b.alert_present),
            "layer_b_alerts": "; ".join(res.layer_b.combined_alerts) if
            getattr(res.layer_b, "combined_alerts", None) else "",
            "layer_c_prob": p,
            "rule_class": res.rule_class,
            "final_class": res.final_class,
            "final_class_label": res.class_label,
            "ai_inferred": bool(res.ai_inferred_liability),
        })

        if i % 100 == 0:
            rate = i / max(time.time() - t0, 0.001)
            eta = (len(df) - i) / max(rate, 0.001)
            print(f"   {i:,}/{len(df):,}   ({rate:.0f}건/초, 남은 시간 약 {eta / 60:.1f}분)")

    elapsed = time.time() - t0
    # 2026-09-03 수정(C-9): 판정 실패 건수와 지표 분모를 일치시킨다.
    #   이전에는 `failed` 를 인쇄하면서 metrics() 분모로 len(y_true) 를 썼다.
    n_input = len(df)
    say(f"완료 — 입력 {n_input:,}건 중 {len(rows):,}건 판정, "
        f"{failed:,}건 판정 실패, {elapsed / 60:.1f}분 소요")
    if len(rows) + failed != n_input:
        say(f"  [경고] 판정 {len(rows):,} + 실패 {failed:,} ≠ 입력 {n_input:,} — 집계를 확인하세요.")
    say()

    if not rows:
        say("[오류] 판정에 성공한 화합물이 없습니다.")
        return 1

    out = pd.DataFrame(rows)
    out["pred_ensemble"] = out["final_class"].isin(POSITIVE_CLASSES).astype(int)
    out["pred_layer_b"] = out["layer_b_alert"].astype(int)

    # 2026-09-03 수정(C-9): Layer C 예측 실패를 «음성»으로 처리하던 fillna(0) 제거.
    #   실패를 0으로 채우면 특이도가 조용히 부풀려진다. 실패 건은 집계만 하고
    #   Layer C 지표 계산에서 제외한다 (Layer B·앙상블 지표에는 영향 없음).
    n_layer_c_failed = int(out["layer_c_prob"].isna().sum())
    c_mask = out["layer_c_prob"].notna()
    out["pred_layer_c"] = None
    out.loc[c_mask, "pred_layer_c"] = (out.loc[c_mask, "layer_c_prob"] >= 0.5).astype(int)

    y = out["actual_ames"].tolist()
    y_c = out.loc[c_mask, "actual_ames"].tolist()
    pred_c = [int(v) for v in out.loc[c_mask, "pred_layer_c"].tolist()]
    prob_c = [float(v) for v in out.loc[c_mask, "layer_c_prob"].tolist()]

    say(BAR)
    say("  결과")
    say(BAR)
    say()
    say(f"  n_layer_c_failed : {n_layer_c_failed:,}건 "
        "(Layer C 예측 실패 — 음성 처리하지 않고 Layer C 지표에서 제외)")
    say()

    m_b = metrics(y, out["pred_layer_b"].tolist())
    m_c = metrics(y_c, pred_c)
    m_e = metrics(y, out["pred_ensemble"].tolist())

    show_metrics("① Layer B 단독 — 구조 경보 (expert rule-based)", m_b, len(out))
    show_metrics(
        "② Layer C 단독 — 통계 모델 (statistical QSAR)", m_c, len(out),
        note=(f"Layer C 예측 실패 {n_layer_c_failed:,}건 제외 후 n = {len(y_c):,}"
              if n_layer_c_failed else ""),
    )
    show_metrics("③ 전체 앙상블 — ICH M7 최종 등급 (Class 1·2·3 = 양성)", m_e, len(out))

    # 2026-09-03 추가(C-9): Layer C 임계값 비의존 지표 + 임계값 스윕
    say("  ④ Layer C 임계값 분석 (ICH M7: NPV 최대화 지점 선택이 핵심)")
    if len(set(y_c)) > 1 and prob_c:
        try:
            from sklearn.metrics import average_precision_score, roc_auc_score

            say(f"      ROC-AUC {float(roc_auc_score(y_c, prob_c)):.3f}   "
                f"PR-AUC {float(average_precision_score(y_c, prob_c)):.3f}   "
                f"(양성 기준선 {sum(y_c) / len(y_c):.3f}, n = {len(y_c):,})")
        except Exception as exc:  # noqa: BLE001
            say(f"      (ROC/PR-AUC 계산 실패: {exc})")
        say(f"      {'임계값':>7} {'Sens':>8} {'Spec':>8} {'NPV':>8}   {'FN':>5} {'FP':>5}")
        sweep = threshold_sweep(y_c, prob_c)
        for r in sweep:
            say(f"      {r['threshold']:>7.2f} {r['sensitivity'] * 100:7.1f}% "
                f"{r['specificity'] * 100:7.1f}% {r['npv'] * 100:7.1f}%   "
                f"{r['FN']:>5,} {r['FP']:>5,}")
        ok_sweep = [r for r in sweep if r["npv"] == r["npv"]]
        if ok_sweep:
            best = max(ok_sweep, key=lambda r: r["npv"])
            say(f"      → NPV 최대: 임계값 {best['threshold']:.2f} 에서 "
                f"NPV {best['npv'] * 100:.1f}% {fmt_ci(best['TN'], best['TN'] + best['FN'])} · "
                f"Sens {best['sensitivity'] * 100:.1f}% · Spec {best['specificity'] * 100:.1f}%")
    else:
        say("      (한 클래스만 남았거나 확률이 없어 임계값 분석을 건너뜁니다)")
    say()

    say(BAR)
    say("  요약 비교표 (논문용)")
    say(BAR)
    say()
    say(f"  {'방법':32s} {'n':>7s} {'Sens':>7s} {'Spec':>7s} {'Conc':>7s} {'NPV':>7s} {'MCC':>7s}")
    say("  " + "-" * 76)
    for name, m in (("Layer B 단독 (expert rule)", m_b),
                    ("Layer C 단독 (statistical)", m_c),
                    ("전체 앙상블 (A+B+C)", m_e)):
        say(f"  {name:32s} {m['n_evaluated']:7,} {m['sensitivity'] * 100:6.1f}% "
            f"{m['specificity'] * 100:6.1f}% "
            f"{m['concordance'] * 100:6.1f}% {m['npv'] * 100:6.1f}% {m['mcc']:7.3f}")
    say("  (Sens/Spec/NPV 의 Wilson 95% CI 는 위 개별 표에 있습니다)")
    say()

    # 등급 분포
    say("  ICH M7 등급 분포")
    dist = out["final_class"].value_counts().sort_index()
    for cls, cnt in dist.items():
        sub = out[out["final_class"] == cls]
        pos = int((sub["actual_ames"] == 1).sum())
        say(f"      Class {cls}: {cnt:5,}건   (실제 양성 {pos:,} / 음성 {cnt - pos:,})")
    say()

    n_ai = int(out["ai_inferred"].sum())
    n_dbhit = int(out["layer_a_db_hit"].sum())
    say(f"  Layer A DB 적중: {n_dbhit:,}건 "
        f"({'검증셋에서 사전 제외되어 정상' if n_dbhit == 0 else '⚠️ 순환논증 확인 필요'})")
    say(f"  AI Inferred Liability 발동: {n_ai:,}건")
    if n_ai:
        sub = out[out["ai_inferred"]]
        say(f"      이 중 실제 Ames 양성: {int((sub['actual_ames'] == 1).sum()):,}건")
    say()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    # 2026-09-03 수정(C-12): --tag 가 있으면 조건별 파일명으로 저장한다.
    out.to_csv(out_csv, index=False, encoding="utf-8-sig")
    out_txt.write_text("\n".join(_log), encoding="utf-8")
    print(f"상세 결과: data/{out_csv.name}")
    print(f"요약 기록: data/{out_txt.name}")
    if not safe_tag:
        print("(팁: --tag ames / --tag tox21 로 실행하면 ablation 조건별로 파일이 분리됩니다)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
