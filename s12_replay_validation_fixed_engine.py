"""고정 검증셋을 현재 엔진으로 다시 태워 3개 층의 성능을 재산출합니다.

왜 필요한가
-----------
`논문자료/02_검증결과/pipeline_validation_results.csv` 는 2026-08-21 시점 엔진이
만든 예측입니다. 2026-09-03 에 구조 경보(p62)·합의 규칙(p65)·M7 라우팅(p66)이
수정되었으므로, 같은 검증셋이라도 예측이 달라집니다. 이 스크립트는 **검증셋과
Layer C 모델은 그대로 두고 엔진만 현재 코드로 바꿔** 다시 판정합니다.

    구 CSV 의 layer_c_prob 와 이 스크립트가 새로 계산한 확률이 1,440/1,440
    일치하는지 검사합니다. 일치하면 «모델은 같고 엔진만 달라졌다» 는 것이
    보장되므로, 성능 차이는 전부 엔진 수정에서 온 것입니다.

s5_run_pipeline_validation.py 와 다른 점
    s5  : 검증셋 → 엔진 → 리포트 (정식 검증. 검증셋 자체를 새로 만들 때 씀)
    s12 : 기존 검증 결과 CSV 를 기준선으로 두고 «수정 전 대비 무엇이 달라졌나»
          를 나란히 보여줍니다. 논문의 수정 전/후 비교표가 여기서 나옵니다.

실행
    set PRIME_LAYER_C_MODEL=ames
    python s12_replay_validation_fixed_engine.py

필요 파일
    논문자료/02_검증결과/pipeline_validation_results.csv   (기준선)
    논문자료/01_데이터셋/ames_train.csv                     (적용범위 계산용)
    models/ames_ecfp_rf.pkl                                (Layer C)
    data/m7_empirical_seed.csv                             (Layer A)

출력
    논문자료/02_검증결과/replay_fixed_engine_report.txt
    논문자료/02_검증결과/replay_fixed_engine_results.csv
    논문자료/02_검증결과/replay_fixed_engine.json

2026-09-04 작성.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

BASELINE = APP_DIR / "논문자료" / "02_검증결과" / "pipeline_validation_results.csv"
TRAIN_CSV = APP_DIR / "논문자료" / "01_데이터셋" / "ames_train.csv"
OUT_DIR = APP_DIR / "논문자료" / "02_검증결과"
OUT_TXT = OUT_DIR / "replay_fixed_engine_report.txt"
OUT_CSV = OUT_DIR / "replay_fixed_engine_results.csv"
OUT_JSON = OUT_DIR / "replay_fixed_engine.json"

AD_THRESHOLD = 0.30
ECFP_RADIUS = 2
ECFP_BITS = 2048
BOOTSTRAP_N = 2000
SEED = 20260904
POSITIVE_CLASSES = {1, 2, 3}

BAR = "=" * 78
_log: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - m) / den * 100, (c + m) / den * 100)


def confusion(y, p):
    TP = int(sum(1 for a, b in zip(y, p) if a == 1 and b == 1))
    FN = int(sum(1 for a, b in zip(y, p) if a == 1 and b == 0))
    FP = int(sum(1 for a, b in zip(y, p) if a == 0 and b == 1))
    TN = int(sum(1 for a, b in zip(y, p) if a == 0 and b == 0))
    return TP, FP, FN, TN


def metrics(y, p, label: str) -> dict:
    TP, FP, FN, TN = confusion(y, p)
    n = TP + FP + FN + TN
    sens = TP / (TP + FN) if TP + FN else float("nan")
    spec = TN / (TN + FP) if TN + FP else float("nan")
    den = math.sqrt(float((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
    return {
        "label": label, "n": n, "TP": TP, "FP": FP, "FN": FN, "TN": TN,
        "sensitivity": sens * 100, "sensitivity_ci": wilson(TP, TP + FN),
        "specificity": spec * 100, "specificity_ci": wilson(TN, TN + FP),
        "concordance": (TP + TN) / n * 100, "concordance_ci": wilson(TP + TN, n),
        "ppv": TP / (TP + FP) * 100 if TP + FP else float("nan"),
        "ppv_ci": wilson(TP, TP + FP),
        "npv": TN / (TN + FN) * 100 if TN + FN else float("nan"),
        "npv_ci": wilson(TN, TN + FN),
        "false_positive_rate": (1 - spec) * 100,
        "false_negative_rate": (1 - sens) * 100,
        "mcc": (TP * TN - FP * FN) / den if den else float("nan"),
        # Cooper 비 = 민감도 ÷ 위양성률. Contrera(2013) 의 "ROC" 와 같은 정의이며
        # ROC 곡선 아래 면적(AUC)이 아닙니다. 비교 가능하도록 같은 정의를 씁니다.
        "cooper_roc": sens / (1 - spec) if spec < 1 else float("inf"),
    }


def mcnemar(y, a, b, name_a: str, name_b: str) -> dict:
    from scipy.stats import binomtest
    ca = [x == t for x, t in zip(a, y)]
    cb = [x == t for x, t in zip(b, y)]
    B = int(sum(1 for i, j in zip(ca, cb) if i and not j))
    C = int(sum(1 for i, j in zip(ca, cb) if (not i) and j))
    p = binomtest(B, B + C, 0.5).pvalue if B + C else 1.0
    return {"pair": f"{name_a} vs {name_b}", "b": B, "c": C, "n_discordant": B + C, "p": float(p)}


def main() -> int:
    import numpy as np
    import pandas as pd
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    from sklearn.metrics import average_precision_score, roc_auc_score

    RDLogger.DisableLog("rdApp.*")

    if not BASELINE.is_file():
        say(f"[오류] 기준선 파일이 없습니다: {BASELINE}")
        say("       먼저 s5_run_pipeline_validation.py 를 실행하십시오.")
        return 1

    say(BAR)
    say("  PRIME-DS — 고정 검증셋 재판정 (수정 전 엔진 대비)")
    say(BAR)
    say()
    base = pd.read_csv(BASELINE)
    say(f"기준선 파일 : {BASELINE.name}")
    say(f"  SHA-256   : {sha256(BASELINE)}")
    say(f"  행 수      : {len(base):,}")

    import p60_m7_empirical_db as p60
    import p63_deepchem_predictor as p63
    import p65_m7_ensemble_matrix as p65

    say(f"Layer C 모델 : {p63.active_model_label()}")
    say()

    conn = p60.get_default_db_connection()
    rows, errors = [], 0
    for i, r in base.iterrows():
        try:
            e = p65.run_m7_ensemble_matrix(str(r["smiles"]), db_connection=conn)
            prob = e.layer_c.gnn.p_mutagenic
            rows.append({
                "smiles": r["smiles"],
                "actual_ames": int(r["actual_ames"]),
                "old_layer_b_alert": int(bool(r["layer_b_alert"])),
                "old_layer_c_prob": r["layer_c_prob"],
                "old_final_class": int(r["final_class"]),
                "old_pred_ensemble": int(r["pred_ensemble"]),
                "new_layer_b_alert": int(e.layer_b.alert_present),
                "new_layer_b_alerts": "; ".join(e.layer_b.combined_alerts),
                "new_layer_c_prob": prob,
                "new_rule_class": int(e.rule_class),
                "new_final_class": int(e.final_class),
                "new_pred_ensemble": int(e.final_class in POSITIVE_CLASSES),
                "new_pred_layer_c": (None if prob is None else int(prob >= 0.5)),
            })
        except Exception as exc:  # noqa: BLE001
            errors += 1
            if errors <= 5:
                say(f"  [경고] {r['smiles'][:48]} — {type(exc).__name__}: {exc}")
        if (i + 1) % 300 == 0:
            print(f"  ... {i + 1}/{len(base)}", file=sys.stderr)

    d = pd.DataFrame(rows)
    say(f"판정 완료 : {len(d):,}건 (실패 {errors}건)")

    # --- 모델 동일성 검사 -------------------------------------------------
    both = d.dropna(subset=["old_layer_c_prob", "new_layer_c_prob"])
    same = int((both["old_layer_c_prob"].round(4) == both["new_layer_c_prob"].round(4)).sum())
    say()
    say(f"Layer C 확률 일치 : {same:,}/{len(both):,}")
    if same == len(both):
        say("  → 모델과 특징 생성은 동일합니다. 아래 성능 차이는 전부 엔진 수정에서 옵니다.")
    else:
        say("  ⚠ 확률이 달라졌습니다. 모델 파일이 기준선과 다를 수 있습니다 —")
        say("    아래 비교를 «엔진 수정 효과»로 해석하지 마십시오.")

    y = d["actual_ames"].tolist()
    say()
    say(f"검증셋 구성 : 양성 {sum(y):,} / 음성 {len(y) - sum(y):,}")

    # --- 성능 ------------------------------------------------------------
    perf = {
        "layer_b_fixed": metrics(y, d["new_layer_b_alert"].tolist(), "Layer B — expert rule (fixed)"),
        "layer_c": metrics(y, d["new_pred_layer_c"].tolist(), "Layer C — statistical (Q)SAR"),
        "ensemble_fixed": metrics(y, d["new_pred_ensemble"].tolist(), "Four-layer ensemble (fixed)"),
        "layer_b_before": metrics(y, d["old_layer_b_alert"].tolist(), "Layer B (before fix)"),
        "ensemble_before": metrics(y, d["old_pred_ensemble"].tolist(), "Ensemble (before fix)"),
    }

    rng = np.random.default_rng(SEED)
    yy = np.array(y)
    for key, col in [("layer_b_fixed", "new_layer_b_alert"),
                     ("layer_c", "new_pred_layer_c"),
                     ("ensemble_fixed", "new_pred_ensemble"),
                     ("layer_b_before", "old_layer_b_alert"),
                     ("ensemble_before", "old_pred_ensemble")]:
        pp = d[col].to_numpy()
        boot = []
        for _ in range(BOOTSTRAP_N):
            idx = rng.integers(0, len(yy), len(yy))
            TP, FP, FN, TN = confusion(yy[idx], pp[idx])
            den = math.sqrt(float((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
            if den:
                boot.append((TP * TN - FP * FN) / den)
        perf[key]["mcc_ci"] = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))

    say()
    say(BAR)
    say("  표 1  외부 검증 성능 (95% CI 는 Wilson score, MCC 는 부트스트랩)")
    say(BAR)
    hdr = f"{'':<30}{'Spec':>7}{'Sens':>7}{'Conc':>7}{'FP':>7}{'FN':>7}{'PPV':>7}{'NPV':>7}{'MCC':>8}{'ROC*':>7}"
    say(hdr)
    say("-" * len(hdr))
    for k in ("layer_b_fixed", "layer_c", "ensemble_fixed", "layer_b_before", "ensemble_before"):
        m = perf[k]
        say(f"{m['label']:<30}{m['specificity']:>6.1f}%{m['sensitivity']:>6.1f}%{m['concordance']:>6.1f}%"
            f"{m['false_positive_rate']:>6.1f}%{m['false_negative_rate']:>6.1f}%{m['ppv']:>6.1f}%"
            f"{m['npv']:>6.1f}%{m['mcc']:>8.3f}{m['cooper_roc']:>7.2f}")
    say()
    say("  * ROC = Cooper 비 (민감도 ÷ 위양성률). Contrera(2013) 와 같은 정의이며 AUC 가 아닙니다.")
    say()
    say("  혼동행렬")
    for k in ("layer_b_fixed", "layer_c", "ensemble_fixed", "layer_b_before", "ensemble_before"):
        m = perf[k]
        say(f"    {m['label']:<30} TP={m['TP']:>4}  FP={m['FP']:>4}  FN={m['FN']:>4}  TN={m['TN']:>4}")

    # --- McNemar ---------------------------------------------------------
    mcs = [
        mcnemar(y, d["new_pred_ensemble"], d["new_layer_b_alert"], "Ensemble", "Layer B"),
        mcnemar(y, d["new_pred_ensemble"], d["new_pred_layer_c"], "Ensemble", "Layer C"),
        mcnemar(y, d["new_layer_b_alert"], d["new_pred_layer_c"], "Layer B", "Layer C"),
        mcnemar(y, d["new_pred_ensemble"], d["old_pred_ensemble"], "Ensemble (fixed)", "Ensemble (before)"),
    ]
    order = sorted(range(len(mcs)), key=lambda i: mcs[i]["p"])
    for rank, i in enumerate(order):
        mcs[i]["p_holm"] = min(1.0, mcs[i]["p"] * (len(mcs) - rank))
    say()
    say(BAR)
    say("  표 2  McNemar 쌍체 비교 (정확검정, Holm 보정)")
    say(BAR)
    say(f"{'비교 (A 대 B)':<40}{'b':>6}{'c':>6}{'p':>13}{'p(Holm)':>13}")
    say("-" * 78)
    for m in mcs:
        say(f"{m['pair']:<40}{m['b']:>6}{m['c']:>6}{m['p']:>13.3e}{m['p_holm']:>13.3e}")
    say("  b = A 만 맞힌 건수, c = B 만 맞힌 건수.")

    # --- ROC / PR --------------------------------------------------------
    ok = d["new_layer_c_prob"].notna()
    auc = float(roc_auc_score(d.loc[ok, "actual_ames"], d.loc[ok, "new_layer_c_prob"]))
    ap = float(average_precision_score(d.loc[ok, "actual_ames"], d.loc[ok, "new_layer_c_prob"]))
    say()
    say(f"Layer C 확률의 ROC-AUC = {auc:.4f}   PR-AUC = {ap:.4f}   (n={int(ok.sum()):,}, 결측 {int((~ok).sum())})")

    # --- 적용범위 ---------------------------------------------------------
    ad = None
    if TRAIN_CSV.is_file():
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=ECFP_RADIUS, fpSize=ECFP_BITS)
        tr = pd.read_csv(TRAIN_CSV)
        trfp = [gen.GetFingerprint(m) for m in (Chem.MolFromSmiles(s) for s in tr["smiles"]) if m]
        dm = []
        for s in d["smiles"]:
            m = Chem.MolFromSmiles(str(s))
            dm.append(float("nan") if m is None
                      else max(DataStructs.BulkTanimotoSimilarity(gen.GetFingerprint(m), trfp)))
        d["max_tanimoto_to_train"] = dm
        d["in_domain"] = d["max_tanimoto_to_train"] >= AD_THRESHOLD
        n_in = int(d["in_domain"].sum())
        n_id = int((d["max_tanimoto_to_train"] >= 0.9999).sum())
        say()
        say(BAR)
        say(f"  표 3  적용범위 (ECFP{ECFP_RADIUS * 2} {ECFP_BITS}bit, 최대 Tanimoto ≥ {AD_THRESHOLD})")
        say(BAR)
        say(f"  in-domain {n_in:,} ({n_in / len(d) * 100:.1f}%)  ·  out-of-domain {len(d) - n_in:,}")
        say(f"  d 분포: p05={d['max_tanimoto_to_train'].quantile(.05):.3f} "
            f"p25={d['max_tanimoto_to_train'].quantile(.25):.3f} "
            f"중앙={d['max_tanimoto_to_train'].median():.3f} "
            f"p75={d['max_tanimoto_to_train'].quantile(.75):.3f}")
        say()
        # 2026-09-15 수정: 이 자리의 문구가 「s1 을 --split scaffold 로 재실행하면
        #   해소된다」고 아직 안 한 것처럼 말하고 있었다. 이미 재실행했고 누출은
        #   0건이다. 실제 건수를 보고 문구를 고르도록 바꾼다.
        say(f"  학습셋과 지문이 완전히 일치(d = 1.000)하는 검증셋 화합물 : "
            f"{n_id:,}건 ({n_id / len(d) * 100:.2f}%)")
        if n_id:
            say("    ⚠ 스캐폴드 분할에서는 0건이어야 한다. s1 [3/6] 이 ECFP 완전 일치")
            say("      화합물을 검증셋에서 제외하므로, 0건이 아니면 s1 을 다시 돌리지")
            say("      않고 이 리포트를 쓴 것이다.")
            say("      ECFP 는 입체화학을 인코딩하지 않아 거울상이성질체가 같은 벡터가")
            say("      되고, full InChIKey 로 중복을 판정하면 학습·검증에 갈린다.")
        else:
            say("    ✅ 0건 — 골격 분할과 ECFP 완전 일치 제외가 정상 작동했다.")
            say("      참고: 이전 랜덤 분할에서는 187건(12.99%)이었다.")
        say()
        say(f"  {'':<22}{'n':>6}{'Sens':>8}{'Spec':>8}")
        for key, col in [("Layer B", "new_layer_b_alert"), ("Layer C", "new_pred_layer_c"),
                         ("Ensemble", "new_pred_ensemble")]:
            for tag, sub in [("in-domain", d[d["in_domain"]]), ("out-of-domain", d[~d["in_domain"]])]:
                TP, FP, FN, TN = confusion(sub["actual_ames"].tolist(), sub[col].tolist())
                se = TP / (TP + FN) * 100 if TP + FN else float("nan")
                sp = TN / (TN + FP) * 100 if TN + FP else float("nan")
                say(f"  {key + ' ' + tag:<22}{len(sub):>6}{se:>7.1f}%{sp:>7.1f}%")
        ad = {"threshold": AD_THRESHOLD, "in_domain": n_in,
              "out_of_domain": len(d) - n_in, "identical_fingerprint": n_id}
    else:
        say(f"\n[건너뜀] 적용범위 — 학습셋 파일이 없습니다: {TRAIN_CSV}")

    # --- 등급 분포 --------------------------------------------------------
    say()
    say("  ICH M7 등급 분포")
    say(f"    수정 전 : {d['old_final_class'].value_counts().sort_index().to_dict()}")
    say(f"    수정 후 : {d['new_final_class'].value_counts().sort_index().to_dict()}")
    changed = int((d["old_pred_ensemble"] != d["new_pred_ensemble"]).sum())
    say(f"  최종 예측이 바뀐 화합물 : {changed}건")

    # 2026-09-15 추가 — 「수정 전/후」 비교가 퇴화했는지 알린다.
    #
    #   이 스크립트는 기준선으로 pipeline_validation_results.csv 를 읽는다.
    #   그런데 s5 를 재실행하면 그 파일이 **이미 수정된 엔진의 산출물**이 되므로
    #   「수정 전」이 「수정 후」와 같아지고 두 줄이 동일하게 나온다.
    #   경고 없이 두면 «엔진 수정이 효과가 없었다»로 읽힌다.
    if changed == 0:
        say()
        say("  " + "!" * 70)
        say("  ⚠ 「수정 전」과 「수정 후」가 완전히 같다 (바뀐 화합물 0건).")
        say("    이 비교는 퇴화했다. 기준선 파일이 이미 수정된 엔진의 산출물이기")
        say("    때문이다 — s5 를 재실행하면 그렇게 된다.")
        say(f"    기준선 : {BASELINE.name}")
        say("    엔진 수정 전후를 비교하려면 수정 전 시점의 기준선 파일이 필요하다")
        say("    (논문자료/_백업_랜덤분할_*/ 의 리포트를 쓸 것).")
        say("    **위 「수정 전」 행을 논문에 인용하지 말 것.**")
        say("  " + "!" * 70)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    payload = {
        "baseline": {"file": str(BASELINE), "sha256": sha256(BASELINE), "rows": len(base)},
        "layer_c_model": p63.active_model_label(),
        "layer_c_prob_identical": {"same": same, "of": len(both)},
        "performance": perf, "mcnemar": mcs,
        "roc_auc": auc, "pr_auc": ap, "applicability_domain": ad,
        "n_prediction_changed": changed,
    }
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=1, default=float),
                        encoding="utf-8")
    say()
    say(BAR)
    say(f"  상세 결과 : {OUT_CSV.relative_to(APP_DIR)}")
    say(f"  기계 판독 : {OUT_JSON.relative_to(APP_DIR)}")
    say(f"  요약 기록 : {OUT_TXT.relative_to(APP_DIR)}")
    say(BAR)
    OUT_TXT.write_text("\n".join(_log) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
