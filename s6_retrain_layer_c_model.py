"""
Layer C **ablation 기준선(Tox21 SR-p53) 학습** 스크립트 (RDKit + scikit-learn).

⚠️ 이 스크립트는 «Ames 재학습»이 아닙니다 (2026-09-03 정정, C-4).
--------------------------------------------------------------------
파일 이름은 문서·배치가 참조하므로 그대로 두었지만, 실제로 하는 일은
**다른 endpoint 모델(Tox21 SR-p53)을 만드는 것**입니다. Ames 모델을 만드는
스크립트는 `s1_setup_layer_c_ames.py` 입니다.

Tox21 SR-p53 은 진핵세포 p53 스트레스 응답이며 세균 변이원성(Ames)이
아닙니다. 이 대체가 Ames 대리 지표로 성립하지 않는다는 증거는 이 프로젝트
자체 산출물에 있습니다 —
    논문자료/02_검증결과/report_tox21.txt
    Layer C 단독 Sensitivity 5.1% · MCC 0.078 (위음성률 94.9%)
따라서 이 모델은 «Ames 대체 예측기»가 아니라, 논문의 ablation 실험에서
«엔드포인트가 어긋난 통계층을 넣으면 어떻게 되는가»를 보여 주는
**기준선(baseline)** 으로만 쓰십시오.

사용법:
    python s6_retrain_layer_c_model.py

학습 조건 (원본 DeepChem 버전과 동일):
    분자 지문 : ECFP / Morgan fingerprint (radius=2, 2048 bits)
    분류기    : RandomForestClassifier(n_estimators=100, random_state=0)
    학습 대상 : Tox21 SR-p53 endpoint

홀드아웃 평가 (2026-09-03 추가, C-4d)
    p63 의 `train_and_save_model` 은 전체 데이터로 fit 한 뒤 저장할 뿐
    홀드아웃이 없어 성능을 알 방법이 없습니다. p63 은 편집 대상이 아니므로,
    이 스크립트가 **자체적으로** 층화 20% 홀드아웃을 떼어 같은 조건으로
    학습·평가하고 sens / spec / MCC 등을 출력·저장합니다.
    (배포 모델은 종전대로 전체 데이터로 학습된 것이며, 홀드아웃 모델은
     성능 보고 전용으로 따로 학습합니다 — 두 모델을 혼동하지 마십시오.)

결과물:
    models/tox21_ecfp_rf_sklearn.pkl        배포용 모델 (전체 데이터 학습)
    data/tox21_baseline_holdout_metrics.txt 홀드아웃 성능 기록
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

# 이 스크립트는 Tox21 기준선 전용이므로, 확인 단계에서도 Tox21 모델을 보도록 고정한다.
# (고정하지 않으면 Ames 모델이 있을 때 그쪽을 불러와 엉뚱한 결과를 보고한다)
# 2026-09-03(C-4c): 이 고정은 **현재 프로세스 전체**에 적용된다. 종료 시 경고로 알린다.
os.environ["PRIME_LAYER_C_MODEL"] = "tox21"

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
HOLDOUT_TXT = DATA_DIR / "tox21_baseline_holdout_metrics.txt"

HOLDOUT_FRACTION = 0.20
HOLDOUT_SEED = 42

BANNER = "=" * 62


def _holdout_evaluation(dp, lines: list[str]) -> None:
    """2026-09-03 추가(C-4d): 층화 20% 홀드아웃으로 기준선 성능을 실측한다.

    p63.train_and_save_model 은 홀드아웃 없이 전체 데이터로 fit 하므로,
    저장된 배포 모델로는 성능을 알 수 없다. p63 은 편집 대상이 아니라
    여기서 같은 조건(ECFP-2048 / RF n=100, seed=0)으로 따로 학습해 평가한다.
    """
    import gzip

    import numpy as np
    import pandas as pd
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split

    def out(msg: str = "") -> None:
        print(msg)
        lines.append(msg)

    out(BANNER)
    out("  ablation 기준선 홀드아웃 평가 (Tox21 SR-p53)")
    out(BANNER)
    out(f"  분할        : 층화 {int(HOLDOUT_FRACTION * 100)}% 홀드아웃 (seed={HOLDOUT_SEED})")
    out("  주의        : 여기서 학습한 모델은 «성능 보고 전용»이며 저장하지 않습니다.")
    out("               배포 모델은 전체 데이터로 학습된 별개의 파일입니다.")

    task = dp._MUTAGEN_TASK
    with gzip.open(dp._TOX21_CSV, "rt") as handle:
        df = pd.read_csv(handle)
    df = df.dropna(subset=["smiles", task])
    df[task] = pd.to_numeric(df[task], errors="coerce")
    df = df.dropna(subset=[task])
    df = df[df[task].isin([0.0, 1.0])]

    X_rows, y_rows = [], []
    skipped = 0
    for _, row in df.iterrows():
        feat = dp.featurize_smiles_one(str(row["smiles"]))
        if feat is None:
            skipped += 1
            continue
        X_rows.append(feat)
        y_rows.append(int(row[task]))

    if len(X_rows) < 100:
        out("  [건너뜀] 홀드아웃 평가에 필요한 데이터가 부족합니다.")
        return

    X = np.vstack(X_rows)
    y = np.asarray(y_rows, dtype=int)
    out(f"  사용 건수   : {len(y):,}건 (지문 실패 {skipped:,}건 제외) · "
        f"양성 {int((y == 1).sum()):,} / 음성 {int((y == 0).sum()):,}")

    Xtr, Xte, ytr, yte = train_test_split(
        X, y, test_size=HOLDOUT_FRACTION, random_state=HOLDOUT_SEED, stratify=y
    )
    model = RandomForestClassifier(n_estimators=100, random_state=0)
    model.fit(Xtr, ytr)

    pos_idx = list(model.classes_).index(1)
    prob = model.predict_proba(Xte)[:, pos_idx]
    pred = (prob >= 0.5).astype(int)

    TP = int(((pred == 1) & (yte == 1)).sum())
    TN = int(((pred == 0) & (yte == 0)).sum())
    FP = int(((pred == 1) & (yte == 0)).sum())
    FN = int(((pred == 0) & (yte == 1)).sum())

    def safe(a, b):
        return a / b if b else float("nan")

    sens = safe(TP, TP + FN)
    spec = safe(TN, TN + FP)
    npv = safe(TN, TN + FN)
    ppv = safe(TP, TP + FP)
    conc = safe(TP + TN, len(yte))
    denom = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
    mcc = ((TP * TN - FP * FN) / denom) if denom else float("nan")

    out(f"  학습 {len(ytr):,}건 / 홀드아웃 {len(yte):,}건")
    out("  혼동행렬 (기준 0.5)      실제 양성   실제 음성")
    out(f"      예측 양성            {TP:8,}   {FP:9,}")
    out(f"      예측 음성            {FN:8,}   {TN:9,}")
    out(f"      Sensitivity          {sens * 100:6.1f} %")
    out(f"      Specificity          {spec * 100:6.1f} %")
    out(f"      Concordance          {conc * 100:6.1f} %")
    out(f"      NPV                  {npv * 100:6.1f} %")
    out(f"      PPV                  {ppv * 100:6.1f} %")
    out(f"      Balanced accuracy    {(sens + spec) / 2 * 100:6.1f} %")
    out(f"      MCC                  {mcc:6.3f}")
    try:
        from sklearn.metrics import average_precision_score, roc_auc_score

        out(f"      ROC-AUC              {float(roc_auc_score(yte, prob)):6.3f}")
        out(f"      PR-AUC               {float(average_precision_score(yte, prob)):6.3f}")
    except Exception as exc:  # noqa: BLE001
        out(f"      (ROC/PR-AUC 계산 실패: {exc})")
    out("  ※ 이 수치는 Tox21 SR-p53 을 맞히는 성능입니다. Ames 변이원성 성능이 아닙니다.")
    out()


def main() -> int:
    log_lines: list[str] = []
    print(BANNER)
    # 2026-09-03 수정(C-4a): «재학습»이 아니라 «ablation 기준선 학습»으로 정정
    print("  PRIME-DS  Layer C — ablation 기준선(Tox21 SR-p53) 학습")
    print("  ECFP-2048 + RandomForest  /  Tox21 SR-p53")
    print("  ⚠️ Ames 재학습이 아닙니다. Ames 모델은 s1_setup_layer_c_ames.py 로 만듭니다.")
    print(BANNER)
    print()

    print(f"Python: {sys.version.split()[0]}  ({sys.executable})")

    try:
        import numpy  # noqa: F401
        import pandas  # noqa: F401
        import rdkit
        import sklearn
    except ImportError as exc:
        print(f"\n[오류] 필요한 라이브러리가 없습니다: {exc}")
        print("       pip install rdkit scikit-learn pandas numpy")
        return 1

    print(f"RDKit: {rdkit.__version__}   scikit-learn: {sklearn.__version__}")
    print()

    import p63_deepchem_predictor as dp

    if dp._TOX21_CSV.exists():
        print(f"학습 데이터: {dp._TOX21_CSV.name} (이미 존재)")
    else:
        print("학습 데이터가 없어 내려받습니다...")
    print()

    started = time.time()
    try:
        dp.train_and_save_model(verbose=True)
    except Exception as exc:  # noqa: BLE001
        print(f"\n[오류] 학습 실패: {type(exc).__name__}: {exc}")
        return 1

    elapsed = time.time() - started
    print(f"\n소요 시간: {elapsed:.1f}초")
    print()

    # 2026-09-03 추가(C-4d): 층화 홀드아웃으로 기준선 성능을 실측한다.
    try:
        _holdout_evaluation(dp, log_lines)
    except Exception as exc:  # noqa: BLE001
        print(f"[경고] 홀드아웃 평가를 건너뜁니다: {type(exc).__name__}: {exc}")
        print()

    # 검증 — 캐시를 비우고 저장된 모델로 실제 예측
    dp._model_ready = False
    dp._model = None
    dp._model_error = None

    print(BANNER)
    print("  동작 확인")
    print(BANNER)

    status = dp.model_status()
    if not status["ready"]:
        print(f"[오류] 모델을 불러오지 못했습니다: {status['error']}")
        return 1

    meta = status.get("meta") or {}
    print(f"사용 모델 : {status.get('model_label')}")
    print(f"모델 파일 : {status.get('model_path')}")

    def _num(*keys):
        """메타 키 이름이 모델마다 다를 수 있어 여러 후보를 시도한다."""
        for k in keys:
            v = meta.get(k)
            if isinstance(v, (int, float)):
                return f"{int(v):,}"
        return "?"

    if meta:
        print(f"학습 건수 : {_num('n_train')}  "
              f"(양성 {_num('n_positive', 'n_train_positive')} / "
              f"음성 {_num('n_negative', 'n_train_negative')})")
    print()

    samples = [
        ("CN(C)N=O", "N-니트로소디메틸아민 (NDMA) — 변이원성 양성 물질"),
        ("Nc1ccccc1", "아닐린 — 방향족 아민"),
        ("CC(=O)O", "아세트산 — 음성 대조"),
        ("CC(=O)Oc1ccccc1C(=O)O", "아스피린 — 음성 대조"),
    ]
    smiles_only = [s for s, _ in samples]
    results = dp.predict_advanced_properties(smiles_only)

    # 2026-09-03 수정(C-4b): «P(변이원성)» 표기는 사실과 다르다.
    #   이 확률은 Tox21 SR-p53 (진핵세포 p53 스트레스 응답) 값이다.
    print(f"{'SMILES':32s} {'P(SR-p53 stress response)':>26s}  판정")
    print("-" * 74)
    ok = 0
    for smi, label in samples:
        row = results.get(smi, {})
        prob = row.get("Ames_Mutagenicity_Prob")
        if prob is None:
            print(f"{smi:32s} {'실패':>26s}  {row.get('Prediction_Class', '')}")
            continue
        ok += 1
        print(f"{smi:32s} {prob:26.3f}  {row.get('Prediction_Class', '')}")
        print(f"{'':32s} {'':26s}  └ {label}")
    print()
    print("  ※ 위 확률은 Tox21 SR-p53 (진핵세포 p53 스트레스 응답)입니다.")
    print("     **Ames 세균 변이원성의 대리 지표가 아닙니다.**")
    print("     자체 검증(논문자료/02_검증결과/report_tox21.txt)에서 이 모델의")
    print("     Layer C 단독 Sensitivity 는 5.1%, MCC 는 0.078 이었습니다.")
    print("     (p63 의 결과 키 이름이 'Ames_Mutagenicity_Prob' 인 것은 하위 호환")
    print("      때문이며, 값의 의미가 Ames 라는 뜻이 아닙니다.)")

    # 홀드아웃 기록 저장
    if log_lines:
        try:
            DATA_DIR.mkdir(parents=True, exist_ok=True)
            HOLDOUT_TXT.write_text("\n".join(log_lines), encoding="utf-8")
            print()
            print(f"홀드아웃 성능 기록: data/{HOLDOUT_TXT.name}")
        except Exception as exc:  # noqa: BLE001
            print(f"(홀드아웃 기록 저장 실패: {exc})")

    # 2026-09-03 추가(C-4c): 환경변수 고정을 종료 시 경고로 알린다.
    print()
    print(BANNER)
    print("  [경고] 이 스크립트는 PRIME_LAYER_C_MODEL='tox21' 을 **현재 프로세스에**")
    print("         고정합니다. 같은 프로세스에서 이어서 실행되는 코드(배치 파일의")
    print("         후속 단계가 아니라, 이 파이썬 프로세스 안에서 import 되는 코드)는")
    print("         Ames 모델이 있어도 Tox21 모델을 보게 됩니다.")
    print("         Ames 조건으로 검증하려면 **새 프로세스**에서 실행하고,")
    print("         필요하면 명시적으로 지정하십시오:")
    print("             set PRIME_LAYER_C_MODEL=ames  &&  python s5_run_pipeline_validation.py --tag ames")
    print("             set PRIME_LAYER_C_MODEL=tox21 &&  python s5_run_pipeline_validation.py --tag tox21")
    print(BANNER)

    if ok == len(samples):
        print()
        print("✅ ablation 기준선 모델이 정상 작동합니다.")
        print("   이 모델은 비교 기준선이며, 실사용 Layer C 는 Ames 모델")
        print("   (models/ames_ecfp_rf.pkl) 입니다.")
        return 0

    print()
    print(f"⚠️  {len(samples) - ok}건 예측 실패. 위 메시지를 확인하세요.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
