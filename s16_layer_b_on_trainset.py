"""
학습셋에 Layer B 구조경보를 적용해 화합물별 판정을 산출한다.

Layer B 보강을 위한 결손 분석에는 «학습셋 5,796건에 대한 Layer B 판정»이
필요하다. s12 는 검증셋만 다루므로 이 스크립트로 따로 만든다.

경보 정의는 실제 판정 경로와 동일하게 가져온다.
    p62_toxicity_alerts.scan_ashby_tennant_alerts   구조경보 18종
    p65_m7_ensemble_matrix._scan_metabolization_fragments  대사활성화 단편 4종

    alert_present = (구조경보 있음) OR (대사 단편 있음)

s5 · s12 의 pred_layer_b 와 같은 정의다. 경보를 추가·수정한 뒤 다시 실행하면
보강 효과를 학습셋에서 먼저 확인할 수 있다.

입력
    논문자료/02_검증결과/external_tools/ames_train_key.csv
        (s13 --stem ames_train 이 만든다. PRIME_ID · smiles · actual_ames)

출력
    논문자료/02_검증결과/external_tools/layer_b_trainset_calls.csv
        PRIME_ID, smiles, actual_ames, prime_rule, alerts, n_alerts

사용법
    py s16_layer_b_on_trainset.py
    py s16_layer_b_on_trainset.py --key 다른키.csv --out 다른출력.csv
"""

from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
EXT_DIR = APP_DIR / "논문자료" / "02_검증결과" / "external_tools"
DEFAULT_KEY = EXT_DIR / "ames_train_key.csv"
DEFAULT_OUT = EXT_DIR / "layer_b_trainset_calls.csv"


def metrics(y_true, y_pred) -> dict:
    tp = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    tn = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 0)
    fp = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)

    def r(a, b):
        return a / b if b else float("nan")

    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "sens": r(tp, tp + fn), "spec": r(tn, tn + fp),
        "conc": r(tp + tn, len(y_true)),
        "ppv": r(tp, tp + fp), "npv": r(tn, tn + fn),
        "mcc": ((tp * tn - fp * fn) / den) if den else float("nan"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--key", default=str(DEFAULT_KEY),
                    help="PRIME_ID · smiles · actual_ames 를 담은 CSV")
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="출력 CSV")
    args = ap.parse_args()

    import pandas as pd
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")

    from p62_toxicity_alerts import scan_ashby_tennant_alerts
    from p65_m7_ensemble_matrix import _scan_metabolization_fragments

    key_path = Path(args.key)
    if not key_path.is_absolute():
        key_path = (APP_DIR / key_path).resolve()
    if not key_path.exists():
        print(f"[오류] 키 파일 없음: {key_path}")
        print("       py s13_export_for_external_tools.py --source \"논문자료\\01_데이터셋\\ames_train.csv\""
              " --prefix TRAIN --stem ames_train  을 먼저 실행하세요.")
        return 1

    df = pd.read_csv(key_path)
    for col in ("PRIME_ID", "smiles", "actual_ames"):
        if col not in df.columns:
            print(f"[오류] '{col}' 열이 없습니다. 있는 열: {list(df.columns)}")
            return 1

    print("=" * 74)
    print("  학습셋에 Layer B 적용")
    print("=" * 74)
    print(f"입력 : {key_path.name}  {len(df):,}건")
    print()

    rows = []
    failed = 0
    counter: Counter[str] = Counter()

    for rec in df.itertuples(index=False):
        smi = str(getattr(rec, "smiles"))
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            failed += 1
            rows.append({
                "PRIME_ID": rec.PRIME_ID, "smiles": smi,
                "actual_ames": int(rec.actual_ames),
                "prime_rule": 0, "alerts": "", "n_alerts": 0,
            })
            continue

        names, _ = scan_ashby_tennant_alerts(mol)
        metab = _scan_metabolization_fragments(mol)
        combined = list(dict.fromkeys(names + metab))

        for name in combined:
            counter[name] += 1

        rows.append({
            "PRIME_ID": rec.PRIME_ID,
            "smiles": smi,
            "actual_ames": int(rec.actual_ames),
            "prime_rule": int(bool(combined)),
            "alerts": "; ".join(combined),
            "n_alerts": len(combined),
        })

    out = pd.DataFrame(rows)
    if failed:
        print(f"⚠ SMILES 파싱 실패 {failed:,}건 — 음성으로 처리했습니다.")
        print()

    y = out["actual_ames"].tolist()
    m = metrics(y, out["prime_rule"].tolist())

    print("-" * 74)
    print("  Layer B 성능 (학습셋)")
    print("-" * 74)
    print(f"  TP {m['TP']:5,}  FP {m['FP']:5,}  FN {m['FN']:5,}  TN {m['TN']:5,}")
    print(f"  민감도 {m['sens']*100:5.1f}%   특이도 {m['spec']*100:5.1f}%   "
          f"일치율 {m['conc']*100:5.1f}%")
    print(f"  PPV    {m['ppv']*100:5.1f}%   NPV    {m['npv']*100:5.1f}%   "
          f"Matthews {m['mcc']:.3f}")
    print()

    # ── 경보별 정밀도 — 어느 경보가 느슨한가 ────────────────────────
    print("-" * 74)
    print("  경보별 정밀도 (학습셋)")
    print("-" * 74)
    print(f"  {'경보':<44} {'건수':>6} {'양성':>6} {'PPV':>7}")
    print("  " + "-" * 66)

    prevalence = out["actual_ames"].mean()
    stats = []
    for name, cnt in counter.most_common():
        hit = out[out["alerts"].str.contains(name, regex=False, na=False)]
        pos = int((hit["actual_ames"] == 1).sum())
        ppv = pos / len(hit) if len(hit) else float("nan")
        stats.append((name, len(hit), pos, ppv))
        print(f"  {name:<44} {len(hit):6,} {pos:6,} {ppv*100:6.1f}%")

    print()
    print(f"  참고 — 학습셋 전체 양성 비율 {prevalence*100:.1f}%")
    print("  PPV 가 이 값에 가깝거나 낮은 경보는 «아무거나 잡는» 경보다.")
    print("  제거하거나 패턴을 좁히면 특이도가 오른다.")
    print()

    weak = [s for s in stats if s[3] == s[3] and s[3] <= prevalence]
    if weak:
        print("  ⚠ 전체 양성 비율 이하인 경보")
        for name, n, pos, ppv in weak:
            print(f"      {name:<44} {n:6,}건  PPV {ppv*100:5.1f}%")
        print()

    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = (APP_DIR / out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False, encoding="utf-8-sig")
    print(f"산출 : {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
