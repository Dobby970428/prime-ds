"""
경보 선택 누출 정량화.

무엇을 묻는가
    2026-09-14 Layer B 경보 보강은 **구 학습셋(랜덤 분할, 5,796건)** 의 결손
    분석에 근거해 이루어졌다. 그 뒤 s1 을 스캐폴드 분할로 재실행하면서
    학습·검증이 다시 갈렸으므로, **새 검증셋(1,288건) 중 일부는 경보를 고를
    때 들여다본 화합물** 이다.

    그 부분집합에서 Layer B 가 유독 잘 맞는다면 경보 집합이 구 학습셋에
    맞춰 조립된 것이고, 두 부분집합의 성능이 비슷하면 선택이 특정 화합물이
    아니라 구조 계열 수준에서 이루어졌다는 근거가 된다.

    Layer B 민감도가 86.1%(구 검증셋) → 74.5%(새 검증셋)로 떨어졌는데,
    이것이 (a) 경보 선택 누출이 걷힌 탓인지 (b) 스캐폴드 분할로 검증셋이
    어려워진 탓인지 가르는 것이 목적이다.

어떻게 가르는가
    구 학습셋과 새 검증셋을 **InChIKey14(연결성)** 로 대조한다. 양쪽 모두
    p74 로 표준화한 뒤 키를 만든다. 구 학습셋은 표준화 전 원본이므로
    이 처리가 없으면 매칭이 샌다.

        노출군   : 새 검증셋 ∩ 구 학습셋   — 경보 선택 때 본 구조
        비노출군 : 새 검증셋 − 구 학습셋   — 처음 보는 구조

입력
    논문자료/02_검증결과/pipeline_validation_results.csv   (s5 산출, pred_layer_b)
    논문자료/_백업_랜덤분할_*/01_데이터셋/ames_train.csv    (구 학습셋 5,796건)

출력
    논문자료/02_검증결과/alert_selection_leakage_report.txt
    논문자료/02_검증결과/alert_selection_leakage.csv

사용법
    py s17_alert_selection_leakage.py
    py s17_alert_selection_leakage.py --old-train "경로\ames_train.csv"
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP / "논문자료"
PIPE = ROOT / "02_검증결과" / "pipeline_validation_results.csv"
OUT_TXT = ROOT / "02_검증결과" / "alert_selection_leakage_report.txt"
OUT_CSV = ROOT / "02_검증결과" / "alert_selection_leakage.csv"

_log: list[str] = []


def w(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


def find_old_train() -> Path | None:
    """가장 최근 백업 폴더의 구 학습셋을 찾는다."""
    cands = sorted(ROOT.glob("_백업_랜덤분할_*/01_데이터셋/ames_train.csv"))
    return cands[-1] if cands else None


def wilson(k: int, n: int) -> tuple[float, float]:
    """Wilson score 95% CI."""
    if n == 0:
        return (float("nan"), float("nan"))
    z = 1.959963984540054
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


def metrics(frame) -> dict:
    y = frame["actual_ames"].astype(int)
    p = frame["pred_layer_b"].astype(int)
    tp = int(((y == 1) & (p == 1)).sum())
    tn = int(((y == 0) & (p == 0)).sum())
    fp = int(((y == 0) & (p == 1)).sum())
    fn = int(((y == 1) & (p == 0)).sum())
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "n": len(frame), "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "sens": tp / (tp + fn) if tp + fn else float("nan"),
        "sens_ci": wilson(tp, tp + fn),
        "spec": tn / (tn + fp) if tn + fp else float("nan"),
        "spec_ci": wilson(tn, tn + fp),
        "prev": (tp + fn) / len(frame) if len(frame) else float("nan"),
        "mcc": ((tp * tn - fp * fn) / den) if den else float("nan"),
    }


def show(label: str, m: dict) -> None:
    w(f"  ── {label}")
    w(f"     n = {m['n']:,}   실제 양성 {m['tp'] + m['fn']:,} "
      f"(유병률 {m['prev'] * 100:.1f}%)")
    w(f"     TP {m['tp']:4,}  FP {m['fp']:4,}  FN {m['fn']:4,}  TN {m['tn']:4,}")
    w(f"     민감도 {m['sens'] * 100:5.1f}% "
      f"[{m['sens_ci'][0] * 100:.1f}–{m['sens_ci'][1] * 100:.1f}]"
      f"   특이도 {m['spec'] * 100:5.1f}% "
      f"[{m['spec_ci'][0] * 100:.1f}–{m['spec_ci'][1] * 100:.1f}]")
    w(f"     Matthews {m['mcc']:.3f}")
    w()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-train", default=None,
                    help="구 학습셋 CSV. 생략하면 최근 백업에서 찾는다.")
    args = ap.parse_args()

    import pandas as pd
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    from p74_structure_standardize import standardized_mol_from_smiles

    old_train = Path(args.old_train) if args.old_train else find_old_train()
    if old_train is None or not old_train.exists():
        w("[오류] 구 학습셋을 찾지 못했습니다.")
        w("       논문자료/_백업_랜덤분할_*/01_데이터셋/ames_train.csv 가 있어야 합니다.")
        w("       --old-train 으로 직접 지정할 수도 있습니다.")
        return 1
    if not PIPE.exists():
        w(f"[오류] s5 산출물이 없습니다: {PIPE}")
        return 1

    w("=" * 78)
    w("  경보 선택 누출 — 새 검증셋 중 구 학습셋 출신은 얼마나 되는가")
    w("=" * 78)
    w()
    w(f"  새 검증셋 : {PIPE.name}")
    w(f"  구 학습셋 : {old_train.parent.parent.name}/{old_train.name}")
    w()

    def key14(smi: str) -> str | None:
        mol = standardized_mol_from_smiles(smi)
        if mol is None:
            return None
        try:
            return Chem.MolToInchiKey(mol)[:14]
        except Exception:  # noqa: BLE001
            return None

    old = pd.read_csv(old_train)
    smi_col = next((c for c in old.columns if "smi" in c.lower()), None)
    if smi_col is None:
        w(f"[오류] 구 학습셋에 SMILES 열이 없습니다: {list(old.columns)}")
        return 1

    w(f"  구 학습셋 {len(old):,}건 표준화 중... (호변이성체 정규화라 느립니다)")
    old_keys = {k for k in (key14(s) for s in old[smi_col].astype(str)) if k}
    w(f"  고유 InChIKey14 {len(old_keys):,}개")
    w()

    val = pd.read_csv(PIPE)
    w(f"  새 검증셋 {len(val):,}건 표준화 중...")
    val["key14"] = [key14(s) for s in val["smiles"].astype(str)]
    n_bad = int(val["key14"].isna().sum())
    if n_bad:
        w(f"  ⚠ 키 생성 실패 {n_bad}건 — 비노출군으로 둡니다.")
    val["seen_in_old_train"] = val["key14"].isin(old_keys)
    w()

    exposed = val[val["seen_in_old_train"]]
    fresh = val[~val["seen_in_old_train"]]

    w("-" * 78)
    w("  노출 분포")
    w("-" * 78)
    w()
    w(f"  노출군   (구 학습셋에 있던 구조) : {len(exposed):,}건 "
      f"({len(exposed) / len(val) * 100:.1f}%)")
    w(f"  비노출군 (처음 보는 구조)        : {len(fresh):,}건 "
      f"({len(fresh) / len(val) * 100:.1f}%)")
    w()
    w("  ※ 스캐폴드 분할은 골격 단위로 가르므로, 구 학습셋 화합물이")
    w("    새 검증셋으로 넘어오는 것은 정상이다. 문제는 그 부분집합에서")
    w("    Layer B 가 유독 잘 맞느냐이다.")
    w()

    if len(exposed) == 0 or len(fresh) == 0:
        w("  한쪽이 비어 비교할 수 없습니다.")
        return 1

    w("-" * 78)
    w("  Layer B 성능 — 노출 여부로 가른 것")
    w("-" * 78)
    w()

    m_all = metrics(val)
    m_exp = metrics(exposed)
    m_fresh = metrics(fresh)

    show("전체", m_all)
    show("노출군 — 경보 선택 때 본 구조", m_exp)
    show("비노출군 — 처음 보는 구조", m_fresh)

    w("-" * 78)
    w("  해석")
    w("-" * 78)
    w()

    d_sens = (m_exp["sens"] - m_fresh["sens"]) * 100
    d_mcc = m_exp["mcc"] - m_fresh["mcc"]
    w(f"  민감도 차 (노출 − 비노출) = {d_sens:+.1f}%p")
    w(f"  Matthews 차               = {d_mcc:+.3f}")
    w()

    # 유병률이 다르면 MCC 비교가 오염된다. 먼저 경고한다.
    d_prev = abs(m_exp["prev"] - m_fresh["prev"]) * 100
    if d_prev > 5:
        w(f"  ⚠ 두 부분집합의 유병률이 {d_prev:.1f}%p 다르다 "
          f"(노출 {m_exp['prev'] * 100:.1f}% · 비노출 {m_fresh['prev'] * 100:.1f}%).")
        w("    Matthews 계수는 유병률에 민감하므로 **민감도·특이도로 판단할 것.**")
        w()

    lo_e, hi_e = m_exp["sens_ci"]
    lo_f, hi_f = m_fresh["sens_ci"]
    overlap = not (hi_e < lo_f or hi_f < lo_e)

    if overlap:
        w("  ▸ 두 부분집합의 민감도 95% CI 가 겹친다.")
        w("    경보 선택이 특정 화합물이 아니라 구조 계열 수준에서 이루어졌다는")
        w("    근거다. 논문에는 이 표를 싣고 «선택 풀에 포함된 N건을 제외해도")
        w("    성능이 유지된다»고 기술하면 된다.")
    else:
        w("  ▸ 두 부분집합의 민감도 95% CI 가 겹치지 않는다.")
        if d_sens > 0:
            w("    노출군에서 유의하게 잘 맞는다 — **경보 집합이 구 학습셋에**")
            w("    **맞춰 조립되었을 가능성이 크다.** 새 학습셋 5,309건으로")
            w("    결손 분석을 다시 수행하고 경보를 재선택할 것을 권한다.")
            w("    그때까지 논문의 Layer B 수치는 비노출군 값을 함께 보고해야 한다.")
        else:
            w("    비노출군에서 오히려 잘 맞는다. 누출 방향이 아니므로")
            w("    유병률·난이도 차이로 설명해야 한다.")
    w()

    w("  ※ 이 분석은 경보 «선택»의 누출만 본다. 지문 누출(ECFP 완전 일치)은")
    w("    s1 [3/6] 이 이미 0건으로 만들었고 s11 (B) 가 확인한다.")
    w()

    cols = ["smiles", "actual_ames", "pred_layer_b", "layer_b_alerts",
            "seen_in_old_train", "key14"]
    val[[c for c in cols if c in val.columns]].to_csv(
        OUT_CSV, index=False, encoding="utf-8-sig")
    OUT_TXT.write_text("\n".join(_log), encoding="utf-8")

    print(f"\n화합물 목록 : {OUT_CSV}")
    print(f"요약 기록   : {OUT_TXT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
