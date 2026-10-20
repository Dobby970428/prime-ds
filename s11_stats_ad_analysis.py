# -*- coding: utf-8 -*-
"""
s11_stats_ad_analysis.py
========================
PRIME-DS — 검증셋 통계 분석 및 적용범위(applicability domain, AD) 산출

이 스크립트는 **읽기 전용**이다. 기존 산출물을 절대 덮어쓰지 않으며,
아래 세 개의 새 파일만 생성한다.

    논문자료/02_검증결과/statistics_and_ad_report.txt   (사람이 읽는 리포트)
    논문자료/02_검증결과/statistics_and_ad.json         (기계 판독용)
    논문자료/02_검증결과/ad_tanimoto_validation.csv     (AD 원자료)

수행하는 계산
-------------
(A) 통계
    1. 세 방식(4층 앙상블 · Layer B 단독 · Layer C 단독)의 혼동행렬과
       민감도·특이도·양성예측도·음성예측도·일치율·균형정확도·Matthews 계수
    2. 비율 지표에 Wilson score 95% CI, Matthews 계수에 부트스트랩 백분위 95% CI
       (화합물 단위 재표집, 2,000회, 고정 seed)
    3. McNemar 정확검정 3쌍 (scipy.stats.binomtest) + Holm 보정
    4. layer_c_prob 기반 ROC-AUC · PR-AUC(average precision)
    5. 임계값 스윕 0.05 간격 — 민감도·특이도·음성예측도·Matthews 계수

(B) 적용범위
    1. ames_train.csv 전 구조의 ECFP4(Morgan r=2, 2048비트) 지문
    2. ames_validation.csv 각 구조의 학습셋 대비 최대 Tanimoto 유사도 d(x)
    3. d < 0.30 → out-of-domain
    4. AD 구간별 성능(세 방식 전부)
    5. 최대 Tanimoto 분포 요약통계
    6. 사례연구 생성 구조(니트로사민 · 부형제 부가물 · 이원 배합 산물)의 AD

주의 — 분할 기준
----------------
2026-09-14 부터 **스캐폴드 분할** 기준이다. `s1_setup_layer_c_ames.py --split
scaffold` 가 구조 표준화(염 제거 · 전하 중화 · 호변이성체 정규화) 후
InChIKey14(연결성)로 중복을 판정하고 Murcko 골격으로 학습·검증을 가른다.
검증셋 n = 1,288 · 양성비 45.7% · ECFP 완전 일치 누출 0건.

이전 랜덤 분할(n = 1,440 · 양성비 54.5% · 누출 187건 12.99%) 기준 산출물은
`논문자료/_백업_랜덤분할_*/` 에 있다. **두 기준의 수치를 섞어 인용하지 말 것.**

지문을 만들기 전에 `p74_structure_standardize` 로 구조를 표준화한다. s1 이
학습셋을 표준화해 저장하므로 질의 구조도 같은 처리를 거쳐야 같은 척도가 된다.
표준화하지 않으면 사례연구 생성 구조의 d 가 낮게 나온다 (2026-09-14 수정).

실행
----
    python s11_stats_ad_analysis.py

재실행 순서상 위치: s1 ~ s5 를 실행한 뒤 s11 을 실행한다.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import BulkTanimotoSimilarity

RDLogger.DisableLog("rdApp.*")

# ─────────────────────────────────────────────────────────────────────────────
# 상수
# ─────────────────────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent

PATH_TRAIN = ROOT / "논문자료" / "01_데이터셋" / "ames_train.csv"
PATH_VALID = ROOT / "논문자료" / "01_데이터셋" / "ames_validation.csv"
PATH_PIPE = ROOT / "논문자료" / "02_검증결과" / "pipeline_validation_results.csv"
PATH_REPORT_AMES = ROOT / "논문자료" / "02_검증결과" / "report_ames.txt"
PATH_REPORT_TOX21 = ROOT / "논문자료" / "02_검증결과" / "report_tox21.txt"
PATH_EXCIP = ROOT / "논문자료" / "03_사례연구" / "excipient_case_study_results.csv"
PATH_BINARY = ROOT / "논문자료" / "03_사례연구" / "binary_compatibility_pairs.csv"
PATH_NDSRI_REPORT = ROOT / "논문자료" / "03_사례연구" / "case_study_report.txt"

OUT_DIR = ROOT / "논문자료" / "02_검증결과"
OUT_TXT = OUT_DIR / "statistics_and_ad_report.txt"
OUT_JSON = OUT_DIR / "statistics_and_ad.json"
OUT_CSV = OUT_DIR / "ad_tanimoto_validation.csv"

AD_THRESHOLD = 0.30            # Tanimoto 기준 in/out-of-domain 경계
ECFP_RADIUS = 2                # Morgan 반경 2 = ECFP4
ECFP_BITS = 2048
BOOTSTRAP_N = 2000
SEED = 20260903
SWEEP_STEP = 0.05
TARGET_SENSITIVITY = 0.95      # ICH M7 맥락 — 위음성 비용이 큼

METHODS = [
    ("pred_ensemble", "4층 통합 앙상블"),
    ("pred_layer_b", "Layer B 단독 (전문가 규칙)"),
    ("pred_layer_c", "Layer C 단독 (통계)"),
]

MCNEMAR_PAIRS = [
    ("pred_ensemble", "pred_layer_b", "앙상블 ↔ Layer B"),
    ("pred_ensemble", "pred_layer_c", "앙상블 ↔ Layer C"),
    ("pred_layer_b", "pred_layer_c", "Layer B ↔ Layer C"),
]


# ─────────────────────────────────────────────────────────────────────────────
# 유틸
# ─────────────────────────────────────────────────────────────────────────────

def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def wilson_ci(successes: int, total: int, z: float = 1.959963984540054):
    """Wilson score 95% 신뢰구간. total = 0 이면 (None, None)."""
    if total <= 0:
        return None, None
    p = successes / total
    denom = 1.0 + z * z / total
    centre = (p + z * z / (2 * total)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total))
    return max(0.0, centre - half), min(1.0, centre + half)


def mcc_from_counts(tp: int, fp: int, fn: int, tn: int):
    num = tp * tn - fp * fn
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if den == 0:
        return None
    return num / den


def confusion(y_true: np.ndarray, y_pred: np.ndarray):
    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    return tp, fp, fn, tn


def ratio_block(successes: int, total: int):
    """비율 + Wilson CI를 하나의 dict 로."""
    if total <= 0:
        return {"value": None, "n_num": successes, "n_den": total,
                "ci_low": None, "ci_high": None}
    lo, hi = wilson_ci(successes, total)
    return {"value": successes / total, "n_num": successes, "n_den": total,
            "ci_low": lo, "ci_high": hi}


def metric_set(y_true: np.ndarray, y_pred: np.ndarray, rng: np.random.Generator | None = None):
    """혼동행렬 + 전 지표 + Wilson CI + (rng 제공 시) Matthews 계수 부트스트랩 CI."""
    tp, fp, fn, tn = confusion(y_true, y_pred)
    n = tp + fp + fn + tn

    out = {
        "n": n,
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "n_actual_pos": tp + fn,
        "n_actual_neg": fp + tn,
        "sensitivity": ratio_block(tp, tp + fn),
        "specificity": ratio_block(tn, tn + fp),
        "ppv": ratio_block(tp, tp + fp),
        "npv": ratio_block(tn, tn + fn),
        "concordance": ratio_block(tp + tn, n),
    }

    sens = out["sensitivity"]["value"]
    spec = out["specificity"]["value"]
    out["balanced_accuracy"] = (
        None if (sens is None or spec is None) else (sens + spec) / 2.0
    )
    out["mcc"] = mcc_from_counts(tp, fp, fn, tn)
    out["mcc_ci_low"] = None
    out["mcc_ci_high"] = None
    out["mcc_bootstrap_n_valid"] = None

    if rng is not None and n > 0:
        idx_all = np.arange(n)
        vals = []
        for _ in range(BOOTSTRAP_N):
            idx = rng.choice(idx_all, size=n, replace=True)
            m = mcc_from_counts(*confusion(y_true[idx], y_pred[idx]))
            if m is not None:
                vals.append(m)
        if vals:
            arr = np.asarray(vals, dtype=float)
            out["mcc_ci_low"] = float(np.percentile(arr, 2.5))
            out["mcc_ci_high"] = float(np.percentile(arr, 97.5))
            out["mcc_bootstrap_n_valid"] = int(arr.size)

    return out


def mcnemar_exact(y_true: np.ndarray, pred_a: np.ndarray, pred_b: np.ndarray):
    """
    쌍체 이분형 비교. 정확검정(이항검정)을 사용한다.
      b = A 만 맞힌 건수, c = B 만 맞힌 건수
    불일치쌍이 적을 수 있으므로 χ² 근사 대신 정확검정을 쓴다.
    """
    correct_a = (pred_a == y_true)
    correct_b = (pred_b == y_true)
    b = int(np.sum(correct_a & ~correct_b))
    c = int(np.sum(~correct_a & correct_b))
    n_disc = b + c
    both_correct = int(np.sum(correct_a & correct_b))
    both_wrong = int(np.sum(~correct_a & ~correct_b))

    if n_disc == 0:
        p = 1.0
    else:
        p = float(binomtest(b, n_disc, 0.5, alternative="two-sided").pvalue)

    return {
        "b_only_a_correct": b,
        "c_only_b_correct": c,
        "n_discordant": n_disc,
        "both_correct": both_correct,
        "both_wrong": both_wrong,
        "p_raw": p,
        "test": "exact binomial (McNemar exact)",
    }


def holm(pvals):
    """Holm-Bonferroni 보정. 입력 순서를 유지한 보정 p 리스트를 돌려준다."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        val = (m - rank) * pvals[i]
        running = max(running, val)
        adj[i] = min(1.0, running)
    return adj


def fmt_pct(block, digits: int = 1) -> str:
    if block is None or block.get("value") is None:
        return "n/a"
    v = block["value"] * 100
    lo, hi = block["ci_low"], block["ci_high"]
    if lo is None:
        return f"{v:.{digits}f}%"
    return f"{v:.{digits}f}% [{lo * 100:.{digits}f}–{hi * 100:.{digits}f}]"


def fmt_num(x, digits: int = 3) -> str:
    return "n/a" if x is None else f"{x:.{digits}f}"


def fmt_p(p: float) -> str:
    if p is None:
        return "n/a"
    if p < 1e-300:
        return "< 1e-300"
    if p < 0.001:
        return f"{p:.3e}"
    return f"{p:.4f}"


# ─────────────────────────────────────────────────────────────────────────────
# 지문 / AD
# ─────────────────────────────────────────────────────────────────────────────

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=ECFP_RADIUS, fpSize=ECFP_BITS)


def fingerprint(smiles: str):
    """
    SMILES → ECFP4 비트벡터. 파싱 실패면 None.

    2026-09-14 수정 — 지문을 만들기 전에 p74 로 구조를 표준화한다.
    s1 이 학습셋·검증셋을 표준화해 저장하는데(s1 395-402행) 여기서는
    하지 않아, 표준화되지 않은 사례연구 생성 구조와 Tanimoto 를 재는
    불일치가 있었다. 표 Y 의 d 값이 이 때문에 낮게 나왔다.
    p73_applicability_domain.fingerprint 도 같은 정의를 쓴다.
    """
    from p74_structure_standardize import standardized_mol_from_smiles

    mol = standardized_mol_from_smiles(smiles)
    if mol is None:
        return None
    return _GEN.GetFingerprint(mol)


def build_train_fps(smiles_list):
    fps, failed = [], []
    for smi in smiles_list:
        fp = fingerprint(smi)
        if fp is None:
            failed.append(smi)
        else:
            fps.append(fp)
    return fps, failed


def max_tanimoto(smiles_list, train_fps):
    """각 구조의 학습셋 대비 최대 Tanimoto. 파싱 실패는 NaN."""
    out = np.full(len(smiles_list), np.nan, dtype=float)
    for i, smi in enumerate(smiles_list):
        fp = fingerprint(smi)
        if fp is None:
            continue
        out[i] = float(max(BulkTanimotoSimilarity(fp, train_fps)))
    return out


def distribution_summary(arr: np.ndarray):
    a = arr[~np.isnan(arr)]
    if a.size == 0:
        return None
    return {
        "n": int(a.size),
        "min": float(np.min(a)),
        "p05": float(np.percentile(a, 5)),
        "p25": float(np.percentile(a, 25)),
        "median": float(np.percentile(a, 50)),
        "p75": float(np.percentile(a, 75)),
        "p95": float(np.percentile(a, 95)),
        "max": float(np.max(a)),
        "mean": float(np.mean(a)),
    }


def parse_pipeline_report(path: Path):
    """
    report_ames.txt / report_tox21.txt 의 세 혼동행렬 블록을 읽는다.
    형식:
        ① Layer B 단독 ...
              혼동행렬        실제양성  실제음성
                예측 양성          556       310
                예측 음성          229       345
    반환: {"layer_b": (tp, fp, fn, tn), "layer_c": ..., "ensemble": ...}
    (읽기 전용. 화합물별 예측이 아니라 집계 혼동행렬이므로 쌍체 검정에는
     사용할 수 없고, Wilson CI 산출에만 사용한다.)
    """
    import re
    if not path.exists():
        return None
    keys = {"①": "layer_b", "②": "layer_c", "③": "ensemble"}
    lines = path.read_text(encoding="utf-8").splitlines()
    out, current = {}, None
    pos_pat = re.compile(r"예측 양성\s+(\d+)\s+(\d+)")
    neg_pat = re.compile(r"예측 음성\s+(\d+)\s+(\d+)")
    pending = {}
    for line in lines:
        for mark, name in keys.items():
            if mark in line:
                current = name
                pending = {}
        if current:
            mp = pos_pat.search(line)
            if mp:
                pending["tp"] = int(mp.group(1))
                pending["fp"] = int(mp.group(2))
            mn = neg_pat.search(line)
            if mn:
                pending["fn"] = int(mn.group(1))
                pending["tn"] = int(mn.group(2))
            if len(pending) == 4:
                out[current] = (pending["tp"], pending["fp"],
                                pending["fn"], pending["tn"])
                pending = {}
                current = None
    return out or None


def metrics_from_counts(tp: int, fp: int, fn: int, tn: int):
    """집계 혼동행렬만으로 산출 가능한 지표 + Wilson CI (부트스트랩 없음)."""
    n = tp + fp + fn + tn
    sens = ratio_block(tp, tp + fn)
    spec = ratio_block(tn, tn + fp)
    return {
        "n": n, "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "sensitivity": sens,
        "specificity": spec,
        "ppv": ratio_block(tp, tp + fp),
        "npv": ratio_block(tn, tn + fn),
        "concordance": ratio_block(tp + tn, n),
        "balanced_accuracy": (None if (sens["value"] is None or spec["value"] is None)
                              else (sens["value"] + spec["value"]) / 2.0),
        "mcc": mcc_from_counts(tp, fp, fn, tn),
    }


def parse_ndsri_smiles(path: Path):
    """
    case_study_report.txt 의 검증 니트로사민 표에서 구조 SMILES 를 추출한다.
    표 행 형식:  "    1   Class 3          0.650    329.38  <SMILES>"
    (사례연구 니트로사민은 CSV 산출물이 없어 리포트에서 읽는다. 읽기 전용.)
    """
    rows = []
    if not path.exists():
        return rows
    import re
    pat = re.compile(
        r"^\s*(\d+)\s+Class\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+(\S+)\s*$"
    )
    current_api = None
    api_pat = re.compile(r"^\s{2}([A-Z][A-Za-z]+)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        m_api = api_pat.match(line)
        if m_api:
            current_api = m_api.group(1)
            continue
        m = pat.match(line)
        if m:
            rows.append({
                "api": current_api,
                "ich_m7": f"Class {m.group(2)}",
                "prob": float(m.group(3)),
                "mw": float(m.group(4)),
                "smiles": m.group(5),
            })
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# 본체
# ─────────────────────────────────────────────────────────────────────────────

def main() -> int:
    for p in (PATH_TRAIN, PATH_VALID, PATH_PIPE):
        if not p.exists():
            print(f"[오류] 필수 입력 파일이 없습니다: {p}", file=sys.stderr)
            return 1

    L: list[str] = []              # 리포트 줄
    J: dict = {}                   # JSON 결과

    def w(line: str = ""):
        L.append(line)

    # ── 헤더 · 입력 파일 지문 ────────────────────────────────────────────
    w("=" * 78)
    w("  PRIME-DS  통계 분석 및 적용범위(applicability domain) 산출")
    w("  s11_stats_ad_analysis.py")
    w("=" * 78)
    w()
    w("[분할 기준]  2026-09-14 이후 — 스캐폴드 분할")
    w("  s1_setup_layer_c_ames.py --split scaffold")
    w("    · 구조 표준화(염 제거 · 전하 중화 · 호변이성체 정규화) 후")
    w("      InChIKey14(연결성) 로 중복 판정")
    w("    · Murcko 골격 분할 — 같은 골격은 학습·검증 한쪽에만")
    w("    · 검증셋 실제 건수·양성비는 아래 (A-1) 에 산출값으로 표시된다")
    w("  지문은 p74_structure_standardize 로 표준화한 뒤 만든다 —")
    w("  학습셋과 질의 구조가 같은 척도여야 하기 때문이다.")
    w()
    w("  ⚠ 이전 랜덤 분할(n = 1,440 · 양성비 54.5% · ECFP 누출 187건)")
    w("    기준 산출물은 논문자료/_백업_랜덤분할_*/ 에 있다.")
    w("    두 기준의 수치를 섞어 인용하지 말 것.")
    w()

    w("-" * 78)
    w("[입력 파일 — SHA-256 및 행 수]")
    w("-" * 78)

    inputs_meta = {}
    input_paths = [
        ("ames_train.csv", PATH_TRAIN),
        ("ames_validation.csv", PATH_VALID),
        ("pipeline_validation_results.csv", PATH_PIPE),
        ("excipient_case_study_results.csv", PATH_EXCIP),
        ("binary_compatibility_pairs.csv", PATH_BINARY),
        ("case_study_report.txt", PATH_NDSRI_REPORT),
    ]
    for name, path in input_paths:
        if not path.exists():
            w(f"  {name:<40s}  (없음)")
            inputs_meta[name] = {"present": False}
            continue
        digest = sha256_of(path)
        if path.suffix == ".csv":
            n_rows = int(len(pd.read_csv(path)))
            rows_desc = f"{n_rows:,} 행"
        else:
            n_rows = len(path.read_text(encoding="utf-8").splitlines())
            rows_desc = f"{n_rows:,} 줄"
        w(f"  {name}")
        w(f"      sha256 = {digest}")
        w(f"      {rows_desc}")
        inputs_meta[name] = {
            "present": True,
            "path": str(path.relative_to(ROOT)),
            "sha256": digest,
            "rows": n_rows,
        }
    w()
    w(f"  난수 seed = {SEED} · 부트스트랩 {BOOTSTRAP_N:,}회 · AD 임계값 = {AD_THRESHOLD}")
    w(f"  지문 = Morgan radius {ECFP_RADIUS} ({ECFP_BITS} bits, ECFP4)")
    w()

    J["meta"] = {
        "script": "s11_stats_ad_analysis.py",
        "split_basis": "scaffold (rdMolStandardize + InChIKey14 dedup + Murcko scaffold split)",
        "split_warning": (
            "2026-09-14 이후 스캐폴드 분할 기준. 이전 랜덤 분할"
            "(n=1,440, 양성비 54.5%, ECFP 완전일치 누출 187건) 산출물과 "
            "수치를 섞어 인용하지 말 것. 지문은 p74_structure_standardize 로 "
            "표준화한 뒤 산출한다."
        ),
        "standardizer_version": __import__(
            "p74_structure_standardize").STANDARDIZER_VERSION,
        "seed": SEED,
        "bootstrap_n": BOOTSTRAP_N,
        "ad_threshold": AD_THRESHOLD,
        "ecfp_radius": ECFP_RADIUS,
        "ecfp_bits": ECFP_BITS,
        "inputs": inputs_meta,
    }

    # ── 데이터 적재 ──────────────────────────────────────────────────────
    pipe = pd.read_csv(PATH_PIPE)
    train = pd.read_csv(PATH_TRAIN)
    valid = pd.read_csv(PATH_VALID)

    y = pipe["actual_ames"].to_numpy(dtype=int)
    n_total = len(pipe)

    # ─────────────────────────────────────────────────────────────────────
    # (A-1,2) 방식별 성능 + CI
    # ─────────────────────────────────────────────────────────────────────
    w("=" * 78)
    w("  (A-1) 방식별 예측 성능 — 검증셋 전체")
    w("=" * 78)
    w(f"  n = {n_total:,}  (실제 양성 {int(y.sum()):,} · 실제 음성 {int((y == 0).sum()):,},"
      f" 양성비 {y.mean() * 100:.1f}%)")
    w("  괄호 안은 Wilson score 95% CI (%). Matthews 계수 CI 는 부트스트랩 백분위.")
    w()

    overall = {}
    for col, label in METHODS:
        rng = np.random.default_rng(SEED)
        m = metric_set(y, pipe[col].to_numpy(dtype=int), rng=rng)
        overall[col] = m
        w(f"  ── {label}  [{col}]")
        w(f"     혼동행렬  TP {m['tp']} / FP {m['fp']} / FN {m['fn']} / TN {m['tn']}"
          f"   (합 {m['tp'] + m['fp'] + m['fn'] + m['tn']})")
        w(f"     민감도       {fmt_pct(m['sensitivity'])}")
        w(f"     특이도       {fmt_pct(m['specificity'])}")
        w(f"     양성예측도   {fmt_pct(m['ppv'])}")
        w(f"     음성예측도   {fmt_pct(m['npv'])}")
        w(f"     일치율       {fmt_pct(m['concordance'])}")
        w(f"     균형정확도   {fmt_num(m['balanced_accuracy'])}")
        w(f"     Matthews 계수 {fmt_num(m['mcc'])} "
          f"[{fmt_num(m['mcc_ci_low'])}–{fmt_num(m['mcc_ci_high'])}]")
        w()

    J["overall_performance"] = overall

    # ─────────────────────────────────────────────────────────────────────
    # (A-3) McNemar
    # ─────────────────────────────────────────────────────────────────────
    w("=" * 78)
    w("  (A-3) McNemar 정확검정 — 쌍체 비교 3쌍")
    w("=" * 78)
    w("  동일 검증셋 · 동일 화합물의 쌍체 관측이므로 McNemar 가 적합하다.")
    w("  불일치쌍이 적을 수 있으므로 χ² 근사가 아닌 정확검정(이항검정)을 사용한다.")
    w("  다중비교는 Holm 보정 (α = 0.05, 3 비교).")
    w()

    mc_results = []
    for a, b, label in MCNEMAR_PAIRS:
        r = mcnemar_exact(y, pipe[a].to_numpy(dtype=int), pipe[b].to_numpy(dtype=int))
        r["label"] = label
        r["method_a"] = a
        r["method_b"] = b
        mc_results.append(r)

    adj = holm([r["p_raw"] for r in mc_results])
    for r, pa in zip(mc_results, adj):
        r["p_holm"] = pa
        r["significant_holm_005"] = bool(pa < 0.05)

    for r in mc_results:
        w(f"  ── {r['label']}   ({r['method_a']} 대 {r['method_b']})")
        w(f"     b (A만 정답) = {r['b_only_a_correct']}    "
          f"c (B만 정답) = {r['c_only_b_correct']}    "
          f"불일치쌍 b+c = {r['n_discordant']}")
        w(f"     둘 다 정답 = {r['both_correct']}    둘 다 오답 = {r['both_wrong']}")
        w(f"     p (정확검정, 보정 전) = {fmt_p(r['p_raw'])}")
        w(f"     p (Holm 보정)         = {fmt_p(r['p_holm'])}"
          f"   → {'유의 (α=0.05)' if r['significant_holm_005'] else '유의하지 않음 (α=0.05)'}")
        w()

    J["mcnemar"] = mc_results

    # ─────────────────────────────────────────────────────────────────────
    # (A-4) ROC-AUC / PR-AUC
    # ─────────────────────────────────────────────────────────────────────
    w("=" * 78)
    w("  (A-4) Layer C 확률 기반 임계값 비의존 성능")
    w("=" * 78)

    prob_raw = pipe["layer_c_prob"]
    n_missing = int(prob_raw.isna().sum())
    mask = prob_raw.notna().to_numpy()
    prob = prob_raw.to_numpy(dtype=float)[mask]
    y_prob = y[mask]

    roc_auc = float(roc_auc_score(y_prob, prob))
    pr_auc = float(average_precision_score(y_prob, prob))
    prevalence = float(y_prob.mean())

    w(f"  layer_c_prob 결측 = {n_missing} 건 (결측은 0 으로 채우지 않고 제외)")
    w(f"  평가에 사용된 n = {int(mask.sum()):,}  "
      f"(양성 {int(y_prob.sum()):,} · 음성 {int((y_prob == 0).sum()):,})")
    w()
    w(f"  ROC-AUC                       = {roc_auc:.4f}")
    w(f"  PR-AUC (average precision)    = {pr_auc:.4f}")
    w(f"  무작위 분류기의 PR-AUC 기준선 = {prevalence:.4f}  (= 양성 유병률)")
    w()
    w(f"  고유 확률값 수 = {int(np.unique(prob).size)}  "
      f"(범위 {prob.min():.3f} ~ {prob.max():.3f})")
    w()

    J["threshold_free"] = {
        "layer_c_prob_missing": n_missing,
        "n_evaluated": int(mask.sum()),
        "roc_auc": roc_auc,
        "pr_auc_average_precision": pr_auc,
        "pr_auc_baseline_prevalence": prevalence,
        "prob_min": float(prob.min()),
        "prob_max": float(prob.max()),
        "prob_unique_values": int(np.unique(prob).size),
    }

    # ROC 곡선 좌표 (Figure 1 작도용)
    fpr, tpr, thr = roc_curve(y_prob, prob)
    J["roc_curve"] = {
        "fpr": [float(x) for x in fpr],
        "tpr": [float(x) for x in tpr],
        "thresholds": [float(x) if np.isfinite(x) else None for x in thr],
    }

    # ─────────────────────────────────────────────────────────────────────
    # (A-5) 임계값 스윕
    # ─────────────────────────────────────────────────────────────────────
    w("=" * 78)
    w("  (A-5) Layer C 임계값 스윕 (0.05 간격)  —  P >= t 이면 양성")
    w("=" * 78)
    w("  ICH M7 맥락에서는 위음성 비용이 크므로 음성예측도(NPV)와 민감도를 중시한다.")
    w()
    w("     t     예측양성   TP   FP   FN   TN    민감도   특이도    NPV    Matthews")
    w("   " + "-" * 74)

    sweep = []
    n_steps = int(round(1.0 / SWEEP_STEP)) + 1
    for k in range(n_steps):
        t = round(k * SWEEP_STEP, 2)
        pred = (prob >= t).astype(int)
        tp, fp, fn, tn = confusion(y_prob, pred)
        sens = tp / (tp + fn) if (tp + fn) else None
        spec = tn / (tn + fp) if (tn + fp) else None
        npv = tn / (tn + fn) if (tn + fn) else None
        mcc = mcc_from_counts(tp, fp, fn, tn)
        row = {
            "threshold": t, "n_pred_pos": int(tp + fp),
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "sensitivity": sens, "specificity": spec, "npv": npv, "mcc": mcc,
            "n_pred_neg": int(tn + fn),
        }
        sweep.append(row)
        w(f"   {t:4.2f}   {tp + fp:6d}  {tp:4d} {fp:4d} {fn:4d} {tn:4d}   "
          f"{'n/a  ' if sens is None else f'{sens * 100:5.1f}'}   "
          f"{'n/a  ' if spec is None else f'{spec * 100:5.1f}'}   "
          f"{'n/a  ' if npv is None else f'{npv * 100:5.1f}'}   "
          f"{'n/a' if mcc is None else f'{mcc:6.3f}'}")
    w()

    # NPV 최대화 임계값 (예측 음성이 존재하는 구간에서만)
    npv_cands = [r for r in sweep if r["npv"] is not None and r["n_pred_neg"] > 0]
    best_npv = max(npv_cands, key=lambda r: (r["npv"], -r["threshold"])) if npv_cands else None

    # 민감도 >= 95% 를 만족하는 최소(=가장 관대하지 않은 쪽) 임계값 중 최댓값
    sens_ok = [r for r in sweep if r["sensitivity"] is not None
               and r["sensitivity"] >= TARGET_SENSITIVITY]
    # 민감도는 t 에 대해 단조 감소하므로, 조건을 만족하는 t 중 가장 큰 값이
    # "민감도 95%를 유지하면서 쓸 수 있는 가장 엄격한 임계값"이다.
    sens_pick = max(sens_ok, key=lambda r: r["threshold"]) if sens_ok else None

    if best_npv:
        if best_npv["n_pred_neg"] < 0.05 * len(y_prob):
            w("  ※ 주의: NPV 를 단독으로 최대화하면 예측 음성 집합이 극도로 작아지는")
            w("     퇴화 해가 선택된다. 아래 값은 그 구간의 예측 음성 건수와 함께 읽을 것.")
        w(f"  ▸ NPV 최대 임계값        t = {best_npv['threshold']:.2f}   "
          f"NPV {best_npv['npv'] * 100:.1f}%  "
          f"(예측 음성 {best_npv['n_pred_neg']}건, 민감도 {best_npv['sensitivity'] * 100:.1f}%, "
          f"특이도 {best_npv['specificity'] * 100:.1f}%)")
    if sens_pick:
        w(f"  ▸ 민감도 ≥ {TARGET_SENSITIVITY * 100:.0f}% 를 만족하는 가장 엄격한 임계값  "
          f"t = {sens_pick['threshold']:.2f}   "
          f"민감도 {sens_pick['sensitivity'] * 100:.1f}%, "
          f"특이도 {sens_pick['specificity'] * 100:.1f}%, "
          f"NPV {sens_pick['npv'] * 100:.1f}%")
    else:
        w(f"  ▸ 민감도 ≥ {TARGET_SENSITIVITY * 100:.0f}% 를 만족하는 임계값이 "
          f"스윕 격자 위에 존재하지 않는다.")

    best_mcc = max((r for r in sweep if r["mcc"] is not None),
                   key=lambda r: r["mcc"], default=None)
    if best_mcc:
        w(f"  ▸ Matthews 계수 최대 임계값  t = {best_mcc['threshold']:.2f}   "
          f"MCC {best_mcc['mcc']:.3f}  "
          f"(민감도 {best_mcc['sensitivity'] * 100:.1f}%, "
          f"특이도 {best_mcc['specificity'] * 100:.1f}%)")
    w()
    w("  참고: 운용 임계값 0.50 은 표 1 의 Layer C 단독 성능과 일치한다.")
    w()

    J["threshold_sweep"] = sweep
    J["threshold_selection"] = {
        "max_npv": best_npv,
        "min_threshold_meeting_sensitivity_target": sens_pick,
        "target_sensitivity": TARGET_SENSITIVITY,
        "max_mcc": best_mcc,
    }

    # ─────────────────────────────────────────────────────────────────────
    # (A-6) 학습 endpoint ablation — Wilson CI (표 2)
    # ─────────────────────────────────────────────────────────────────────
    w("=" * 78)
    w("  (A-6) 학습 endpoint ablation (표 2) — 집계 혼동행렬 기반 Wilson CI")
    w("=" * 78)
    w("  근거: report_ames.txt · report_tox21.txt")
    w("  ※ 두 파일은 **집계 혼동행렬만** 담고 있고 화합물별 예측이 없다.")
    w("    따라서 여기서는 비율 지표의 Wilson CI 까지만 낸다.")
    w("    2026-09-15 — 쌍체 McNemar 검정과 쌍체 차이의 부트스트랩 CI 는")
    w("    s18_ablation_paired_tests.py 가 화합물별 예측 CSV 두 개를 직접")
    w("    읽어 산출한다. 산출물은 아래에 있다.")
    w("        논문자료/02_검증결과/ablation_paired_tests_report.txt")
    w()

    ablation = {}
    for tag, rpath in (("ames", PATH_REPORT_AMES), ("tox21", PATH_REPORT_TOX21)):
        parsed = parse_pipeline_report(rpath)
        if not parsed:
            w(f"  {tag}: {rpath.name} 를 읽을 수 없어 생략")
            ablation[tag] = None
            continue
        ablation[tag] = {}
        w(f"  ── 학습 endpoint = {tag}   ({rpath.name})")
        for key, label in (("layer_c", "Layer C 단독"), ("ensemble", "4층 앙상블"),
                           ("layer_b", "Layer B 단독")):
            if key not in parsed:
                continue
            m = metrics_from_counts(*parsed[key])
            ablation[tag][key] = m
            w(f"     {label}: TP {m['tp']} / FP {m['fp']} / FN {m['fn']} / TN {m['tn']}")
            w(f"        민감도 {fmt_pct(m['sensitivity'])} · "
              f"특이도 {fmt_pct(m['specificity'])} · NPV {fmt_pct(m['npv'])}")
            w(f"        일치율 {fmt_pct(m['concordance'])} · "
              f"Matthews 계수 {fmt_num(m['mcc'])}  (CI 없음 — 화합물별 예측 부재)")
        w()

    # 2026-09-14 추가 — 두 리포트가 같은 검증셋에서 나왔는지 확인한다.
    # s1 을 재분할한 뒤 s6(Tox21 판)을 다시 돌리지 않으면 여기서 n 이 어긋난다.
    # 예전에는 "동일 검증셋 n = 1,440" 이라고 적어두기만 해서 이 사고가
    # 조용히 지나갔다.
    ab_n = {}
    for tag in ("ames", "tox21"):
        block = ablation.get(tag)
        if block and "layer_c" in block:
            m = block["layer_c"]
            ab_n[tag] = m["tp"] + m["fp"] + m["fn"] + m["tn"]

    split_mismatch = len(ab_n) == 2 and ab_n["ames"] != ab_n["tox21"]
    if split_mismatch:
        w("  " + "!" * 70)
        w(f"  ⚠ 두 리포트의 검증셋이 다르다 — ames n = {ab_n['ames']:,} · "
          f"tox21 n = {ab_n['tox21']:,}")
        w("    분할이 다른 두 산출물을 나란히 놓은 것이므로 **표 2 를 이대로")
        w("    인용하면 안 된다.** s6_retrain_layer_c_model.py 를 현재 분할로")
        w("    재실행한 뒤 report_tox21.txt 를 갱신할 것.")
        w("  " + "!" * 70)
        w()
    elif len(ab_n) == 2:
        w(f"  두 리포트 모두 n = {ab_n['ames']:,} — 같은 검증셋이다.")
        w()

    J["ablation_split_consistent"] = (not split_mismatch) if len(ab_n) == 2 else None
    J["ablation_n"] = ab_n

    if ablation.get("ames") and ablation.get("tox21"):
        sa = ablation["ames"]["layer_c"]["sensitivity"]["value"]
        st = ablation["tox21"]["layer_c"]["sensitivity"]["value"]
        w(f"  ▸ Layer C 민감도 차 (Ames − Tox21) = {(sa - st) * 100:.1f}%p")
        if split_mismatch:
            w("    (검증셋이 달라 이 차이는 학습 endpoint 효과가 아니다 — 인용 금지)")
        else:
            w("    이 차이의 95% CI 와 쌍체 McNemar 검정은 s18 이 산출한다")
            w("    (화합물 단위 부트스트랩 · 같은 seed·횟수).")
        w("    (두 비율의 CI 는 위에 각각 표시. 쌍체 차이의 CI 는 화합물별 예측이")
        w("     없어 산출하지 않았다 — 지어내지 않고 자리표시자로 둔다.)")
        w()

    J["ablation_endpoint"] = ablation
    J["ablation_note"] = (
        "report_ames.txt / report_tox21.txt 는 집계 혼동행렬만 담고 있어 "
        "Wilson CI 는 산출 가능하나 쌍체 McNemar(비교 ③)는 산출 불가."
    )

    # ─────────────────────────────────────────────────────────────────────
    # (B) 적용범위
    # ─────────────────────────────────────────────────────────────────────
    w("=" * 78)
    w("  (B) 적용범위 (applicability domain)")
    w("=" * 78)
    w(f"  d(x) = max over t in TRAIN  Tanimoto( ECFP4(x), ECFP4(t) )")
    w(f"  TRAIN = ames_train.csv ({len(train):,}건),  "
      f"ECFP4 = Morgan 반경 {ECFP_RADIUS} · {ECFP_BITS}비트")
    w(f"  d(x) < {AD_THRESHOLD:.2f}  →  out-of-domain")
    w()

    train_fps, train_failed = build_train_fps(train["smiles"].tolist())
    w(f"  학습셋 지문 생성 : {len(train_fps):,}건 성공 · {len(train_failed)}건 파싱 실패")
    w()

    d_valid = max_tanimoto(valid["smiles"].tolist(), train_fps)
    valid_ad = valid[["smiles", "label"]].copy()
    valid_ad["max_tanimoto"] = d_valid

    # pipeline 결과와 SMILES 로 결합
    merged = pipe.merge(valid_ad[["smiles", "max_tanimoto"]], on="smiles", how="left")
    n_unmatched = int(merged["max_tanimoto"].isna().sum())
    w(f"  검증셋 ↔ 파이프라인 결과 결합 : {len(merged):,}행 중 "
      f"미결합·파싱실패 {n_unmatched}건")
    w()

    dist = distribution_summary(merged["max_tanimoto"].to_numpy(dtype=float))
    w("  ── 최대 Tanimoto 분포 (검증셋)")
    if dist:
        w(f"     n = {dist['n']:,}")
        w(f"     min {dist['min']:.4f} · p05 {dist['p05']:.4f} · p25 {dist['p25']:.4f} · "
          f"중앙값 {dist['median']:.4f}")
        w(f"     p75 {dist['p75']:.4f} · p95 {dist['p95']:.4f} · max {dist['max']:.4f} · "
          f"평균 {dist['mean']:.4f}")
    w()

    merged["in_domain"] = merged["max_tanimoto"] >= AD_THRESHOLD
    ok = merged["max_tanimoto"].notna()

    n_in = int((ok & merged["in_domain"]).sum())
    n_out = int((ok & ~merged["in_domain"]).sum())
    pos_in = float(merged.loc[ok & merged["in_domain"], "actual_ames"].mean()) if n_in else None
    pos_out = float(merged.loc[ok & ~merged["in_domain"], "actual_ames"].mean()) if n_out else None

    # 학습셋과 ECFP 지문이 완전히 일치하는 검증셋 화합물 — 누출 진단 지표
    dvals = merged["max_tanimoto"].to_numpy(dtype=float)
    n_exact = int(np.sum(dvals >= 0.999999))
    n_near = int(np.sum(dvals >= 0.90))
    w("  ── 학습셋과의 근접도 (누출 진단)")
    w(f"     d = 1.000 (ECFP 완전 일치) : {n_exact}건 "
      f"({n_exact / len(merged) * 100:.2f}%)")
    w(f"     d ≥ 0.90                   : {n_near}건 "
      f"({n_near / len(merged) * 100:.2f}%)")
    w("     ※ 스캐폴드 분할에서는 d = 1.000 이 0건이어야 정상이다. s1 이")
    w("       [3/6] 단계에서 ECFP 완전 일치 화합물을 검증셋에서 제외한다.")
    w("       0건이 아니면 s1 을 다시 돌리지 않고 이 리포트를 쓴 것이다.")
    w("       참고 — 이전 랜덤 분할에서는 187건(12.99%)이었다.")
    w()

    w(f"  ── AD 구간 분포 (임계값 {AD_THRESHOLD:.2f})")
    w(f"     in-domain  (d ≥ {AD_THRESHOLD:.2f}) : {n_in:,}건 "
      f"({n_in / len(merged) * 100:.1f}%)  ·  실제 양성비 "
      f"{'n/a' if pos_in is None else f'{pos_in * 100:.1f}%'}")
    w(f"     out-of-domain (d < {AD_THRESHOLD:.2f}) : {n_out:,}건 "
      f"({n_out / len(merged) * 100:.1f}%)  ·  실제 양성비 "
      f"{'n/a' if pos_out is None else f'{pos_out * 100:.1f}%'}")
    w()

    J["ad_validation"] = {
        "train_n": int(len(train)),
        "train_fp_ok": len(train_fps),
        "train_fp_failed": len(train_failed),
        "validation_n": int(len(merged)),
        "unmatched_or_unparsed": n_unmatched,
        "distribution": dist,
        "n_exact_ecfp_match_to_train": n_exact,
        "n_near_match_d_ge_090": n_near,
        "n_in_domain": n_in,
        "n_out_of_domain": n_out,
        "positive_rate_in_domain": pos_in,
        "positive_rate_out_of_domain": pos_out,
    }

    # AD 구간별 성능 — 세 방식 전부 (표 X)
    w("=" * 78)
    w("  (B-3) 표 X — 적용범위 구간별 성능 (세 방식)")
    w("=" * 78)

    strata = {
        "in_domain": merged[ok & merged["in_domain"]],
        "out_of_domain": merged[ok & ~merged["in_domain"]],
    }
    ad_perf = {}
    for sname, sdf in strata.items():
        ad_perf[sname] = {}
        y_s = sdf["actual_ames"].to_numpy(dtype=int)
        w(f"  ── {sname}   n = {len(sdf):,}  "
          f"(실제 양성 {int(y_s.sum())} · 음성 {int((y_s == 0).sum())})")
        for col, label in METHODS:
            rng = np.random.default_rng(SEED)
            m = metric_set(y_s, sdf[col].to_numpy(dtype=int), rng=rng)
            ad_perf[sname][col] = m
            w(f"     {label}")
            w(f"        TP {m['tp']} / FP {m['fp']} / FN {m['fn']} / TN {m['tn']}")
            w(f"        민감도 {fmt_pct(m['sensitivity'])} · "
              f"특이도 {fmt_pct(m['specificity'])}")
            w(f"        NPV {fmt_pct(m['npv'])} · PPV {fmt_pct(m['ppv'])}")
            w(f"        일치율 {fmt_pct(m['concordance'])} · "
              f"균형정확도 {fmt_num(m['balanced_accuracy'])}")
            w(f"        Matthews 계수 {fmt_num(m['mcc'])} "
              f"[{fmt_num(m['mcc_ci_low'])}–{fmt_num(m['mcc_ci_high'])}]")
        w()

    # AD 구간별 ROC/PR (Layer C 확률)
    ad_auc = {}
    for sname, sdf in strata.items():
        sub = sdf[sdf["layer_c_prob"].notna()]
        ys = sub["actual_ames"].to_numpy(dtype=int)
        if len(np.unique(ys)) < 2:
            ad_auc[sname] = {"n": int(len(sub)), "roc_auc": None, "pr_auc": None,
                             "note": "한 클래스만 존재하여 산출 불가"}
            continue
        ps = sub["layer_c_prob"].to_numpy(dtype=float)
        ad_auc[sname] = {
            "n": int(len(sub)),
            "roc_auc": float(roc_auc_score(ys, ps)),
            "pr_auc": float(average_precision_score(ys, ps)),
            "prevalence": float(ys.mean()),
        }
    w("  ── AD 구간별 Layer C 임계값 비의존 성능")
    for sname, v in ad_auc.items():
        if v.get("roc_auc") is None:
            w(f"     {sname}: {v.get('note')}")
        else:
            w(f"     {sname}: n = {v['n']:,} · ROC-AUC {v['roc_auc']:.4f} · "
              f"PR-AUC {v['pr_auc']:.4f} (기준선 {v['prevalence']:.4f})")
    w()

    J["ad_stratified_performance"] = ad_perf
    J["ad_stratified_auc"] = ad_auc

    # ─────────────────────────────────────────────────────────────────────
    # (B-5) 사례연구 구조의 AD
    # ─────────────────────────────────────────────────────────────────────
    w("=" * 78)
    w("  (B-5) 표 Y — 사례연구 생성 구조의 적용범위")
    w("=" * 78)
    w("  학습은 Ames 실측 화합물이지만 적용 대상은 SMIRKS 로 갓 생성된")
    w("  NDSRI · 부형제 부가물이다. 두 화학공간의 거리가 본 연구의 핵심 한계다.")
    w()

    case_groups = []

    ndsri_rows = parse_ndsri_smiles(PATH_NDSRI_REPORT)
    if ndsri_rows:
        case_groups.append((
            "3.5절 니트로사민 (case_study_report.txt)",
            [r["smiles"] for r in ndsri_rows],
        ))

    if PATH_EXCIP.exists():
        exc = pd.read_csv(PATH_EXCIP)
        # SMILES 컬럼 확인: 'smiles' = 생성된 불순물, 'api_smiles' = 원료의약품
        case_groups.append((
            "3.6절 부형제 부가물 (excipient_case_study_results.csv, 컬럼 'smiles')",
            exc["smiles"].tolist(),
        ))

    if PATH_BINARY.exists():
        binp = pd.read_csv(PATH_BINARY)
        smi_col = "Impurity_SMILES"
        if smi_col in binp.columns:
            vals = binp[smi_col].dropna().astype(str).tolist()
            vals = [v for v in vals if v.strip() and v.strip() != "-"]
            if vals:
                case_groups.append((
                    f"3.7절 이원 배합 산물 (binary_compatibility_pairs.csv, 컬럼 '{smi_col}')",
                    vals,
                ))

    case_results = {}
    for gname, smi_list in case_groups:
        d = max_tanimoto(smi_list, train_fps)
        n_parsed = int(np.sum(~np.isnan(d)))
        n_fail = int(np.sum(np.isnan(d)))
        dv = d[~np.isnan(d)]
        n_in_c = int(np.sum(dv >= AD_THRESHOLD))
        n_out_c = int(np.sum(dv < AD_THRESHOLD))
        summ = distribution_summary(d)
        uniq = len(set(smi_list))
        case_results[gname] = {
            "n_structures": len(smi_list),
            "n_unique_smiles": uniq,
            "n_parsed": n_parsed,
            "n_parse_failed": n_fail,
            "n_in_domain": n_in_c,
            "n_out_of_domain": n_out_c,
            "fraction_out_of_domain": (n_out_c / n_parsed) if n_parsed else None,
            "distribution": summ,
        }
        w(f"  ── {gname}")
        w(f"     구조 {len(smi_list)}건 (고유 SMILES {uniq}건) · 파싱 성공 {n_parsed} · 실패 {n_fail}")
        w(f"     in-domain (d ≥ {AD_THRESHOLD:.2f})    {n_in_c}건")
        w(f"     out-of-domain (d < {AD_THRESHOLD:.2f}) {n_out_c}건  "
          f"({'n/a' if not n_parsed else f'{n_out_c / n_parsed * 100:.1f}%'})")
        if summ:
            w(f"     d  min {summ['min']:.4f} · 중앙값 {summ['median']:.4f} · "
              f"max {summ['max']:.4f}")
        w()

    J["ad_case_study"] = case_results

    # ─────────────────────────────────────────────────────────────────────
    # 출력 파일
    # ─────────────────────────────────────────────────────────────────────
    out_cols = ["smiles", "max_tanimoto", "in_domain", "actual_ames",
                "pred_ensemble", "pred_layer_b", "pred_layer_c", "layer_c_prob"]
    ad_csv = merged[out_cols].copy()
    ad_csv["in_domain"] = ad_csv["in_domain"].astype(int)
    ad_csv.to_csv(OUT_CSV, index=False, encoding="utf-8")

    w("=" * 78)
    w("  산출 파일")
    w("=" * 78)
    w(f"  {OUT_TXT.relative_to(ROOT)}")
    w(f"  {OUT_JSON.relative_to(ROOT)}")
    w(f"  {OUT_CSV.relative_to(ROOT)}   ({len(ad_csv):,}행)")
    w()
    w("  ※ 본 스크립트는 읽기 전용이며 기존 산출물을 덮어쓰지 않는다.")
    w("  ※ 재실행 순서: s1 → p61 → s2 → s3 → s4 → s5 → **s11**")
    w()

    OUT_TXT.write_text("\n".join(L) + "\n", encoding="utf-8")
    OUT_JSON.write_text(
        json.dumps(J, ensure_ascii=False, indent=2, allow_nan=True),
        encoding="utf-8",
    )

    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
