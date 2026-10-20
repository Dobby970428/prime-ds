"""
외부 도구(Derek Nexus · Sarah Nexus · Toxtree) 결과 채점 및 PRIME-DS 대조.

논문 초안 2.12절의 채점 관례를 그대로 구현한다.

    Derek   : 추론 수준이 EQUIVOCAL 이상이면 양성.
              CERTAIN · PROBABLE · PLAUSIBLE · EQUIVOCAL → 양성
              DOUBTED · IMPROBABLE · IMPOSSIBLE · INACTIVE · (판정 없음) → 음성
    Sarah   : 전체 정성 판정을 사용. equivocal 및 적용범위 밖은
              음성으로 편입하지 않고 별도로 집계한 뒤, 결과를 두 번 보고한다.
                (a) equivocal 을 음성으로 처리
                (b) equivocal 을 제외하고 커버리지를 명시
    결합    : 어느 한쪽이라도 양성이면 양성 (ICH M7 요건에 대응)

입력
    논문자료/02_검증결과/external_tools/ames_validation_key.csv      정답 (PRIME_ID ↔ actual_ames)
    논문자료/02_검증결과/replay_fixed_engine_results.csv              PRIME-DS 화합물별 판정
    외부 도구 결과 파일 (TSV/CSV) — 아래 --derek / --sarah 로 지정

출력
    논문자료/02_검증결과/external_tools/external_comparison_report.txt
    논문자료/02_검증결과/external_tools/external_comparison_merged.csv

사용법
    python s14_score_external_tools.py --derek "Data Export (Delimited) ....tsv"
    python s14_score_external_tools.py --derek derek.tsv --sarah sarah.tsv

지표 정의는 s11_stats_ad_analysis.py 와 동일하다 (Cooper et al. 1979, Wilson CI).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

APP_DIR = Path(__file__).resolve().parent
EXT_DIR = APP_DIR / "논문자료" / "02_검증결과" / "external_tools"
KEY_CSV = EXT_DIR / "ames_validation_key.csv"
PRIME_CSV = APP_DIR / "논문자료" / "02_검증결과" / "replay_fixed_engine_results.csv"

OUT_TXT = EXT_DIR / "external_comparison_report.txt"
OUT_CSV = EXT_DIR / "external_comparison_merged.csv"

# 2.12절 — Derek 채점 경계
DEREK_POSITIVE = {"CERTAIN", "PROBABLE", "PLAUSIBLE", "EQUIVOCAL"}
DEREK_NEGATIVE = {"DOUBTED", "IMPROBABLE", "IMPOSSIBLE", "INACTIVE", "CONTRADICTED", ""}

_log: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


# ---------------------------------------------------------------------------
# 지표
# ---------------------------------------------------------------------------
def wilson(k: int, n: int, z: float = 1.959963985):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, centre - half), min(1.0, centre + half)


def metrics(y_true, y_pred) -> dict:
    tp = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    tn = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 0)
    fp = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)

    def ratio(a, b):
        return a / b if b else float("nan")

    sens = ratio(tp, tp + fn)
    spec = ratio(tn, tn + fp)
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "n": len(y_true),
        "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "sens": sens, "sens_ci": wilson(tp, tp + fn),
        "spec": spec, "spec_ci": wilson(tn, tn + fp),
        "conc": ratio(tp + tn, len(y_true)), "conc_ci": wilson(tp + tn, len(y_true)),
        "ppv": ratio(tp, tp + fp), "npv": ratio(tn, tn + fn),
        "fn_rate": 1 - sens if sens == sens else float("nan"),
        "mcc": ((tp * tn - fp * fn) / den) if den else float("nan"),
        "cooper": ratio(sens, 1 - spec) if spec == spec and spec != 1 else float("nan"),
    }


def show(title: str, m: dict) -> None:
    say(f"  {title}")
    say(f"      TP {m['TP']:4d}  FP {m['FP']:4d}  FN {m['FN']:4d}  TN {m['TN']:4d}   (n = {m['n']:,})")
    say(f"      민감도 {m['sens']*100:5.1f}% [{m['sens_ci'][0]*100:.1f}–{m['sens_ci'][1]*100:.1f}]"
        f"   특이도 {m['spec']*100:5.1f}% [{m['spec_ci'][0]*100:.1f}–{m['spec_ci'][1]*100:.1f}]")
    say(f"      일치율 {m['conc']*100:5.1f}% [{m['conc_ci'][0]*100:.1f}–{m['conc_ci'][1]*100:.1f}]"
        f"   위음성률 {m['fn_rate']*100:5.1f}%")
    say(f"      PPV {m['ppv']*100:5.1f}%   NPV {m['npv']*100:5.1f}%"
        f"   Matthews {m['mcc']:.3f}   Cooper ROC {m['cooper']:.2f}")
    say()


def mcnemar(y, a, b):
    """정확검정(이항, 양측). b = A만 정답, c = B만 정답."""
    from scipy.stats import binomtest

    only_a = sum(1 for t, x, z in zip(y, a, b) if x == t and z != t)
    only_b = sum(1 for t, x, z in zip(y, a, b) if x != t and z == t)
    n = only_a + only_b
    p = binomtest(only_a, n, 0.5).pvalue if n else float("nan")
    return only_a, only_b, p


def holm(pvals: list[float]) -> list[float]:
    idx = sorted(range(len(pvals)), key=lambda i: pvals[i])
    out = [float("nan")] * len(pvals)
    running = 0.0
    for rank, i in enumerate(idx):
        adj = min(1.0, pvals[i] * (len(pvals) - rank))
        running = max(running, adj)
        out[i] = running
    return out


# ---------------------------------------------------------------------------
# 도구별 파서
# ---------------------------------------------------------------------------
def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def exact_column(frame: pd.DataFrame, *names: str) -> str | None:
    """정확한 열 이름을 우선 찾는다 (대소문자·공백 무시)."""
    lookup = {c.strip().lower(): c for c in frame.columns}
    for name in names:
        hit = lookup.get(name.strip().lower())
        if hit is not None:
            return hit
    return None


def find_column(frame: pd.DataFrame, *needles: str) -> str | None:
    for col in frame.columns:
        low = col.lower()
        if all(n.lower() in low for n in needles):
            return col
    return None


def id_column(frame: pd.DataFrame) -> str | None:
    """PRIME_ID 열을 찾는다. 없으면 Compound Name / Structure ID 순."""
    return (exact_column(frame, "PRIME_ID", "PRIME ID")
            or exact_column(frame, "Compound Name", "Structure ID")
            or find_column(frame, "compound", "name"))


def parse_derek(path: Path) -> pd.DataFrame:
    frame = read_table(path)

    id_col = id_column(frame)
    reason_col = exact_column(frame, "Mutagenicity in vitro bacterium reasoning")
    if reason_col is None:
        reason_col = find_column(frame, "mutagenicity", "bacterium", "reasoning")
    if reason_col and "overview" in reason_col.lower():
        reason_col = None
    if reason_col is None:
        for col in frame.columns:
            low = col.lower()
            if "bacterium" in low and low.endswith("reasoning"):
                reason_col = col
                break
    if id_col is None or reason_col is None:
        raise SystemExit(
            f"[오류] Derek 파일에서 열을 찾지 못했습니다.\n"
            f"       식별자 열: {id_col}\n       판정 열: {reason_col}\n"
            f"       열 목록 앞부분: {list(frame.columns)[:8]}"
        )

    level = frame[reason_col].str.strip().str.upper()
    unknown = sorted(set(level) - DEREK_POSITIVE - DEREK_NEGATIVE)
    if unknown:
        say(f"  ⚠ Derek 판정값 중 분류되지 않은 값: {unknown}")
        say(f"    (음성으로 처리했습니다. 2.12절 규칙에 없는 값이면 확인 필요)")
        say()

    out = pd.DataFrame({
        "PRIME_ID": frame[id_col].str.strip(),
        "derek_level": level,
        "derek_pos": level.isin(DEREK_POSITIVE).astype(int),
    })
    return out.drop_duplicates("PRIME_ID")


def parse_toxtree(path: Path, force_col: str | None = None) -> pd.DataFrame:
    """
    Toxtree (Benigni/Bossa rulebase) 결과 파서.

    채점 (2.12절): 세균 변이원성 관련 구조 경보가 **하나라도** 있으면 양성.
    Toxtree 판본마다 열 이름이 달라서, 아래 후보를 순서대로 찾고
    실패하면 열 목록을 출력한다. --toxtree-col 로 직접 지정할 수도 있다.
    """
    frame = read_table(path)

    id_col = id_column(frame)
    if id_col is None:
        id_col = exact_column(frame, "Name", "Compound", "ID", "Molecule")

    if force_col:
        alert_col = exact_column(frame, force_col) or force_col
        if alert_col not in frame.columns:
            raise SystemExit(f"[오류] 지정한 열이 없습니다: {force_col}\n"
                             f"       열 목록: {list(frame.columns)}")
        alert_cols = [alert_col]
    else:
        # 양성 방향 열만 쓴다.
        #   · genotoxic carcinogenicity 구조 경보
        #   · "Potential ... TA100 mutagen" QSAR 양성 플래그
        # 제외해야 하는 것 —
        #   · "Unlikely to be a ... mutagen"  (음성 증거)
        #   · "Negative for ..."               (음성 증거)
        #   · "nongenotoxic"                   (비유전독성 발암성, Ames 와 무관)
        exclude = ("unlikely", "negative", "nongenotox", "applicable", "error")

        alert_cols = []
        for col in frame.columns:
            low = col.lower()
            if any(x in low for x in exclude):
                continue
            if "alert" in low and "genotox" in low:
                alert_cols.append(col)
            elif "potential" in low and "mutagen" in low:
                alert_cols.append(col)

        alert_cols = list(dict.fromkeys(alert_cols))  # 중복 제거

    if id_col is None or not alert_cols:
        raise SystemExit(
            f"[오류] Toxtree 파일에서 열을 찾지 못했습니다.\n"
            f"       식별자 열: {id_col}\n       경보 열 후보: {alert_cols}\n"
            f"       열 목록: {list(frame.columns)}\n"
            f"       → --toxtree-col \"열이름\" 으로 직접 지정하거나 열 목록을 알려주세요."
        )

    say(f"  Toxtree 경보 열: {alert_cols}")

    def is_positive(row) -> int:
        for col in alert_cols:
            val = str(row[col]).strip().lower()
            if val in ("", "no", "none", "nan", "0", "false", "-"):
                continue
            return 1
        return 0

    out = pd.DataFrame({
        "PRIME_ID": frame[id_col].astype(str).str.strip(),
        "toxtree_pos": frame.apply(is_positive, axis=1).astype(int),
    })
    return out.drop_duplicates("PRIME_ID")


def parse_sarah(path: Path) -> pd.DataFrame:
    frame = read_table(path)

    id_col = id_column(frame)
    call_col = exact_column(frame, "Sarah Result", "Result")
    if call_col is None:
        call_col = find_column(frame, "sarah", "result")
    conf_col = exact_column(frame, "Probability", "Confidence")

    if id_col is None or call_col is None:
        raise SystemExit(
            f"[오류] Sarah 파일에서 열을 찾지 못했습니다.\n"
            f"       식별자 열: {id_col}\n       판정 열: {call_col}\n"
            f"       열 목록: {list(frame.columns)}\n"
            f"       → 열 이름을 알려주시면 파서를 맞추겠습니다."
        )

    call = frame[call_col].str.strip().str.lower()
    known = {"positive", "negative", "equivocal", "outside domain", "out of domain", ""}
    unknown = sorted(set(call) - known)
    if unknown:
        say(f"  ⚠ Sarah 판정값 중 분류되지 않은 값: {unknown}")
        say()

    out = pd.DataFrame({
        "PRIME_ID": frame[id_col].str.strip(),
        "sarah_call": call,
        "sarah_conf": frame[conf_col].str.strip() if conf_col else "",
    })
    out["sarah_pos"] = call.str.startswith("positive").astype(int)
    out["sarah_equivocal"] = call.str.contains("equivocal").astype(int)
    out["sarah_ood"] = call.str.contains("domain").astype(int)

    # 학습셋 예시가 가설을 덮어쓴 건 = Sarah 가 그 구조를 이미 알고 있었다는 표시.
    pos_ov = exact_column(frame, "Positive Hypothesis Count (Overruled by Training Set Example)")
    neg_ov = exact_column(frame, "Negative Hypothesis Count (Overruled by Training Set Example)")
    if pos_ov and neg_ov:
        counts = (pd.to_numeric(frame[pos_ov], errors="coerce").fillna(0)
                  + pd.to_numeric(frame[neg_ov], errors="coerce").fillna(0))
        out["sarah_overruled"] = (counts > 0).astype(int)
    else:
        out["sarah_overruled"] = 0

    return out.drop_duplicates("PRIME_ID")


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--derek", type=str, default=None, help="Derek 결과 TSV/CSV")
    ap.add_argument("--sarah", type=str, default=None, help="Sarah 결과 TSV/CSV")
    ap.add_argument("--toxtree", type=str, default=None, help="Toxtree 결과 TSV/CSV")
    ap.add_argument("--toxtree-col", type=str, default=None,
                    help="Toxtree 경보 열 이름을 직접 지정 (자동 탐지 실패 시)")
    ap.add_argument("--key", type=str, default=None,
                    help="정답 키 CSV. 기본은 검증셋. "
                         "학습셋 결손 분석 시 ames_train_key.csv 를 지정할 것.")
    ap.add_argument("--no-prime", action="store_true",
                    help="PRIME-DS 화합물별 결과와 결합하지 않는다. "
                         "학습셋처럼 s12 산출물이 없는 경우에 쓴다.")
    ap.add_argument("--allow-missing-prime", action="store_true",
                    help="PRIME-DS 결과 결측이 5%% 를 넘어도 강행한다. "
                         "검증셋이 어긋난 채 지표가 산출될 수 있으므로 "
                         "논문용으로는 쓰지 말 것.")
    ap.add_argument("--out-prefix", type=str, default="external_comparison",
                    help="출력 파일 이름 앞부분")
    args = ap.parse_args()

    if not (args.derek or args.sarah or args.toxtree):
        ap.error("--derek, --sarah, --toxtree 중 최소 하나는 필요합니다.")

    say("=" * 78)
    say("  외부 도구 대조 — 동일 검증셋")
    say("=" * 78)
    say()

    key_path = Path(args.key) if args.key else KEY_CSV
    if not key_path.is_absolute():
        key_path = (APP_DIR / key_path).resolve()
    if not key_path.exists():
        say(f"[오류] 키 파일 없음: {key_path}")
        return 1

    key = pd.read_csv(key_path)
    say(f"정답 키 : {key_path.name}  {len(key):,}건")

    if args.no_prime:
        merged = key.copy()
        say("PRIME-DS 결합 : 건너뜀 (--no-prime)")
    else:
        if not PRIME_CSV.exists():
            say(f"[오류] PRIME-DS 결과 없음: {PRIME_CSV}")
            say("       먼저 python s12_replay_validation_fixed_engine.py 를 실행하거나,")
            say("       학습셋 분석이면 --no-prime 을 붙이세요.")
            return 1

        prime = pd.read_csv(PRIME_CSV)
        prime_cols = ["smiles", "new_layer_b_alert", "new_pred_layer_c", "new_pred_ensemble"]
        if "max_tanimoto_to_train" in prime.columns:
            prime_cols.append("max_tanimoto_to_train")
        merged = key.merge(
            prime[prime_cols],
            on="smiles", how="left", validate="one_to_one",
        )
        merged = merged.rename(columns={
            "new_layer_b_alert": "prime_rule",
            "new_pred_layer_c": "prime_stat",
            "new_pred_ensemble": "prime_ensemble",
        })
        n_missing = int(merged["prime_ensemble"].isna().sum())
        say(f"정답 + PRIME-DS 결합 : {len(merged):,}건 (결측 {n_missing:,}건)")

        # 2026-09-14 추가 — 결측이 많으면 멈춘다.
        #
        # s1 을 스캐폴드 분할로 재실행한 뒤 s12 를 다시 돌리지 않아,
        # replay_fixed_engine_results.csv 가 구 검증셋(1,440건)인 채로
        # 새 정답 키(1,288건)와 병합된 일이 있었다. SMILES 가 우연히 겹친
        # 229건만 붙었는데도 그대로 지표를 산출해 버렸다.
        # 결측을 경고만 하고 넘어가면 «틀린 수치»가 조용히 리포트에 실린다.
        if n_missing:
            ratio = n_missing / len(merged)
            say()
            say("  " + "!" * 74)
            say(f"  ⚠ PRIME-DS 결과가 {n_missing:,}건 비어 있다 "
                f"({ratio * 100:.1f}%).")
            say(f"    근거 파일 : {PRIME_CSV.name}")
            say("    정답 키와 다른 검증셋에서 나온 파일일 가능성이 크다.")
            say("    s1 을 재분할했다면 s5 → 동기화 → s12 를 다시 돌린 뒤")
            say("    이 스크립트를 실행할 것.")
            say("  " + "!" * 74)
            say()
            if ratio > 0.05 and not args.allow_missing_prime:
                say("  결측이 5% 를 넘어 중단한다. PRIME-DS 지표를 남은 일부로")
                say("  계산하면 검증셋이 다른 수치가 만들어진다.")
                say("  외부 도구만 채점하려면 --no-prime,")
                say("  그래도 강행하려면 --allow-missing-prime 을 붙일 것.")
                return 1

    if args.derek:
        d = parse_derek(Path(args.derek))
        merged = merged.merge(d, on="PRIME_ID", how="left")
        miss = int(merged["derek_pos"].isna().sum())
        say(f"Derek 결합           : {len(d):,}건 매칭, 누락 {miss:,}건")
        if miss:
            say("  ⚠ 누락분은 채점에서 제외됩니다. PRIME_ID 가 어긋났는지 확인하세요.")

    if args.sarah:
        s = parse_sarah(Path(args.sarah))
        merged = merged.merge(s, on="PRIME_ID", how="left")
        miss = int(merged["sarah_pos"].isna().sum())
        say(f"Sarah 결합           : {len(s):,}건 매칭, 누락 {miss:,}건")

    if args.toxtree:
        t = parse_toxtree(Path(args.toxtree), args.toxtree_col)
        merged = merged.merge(t, on="PRIME_ID", how="left")
        miss = int(merged["toxtree_pos"].isna().sum())
        say(f"Toxtree 결합         : {len(t):,}건 매칭, 누락 {miss:,}건")
    say()

    y_all = merged["actual_ames"].tolist()

    has_prime = not args.no_prime

    if has_prime:
        say("=" * 78)
        say("  PRIME-DS (기준)")
        say("=" * 78)
        say()
        base = merged.dropna(subset=["prime_ensemble"])
        show("규칙 계층 (Layer B)", metrics(base["actual_ames"].tolist(), base["prime_rule"].astype(int).tolist()))
        show("통계 계층 (Layer C)", metrics(base["actual_ames"].tolist(), base["prime_stat"].astype(int).tolist()))
        show("4계층 앙상블", metrics(base["actual_ames"].tolist(), base["prime_ensemble"].astype(int).tolist()))

    comparisons: list[tuple[str, list, list, list]] = []

    if args.derek:
        say("=" * 78)
        say("  Derek Nexus  (EQUIVOCAL 이상 = 양성)")
        say("=" * 78)
        say()
        sub = merged.dropna(subset=["derek_pos"])
        show("Derek Nexus", metrics(sub["actual_ames"].tolist(), sub["derek_pos"].astype(int).tolist()))

        say("  추론 수준 분포")
        for level, cnt in sub["derek_level"].value_counts().items():
            part = sub[sub["derek_level"] == level]
            pos = int((part["actual_ames"] == 1).sum())
            say(f"      {level:<14} {cnt:5,}건   (실제 양성 {pos:,} / 음성 {cnt-pos:,})")
        say()

        if has_prime:
            pair = merged.dropna(subset=["derek_pos", "prime_rule"])
            comparisons.append((
                "PRIME-DS 규칙 계층 vs Derek Nexus",
                pair["actual_ames"].tolist(),
                pair["prime_rule"].astype(int).tolist(),
                pair["derek_pos"].astype(int).tolist(),
            ))

    if args.sarah:
        say("=" * 78)
        say("  Sarah Nexus")
        say("=" * 78)
        say()
        sub = merged.dropna(subset=["sarah_pos"])
        n_eq = int(sub["sarah_equivocal"].sum())
        n_ood = int(sub["sarah_ood"].sum())
        say(f"  equivocal {n_eq:,}건 · 적용범위 밖 {n_ood:,}건 "
            f"(커버리지 {(len(sub)-n_eq-n_ood)/len(sub)*100:.1f}%)")
        say()

        show("(a) equivocal 을 음성으로 처리",
             metrics(sub["actual_ames"].tolist(), sub["sarah_pos"].astype(int).tolist()))

        keep = sub[(sub["sarah_equivocal"] == 0) & (sub["sarah_ood"] == 0)]
        show(f"(b) equivocal · 적용범위 밖 제외  (n = {len(keep):,})",
             metrics(keep["actual_ames"].tolist(), keep["sarah_pos"].astype(int).tolist()))

        if has_prime:
            pair = merged.dropna(subset=["sarah_pos", "prime_stat"])
            comparisons.append((
                "PRIME-DS 통계 계층 vs Sarah Nexus",
                pair["actual_ames"].tolist(),
                pair["prime_stat"].astype(int).tolist(),
                pair["sarah_pos"].astype(int).tolist(),
            ))

    if args.toxtree:
        say("=" * 78)
        say("  Toxtree  (Benigni/Bossa — 구조 경보가 하나라도 있으면 양성)")
        say("=" * 78)
        say()
        sub = merged.dropna(subset=["toxtree_pos"])
        show("Toxtree", metrics(sub["actual_ames"].tolist(), sub["toxtree_pos"].astype(int).tolist()))

        if has_prime:
            pair = merged.dropna(subset=["toxtree_pos", "prime_rule"])
            comparisons.append((
                "PRIME-DS 규칙 계층 vs Toxtree",
                pair["actual_ames"].tolist(),
                pair["prime_rule"].astype(int).tolist(),
                pair["toxtree_pos"].astype(int).tolist(),
            ))

    if args.derek and args.sarah:
        say("=" * 78)
        say("  Derek + Sarah 결합  (어느 한쪽이라도 양성이면 양성)")
        say("=" * 78)
        say()
        sub = merged.dropna(subset=["derek_pos", "sarah_pos"])
        combo = ((sub["derek_pos"].astype(int) + sub["sarah_pos"].astype(int)) > 0).astype(int)
        show("Derek + Sarah", metrics(sub["actual_ames"].tolist(), combo.tolist()))

        if has_prime:
            # 2026-09-14 수정 — 쌍체 비교에는 PRIME 쪽 결측도 함께 떨어뜨려야
            # 한다. 458·475행은 그렇게 하는데 여기만 빠져 있어, PRIME 결과가
            # 하나라도 비면 astype(int) 에서 IntCastingNaNError 로 죽었다.
            pair = sub.dropna(subset=["prime_ensemble"])
            if len(pair) < len(sub):
                say(f"  ※ PRIME-DS 결과가 없는 {len(sub) - len(pair):,}건은 "
                    f"쌍체 비교에서 제외했다 (n = {len(pair):,}).")
                say()
            combo_pair = ((pair["derek_pos"].astype(int)
                           + pair["sarah_pos"].astype(int)) > 0).astype(int)
            comparisons.append((
                "PRIME-DS 앙상블 vs Derek + Sarah",
                pair["actual_ames"].tolist(),
                pair["prime_ensemble"].astype(int).tolist(),
                combo_pair.tolist(),
            ))

    # ------------------------------------------------------------------
    # 학습셋 노출 보정 — 세 시스템 모두 이 검증셋을 어느 정도 «본» 상태다.
    # ------------------------------------------------------------------
    if args.sarah and "sarah_overruled" in merged.columns:
        say("=" * 78)
        say("  학습셋 노출 — Sarah 의 가설이 학습셋 예시로 덮어써진 경우")
        say("=" * 78)
        say()
        sub = merged.dropna(subset=["sarah_pos"])
        n_ov = int(sub["sarah_overruled"].sum())
        say(f"  덮어써진 화합물 : {n_ov:,} / {len(sub):,}건 ({n_ov/len(sub)*100:.1f}%)")
        say("  이 구조들에 대해 Sarah 는 예측이 아니라 학습셋 조회로 답한 것에 가깝다.")
        say()

        for flag, label in ((1, "학습셋 예시 있음 (기억)"), (0, "학습셋 예시 없음 (예측)")):
            part = sub[sub["sarah_overruled"] == flag]
            if len(part):
                show(f"Sarah — {label}  (n = {len(part):,})",
                     metrics(part["actual_ames"].tolist(),
                             part["sarah_pos"].astype(int).tolist()))

    if has_prime and args.sarah and args.derek and "sarah_overruled" in merged.columns:
        say("=" * 78)
        say("  참고 — 양쪽 모두 «처음 보는» 부분집합")
        say("=" * 78)
        say()
        say("  Sarah 학습셋 예시가 없고, PRIME-DS 학습셋과 지문이 완전히 일치하지도")
        say("  않는 화합물만 남긴 것이다. 세 시스템의 노출을 동시에 줄인 조건이므로")
        say("  가장 공정에 가까운 비교이나, 표본이 줄고 난이도가 높은 쪽으로 치우친다.")
        say()

        # 2026-09-14 수정 — 아래에서 prime_* 세 열을 astype(int) 하므로
        # PRIME 쪽 결측도 함께 떨어뜨린다. 그러지 않으면 여기서도 죽는다.
        drop_cols = ["sarah_pos", "derek_pos"]
        if has_prime:
            drop_cols += ["prime_rule", "prime_stat", "prime_ensemble"]
        clean = merged.dropna(subset=drop_cols)
        clean = clean[clean["sarah_overruled"] == 0]
        if "max_tanimoto_to_train" in clean.columns:
            clean = clean[clean["max_tanimoto_to_train"] < 1.0]
        say(f"  부분집합 n = {len(clean):,} "
            f"(양성 {int((clean['actual_ames']==1).sum()):,} / "
            f"음성 {int((clean['actual_ames']==0).sum()):,})")
        say()

        if len(clean) >= 30:
            y = clean["actual_ames"].tolist()
            show("PRIME-DS 규칙 계층", metrics(y, clean["prime_rule"].astype(int).tolist()))
            show("Derek Nexus", metrics(y, clean["derek_pos"].astype(int).tolist()))
            show("PRIME-DS 통계 계층", metrics(y, clean["prime_stat"].astype(int).tolist()))
            show("Sarah Nexus", metrics(y, clean["sarah_pos"].astype(int).tolist()))
            show("PRIME-DS 앙상블", metrics(y, clean["prime_ensemble"].astype(int).tolist()))
            combo = ((clean["derek_pos"].astype(int) + clean["sarah_pos"].astype(int)) > 0).astype(int)
            show("Derek + Sarah", metrics(y, combo.tolist()))
        else:
            say("  ⚠ 표본이 30건 미만이라 지표를 산출하지 않았다.")
            say()

    if comparisons:
        say("=" * 78)
        say("  McNemar 쌍체 비교 (정확검정, Holm 보정) — 동종 대 동종")
        say("=" * 78)
        say()
        rows = []
        for name, y, a, b in comparisons:
            only_a, only_b, p = mcnemar(y, a, b)
            rows.append((name, only_a, only_b, p))
        adj = holm([r[3] for r in rows])
        say(f"  {'비교 (A 대 B)':<40} {'b':>5} {'c':>5} {'p':>12} {'p(Holm)':>12}")
        say("  " + "-" * 76)
        for (name, b_, c_, p), pa in zip(rows, adj):
            say(f"  {name:<40} {b_:5d} {c_:5d} {p:12.3e} {pa:12.3e}")
        say()
        say("  b = A 만 맞힌 건수, c = B 만 맞힌 건수.")
        say()

    EXT_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = EXT_DIR / f"{args.out_prefix}_merged.csv"
    out_txt = EXT_DIR / f"{args.out_prefix}_report.txt"
    merged.to_csv(out_csv, index=False, encoding="utf-8-sig")
    out_txt.write_text("\n".join(_log), encoding="utf-8")
    print(f"\n상세 결과 : {out_csv}")
    print(f"요약 기록 : {out_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
