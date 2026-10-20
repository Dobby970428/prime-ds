"""
Layer B 구조경보 보강을 위한 결손 분석.

목적
    Layer B(Ashby–Tennant 18종)가 놓쳤는데 Toxtree(Benigni/Bossa) 또는
    Derek 이 잡아낸 실제 Ames 양성 화합물을 모아, 어떤 구조 유형이
    비어 있는지 드러낸다. 경보를 «아무거나 늘리는» 대신 실측 결손을
    근거로 추가하기 위한 자료다.

입력
    논문자료/02_검증결과/external_tools/external_comparison_merged.csv
        (s14_score_external_tools.py 가 만든다. 먼저 실행할 것)
    Derek 결과 TSV  — 경보 이름을 가져오기 위해 선택적으로 지정

출력
    논문자료/02_검증결과/external_tools/layer_b_gap_report.txt
    논문자료/02_검증결과/external_tools/layer_b_gap_compounds.csv
        Layer B 위음성 + 타 도구 양성 화합물 목록 (SMILES 포함)

사용법
    py s15_layer_b_gap_analysis.py
    py s15_layer_b_gap_analysis.py --derek "Data Export ....tsv"
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd

APP_DIR = Path(__file__).resolve().parent
EXT_DIR = APP_DIR / "논문자료" / "02_검증결과" / "external_tools"
MERGED = EXT_DIR / "external_comparison_merged.csv"
OUT_TXT = EXT_DIR / "layer_b_gap_report.txt"
OUT_CSV = EXT_DIR / "layer_b_gap_compounds.csv"

_log: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


def read_derek_alerts(path: Path) -> pd.DataFrame:
    """Derek 결과에서 PRIME_ID 별로 걸린 경보 이름을 모은다."""
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    id_col = None
    for cand in ("PRIME_ID", "Compound Name"):
        for col in frame.columns:
            if col.strip().lower() == cand.lower():
                id_col = col
                break
        if id_col:
            break
    if id_col is None:
        raise SystemExit(f"[오류] Derek 파일에서 식별자 열을 찾지 못했습니다: {list(frame.columns)[:8]}")

    alert_cols = [c for c in frame.columns
                  if "bacterium" in c.lower() and "alert" in c.lower() and c.lower().endswith("name")]
    if not alert_cols:
        raise SystemExit("[오류] Derek 경보 이름 열을 찾지 못했습니다.")

    def collect(row) -> str:
        names = [str(row[c]).strip() for c in alert_cols if str(row[c]).strip()]
        return "; ".join(dict.fromkeys(names))

    return pd.DataFrame({
        "PRIME_ID": frame[id_col].str.strip(),
        "derek_alerts": frame.apply(collect, axis=1),
    }).drop_duplicates("PRIME_ID")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--derek", type=str, default=None,
                    help="Derek 결과 TSV (경보 이름을 함께 보려면 지정)")
    ap.add_argument("--merged", type=str, default=None,
                    help="s14 산출 merged CSV. 기본은 검증셋. "
                         "학습셋이면 train_comparison_merged.csv 를 지정.")
    ap.add_argument("--layer-b", type=str, default=None,
                    help="Layer B 판정 CSV (s16 산출). merged 에 prime_rule 이 "
                         "없을 때 여기서 가져온다.")
    ap.add_argument("--out-prefix", type=str, default="layer_b_gap",
                    help="출력 파일 이름 앞부분")
    args = ap.parse_args()

    merged_path = Path(args.merged) if args.merged else MERGED
    if not merged_path.is_absolute():
        merged_path = (APP_DIR / merged_path).resolve()
    if not merged_path.exists():
        say(f"[오류] 먼저 s14 를 실행하세요. 없는 파일: {merged_path}")
        return 1

    df = pd.read_csv(merged_path)

    # Layer B 판정이 없으면 s16 산출물에서 붙인다 (학습셋 경로).
    if "prime_rule" not in df.columns:
        lb_path = Path(args.layer_b) if args.layer_b else (
            EXT_DIR / "layer_b_trainset_calls.csv")
        if not lb_path.is_absolute():
            lb_path = (APP_DIR / lb_path).resolve()
        if not lb_path.exists():
            say(f"[오류] Layer B 판정 파일이 없습니다: {lb_path}")
            say("       py s16_layer_b_on_trainset.py 를 먼저 실행하세요.")
            return 1
        lb = pd.read_csv(lb_path)[["PRIME_ID", "prime_rule", "alerts"]]
        lb = lb.rename(columns={"alerts": "layer_b_alerts"})
        df = df.merge(lb, on="PRIME_ID", how="left")
        say(f"Layer B 판정 결합 : {lb_path.name}  {len(lb):,}건")
        say()

    needed = {"actual_ames", "prime_rule"}
    if not needed.issubset(df.columns):
        say(f"[오류] 필요한 열이 없습니다: {needed - set(df.columns)}")
        return 1

    if args.derek:
        df = df.merge(read_derek_alerts(Path(args.derek)), on="PRIME_ID", how="left")

    say("=" * 78)
    say("  Layer B 결손 분석 — 어떤 구조를 놓치고 있는가")
    say("=" * 78)
    say()

    total_pos = int((df["actual_ames"] == 1).sum())
    fn = df[(df["actual_ames"] == 1) & (df["prime_rule"] == 0)]
    say(f"검증셋 실제 양성        : {total_pos:,}건")
    say(f"Layer B 위음성 (놓침)   : {len(fn):,}건  ({len(fn)/total_pos*100:.1f}%)")
    say()

    # ── 다른 도구는 잡았는가 ────────────────────────────────────────
    say("-" * 78)
    say("  놓친 것 중 다른 도구가 잡아낸 비율")
    say("-" * 78)
    say()

    for col, label in (("toxtree_pos", "Toxtree"), ("derek_pos", "Derek"),
                       ("sarah_pos", "Sarah")):
        if col not in df.columns:
            continue
        caught = int((fn[col] == 1).sum())
        say(f"  {label:<10} {caught:4,} / {len(fn):,}건 잡음 "
            f"({caught/len(fn)*100:.1f}%)")
    say()

    rule_tools = [c for c in ("toxtree_pos", "derek_pos") if c in df.columns]
    if rule_tools:
        any_rule = fn[rule_tools].max(axis=1)
        n_any = int((any_rule == 1).sum())
        say(f"  규칙 기반 도구 중 하나라도 잡은 것 : {n_any:,}건 "
            f"({n_any/len(fn)*100:.1f}%)")
        say(f"  세 규칙 시스템 모두 놓친 것        : {len(fn)-n_any:,}건")
        say()
        say("  → 앞의 숫자가 «경보 보강으로 회수 가능한 상한»이다.")
        say("    뒤의 숫자는 구조경보 방식 자체의 한계에 가깝다.")
        say()

    # ── 회수 가능분의 구조 유형 ─────────────────────────────────────
    if "derek_alerts" in df.columns and rule_tools:
        say("-" * 78)
        say("  회수 가능한 위음성에서 Derek 이 지목한 경보 (빈도순)")
        say("-" * 78)
        say()
        target = fn[any_rule == 1]
        counter: Counter[str] = Counter()
        for cell in target["derek_alerts"].fillna(""):
            for name in str(cell).split(";"):
                name = name.strip()
                if name:
                    counter[name] += 1

        if counter:
            say(f"  {'경보 이름':<52} {'건수':>6}")
            say("  " + "-" * 60)
            for name, cnt in counter.most_common(30):
                say(f"  {name:<52} {cnt:6,}")
            say()
            say("  → 상위 항목이 Layer B 에 추가할 후보다.")
            say("    기존 18종과 겹치는지 p62_toxicity_alerts.py 와 대조할 것.")
        else:
            say("  (Derek 경보 이름이 비어 있다)")
        say()

    # ── 반대 방향 — Layer B 만 잡은 것 ──────────────────────────────
    if "toxtree_pos" in df.columns:
        say("-" * 78)
        say("  참고 — Layer B 만 잡은 진양성 (제거하면 손해 보는 부분)")
        say("-" * 78)
        say()
        only_b = df[(df["actual_ames"] == 1) & (df["prime_rule"] == 1)
                    & (df["toxtree_pos"] == 0)]
        say(f"  Toxtree 는 놓쳤는데 Layer B 가 잡은 양성 : {len(only_b):,}건")
        say("  → 기존 경보를 버리지 말고 «보강»해야 하는 근거다.")
        say()

    # ── 위양성 쪽 ───────────────────────────────────────────────────
    say("-" * 78)
    say("  위양성 — 특이도를 깎는 쪽")
    say("-" * 78)
    say()
    fp = df[(df["actual_ames"] == 0) & (df["prime_rule"] == 1)]
    say(f"  Layer B 위양성 : {len(fp):,}건")
    if "toxtree_pos" in df.columns:
        both = int((fp["toxtree_pos"] == 1).sum())
        say(f"    이 중 Toxtree 도 양성으로 본 것 : {both:,}건 ({both/len(fp)*100:.1f}%)")
        say(f"    Layer B 만 잘못 잡은 것          : {len(fp)-both:,}건")
        say()
        say("  → 뒤의 숫자가 «느슨한 경보»로 생긴 고유 위양성이다.")
        say("    경보별 정밀도를 재서 이 부분을 줄이면 특이도가 오른다.")
    say()

    # ── 산출 ────────────────────────────────────────────────────────
    cols = [c for c in ("PRIME_ID", "smiles", "actual_ames", "prime_rule",
                        "toxtree_pos", "derek_pos", "sarah_pos", "derek_alerts")
            if c in df.columns]
    out = fn[cols].copy()
    out.insert(0, "gap_type", "layer_b_false_negative")
    if rule_tools:
        out["recoverable"] = any_rule.values

    out_csv = EXT_DIR / f"{args.out_prefix}_compounds.csv"
    out_txt = EXT_DIR / f"{args.out_prefix}_report.txt"
    out.to_csv(out_csv, index=False, encoding="utf-8-sig")
    out_txt.write_text("\n".join(_log), encoding="utf-8")

    print(f"\n화합물 목록 : {out_csv}")
    print(f"요약 기록   : {out_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
