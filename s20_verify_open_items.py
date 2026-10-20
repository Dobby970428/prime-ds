"""
투고 전 미확인 항목 점검.

원고에 숫자로 들어가는데 아직 근거를 확인하지 않은 것들을 한 번에 센다.
추측으로 쓰지 않기 위한 스크립트다.

① 발암성 양성 레코드 — 20건인가 19종인가
    s3_fix_layer_a_carcinogens 가 Azobenzene 의 잘못된 구조를 4-Aminobiphenyl
    로 교체하는데, 4-Aminobiphenyl 레코드가 이미 있어 **같은 구조가 두 줄**이
    된다. 원고가 「발암성 양성 20건」 이라고 쓰면 고유 화합물 수와 다르다.

② Derek 추론 수준 빈칸 — 몇 건이고 무엇인가
    스캐폴드 검증셋에서 15건(실제 양성 4 / 음성 11)이 비어 있고 음성으로
    채점된다. 처리 실패인지 「보고할 것 없음」 인지 확인한다. Derek 은
    Show Negative Predictions 를 켜면 INACTIVE 를 명시하므로, 빈칸은
    INACTIVE 와 다른 상태다.

③ 구조 경보 개수 — 23종이 맞는가
    ASHBY_TENNANT_ALERTS 의 SMARTS 수 + 알고리즘 PAH 2종.

④ 부형제 라이브러리 구성 — 143종의 출처 내역
    원고 §2.3 이 USP/NF · Handbook · FDA IID 라고만 적고 있다.

출력
    논문자료/02_검증결과/open_items_check.txt

사용법
    py s20_verify_open_items.py
    py s20_verify_open_items.py --derek "논문자료\\02_검증결과\\external_tools\\scaffold_derek.tsv"
"""

from __future__ import annotations

import argparse
import sqlite3
from collections import Counter
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP / "논문자료"
EXT = ROOT / "02_검증결과" / "external_tools"
DB = APP / "cache" / "m7_empirical.db"
KEY = EXT / "ames_validation_key.csv"
OUT_TXT = ROOT / "02_검증결과" / "open_items_check.txt"

_log: list[str] = []


def w(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


# ---------------------------------------------------------------------------
def check_carcinogens() -> None:
    w("=" * 78)
    w("  ① 발암성 양성 레코드 — 건수 대 고유 구조 수")
    w("=" * 78)
    w()

    if not DB.exists():
        w(f"  [건너뜀] 실측 DB 없음: {DB}")
        w()
        return

    conn = sqlite3.connect(DB)
    try:
        cur = conn.execute("PRAGMA table_info(m7_empirical)")
        cols = {r[1] for r in cur.fetchall()}
        if not cols:
            for cand in ("records", "empirical", "compounds"):
                cur = conn.execute(f"PRAGMA table_info({cand})")
                cols = {r[1] for r in cur.fetchall()}
                if cols:
                    table = cand
                    break
            else:
                cur = conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'")
                w(f"  [오류] 표를 찾지 못했다. 있는 표: {[r[0] for r in cur.fetchall()]}")
                w()
                return
        else:
            table = "m7_empirical"

        carc_col = next((c for c in cols if "carcinog" in c.lower()), None)
        smi_col = next((c for c in cols if "smiles" in c.lower()), None)
        name_col = next((c for c in cols if "name" in c.lower()), None)
        if not (carc_col and smi_col):
            w(f"  [오류] 필요한 열이 없다. 있는 열: {sorted(cols)}")
            w()
            return

        rows = conn.execute(
            f"SELECT {name_col or 'rowid'}, {smi_col}, {carc_col} FROM {table} "
            f"WHERE LOWER({carc_col}) LIKE '%pos%'"
        ).fetchall()
    finally:
        conn.close()

    w(f"  발암성 양성 레코드 : {len(rows)}건")
    smis = [r[1] for r in rows]
    uniq = sorted(set(smis))
    w(f"  고유 구조          : {len(uniq)}종")
    w()

    dupes = {s: c for s, c in Counter(smis).items() if c > 1}
    if dupes:
        w("  ⚠ 같은 구조가 여러 줄인 것")
        for smi, n in dupes.items():
            names = [str(r[0]) for r in rows if r[1] == smi]
            w(f"      {smi}   {n}줄 — {' / '.join(names)}")
        w()
        w("  → 원고에 「발암성 양성 N건」 을 쓸 때 레코드 수인지 고유 구조 수인지")
        w("    명시할 것. s3_fix_layer_a_carcinogens 가 Azobenzene 의 잘못된")
        w("    구조를 4-Aminobiphenyl 로 교체하면서 생긴 중복이다.")
    else:
        w("  ✅ 중복 없음 — 레코드 수와 고유 구조 수가 같다.")
    w()


def check_derek(path: Path | None) -> None:
    w("=" * 78)
    w("  ② Derek 추론 수준 빈칸")
    w("=" * 78)
    w()

    if path is None or not path.exists():
        w(f"  [건너뜀] Derek 파일이 없다: {path}")
        w("           --derek 로 지정할 것.")
        w()
        return

    import pandas as pd

    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    def find(substr: str, *, ends: str | None = None) -> str | None:
        for c in frame.columns:
            lc = c.strip().lower()
            if substr in lc and (ends is None or lc.endswith(ends)):
                return c
        return None

    reason_col = find("bacterium reasoning")
    id_col = next((c for c in frame.columns
                   if c.strip().lower() in ("prime_id", "compound name")), None)
    alert_cnt = find("bacterium alert count")
    unfired = find("unfired endpoints")

    if reason_col is None:
        w(f"  [오류] 추론 수준 열을 찾지 못했다: {list(frame.columns)[:10]}")
        w()
        return

    blank = frame[frame[reason_col].str.strip() == ""]
    w(f"  전체 {len(frame):,}건 중 추론 수준이 빈 것 : {len(blank):,}건")
    w()

    if blank.empty:
        w("  ✅ 빈칸 없음.")
        w()
        return

    w(f"  {'PRIME_ID':<14}{'Alert count':>12}  Unfired endpoints (앞부분)")
    w("  " + "-" * 70)
    for _, r in blank.iterrows():
        pid = str(r[id_col]).strip() if id_col else "?"
        ac = str(r[alert_cnt]).strip() if alert_cnt else ""
        uf = str(r[unfired]).strip()[:40] if unfired else ""
        w(f"  {pid:<14}{ac:>12}  {uf}")
    w()

    # Alert count 가 비어 있으면 «처리 자체가 안 됨», 0 이면 «경보 없음».
    if alert_cnt:
        empty_ac = int((blank[alert_cnt].str.strip() == "").sum())
        zero_ac = int((blank[alert_cnt].str.strip() == "0").sum())
        w(f"  Alert count 가 빈 것 : {empty_ac}건   0 인 것 : {zero_ac}건")
        w()
        if empty_ac == len(blank):
            w("  → Alert count 까지 비어 있다. **구조 처리 자체가 실패한 것**으로")
            w("    보인다. 원고에 «Derek 이 N건을 처리하지 못했다» 고 적고")
            w("    채점에서 제외하거나 음성 처리했음을 밝힐 것.")
        elif zero_ac == len(blank):
            w("  → Alert count 가 0 이다. 경보가 없어 추론 수준이 비었을 뿐이며")
            w("    사실상 INACTIVE 와 같다. 음성 채점이 타당하다.")
        else:
            w("  → 두 경우가 섞여 있다. 개별 확인이 필요하다.")
    w()

    if KEY.exists():
        import pandas as pd2  # noqa: F401
        key = pd.read_csv(KEY)
        if id_col and "PRIME_ID" in key.columns:
            merged = blank[[id_col]].rename(columns={id_col: "PRIME_ID"}).merge(
                key[["PRIME_ID", "actual_ames"]], on="PRIME_ID", how="left")
            pos = int((merged["actual_ames"] == 1).sum())
            neg = int((merged["actual_ames"] == 0).sum())
            w(f"  실제 라벨 — 양성 {pos}건 · 음성 {neg}건")
            w(f"  음성으로 채점하면 위음성 {pos}건이 추가된다"
              f" (민감도에 {pos}/588 = {pos / 588 * 100:.1f}%p 영향).")
            w()


def check_alert_count() -> None:
    w("=" * 78)
    w("  ③ 구조 경보 개수")
    w("=" * 78)
    w()
    try:
        from p62_toxicity_alerts import ASHBY_TENNANT_ALERTS
        from p65_m7_ensemble_matrix import METABOLIZATION_FRAGMENT_SMARTS
    except Exception as exc:  # noqa: BLE001
        w(f"  [오류] import 실패: {exc}")
        w()
        return

    n_smarts = len(ASHBY_TENNANT_ALERTS)
    w(f"  ASHBY_TENNANT_ALERTS (SMARTS)      : {n_smarts}종")
    w("  알고리즘 PAH (_scan_pah, SA18·SA19) : 2종")
    w(f"  → Layer B 구조 경보 합계            : {n_smarts + 2}종")
    w()
    w(f"  METABOLIZATION_FRAGMENT_SMARTS     : {len(METABOLIZATION_FRAGMENT_SMARTS)}종")
    w()

    # 대사 단편 중 구조 경보와 SMARTS 가 같은 것
    alert_map = {smarts: name for name, smarts in ASHBY_TENNANT_ALERTS}
    dup = [(mname, alert_map[ms]) for mname, ms in METABOLIZATION_FRAGMENT_SMARTS
           if ms in alert_map]
    if dup:
        w("  ⚠ 대사 단편과 구조 경보가 **같은 SMARTS** 인 쌍")
        for mname, aname in dup:
            w(f"      {mname}")
            w(f"        = {aname}")
        w()
        w("  → 판정에는 영향이 없다(둘 다 양성으로 셈). 다만 표 S1 에 같은 패턴이")
        w("    두 줄 실리므로 각주로 밝힐 것. 건수와 PPV 가 소수점까지 일치하는")
        w("    것이 증거다.")
    else:
        w("  ✅ 중복 SMARTS 없음.")
    w()


def check_excipients() -> None:
    w("=" * 78)
    w("  ④ 부형제 라이브러리 구성")
    w("=" * 78)
    w()
    try:
        from p72_excipient_interaction_engine import EXCIPIENT_LIBRARY
    except Exception as exc:  # noqa: BLE001
        w(f"  [오류] import 실패: {exc}")
        w()
        return

    w(f"  총 {len(EXCIPIENT_LIBRARY):,}종")
    w()

    for field in ("source", "origin", "reference", "category"):
        vals = Counter(str((v or {}).get(field, "")).strip()
                       for v in EXCIPIENT_LIBRARY.values())
        vals.pop("", None)
        if vals:
            w(f"  '{field}' 별 분포")
            for k, n in vals.most_common(20):
                w(f"      {k[:52]:<54} {n:4,}종")
            w()

    sample = next(iter(EXCIPIENT_LIBRARY.values()), None)
    if isinstance(sample, dict):
        w(f"  항목이 가진 필드 : {sorted(sample.keys())}")
        w()
        w("  → 출처를 나타내는 필드가 없으면 원고 §2.3 의 «USP/NF · Handbook ·")
        w("    FDA IID» 내역을 코드에서 셀 수 없다. 라이브러리를 만든 근거 문서를")
        w("    찾거나, 출처 필드를 추가한 뒤 다시 셀 것.")
    w()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--derek", default=str(EXT / "scaffold_derek.tsv"))
    args = ap.parse_args()

    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    w()
    w("투고 전 미확인 항목 점검 — 추측 대신 코드와 데이터에서 센다")
    w()

    check_carcinogens()
    check_derek(Path(args.derek))
    check_alert_count()
    check_excipients()

    OUT_TXT.parent.mkdir(parents=True, exist_ok=True)
    OUT_TXT.write_text("\n".join(_log), encoding="utf-8")
    print(f"\n기록 : {OUT_TXT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
