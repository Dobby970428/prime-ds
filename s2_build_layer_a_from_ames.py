"""
Layer A 실측 DB(cache/m7_empirical.db)를 Ames 학습셋으로 채우는 스크립트.

배경
----
Layer A는 "이 구조, 실제로 Ames 시험한 기록이 있나?"를 조회하는 층입니다.
그런데 현재 DB에는 화합물 50종만 들어 있어 사실상 작동하지 않습니다.
(drug_index.db는 13,040건이지만 이름·CID·구조만 있고 독성 결과가 없어 사용 불가)

이 스크립트는 data/ames_train.csv(Ames 라벨 보유)를 적재합니다.
검증셋(ames_validation.csv)과는 설계상 겹치지 않으므로 순환논증이 발생하지 않으며,
실행 전에 겹침 여부를 직접 확인하고 하나라도 겹치면 중단합니다.

기존 50건은 그대로 보존합니다. 이 중 일부는 발암성(carcinogenicity) 기록을 갖고 있어
ICH M7 Class 1 판정에 필요하기 때문입니다.

⚠️ 근거 계층의 성격에 대한 주의 (2026-09-03)
-------------------------------------------
여기서 적재하는 것은 **TDC가 컴파일한 Ames 데이터셋**이지, Hansen/ISSCAN 같은
1차 실험 기록이 아닙니다. 그런데 p60의 보고서 라벨은
`"Empirical DB Source Match (Hansen/ISSCAN)"` 로 **고정**되어 있어, 적재분이
Hansen/ISSCAN 실측 기록처럼 표시됩니다.

이 스크립트에서 할 수 있는 대응은 다음 두 가지뿐입니다.
  · `compound_name` 에 «1차 실험 기록이 아님» 을 명시적으로 적어 둔다
  · `source` 태그를 큐레이션 50건과 확실히 구분되게 유지한다
p60 은 편집 대상이 아니므로, 보고서 라벨 자체의 수정은 아래 권고로 남깁니다.

  [p60 권고] `EmpiricalM7Hit.to_report_dict()` 의 고정 키
             "Empirical DB Source Match (Hansen/ISSCAN)" 를
             "Empirical DB Source Match" 로 바꾸고, 출처 문자열을 그대로
             노출해 1차 기록과 컴파일 데이터셋을 구분할 것.

사용법
------
    python s2_build_layer_a_from_ames.py            # 적재
    python s2_build_layer_a_from_ames.py --restore  # 백업본으로 되돌리기

되돌리기
--------
    cache/m7_empirical.db.backup 파일을 m7_empirical.db 로 이름 바꾸면 원상복구됩니다.
"""

from __future__ import annotations

import shutil
import sqlite3
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "cache" / "m7_empirical.db"
BACKUP_PATH = BASE_DIR / "cache" / "m7_empirical.db.backup"
TRAIN_CSV = BASE_DIR / "data" / "ames_train.csv"
# 2026-09-03 수정(C-3): 동결본이 있으면 그쪽을 읽는다 (읽기 전용).
FROZEN_VALID_CSV = BASE_DIR / "data" / "frozen" / "ames_validation_v1.csv"
WORKING_VALID_CSV = BASE_DIR / "data" / "ames_validation.csv"
VALID_CSV = FROZEN_VALID_CSV if FROZEN_VALID_CSV.exists() else WORKING_VALID_CSV

SOURCE_TAG = "TDC-AMES (training split)"
# 2026-09-03 추가(C-11): 적재분이 1차 실험 기록이 아님을 레코드 자체에 남긴다.
#   p60의 보고서 라벨이 "Hansen/ISSCAN"으로 고정되어 있어, compound_name 이
#   비어 있으면 컴파일 데이터셋이 실측 기록처럼 표시된다.
COMPOUND_NAME_TAG = "TDC-AMES compiled dataset (not a primary experimental record)"
BAR = "=" * 66


def restore() -> int:
    if not BACKUP_PATH.exists():
        print(f"[오류] 백업 파일이 없습니다: {BACKUP_PATH}")
        return 1
    for suffix in ("-wal", "-shm"):
        p = Path(str(DB_PATH) + suffix)
        if p.exists():
            p.unlink()
    shutil.copy2(BACKUP_PATH, DB_PATH)
    print(f"복구 완료: {BACKUP_PATH.name} → {DB_PATH.name}")
    return 0


def db_stats(path: Path) -> tuple[int, int]:
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        n = conn.execute("SELECT COUNT(*) FROM m7_empirical").fetchone()[0]
        u = conn.execute("SELECT COUNT(DISTINCT inchikey) FROM m7_empirical").fetchone()[0]
        return n, u
    finally:
        conn.close()


def main() -> int:
    if "--restore" in sys.argv:
        return restore()

    print(BAR)
    print("  PRIME-DS  Layer A 실측 DB 구축")
    print(BAR)
    print()

    for p in (DB_PATH, TRAIN_CSV, VALID_CSV):
        if not p.exists():
            print(f"[오류] 파일이 없습니다: {p}")
            if p is TRAIN_CSV:
                # 2026-09-03 수정(F-2): 실제 파일명으로 정정
                print("       먼저 python s1_setup_layer_c_ames.py 를 실행하세요.")
            return 1

    import pandas as pd

    train = pd.read_csv(TRAIN_CSV)
    valid = pd.read_csv(VALID_CSV)
    need = {"inchikey", "smiles", "label"}
    if not need.issubset(train.columns):
        print(f"[오류] ames_train.csv 에 {need} 열이 필요합니다.")
        return 1

    print(f"학습셋 : {len(train):,}건  (양성 {int((train.label == 1).sum()):,} / "
          f"음성 {int((train.label == 0).sum()):,})")
    print(f"검증셋 : {len(valid):,}건   ({VALID_CSV.relative_to(BASE_DIR)}"
          f"{' · 동결본' if VALID_CSV == FROZEN_VALID_CSV else ' · 작업본(동결본 없음)'})")
    print()

    # --- 안전장치: 검증셋과 겹치면 중단 ---
    # 2026-09-03 수정(C-2): p60의 실제 조회가 inchikey OR canonical_smiles 이므로
    #                       겹침 확인도 두 키 모두로 한다.
    print("[1/4] 검증셋 겹침 확인 (InChIKey + canonical SMILES · 순환논증 방지)")
    overlap = set(train["inchikey"]) & set(valid["inchikey"])
    overlap_smi = set(train["smiles"]) & set(valid["smiles"])
    if overlap or overlap_smi:
        print(f"      [중단] 검증셋과 InChIKey {len(overlap):,}건 · "
              f"canonical SMILES {len(overlap_smi):,}건이 겹칩니다. 적재하지 않습니다.")
        return 1
    print("      겹침 0건 (InChIKey · canonical SMILES 모두) — 안전합니다.")
    print()

    before_n, before_u = db_stats(DB_PATH)
    print(f"[2/4] 백업 생성 (현재 {before_n:,}건)")
    shutil.copy2(DB_PATH, BACKUP_PATH)
    print(f"      cache/{BACKUP_PATH.name}")
    print()

    conn = sqlite3.connect(DB_PATH)
    try:
        conn.execute("PRAGMA journal_mode=WAL")

        existing_keys = {r[0] for r in conn.execute("SELECT inchikey FROM m7_empirical")}
        curated = conn.execute(
            "SELECT COUNT(*) FROM m7_empirical WHERE source <> ?", (SOURCE_TAG,)
        ).fetchone()[0]

        print(f"[3/4] 기존 TDC 적재분 제거 후 재적재 (기존 큐레이션 {curated:,}건은 보존)")
        removed = conn.execute(
            "DELETE FROM m7_empirical WHERE source = ?", (SOURCE_TAG,)
        ).rowcount
        if removed > 0:
            print(f"      이전 적재분 {removed:,}건 제거")
            existing_keys = {r[0] for r in conn.execute("SELECT inchikey FROM m7_empirical")}

        rows = []
        skipped_dup = 0
        for rec in train.itertuples(index=False):
            key = str(rec.inchikey)
            if key in existing_keys:
                skipped_dup += 1  # 기존 큐레이션 레코드 우선
                continue
            rows.append((
                key,
                str(rec.smiles),
                # 2026-09-03 수정(C-11): compound_name 을 비워 두지 않고 출처를 명시한다.
                COMPOUND_NAME_TAG,
                SOURCE_TAG,
                "positive" if int(rec.label) == 1 else "negative",
                "unknown",  # 발암성 정보는 이 데이터셋에 없음
            ))

        conn.executemany(
            "INSERT INTO m7_empirical "
            "(inchikey, canonical_smiles, compound_name, source, ames_result, carcinogenicity_result) "
            "VALUES (?,?,?,?,?,?)",
            rows,
        )
        conn.commit()
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        conn.commit()
        print(f"      신규 적재 {len(rows):,}건 · 기존과 중복되어 건너뜀 {skipped_dup:,}건")
        print()
    finally:
        conn.close()

    after_n, after_u = db_stats(DB_PATH)
    print("[4/4] 검증")
    print(f"      레코드 {before_n:,} → {after_n:,}건   (고유 구조 {before_u:,} → {after_u:,})")

    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    try:
        print()
        print("      출처별 구성")
        for s, n in conn.execute(
            "SELECT source, COUNT(*) FROM m7_empirical GROUP BY source ORDER BY 2 DESC LIMIT 6"
        ):
            print(f"        {n:6,}   {s}")

        print()
        print("      Ames 결과 구성")
        for s, n in conn.execute(
            "SELECT ames_result, COUNT(*) FROM m7_empirical GROUP BY ames_result ORDER BY 2 DESC"
        ):
            print(f"        {n:6,}   {s}")

        # 검증셋이 실제로 조회되지 않는지 최종 확인
        print()
        print("      검증셋 조회 테스트 (0건이어야 정상)")
        hits = 0
        for smi in valid["smiles"].head(400):
            r = conn.execute(
                "SELECT 1 FROM empirical_toxicology_db WHERE smiles = ? LIMIT 1", (str(smi),)
            ).fetchone()
            if r:
                hits += 1
        print(f"        검증셋 앞 400건 중 DB 적중: {hits}건 "
              f"{'✅ 정상' if hits == 0 else '⚠️ 순환논증 위험'}")

        # 학습셋은 조회돼야 정상
        found = 0
        for smi in train["smiles"].head(200):
            r = conn.execute(
                "SELECT 1 FROM empirical_toxicology_db WHERE smiles = ? LIMIT 1", (str(smi),)
            ).fetchone()
            if r:
                found += 1
        print(f"        학습셋 앞 200건 중 DB 적중: {found}건 "
              f"{'✅ 정상' if found > 180 else '⚠️ 확인 필요'}")
    finally:
        conn.close()

    print()
    print(BAR)
    print("완료. PRIME-DS를 재시작하면 Layer A가 새 DB를 사용합니다.")
    # 2026-09-03 수정(F-2): 실제 파일명으로 정정
    print("되돌리려면:  python s2_build_layer_a_from_ames.py --restore")
    print()
    # 2026-09-03 추가(C-11): 근거 계층 성격을 실행 시점에도 명시한다.
    print("  [주의] 방금 적재한 레코드는 TDC가 컴파일한 Ames 데이터셋이며")
    print("         Hansen/ISSCAN 같은 1차 실험 기록이 아닙니다.")
    print(f"         source        = {SOURCE_TAG}")
    print(f"         compound_name = {COMPOUND_NAME_TAG}")
    print("         p60의 보고서 라벨은 \"Empirical DB Source Match (Hansen/ISSCAN)\"로")
    print("         고정되어 있어 이 구분이 보고서에 반영되지 않습니다.")
    print("         → p60 권고: 해당 고정 라벨에서 (Hansen/ISSCAN)을 떼고")
    print("           source 문자열을 그대로 노출하도록 변경할 것.")
    print(BAR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
