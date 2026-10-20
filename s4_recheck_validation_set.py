"""
외부검증셋과 Layer A 실측 DB의 중복을 다시 점검하고 «제외 사유»를 표기한다.

왜 필요한가
-----------
검증셋은 s1_setup_layer_c_ames.py 실행 시점의 DB 내용을 기준으로 걸러졌다.
이후 DB가 변경되면(레코드 추가, 구조 수정 등) 새로운 중복이 생길 수 있다.

실제 사례: s3_fix_layer_a_carcinogens.py 로 2-Aminoanthraquinone 의 구조를
잘못된 값에서 정답으로 교체한 결과, 그 정답 구조가 검증셋에 존재하여
Layer A 조회가 1건 발생하였다. 예측이 아니라 조회이므로 제외해야 한다.

2026-09-03 수정(C-3) — 검증셋을 잠근다
--------------------------------------
이전 버전은 검증셋 CSV에서 행을 **삭제하고 원본을 제자리에서 덮어썼다**.
«DB 수정 → 지표 확인 → 검증셋 축소 → 재측정» 이 반복되면 심사자 관점에서
검증셋이 locked 상태가 아니다. 그래서 다음과 같이 바꿨다.

  · 동결본 data/frozen/ames_validation_v1.csv 를 **읽기 전용**으로 읽는다.
  · 행을 삭제하지 않는다. 대신 작업본에 `excluded_reason` 컬럼을 채운다.
  · s5 가 그 컬럼으로 필터링한다. 그러면 «몇 건을 왜 뺐는지»가 그대로
    Methods 에 들어간다.
  · s1 이 동결본을 만들어 두었다면, 이 스크립트는 사후 수습이 아니라
    무결성 점검(assert)으로 남는다 — 새 중복이 나오면 그것 자체가
    보고 대상이다.

사용법
------
    python s4_recheck_validation_set.py --check    # 확인만
    python s4_recheck_validation_set.py            # excluded_reason 표기
    python s4_recheck_validation_set.py --restore  # 되돌리기

실행 후에는 반드시 재검증하십시오.
    python s5_run_pipeline_validation.py
"""

from __future__ import annotations

import hashlib
import shutil
import sqlite3
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
FROZEN_DIR = DATA_DIR / "frozen"
DB_PATH = BASE_DIR / "cache" / "m7_empirical.db"
# 2026-09-03 수정(C-3a): 동결본을 원본(읽기 전용)으로, 작업본에만 표기한다.
FROZEN_VALID_CSV = FROZEN_DIR / "ames_validation_v1.csv"
FROZEN_VALID_SHA = FROZEN_DIR / "ames_validation_v1.sha256"
VALID_CSV = DATA_DIR / "ames_validation.csv"
BACKUP_CSV = DATA_DIR / "ames_validation.csv.before_recheck"
LOG_TXT = DATA_DIR / "validation_recheck_log.txt"

EXCLUDE_TAG = "layer_a_db_match"

BAR = "=" * 74
_log: list[str] = []


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def say(m: str = "") -> None:
    print(m)
    _log.append(m)


def restore() -> int:
    if not BACKUP_CSV.exists():
        print(f"[오류] 백업이 없습니다: {BACKUP_CSV}")
        return 1
    shutil.copy2(BACKUP_CSV, VALID_CSV)
    print(f"복구 완료: {BACKUP_CSV.name} → {VALID_CSV.name}")
    return 0


def main() -> int:
    if "--restore" in sys.argv:
        return restore()
    check_only = "--check" in sys.argv

    say(BAR)
    say("  외부검증셋 · 실측 DB 중복 재점검" + (" (확인 전용)" if check_only else ""))
    say(BAR)
    say()

    if not DB_PATH.exists():
        say(f"[오류] 파일이 없습니다: {DB_PATH}")
        return 1

    # 2026-09-03 수정(C-3a): 동결본이 있으면 그것을 정본으로 삼는다.
    if FROZEN_VALID_CSV.exists():
        source_csv = FROZEN_VALID_CSV
        digest = sha256_of(FROZEN_VALID_CSV)
        say(f"검증셋 정본 : data/frozen/{FROZEN_VALID_CSV.name}  (동결본 · 읽기 전용)")
        say(f"SHA-256     : {digest}")
        if FROZEN_VALID_SHA.exists():
            recorded = FROZEN_VALID_SHA.read_text(encoding="utf-8").split()
            if recorded and recorded[0] != digest:
                say("[중단] 동결본 해시가 기록과 다릅니다. 검증셋이 변경되었습니다.")
                say(f"        기록 {recorded[0]}")
                say(f"        실제 {digest}")
                return 1
            say("            해시 일치 — 검증셋이 동결 이후 변경되지 않았습니다.")
    elif VALID_CSV.exists():
        source_csv = VALID_CSV
        say(f"검증셋 정본 : data/{VALID_CSV.name}  (동결본 없음 — 작업본 사용)")
        say("[경고] data/frozen/ames_validation_v1.csv 가 없습니다.")
        say("       s1_setup_layer_c_ames.py 를 다시 실행해 검증셋을 동결하십시오.")
    else:
        say(f"[오류] 검증셋을 찾을 수 없습니다: {FROZEN_VALID_CSV} / {VALID_CSV}")
        return 1

    import pandas as pd
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")

    def canon(s: str) -> str | None:
        m = Chem.MolFromSmiles(str(s)) if s else None
        return Chem.MolToSmiles(m) if m else None

    def ikey(s: str) -> str | None:
        m = Chem.MolFromSmiles(str(s)) if s else None
        return Chem.MolToInchiKey(m) if m else None

    v = pd.read_csv(source_csv)
    # 2026-09-03 수정(C-3b): 행을 지우지 않고 사유 컬럼으로 관리한다.
    if "excluded_reason" not in v.columns:
        v["excluded_reason"] = ""
    v["excluded_reason"] = v["excluded_reason"].fillna("").astype(str)
    say(f"현재 검증셋 : {len(v):,}건 "
        f"(양성 {int((v['label'] == 1).sum()):,} · 음성 {int((v['label'] == 0).sum()):,})")
    already = int((v["excluded_reason"] != "").sum())
    if already:
        say(f"              이미 제외 표기된 행 {already:,}건")

    # --- DB 색인 ---
    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    db_rows = list(conn.execute(
        "SELECT compound_name, canonical_smiles, inchikey, ames_result, "
        "carcinogenicity_result, source FROM m7_empirical"))
    conn.close()

    db_by_canon: dict[str, tuple] = {}
    db_by_ikey: dict[str, tuple] = {}
    for nm, smi, ik, ames, carc, src in db_rows:
        c = canon(smi)
        if c:
            db_by_canon.setdefault(c, (nm, ames, carc, src))
        k = ik or ikey(smi)
        if k:
            db_by_ikey.setdefault(k, (nm, ames, carc, src))

    say(f"실측 DB     : {len(db_rows):,}건 "
        f"(고유 구조 {len(db_by_canon):,})")
    say()

    # --- 대조 ---
    say("[1/3] canonical SMILES · InChIKey 양쪽으로 대조")
    hits: list[dict] = []
    for i, r in enumerate(v.itertuples()):
        smi = str(getattr(r, "smiles", ""))
        c = canon(smi)
        k = getattr(r, "inchikey", None) or ikey(smi)
        found = None
        how = ""
        if c and c in db_by_canon:
            found, how = db_by_canon[c], "SMILES"
        elif k and k in db_by_ikey:
            found, how = db_by_ikey[k], "InChIKey"
        if found:
            hits.append({"row": i, "smiles": c or smi, "label": int(getattr(r, "label")),
                         "db_name": found[0], "db_ames": found[1], "db_carc": found[2],
                         "db_source": found[3], "matched_by": how})

    say(f"      중복 {len(hits)}건")
    say()

    if hits:
        say("[2/3] 중복 상세  ← Methods에 그대로 기재할 «무엇을 왜 뺐는가» 목록")
        for h in hits:
            say(f"      · {h['db_name'] or '(이름 없음)'}")
            say(f"        구조      {h['smiles']}")
            say(f"        검증셋    정답 {'양성' if h['label'] == 1 else '음성'}")
            say(f"        DB 기록   Ames={h['db_ames']} · 발암성={h['db_carc']}")
            say(f"        출처      {h['db_source']}   (일치: {h['matched_by']})")
        say()
    else:
        say("[2/3] 중복 없음 — 검증셋이 실측 DB와 완전히 분리되어 있습니다.")
        say()

    if check_only:
        say("확인 전용 모드입니다. 표기하려면 --check 없이 실행하세요.")
        say("    python s4_recheck_validation_set.py")
        LOG_TXT.write_text("\n".join(_log), encoding="utf-8")
        return 0

    if not hits:
        # 2026-09-03 수정(C-3): 사후 수습이 아니라 무결성 점검으로 남는다.
        say("[3/3] 변경할 것이 없습니다 — 검증셋 무결성 확인 통과.")
        if source_csv == FROZEN_VALID_CSV or not VALID_CSV.exists():
            # 작업본이 없거나 낡았을 수 있으므로 동결본에서 다시 만들어 둔다.
            v.to_csv(VALID_CSV, index=False, encoding="utf-8-sig")
            say(f"      작업본 갱신: data/{VALID_CSV.name} "
                f"({len(v):,}건 · 제외 표기 {int((v['excluded_reason'] != '').sum()):,}건)")
        LOG_TXT.write_text("\n".join(_log), encoding="utf-8")
        return 0

    # 2026-09-03 수정(C-3b): 행 삭제 → excluded_reason 표기
    say("[3/3] 검증셋에 제외 사유 표기 (행은 삭제하지 않습니다)")
    if VALID_CSV.exists():
        shutil.copy2(VALID_CSV, BACKUP_CSV)
        say(f"      작업본 백업: data/{BACKUP_CSV.name}")

    for h in hits:
        prev = str(v.at[v.index[h["row"]], "excluded_reason"] or "")
        reason = (f"{EXCLUDE_TAG}:{h['matched_by']}"
                  f":{(h['db_name'] or 'unnamed').replace(',', ' ')}"
                  f":{(h['db_source'] or 'unknown').replace(',', ' ')}")
        v.at[v.index[h["row"]], "excluded_reason"] = (
            reason if not prev else f"{prev}|{reason}")

    v.to_csv(VALID_CSV, index=False, encoding="utf-8-sig")
    n_excl = int((v["excluded_reason"] != "").sum())
    kept = v[v["excluded_reason"] == ""]
    say(f"      작업본 저장: data/{VALID_CSV.name}")
    say(f"      전체 {len(v):,}건 유지 · 제외 표기 {n_excl:,}건 · 분석 대상 {len(kept):,}건")
    say(f"      분석 대상 구성: 양성 {int((kept['label'] == 1).sum()):,} · "
        f"음성 {int((kept['label'] == 0).sum()):,}")
    say(f"      (동결본 data/frozen/{FROZEN_VALID_CSV.name} 는 변경하지 않았습니다)")
    say()
    say(BAR)
    say("  완료 — 이제 재검증을 실행하십시오")
    say(BAR)
    # 2026-09-03 수정(F-2): 실제 파일명으로 정정
    say("    python s5_run_pipeline_validation.py")
    say()
    say("  되돌리려면:  python s4_recheck_validation_set.py --restore")

    LOG_TXT.write_text("\n".join(_log), encoding="utf-8")
    print(f"\n기록 저장: data/{LOG_TXT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
