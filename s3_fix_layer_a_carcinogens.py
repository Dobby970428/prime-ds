"""
Layer A 실측 DB의 발암성 양성 레코드 구조 오류 수정.

배경
----
cache/m7_empirical.db 에서 carcinogenicity_result='positive' 인 레코드는
ICH M7 Class 1(최고 위험 등급) 판정에 직접 사용되며, 이 등급은 AI 예측으로
변경되지 않는다. 그런데 검토 결과 20건 중 10건에서 화합물명과 구조가
일치하지 않는 것이 확인되었다.

예) Benz[a]pyrene 로 등록된 구조의 분자식이 C14H10 (정답 C20H12)
    N-Nitroso-N-methylurethane 에 니트로소기가 없음

또한 5건은 원본 CSV의 쉼표 파싱 오류로 이름과 출처가 뒤바뀌어 있다.
예) 이름='1', 출처='2-Dichloroethane'  ←  원래 '1,2-Dichloroethane'

이 스크립트는 구조를 검증된 값으로 교체하고, 중복 레코드를 제거한다.

사용법
------
    python s3_fix_layer_a_carcinogens.py            # 수정
    python s3_fix_layer_a_carcinogens.py --check    # 확인만 (수정하지 않음)
    python s3_fix_layer_a_carcinogens.py --restore  # 되돌리기
"""

from __future__ import annotations

import shutil
import sqlite3
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "cache" / "m7_empirical.db"
BACKUP = BASE_DIR / "cache" / "m7_empirical.db.before_carc_fix"

BAR = "=" * 74

# (현재 DB의 구조 또는 이름) → (올바른 이름, 올바른 SMILES, 예상 분자식)
# 모든 SMILES는 RDKit으로 분자식을 검증하였다.
CORRECTIONS: list[tuple[str, str, str, str]] = [
    # match_key(기존 canonical_smiles), 올바른 이름, 올바른 SMILES, 예상 분자식
    ("O=C1c2ccccc2C(=O)c2c(O)cccc21", "2-Aminoanthraquinone",
     "Nc1ccc2c(c1)C(=O)c1ccccc1C2=O", "C14H9NO2"),
    ("Nc1ccc(N)c(N)c1", "4,4'-Diaminobiphenyl (Benzidine)",
     "Nc1ccc(-c2ccc(N)cc2)cc1", "C12H12N2"),
    ("c1ccc(N=Nc2ccccc2)cc1", "4-Aminobiphenyl",
     "Nc1ccc(-c2ccccc2)cc1", "C12H11N"),
    ("N=Nc1ccc(-c2ccccc2)cc1", "Azobenzene",
     "c1ccc(N=Nc2ccccc2)cc1", "C12H10N2"),
    ("O=C1C=Cc2cc3ccccc3cc21", "Benz[a]anthracene",
     "c1ccc2cc3c(ccc4ccccc43)cc2c1", "C18H12"),
    ("c1ccc2cc3ccccc3cc2c1", "Benzo[a]pyrene",
     "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34", "C20H12"),
    ("CNC(=O)N=O", "N-Nitroso-N-methylurea (MNU)",
     "CN(N=O)C(N)=O", "C2H5N3O2"),
    ("CNC(=O)OC", "N-Nitroso-N-methylurethane",
     "CCOC(=O)N(C)N=O", "C4H8N2O3"),
]

# 이름만 깨진 것 (구조는 정상) — 쉼표 파싱 오류 복구
NAME_FIXES: list[tuple[str, str]] = [
    ("Nc1ccccc1N", "1,2-Phenylenediamine"),
    ("Nc1ccc(N)cc1", "1,4-Phenylenediamine"),
    ("ClCCCl", "1,2-Dichloroethane"),
    ("BrCCBr", "1,2-Dibromoethane"),
    ("O=[N+]([O-])c1ccc([N+](=O)[O-])c([N+](=O)[O-])c1", "1,2,4-Trinitrobenzene"),
]

# 중복/불완전 레코드 — 삭제 대상
DELETE_SMILES: list[tuple[str, str]] = [
    ("N=Nc1ccc(N)cc1", "4-Aminobiphenyl (alt) — 4-Aminobiphenyl 와 중복"),
    ("N=Nc1ccccc1", "Azobenzene (core) — 구조 불완전, Azobenzene 와 중복"),
]

LEGACY_SOURCE = "Legacy curated (출처 미확인)"


def restore() -> int:
    if not BACKUP.exists():
        print(f"[오류] 백업이 없습니다: {BACKUP}")
        return 1
    for sfx in ("-wal", "-shm"):
        p = Path(str(DB_PATH) + sfx)
        if p.exists():
            p.unlink()
    shutil.copy2(BACKUP, DB_PATH)
    print(f"복구 완료: {BACKUP.name} → {DB_PATH.name}")
    return 0


def main() -> int:
    if "--restore" in sys.argv:
        return restore()
    check_only = "--check" in sys.argv

    print(BAR)
    print("  Layer A 실측 DB — 발암성 양성 레코드 구조 검증" + (" (확인 전용)" if check_only else " 및 수정"))
    print(BAR)
    print()

    if not DB_PATH.exists():
        print(f"[오류] DB가 없습니다: {DB_PATH}")
        return 1

    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors as D

    RDLogger.DisableLog("rdApp.*")

    def canon(smi: str) -> str | None:
        m = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(m) if m else None

    def formula(smi: str) -> str:
        m = Chem.MolFromSmiles(smi)
        return D.CalcMolFormula(m) if m else "?"

    # --- 1단계: 교체할 구조를 먼저 검증 ---
    print("[1/4] 교체할 구조의 분자식 검증")
    bad = 0
    for _, name, new_smi, want in CORRECTIONS:
        got = formula(new_smi)
        ok = got == want
        if not ok:
            bad += 1
        print(f"      {'O' if ok else 'X'} {name[:38]:<40} {got:<12} (예상 {want})")
    if bad:
        print(f"\n[중단] 교체 구조 {bad}건이 예상 분자식과 다릅니다. 스크립트를 확인하세요.")
        return 1
    print("      모두 일치")
    print()

    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    rows = list(conn.execute(
        "SELECT id, compound_name, source, canonical_smiles, ames_result "
        "FROM m7_empirical WHERE carcinogenicity_result='positive'"))
    # 이름 손상 복구는 발암성 여부와 무관하게 적용해야 하므로 별도로 조회한다.
    # (예: 트리니트로벤젠은 carcinogenicity_result='unknown' 이지만 이름이 깨져 있다)
    broken = list(conn.execute(
        "SELECT id, compound_name, source, canonical_smiles, ames_result "
        "FROM m7_empirical WHERE compound_name IN ('1','2','3') OR compound_name IS NULL"))
    conn.close()
    print(f"[2/4] 발암성 양성 레코드: {len(rows)}건 · 이름 손상 레코드: {len(broken)}건")
    print()

    by_canon: dict[str, list] = {}
    for r in rows:
        c = canon(r[3]) or r[3]
        by_canon.setdefault(c, []).append(r)

    broken_canon: dict[str, list] = {}
    for r in broken:
        c = canon(r[3]) or r[3]
        broken_canon.setdefault(c, []).append(r)

    plan_struct, plan_name, plan_del, missing = [], [], [], []

    for match_smi, name, new_smi, _ in CORRECTIONS:
        c = canon(match_smi) or match_smi
        if c in by_canon:
            for r in by_canon[c]:
                plan_struct.append((r[0], r[1], name, new_smi))
        else:
            missing.append(f"구조 교체 대상 없음: {name}")

    for match_smi, name in NAME_FIXES:
        c = canon(match_smi) or match_smi
        targets = {r[0]: r for r in by_canon.get(c, [])}
        targets.update({r[0]: r for r in broken_canon.get(c, [])})
        if targets:
            for r in targets.values():
                if r[1] != name:
                    plan_name.append((r[0], r[1], name))
        else:
            missing.append(f"이름 복구 대상 없음: {name}")

    for match_smi, why in DELETE_SMILES:
        c = canon(match_smi) or match_smi
        if c in by_canon:
            for r in by_canon[c]:
                plan_del.append((r[0], r[1], why))

    print("[3/4] 수행할 작업")
    print(f"      구조 교체 {len(plan_struct)}건 · 이름 복구 {len(plan_name)}건 · 삭제 {len(plan_del)}건")
    print()
    if plan_struct:
        print("      [구조 교체]")
        for _id, old, new, smi in plan_struct:
            print(f"        · {str(old)[:34]:<36} → {new[:34]}")
            print(f"          {smi}")
    if plan_name:
        print("      [이름 복구]")
        for _id, old, new in plan_name:
            print(f"        · {str(old)!r} → {new}")
    if plan_del:
        print("      [삭제]")
        for _id, old, why in plan_del:
            print(f"        · {str(old)[:34]:<36} {why}")
    if missing:
        print("      [경고]")
        for m in missing:
            print(f"        · {m}")
    print()

    if check_only:
        # 2026-09-03 수정(F-2): 실제 파일명 기준으로 안내
        print("확인 전용 모드입니다. 수정하려면 --check 없이 다시 실행하세요.")
        print("    python s3_fix_layer_a_carcinogens.py")
        return 0

    if not (plan_struct or plan_name or plan_del):
        print("수정할 것이 없습니다.")
        return 0

    shutil.copy2(DB_PATH, BACKUP)
    print(f"[4/4] 백업 생성: cache/{BACKUP.name}")

    conn = sqlite3.connect(DB_PATH)
    try:
        conn.execute("PRAGMA journal_mode=WAL")
        for _id, _old, new_name, new_smi in plan_struct:
            m = Chem.MolFromSmiles(new_smi)
            conn.execute(
                "UPDATE m7_empirical SET compound_name=?, canonical_smiles=?, inchikey=?, source=? "
                "WHERE id=?",
                (new_name, Chem.MolToSmiles(m), Chem.MolToInchiKey(m), LEGACY_SOURCE, _id))
        for _id, _old, new_name in plan_name:
            conn.execute("UPDATE m7_empirical SET compound_name=?, source=? WHERE id=?",
                         (new_name, LEGACY_SOURCE, _id))
        for _id, _old, _why in plan_del:
            conn.execute("DELETE FROM m7_empirical WHERE id=?", (_id,))
        conn.commit()
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        conn.commit()
    finally:
        conn.close()

    # --- 검증 ---
    print()
    print(BAR)
    print("  수정 후 발암성 양성 레코드")
    print(BAR)
    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    final = list(conn.execute(
        "SELECT compound_name, canonical_smiles, ames_result, source FROM m7_empirical "
        "WHERE carcinogenicity_result='positive' ORDER BY compound_name"))
    total = conn.execute("SELECT COUNT(*) FROM m7_empirical").fetchone()[0]
    conn.close()
    print(f"{'이름':<38} {'분자식':<12} {'Ames':<9} 구조")
    print("-" * 100)
    for nm, smi, ames, _src in final:
        print(f"{str(nm)[:36]:<38} {formula(smi):<12} {ames:<9} {smi[:34]}")
    print()
    print(f"발암성 양성 {len(final)}건 · 전체 레코드 {total:,}건")
    print()
    # 2026-09-03 수정(F-2): 실제 파일명으로 정정
    print("되돌리려면:  python s3_fix_layer_a_carcinogens.py --restore")
    print()
    # 2026-09-03 추가(C-3): DB를 바꿨으므로 검증셋 무결성 재확인이 필요하다.
    print("DB 내용이 바뀌었습니다. 검증셋 무결성을 다시 확인하십시오:")
    print("    python s4_recheck_validation_set.py --check")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
