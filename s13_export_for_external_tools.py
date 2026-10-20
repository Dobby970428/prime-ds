"""외부 (Q)SAR 도구(Derek Nexus · Sarah Nexus · Toxtree)에 넣을 입력 파일을 만듭니다.

왜 별도 스크립트인가
--------------------
논문 표 8(head-to-head)은 «모든 도구가 동일한 구조를 본다» 는 것이 성립해야 의미가
있습니다. 각자 CSV 를 열어 각자 방식으로 구조를 만들면 그 전제가 깨집니다. 이 스크립트는
검증셋을 RDKit 으로 한 번 표준화해 SD file 로 내보내고, 같은 구조에서 SMILES 목록도
함께 만들어, 어떤 도구를 쓰든 입력 구조가 동일하도록 고정합니다.

정답(실측 Ames)은 예측 파일에 넣지 않습니다
--------------------------------------------
예측을 돌리는 사람이 정답을 볼 수 있으면 «무의식적 조정» 을 배제할 수 없고, 심사에서
지적받습니다. 그래서 두 파일로 나눕니다.

    ames_validation_blinded.sdf   ← 도구에 넣는 파일. PRIME_ID 와 구조만 있음
    ames_validation_key.csv       ← 채점용. PRIME_ID ↔ 실측 Ames. 예측이 끝난 뒤에만 사용

Derek / Sarah 결과를 받아 채점할 때는 PRIME_ID 로 조인하십시오.

실행
    python s13_export_for_external_tools.py
    python s13_export_for_external_tools.py --source 논문자료/01_데이터셋/ames_validation.csv

출력 (논문자료/02_검증결과/external_tools/)
    ames_validation_blinded.sdf     Derek Nexus · Sarah Nexus · Toxtree 입력
    ames_validation_blinded.smi     SMILES 입력을 받는 도구용 (탭 구분: SMILES<TAB>PRIME_ID)
    ames_validation_key.csv         PRIME_ID, smiles, inchikey, actual_ames
    export_manifest.txt             건수·해시·표준화 실패 목록

2026-09-04 작성.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

DEFAULT_SOURCE = APP_DIR / "논문자료" / "01_데이터셋" / "ames_validation.csv"
OUT_DIR = APP_DIR / "논문자료" / "02_검증결과" / "external_tools"

BAR = "=" * 78


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(DEFAULT_SOURCE),
                    help="검증셋 CSV (smiles, label 컬럼 필요)")
    ap.add_argument("--prefix", default="PRIME", help="식별자 접두사")
    ap.add_argument("--stem", default="ames_validation",
                    help="출력 파일 이름의 앞부분 "
                         "(기본 ames_validation → ames_validation_blinded.sdf 등). "
                         "학습셋을 내보낼 때는 --stem ames_train 으로 분리할 것.")
    args = ap.parse_args()

    import pandas as pd
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")

    src = Path(args.source)
    if not src.is_file():
        print(f"[오류] 검증셋을 찾을 수 없습니다: {src}")
        return 1

    df = pd.read_csv(src)
    if "smiles" not in df.columns:
        print(f"[오류] 'smiles' 컬럼이 없습니다. 있는 컬럼: {list(df.columns)}")
        return 1
    label_col = "label" if "label" in df.columns else None

    print(BAR)
    print("  외부 (Q)SAR 도구용 입력 파일 생성")
    print(BAR)
    print(f"원본     : {src}")
    print(f"  SHA-256: {sha256(src)}")
    print(f"  행 수   : {len(df):,}")
    print()

    # 표준화 — 모든 도구가 같은 구조를 보도록 한 번만 수행합니다.
    chooser = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sdf_path = OUT_DIR / f"{args.stem}_blinded.sdf"
    smi_path = OUT_DIR / f"{args.stem}_blinded.smi"
    key_path = OUT_DIR / f"{args.stem}_key.csv"

    # RDKit SDWriter 는 경로 문자열을 C 계층으로 넘기기 때문에 한글이 섞인
    # 경로에서 실패한다. 파이썬이 먼저 연 파일 객체를 넘기면 우회된다.
    sdf_handle = open(sdf_path, "w", encoding="utf-8", newline="")
    writer = Chem.SDWriter(sdf_handle)
    smi_lines: list[str] = []
    key_rows: list[dict] = []
    failed: list[tuple[int, str]] = []

    for i, row in df.iterrows():
        raw = str(row["smiles"])
        mol = Chem.MolFromSmiles(raw)
        if mol is None:
            failed.append((i, raw))
            continue
        try:
            mol = chooser.choose(mol)
            mol = uncharger.uncharge(mol)
            Chem.SanitizeMol(mol)
        except Exception:  # noqa: BLE001
            failed.append((i, raw))
            continue

        pid = f"{args.prefix}-{i + 1:05d}"
        can = Chem.MolToSmiles(mol, canonical=True)
        try:
            ik = Chem.MolToInchiKey(mol)
        except Exception:  # noqa: BLE001
            ik = ""

        # 예측 파일에는 식별자와 구조만. 정답은 넣지 않습니다.
        mol.SetProp("_Name", pid)
        mol.SetProp("PRIME_ID", pid)
        mol.SetProp("Canonical_SMILES", can)
        mol.SetProp("InChIKey", ik)
        mol.SetProp("Molecular_Formula", rdMolDescriptors.CalcMolFormula(mol))
        writer.write(mol)

        smi_lines.append(f"{can}\t{pid}")
        key_rows.append({
            "PRIME_ID": pid,
            "smiles": can,
            "inchikey": ik,
            "actual_ames": (int(row[label_col]) if label_col is not None else ""),
        })

    writer.close()
    sdf_handle.close()
    smi_path.write_text("\n".join(smi_lines) + "\n", encoding="utf-8")
    pd.DataFrame(key_rows).to_csv(key_path, index=False, encoding="utf-8-sig")

    n = len(key_rows)
    pos = sum(1 for r in key_rows if r["actual_ames"] == 1) if label_col else 0
    print(f"내보낸 구조 : {n:,}건 (표준화 실패 {len(failed)}건)")
    if label_col:
        print(f"  실측 양성 : {pos:,} / 음성 : {n - pos:,}   ← key 파일에만 있습니다")
    print()
    print(f"  {sdf_path.relative_to(APP_DIR)}   ← Derek Nexus · Sarah Nexus · Toxtree 입력")
    print(f"  {smi_path.relative_to(APP_DIR)}   ← SMILES 입력을 받는 도구용")
    print(f"  {key_path.relative_to(APP_DIR)}   ← 채점용. 예측이 끝난 뒤에만 여십시오")

    manifest = [
        BAR,
        "  외부 도구 입력 파일 manifest",
        BAR,
        f"원본            : {src}",
        f"원본 SHA-256    : {sha256(src)}",
        f"원본 행 수      : {len(df):,}",
        f"내보낸 구조     : {n:,}",
        f"표준화 실패     : {len(failed)}",
        "",
        "표준화 절차 (모든 도구가 동일 구조를 보도록 1회만 수행)",
        "  1. RDKit MolFromSmiles",
        "  2. LargestFragmentChooser  — 염·용매 제거",
        "  3. Uncharger               — 중성화",
        "  4. SanitizeMol",
        "  5. canonical SMILES 재생성",
        "",
        f"blinded SDF SHA-256 : {sha256(sdf_path)}",
        f"blinded SMI SHA-256 : {sha256(smi_path)}",
        f"key CSV SHA-256     : {sha256(key_path)}",
        "",
        "SD file 필드: PRIME_ID, Canonical_SMILES, InChIKey, Molecular_Formula",
        "실측 Ames 결과는 SD file 에 포함되어 있지 않습니다 (의도된 설계).",
    ]
    if failed:
        manifest += ["", "표준화 실패 목록 (행 번호, 원본 SMILES)"]
        manifest += [f"  {i}\t{s}" for i, s in failed[:50]]
        if len(failed) > 50:
            manifest.append(f"  ... 외 {len(failed) - 50}건")
    manifest_path = OUT_DIR / f"{args.stem}_export_manifest.txt"
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(f"  {manifest_path.relative_to(APP_DIR)}")
    print(BAR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
