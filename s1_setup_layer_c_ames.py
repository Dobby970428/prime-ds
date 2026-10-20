"""
Layer C를 '진짜 Ames 데이터'로 재학습하고 외부검증셋을 만드는 스크립트.

하는 일
--------
 1. 내려받은 Ames 데이터 파일을 찾아서 읽는다 (csv / tab / tsv / xlsx)
 2. RDKit `rdMolStandardize` 로 구조를 표준화한다
      · LargestFragmentChooser  — 염/부성분을 떼고 최대 조각만 남김
      · Uncharger               — 형식전하를 중화
      · TautomerEnumerator.Canonicalize — 호변이성체를 대표형으로 통일
    표준화한 구조의 **InChIKey 앞 14자(연결성 블록)** 로 중복을 제거한다.
    ※ 전체 InChIKey는 입체층까지 포함하므로 거울상이성질체를 서로 다른
      화합물로 통과시킨다. 그런데 Layer C의 지문
      `GetMorganGenerator(radius=2, fpSize=2048)` 은 useChirality 기본값이
      False라 L-/D- 쌍의 ECFP가 비트 단위로 동일하다. 즉 전체 InChIKey
      기준으로 중복제거하면 학습·검증에 «같은 특징 벡터»가 양쪽에 남는다.
 3. 학습 80% / 검증 20%로 분할한다
      · scaffold (기본) — Murcko 스캐폴드 분할. 같은 골격은 한쪽에만 간다.
      · random          — 기존 층화 랜덤 분할 (재현·비교용, --split random)
    분할 직후 학습셋과 검증셋 ECFP 행렬에 완전 일치 행이 없는지 검사하고,
    남아 있으면 건수·예시를 출력한 뒤 검증셋에서 제외한다(0건 보장).
 4. 검증셋에서 cache/m7_empirical.db 에 이미 있는 화합물을 제외한다
    → Layer A가 정답을 조회해 버리는 순환논증 방지
    ※ 실제 런타임 조회(p60)는 `WHERE inchikey = ? OR canonical_smiles = ?`
      이므로 **두 키 모두** 로 제외한다.
 5. 학습셋으로 RandomForest를 학습해 models/ames_ecfp_rf.pkl 로 저장한다
 6. 검증셋으로 Layer C 단독 성능을 계산해 보여준다
    (ROC-AUC / PR-AUC / 임계값 스윕 / Wilson 95% CI 포함)
 7. 검증셋을 SHA-256과 함께 data/frozen/ 에 동결한다 (이후 읽기 전용)

사용법
------
    python s1_setup_layer_c_ames.py                      # data/ 폴더에서 자동 탐색
    python s1_setup_layer_c_ames.py 파일경로                # 파일을 직접 지정
    python s1_setup_layer_c_ames.py --split random       # 기존 랜덤 분할로 재현
    python s1_setup_layer_c_ames.py --split scaffold     # (기본) Murcko 스캐폴드 분할

결과물
------
    models/ames_ecfp_rf.pkl              학습된 Layer C 모델 (자동으로 사용됨)
    data/ames_train.csv                  학습에 쓴 화합물
    data/ames_validation.csv             검증에 쓴 화합물 (작업본)
    data/frozen/ames_validation_v1.csv   동결된 external validation set (읽기 전용)
    data/frozen/ames_validation_v1.sha256
    data/ames_setup_report.txt           위 과정의 요약 기록
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
FROZEN_DIR = DATA_DIR / "frozen"          # 2026-09-03 수정(C-3): 검증셋 동결 폴더
MODEL_DIR = BASE_DIR / "models"
EMPIRICAL_DB = BASE_DIR / "cache" / "m7_empirical.db"

TRAIN_CSV = DATA_DIR / "ames_train.csv"
VALID_CSV = DATA_DIR / "ames_validation.csv"
FROZEN_VALID_CSV = FROZEN_DIR / "ames_validation_v1.csv"
FROZEN_VALID_SHA = FROZEN_DIR / "ames_validation_v1.sha256"
MODEL_PKL = MODEL_DIR / "ames_ecfp_rf.pkl"
REPORT_TXT = DATA_DIR / "ames_setup_report.txt"

TEST_FRACTION = 0.20
SPLIT_SEED = 42
RF_TREES = 100
RF_SEED = 0
DEFAULT_SPLIT = "scaffold"                # 2026-09-03 수정(C-1): 기본 분할 방식

BAR = "=" * 66
_log_lines: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    _log_lines.append(msg)


# ---------------------------------------------------------------- 파일 찾기
SMILES_HINTS = ("smiles", "canonical_smiles", "drug", "structure", "mol")
LABEL_HINTS = ("y", "label", "activity", "ames", "mutagen", "class", "outcome", "target")

POSITIVE_WORDS = {"1", "1.0", "pos", "positive", "mutagen", "mutagenic", "active",
                  "yes", "true", "p", "+"}
NEGATIVE_WORDS = {"0", "0.0", "neg", "negative", "nonmutagen", "non-mutagen",
                  "nonmutagenic", "inactive", "no", "false", "n", "-"}


def find_input_file(explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit)
        if not p.is_absolute():
            p = (BASE_DIR / p).resolve()
        if not p.exists():
            raise SystemExit(f"[오류] 파일을 찾을 수 없습니다: {p}")
        return p

    patterns = ("*ames*.csv", "*ames*.tab", "*ames*.tsv", "*ames*.txt", "*ames*.xlsx",
                "*AMES*.csv", "*AMES*.tab", "*mutagen*.csv")
    seen: list[Path] = []
    for folder in (DATA_DIR, BASE_DIR):
        for pat in patterns:
            for hit in folder.glob(pat):
                if hit.name in {TRAIN_CSV.name, VALID_CSV.name} or hit in seen:
                    continue
                seen.append(hit)

    if not seen:
        raise SystemExit(
            "[오류] Ames 데이터 파일을 찾지 못했습니다.\n"
            f"       내려받은 파일을 다음 폴더에 넣어주세요:\n       {DATA_DIR}\n"
            "       (파일명에 'ames'가 들어가야 자동으로 찾습니다)\n"
            # 2026-09-03 수정(F-2): 실제 파일명으로 정정
            "       또는:  python s1_setup_layer_c_ames.py 파일경로"
        )
    seen.sort(key=lambda p: p.stat().st_size, reverse=True)
    return seen[0]


def read_table(path: Path):
    import pandas as pd

    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xlsm"):
        return pd.read_excel(path)
    for sep in ("\t", ",", ";"):
        try:
            df = pd.read_csv(path, sep=sep, engine="python")
            if df.shape[1] >= 2:
                return df
        except Exception:  # noqa: BLE001, S110
            continue
    raise SystemExit(f"[오류] 파일을 읽을 수 없습니다: {path.name}")


def pick_column(columns: list[str], hints: tuple[str, ...]) -> str | None:
    lowered = {c: str(c).strip().lower() for c in columns}
    for hint in hints:  # 완전 일치 우선
        for col, low in lowered.items():
            if low == hint:
                return col
    for hint in hints:  # 부분 일치
        for col, low in lowered.items():
            if hint in low:
                return col
    return None


def normalize_label(value) -> int | None:
    import pandas as pd

    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip().lower()
    if text in POSITIVE_WORDS:
        return 1
    if text in NEGATIVE_WORDS:
        return 0
    try:
        num = float(text)
    except ValueError:
        return None
    if num == 1:
        return 1
    if num == 0:
        return 0
    return None


# ------------------------------------------------- 2026-09-03 추가(C-1): 구조 표준화
_STD_CACHE: dict[str, object] = {}


def _standardizers():
    """rdMolStandardize 도구를 한 번만 만들어 재사용한다."""
    if "ready" not in _STD_CACHE:
        from rdkit.Chem.MolStandardize import rdMolStandardize

        _STD_CACHE["lfc"] = rdMolStandardize.LargestFragmentChooser()
        _STD_CACHE["unc"] = rdMolStandardize.Uncharger()
        _STD_CACHE["taut"] = rdMolStandardize.TautomerEnumerator()
        _STD_CACHE["ready"] = True
    return _STD_CACHE["lfc"], _STD_CACHE["unc"], _STD_CACHE["taut"]


def standardize_mol(mol):
    """염 제거 → 전하 중화 → 호변이성체 정규화. 실패하면 None."""
    from rdkit import Chem

    lfc, unc, taut = _standardizers()
    try:
        m = lfc.choose(mol)
        m = unc.uncharge(m)
        m = taut.Canonicalize(m)
        Chem.SanitizeMol(m)
        return m
    except Exception:  # noqa: BLE001
        return None


def skeleton_key(inchikey: str) -> str:
    """InChIKey 앞 14자 = 연결성(골격) 블록. 입체층·전하층을 뺀 키."""
    return str(inchikey or "")[:14]


def murcko_scaffold(smiles: str, idx: int) -> str:
    """Murcko 스캐폴드 SMILES. 고리가 없는 분자는 각각 독립 그룹으로 둔다.

    무고리 분자를 하나의 '' 그룹으로 묶으면 (Ames에서 중요한 니트로사민 등)
    사슬형 화합물이 통째로 한쪽에만 몰려 검증셋 구성이 크게 왜곡된다.
    """
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"__invalid__{idx}"
    try:
        sc = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:  # noqa: BLE001
        sc = ""
    return sc if sc else f"__acyclic__{idx}"


def scaffold_split(frame, test_fraction: float):
    """Murcko 스캐폴드 분할 — 같은 골격은 학습/검증 한쪽에만 배정한다.

    큰 스캐폴드 그룹부터 학습셋을 채우고(DeepChem ScaffoldSplitter와 동일한
    순서), 학습 목표치를 넘기는 그룹부터 검증셋으로 보낸다. 결정적이며
    난수 seed에 의존하지 않는다.
    """
    scaffolds = [murcko_scaffold(s, i) for i, s in enumerate(frame["smiles"])]
    groups: dict[str, list[int]] = {}
    for i, sc in enumerate(scaffolds):
        groups.setdefault(sc, []).append(i)

    ordered = sorted(groups.values(), key=lambda g: (-len(g), g[0]))
    train_target = len(frame) - int(round(len(frame) * test_fraction))
    tr_idx: list[int] = []
    te_idx: list[int] = []
    for g in ordered:
        if len(tr_idx) + len(g) <= train_target:
            tr_idx.extend(g)
        else:
            te_idx.extend(g)

    train_df = frame.iloc[sorted(tr_idx)].copy()
    test_df = frame.iloc[sorted(te_idx)].copy()
    # 스캐폴드 라벨을 함께 돌려준다 — 교집합을 실제로 확인하기 위해서.
    # (무고리 분자의 대체 라벨은 행 위치에 의존하므로 여기서 만든 값을 써야 한다)
    train_df["_scaffold"] = [scaffolds[i] for i in sorted(tr_idx)]
    test_df["_scaffold"] = [scaffolds[i] for i in sorted(te_idx)]
    return train_df, test_df, len(groups)


def ecfp_row_bytes(frame, featurize_one) -> list[bytes | None]:
    """행별 ECFP 비트열(bytes). 지문 생성 실패는 None."""
    out: list[bytes | None] = []
    for smi in frame["smiles"]:
        f = featurize_one(smi)
        out.append(None if f is None else bytes(memoryview(f.astype("uint8"))))
    return out


# ------------------------------------------------- 2026-09-03 추가(C-9): 통계 헬퍼
def wilson_ci(successes: int, total: int, z: float = 1.959963984540054):
    """이항비율의 Wilson score 95% 신뢰구간 (하한, 상한). total=0이면 NaN."""
    if total <= 0:
        return float("nan"), float("nan")
    p = successes / total
    denom = 1.0 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    half = (z * ((p * (1 - p) / total + z * z / (4 * total * total)) ** 0.5)) / denom
    return max(0.0, center - half), min(1.0, center + half)


def fmt_ci(successes: int, total: int) -> str:
    lo, hi = wilson_ci(successes, total)
    if lo != lo:  # NaN
        return "(CI 계산 불가)"
    return f"[{lo * 100:.1f}–{hi * 100:.1f}]"


def threshold_sweep(y_true, prob, step: float = 0.05):
    """임계값 스윕 — ICH M7에서는 NPV를 최대화하는 지점 선택이 핵심이다."""
    rows = []
    t = 0.05
    while t <= 0.9501:
        pred = [1 if p >= t else 0 for p in prob]
        TP = sum(1 for a, b in zip(y_true, pred) if a == 1 and b == 1)
        TN = sum(1 for a, b in zip(y_true, pred) if a == 0 and b == 0)
        FP = sum(1 for a, b in zip(y_true, pred) if a == 0 and b == 1)
        FN = sum(1 for a, b in zip(y_true, pred) if a == 1 and b == 0)
        sens = TP / (TP + FN) if (TP + FN) else float("nan")
        spec = TN / (TN + FP) if (TN + FP) else float("nan")
        npv = TN / (TN + FN) if (TN + FN) else float("nan")
        rows.append({"threshold": round(t, 2), "TP": TP, "TN": TN, "FP": FP, "FN": FN,
                     "sensitivity": sens, "specificity": spec, "npv": npv})
        t += step
    return rows


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- 메인
def main() -> int:
    say(BAR)
    say("  PRIME-DS  Layer C — Ames 데이터 기반 재학습 및 검증셋 구축")
    say(BAR)
    say()

    try:
        import numpy as np
        import pandas as pd
        from rdkit import Chem, RDLogger
        from sklearn.ensemble import RandomForestClassifier
    except ImportError as exc:
        say(f"[오류] 라이브러리 없음: {exc}")
        say("       pip install rdkit scikit-learn pandas numpy")
        return 1

    RDLogger.DisableLog("rdApp.*")

    # 2026-09-03 수정(C-1): --split random|scaffold 인자 추가 (기본 scaffold)
    argv = sys.argv[1:]
    split_mode = DEFAULT_SPLIT
    positional: list[str] = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--split" and i + 1 < len(argv):
            split_mode = argv[i + 1].strip().lower()
            i += 2
            continue
        if a.startswith("--split="):
            split_mode = a.split("=", 1)[1].strip().lower()
            i += 1
            continue
        if a.startswith("--"):
            say(f"[오류] 알 수 없는 옵션: {a}")
            say("       사용 가능한 옵션: --split random|scaffold")
            return 1
        positional.append(a)
        i += 1

    if split_mode not in ("random", "scaffold"):
        say(f"[오류] --split 값은 random 또는 scaffold 여야 합니다 (받은 값: {split_mode!r})")
        return 1

    src = find_input_file(positional[0] if positional else None)
    say(f"입력 파일: {src.name}  ({src.stat().st_size / 1024:.0f} KB)")
    say(f"분할 방식: {split_mode}"
        f"{'  (Murcko 스캐폴드 — 같은 골격은 한쪽에만)' if split_mode == 'scaffold' else '  (층화 랜덤 — 기존 재현용)'}")

    df = read_table(src)
    say(f"원본 행 수: {len(df):,}   열: {list(df.columns)[:8]}{' ...' if df.shape[1] > 8 else ''}")

    smi_col = pick_column(list(df.columns), SMILES_HINTS)
    lab_col = pick_column([c for c in df.columns if c != smi_col], LABEL_HINTS)
    if smi_col is None or lab_col is None:
        say(f"[오류] 구조/라벨 열을 찾지 못했습니다. 열 목록: {list(df.columns)}")
        return 1
    say(f"구조 열: '{smi_col}'   라벨 열: '{lab_col}'")
    say()

    # --- 1. 표준화 및 중복 제거 -------------------------------------
    # 2026-09-03 수정(C-1a/C-1b/C-14): rdMolStandardize 표준화 + InChIKey14 중복제거,
    #                                  라벨 동일 중복(dup_concordant) 카운터 추가
    say("[1/6] 구조 표준화 (염 제거 · 전하 중화 · 호변이성체 정규화) 및 중복 제거")
    records: dict[str, dict] = {}
    bad_struct = bad_label = 0
    std_failed = 0
    dup_concordant = 0          # InChIKey14 기준 — 라벨이 같아 흡수된 중복 (행 수)
    dup_conflict_rows = 0       # 라벨이 달라 상충으로 표시된 중복 (행 수)
    dup_concordant_full = 0     # 참고: 기존 기준(전체 InChIKey)이었다면 흡수됐을 중복
    seen_full: dict[str, int] = {}
    for _, row in df.iterrows():
        label = normalize_label(row[lab_col])
        if label is None:
            bad_label += 1
            continue
        mol = Chem.MolFromSmiles(str(row[smi_col]).strip())
        if mol is None:
            bad_struct += 1
            continue
        std = standardize_mol(mol)
        if std is None:
            std_failed += 1
            std = mol  # 표준화 실패 시 원본 구조로 진행 (건수는 리포트에 남긴다)

        full_key = Chem.MolToInchiKey(std)
        key = skeleton_key(full_key)          # InChIKey 앞 14자 = 연결성 블록
        canon = Chem.MolToSmiles(std)

        if full_key in seen_full and seen_full[full_key] == label:
            dup_concordant_full += 1
        seen_full.setdefault(full_key, label)

        if key in records:
            if records[key]["label"] != label:
                records[key]["conflict"] = True
                dup_conflict_rows += 1
            else:
                dup_concordant += 1
            continue
        records[key] = {"skeleton_key": key, "inchikey": full_key, "smiles": canon,
                        "label": label, "conflict": False}

    conflicts = [k for k, v in records.items() if v["conflict"]]
    for k in conflicts:
        del records[k]

    clean = pd.DataFrame(list(records.values()))[["skeleton_key", "inchikey", "smiles", "label"]]
    n_pos = int((clean["label"] == 1).sum())
    n_neg = int((clean["label"] == 0).sum())
    say(f"      구조 오류 제외 {bad_struct:,}건 · 라벨 불명 제외 {bad_label:,}건 "
        f"· 표준화 실패(원본 사용) {std_failed:,}건")
    say(f"      라벨 동일 중복 흡수 {dup_concordant:,}건   ← InChIKey14(연결성) 기준")
    say(f"      (참고) 전체 InChIKey 기준이었다면 흡수됐을 중복 {dup_concordant_full:,}건")
    say(f"      라벨 상충 제외 — 골격 {len(conflicts):,}개 "
        f"(대표행 {len(conflicts):,} + 중복행 {dup_conflict_rows:,} = "
        f"{len(conflicts) + dup_conflict_rows:,}행)")
    say(f"      → 사용 가능 {len(clean):,}건  (양성 {n_pos:,} / 음성 {n_neg:,}, "
        f"양성 비율 {n_pos / max(len(clean), 1) * 100:.1f}%)")
    # 2026-09-03 추가(C-14): 원본 행 수와 결과 건수의 차이를 남김없이 설명한다.
    accounted = (bad_struct + bad_label + dup_concordant + dup_conflict_rows
                 + len(conflicts) + len(clean))
    say(f"      행 수 대조: 원본 {len(df):,} = 구조오류 {bad_struct:,} + 라벨불명 {bad_label:,}"
        f" + 라벨동일중복 {dup_concordant:,} + 상충중복행 {dup_conflict_rows:,}"
        f" + 상충대표행 {len(conflicts):,} + 사용 {len(clean):,}"
        f"   → 합계 {accounted:,}"
        f"{'  ✅ 일치' if accounted == len(df) else '  ⚠️ 불일치 — 확인 필요'}")
    say()

    if len(clean) < 200:
        say("[오류] 데이터가 너무 적습니다. 파일을 확인해주세요.")
        return 1

    # --- 2. 분할 ---------------------------------------------------
    # 2026-09-03 수정(C-1c): Murcko 스캐폴드 분할을 기본값으로, 랜덤은 --split random
    n_scaffolds = 0
    if split_mode == "scaffold":
        say(f"[2/6] 학습 {100 - int(TEST_FRACTION * 100)}% / 검증 {int(TEST_FRACTION * 100)}% "
            "Murcko 스캐폴드 분할 (결정적 · seed 무관)")
        train_df, test_df, n_scaffolds = scaffold_split(clean, TEST_FRACTION)
        # 스캐폴드 교집합을 실제로 계산해서 보고한다 (설계상 0이어야 한다)
        tr_sc = set(train_df["_scaffold"])
        te_sc = set(test_df["_scaffold"])
        overlap_sc = len(tr_sc & te_sc)
        train_df = train_df.drop(columns=["_scaffold"])
        test_df = test_df.drop(columns=["_scaffold"])
        assert overlap_sc == 0, f"스캐폴드가 학습·검증에 걸쳐 있습니다: {overlap_sc}개"
        say(f"      고유 스캐폴드 그룹 {n_scaffolds:,}개")
        say(f"      학습 {len(train_df):,}건 · 검증 {len(test_df):,}건 "
            f"(스캐폴드 교집합 {overlap_sc}개 — 실측, 학습 {len(tr_sc):,} / 검증 {len(te_sc):,}개)")
        say(f"      학습 양성비 {(train_df['label'] == 1).mean() * 100:.1f}% · "
            f"검증 양성비 {(test_df['label'] == 1).mean() * 100:.1f}%")
        say("      ※ 스캐폴드 분할은 층화하지 않으므로 양쪽 양성비가 다를 수 있습니다"
            " (의도된 동작).")
    else:
        say(f"[2/6] 학습 {100 - int(TEST_FRACTION * 100)}% / 검증 {int(TEST_FRACTION * 100)}% "
            f"층화 랜덤 분할 (seed={SPLIT_SEED})")
        from sklearn.model_selection import train_test_split

        train_df, test_df = train_test_split(
            clean, test_size=TEST_FRACTION, random_state=SPLIT_SEED, stratify=clean["label"]
        )
        say(f"      학습 {len(train_df):,}건 · 검증 {len(test_df):,}건")
        say("      ※ 랜덤 분할은 스캐폴드가 학습·검증에 걸쳐 있어 낙관 편향이 있습니다.")
    say()

    # --- 3. ECFP 완전 일치 누출 검사 --------------------------------
    # 2026-09-03 추가(C-1d): 학습/검증 ECFP 행렬에 완전 일치 행이 없어야 한다.
    say("[3/6] 학습·검증 ECFP 완전 일치(특징벡터 누출) 검사")
    import p63_deepchem_predictor as dp

    tr_bytes = ecfp_row_bytes(train_df, dp.featurize_smiles_one)
    te_bytes = ecfp_row_bytes(test_df, dp.featurize_smiles_one)
    tr_index: dict[bytes, int] = {}
    for i, b in enumerate(tr_bytes):
        if b is not None:
            tr_index.setdefault(b, i)

    leak_rows = [i for i, b in enumerate(te_bytes) if b is not None and b in tr_index]
    n_leak = len(leak_rows)
    say(f"      지문 radius={dp.FP_RADIUS} · bits={dp.FP_BITS} (useChirality=False)")
    if n_leak:
        say(f"      ⚠️ 완전 일치 {n_leak:,}건 발견 — 검증셋에서 제외합니다.")
        for i in leak_rows[:5]:
            j = tr_index[te_bytes[i]]
            say(f"         · 검증 {str(test_df.iloc[i]['smiles'])[:52]:54} (y={int(test_df.iloc[i]['label'])})")
            say(f"           학습 {str(train_df.iloc[j]['smiles'])[:52]:54} (y={int(train_df.iloc[j]['label'])})")
        if n_leak > 5:
            say(f"         ... 외 {n_leak - 5:,}건")
        keep = [i for i in range(len(test_df)) if i not in set(leak_rows)]
        test_df = test_df.iloc[keep].copy()
        te_bytes = [te_bytes[i] for i in keep]
    else:
        say("      완전 일치 0건")

    residual = sum(1 for b in te_bytes if b is not None and b in tr_index)
    assert residual == 0, f"ECFP 완전 일치 행이 {residual}건 남아 있습니다"
    say(f"      검사 통과 — 검증셋 {len(test_df):,}건 (제외 {n_leak:,}건)")
    say()

    # --- 4. Layer A 중복 제거 --------------------------------------
    # 2026-09-03 수정(C-2): p60의 실제 조회는
    #   `WHERE inchikey = ? OR canonical_smiles = ?` 이므로 두 키 모두로 제외한다.
    #   (기존에는 inchikey만 봐서 canonical SMILES 경로로 누출이 발생했다)
    say("[4/6] 검증셋에서 실측 DB(m7_empirical.db) 중복 화합물 제외 (InChIKey + canonical SMILES)")
    db_keys: set[str] = set()
    db_smiles: set[str] = set()
    if EMPIRICAL_DB.exists():
        import sqlite3

        try:
            conn = sqlite3.connect(f"file:{EMPIRICAL_DB}?mode=ro", uri=True)
            rows = conn.execute(
                "SELECT inchikey, canonical_smiles FROM m7_empirical").fetchall()
            conn.close()
            for ik, smi in rows:
                if ik:
                    db_keys.add(str(ik))
                if smi:
                    db_smiles.add(str(smi))
                    m = Chem.MolFromSmiles(str(smi))
                    if m is not None:
                        db_smiles.add(Chem.MolToSmiles(m))
                        if not ik:
                            db_keys.add(Chem.MolToInchiKey(m))
            say(f"      실측 DB 보유 구조: InChIKey {len(db_keys):,}개 · "
                f"canonical SMILES {len(db_smiles):,}개")
        except Exception as exc:  # noqa: BLE001
            say(f"      (DB를 읽지 못했습니다: {exc} — 중복 제거를 건너뜁니다)")
    else:
        say("      (m7_empirical.db 없음 — 중복 제거를 건너뜁니다)")

    before = len(test_df)
    removed_by_ikey = removed_by_smiles = 0
    if db_keys or db_smiles:
        hit_ikey = test_df["inchikey"].isin(db_keys)
        hit_smi = test_df["smiles"].isin(db_smiles)
        removed_by_ikey = int((hit_ikey & ~hit_smi).sum())
        removed_by_smiles = int((hit_smi & ~hit_ikey).sum())
        both = int((hit_ikey & hit_smi).sum())
        test_df = test_df[~(hit_ikey | hit_smi)].copy()
        say(f"      InChIKey 단독 일치 {removed_by_ikey:,}건 · "
            f"canonical SMILES 단독 일치 {removed_by_smiles:,}건 · 양쪽 일치 {both:,}건")
        if removed_by_smiles:
            say("      ※ SMILES 단독 일치분은 기존(InChIKey 전용) 기준으로는 걸러지지 않던"
                " 누출 경로입니다.")
    removed = before - len(test_df)
    say(f"      검증셋에서 {removed:,}건 제외 → 최종 검증셋 {len(test_df):,}건")
    t_pos = int((test_df["label"] == 1).sum())
    t_neg = int((test_df["label"] == 0).sum())
    say(f"      검증셋 구성: 양성 {t_pos:,} / 음성 {t_neg:,}")
    say()

    if t_pos < 10 or t_neg < 10:
        say("[경고] 검증셋의 한쪽 클래스가 너무 적습니다. 지표 해석에 주의하세요.")
        say()

    # --- 5. 학습 ---------------------------------------------------
    say(f"[5/6] RandomForest 학습 (트리 {RF_TREES}개, seed={RF_SEED})")

    def featurize(frame):
        X, keep = [], []
        for i, smi in enumerate(frame["smiles"]):
            f = dp.featurize_smiles_one(smi)
            if f is not None:
                X.append(f)
                keep.append(i)
        return np.vstack(X), keep

    Xtr, keep_tr = featurize(train_df)
    ytr = train_df["label"].to_numpy()[keep_tr]

    model = RandomForestClassifier(n_estimators=RF_TREES, random_state=RF_SEED)
    model.fit(Xtr, ytr)
    say(f"      학습 완료 — {len(ytr):,}건 사용")

    import pickle

    meta = {
        "model_label": dp.AMES_MODEL_LABEL,
        "endpoint": "Ames bacterial reverse mutation",
        "source_file": src.name,
        "fp_radius": dp.FP_RADIUS,
        "fp_bits": dp.FP_BITS,
        "n_estimators": RF_TREES,
        "random_state": RF_SEED,
        "n_train": int(len(ytr)),
        "n_train_positive": int((ytr == 1).sum()),
        "n_train_negative": int((ytr == 0).sum()),
        "split_seed": SPLIT_SEED,
        "test_fraction": TEST_FRACTION,
        "n_validation": int(len(test_df)),
        "excluded_by_empirical_db": int(removed),
        # 2026-09-03 추가: 재현에 필요한 분할·중복제거 조건을 메타에 남긴다
        "split_mode": split_mode,
        "dedup_key": "InChIKey[:14] (connectivity block)",
        "standardization": "LargestFragmentChooser + Uncharger + TautomerEnumerator.Canonicalize",
        "n_scaffold_groups": int(n_scaffolds),
        "dup_concordant": int(dup_concordant),
        "dup_concordant_full_inchikey": int(dup_concordant_full),
        "n_label_conflict": int(len(conflicts)),
        "n_label_conflict_rows": int(dup_conflict_rows),
        "ecfp_leak_removed": int(n_leak),
        "excluded_by_db_inchikey_only": int(removed_by_ikey),
        "excluded_by_db_smiles_only": int(removed_by_smiles),
    }
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    with open(MODEL_PKL, "wb") as fh:
        pickle.dump({"model": model, "meta": meta}, fh)
    say(f"      저장: models/{MODEL_PKL.name} ({MODEL_PKL.stat().st_size / 1024 / 1024:.1f} MB)")
    say()

    # --- 6. 검증 ---------------------------------------------------
    say("[6/6] 검증셋으로 Layer C 단독 성능 평가")
    Xte, keep_te = featurize(test_df)
    yte = test_df["label"].to_numpy()[keep_te]
    n_featurize_failed = len(test_df) - len(keep_te)
    if n_featurize_failed:
        say(f"      지문 생성 실패 {n_featurize_failed:,}건 — 지표 계산에서 제외 "
            "(음성으로 처리하지 않습니다)")
    pos_idx = list(model.classes_).index(1)
    prob = model.predict_proba(Xte)[:, pos_idx]
    pred = (prob >= 0.5).astype(int)

    TP = int(((pred == 1) & (yte == 1)).sum())
    TN = int(((pred == 0) & (yte == 0)).sum())
    FP = int(((pred == 1) & (yte == 0)).sum())
    FN = int(((pred == 0) & (yte == 1)).sum())

    def safe(a, b):
        return a / b if b else float("nan")

    sens = safe(TP, TP + FN)
    spec = safe(TN, TN + FP)
    conc = safe(TP + TN, len(yte))
    npv = safe(TN, TN + FN)
    ppv = safe(TP, TP + FP)
    bacc = (sens + spec) / 2
    denom = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
    mcc = ((TP * TN - FP * FN) / denom) if denom else float("nan")

    say()
    say(f"      평가 대상 {len(yte):,}건 (검증셋 {len(test_df):,}건 − 지문 실패 {n_featurize_failed:,}건)")
    say("      혼동행렬 (기준 0.5)")
    say("                        실제 양성    실제 음성")
    say(f"        예측 양성       {TP:8,}    {FP:9,}")
    say(f"        예측 음성       {FN:8,}    {TN:9,}")
    say()
    # 2026-09-03 수정(C-9): Wilson score 95% CI 병기
    say(f"        Sensitivity  (민감도)      {sens * 100:6.1f} %  {fmt_ci(TP, TP + FN)}")
    say(f"        Specificity  (특이도)      {spec * 100:6.1f} %  {fmt_ci(TN, TN + FP)}")
    say(f"        Concordance  (일치율)      {conc * 100:6.1f} %  {fmt_ci(TP + TN, len(yte))}")
    say(f"        NPV          (음성예측도)   {npv * 100:6.1f} %  {fmt_ci(TN, TN + FN)}   ← ICH M7 핵심 지표")
    say(f"        PPV          (양성예측도)   {ppv * 100:6.1f} %  {fmt_ci(TP, TP + FP)}")
    say(f"        Balanced accuracy          {bacc * 100:6.1f} %")
    say(f"        MCC                        {mcc:6.3f}")
    say(f"        False negative rate        {(1 - sens) * 100:6.1f} %")
    say()

    # 2026-09-03 추가(C-9): 임계값 비의존 지표 (ROC-AUC / PR-AUC)
    try:
        from sklearn.metrics import average_precision_score, roc_auc_score

        if len(set(yte.tolist())) > 1:
            roc = float(roc_auc_score(yte, prob))
            prauc = float(average_precision_score(yte, prob))
            say(f"        ROC-AUC                    {roc:6.3f}   (임계값 비의존)")
            say(f"        PR-AUC (average precision) {prauc:6.3f}   "
                f"기준선(양성비) {float((yte == 1).mean()):.3f}")
        else:
            roc = prauc = float("nan")
            say("        (검증셋에 한 클래스만 있어 ROC/PR-AUC를 계산할 수 없습니다)")
    except Exception as exc:  # noqa: BLE001
        roc = prauc = float("nan")
        say(f"        (ROC/PR-AUC 계산 실패: {exc})")
    say()

    # 2026-09-03 추가(C-9): 임계값 스윕 — ICH M7은 NPV 최대화 지점이 핵심
    say("      임계값 스윕 (0.05 간격) — ICH M7에서는 NPV를 최대화하는 지점을 씁니다")
    say(f"        {'임계값':>7} {'Sens':>8} {'Spec':>8} {'NPV':>8}   {'FN':>5} {'FP':>5}")
    sweep = threshold_sweep(yte.tolist(), prob.tolist())
    for r in sweep:
        say(f"        {r['threshold']:>7.2f} {r['sensitivity'] * 100:7.1f}% "
            f"{r['specificity'] * 100:7.1f}% {r['npv'] * 100:7.1f}%   "
            f"{r['FN']:>5,} {r['FP']:>5,}")
    valid_sweep = [r for r in sweep if r["npv"] == r["npv"]]
    if valid_sweep:
        best = max(valid_sweep, key=lambda r: r["npv"])
        say(f"        → NPV 최대: 임계값 {best['threshold']:.2f} 에서 "
            f"NPV {best['npv'] * 100:.1f}% {fmt_ci(best['TN'], best['TN'] + best['FN'])} · "
            f"Sens {best['sensitivity'] * 100:.1f}% · Spec {best['specificity'] * 100:.1f}%")
    say()

    # 저장
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    # 2026-09-03 수정(C-3b): 이후 단계가 «행 삭제»가 아니라 «사유 표기»로 제외하도록
    #                        excluded_reason 컬럼을 처음부터 만들어 둔다.
    train_df = train_df.copy()
    test_df = test_df.copy()
    if "excluded_reason" not in test_df.columns:
        test_df["excluded_reason"] = ""
    train_df.to_csv(TRAIN_CSV, index=False, encoding="utf-8-sig")
    test_df.to_csv(VALID_CSV, index=False, encoding="utf-8-sig")
    say(f"      학습셋 저장: data/{TRAIN_CSV.name}")
    say(f"      검증셋 저장: data/{VALID_CSV.name}   (작업본)")

    # 2026-09-03 추가(C-3a): 검증셋 동결 — 이후 스크립트는 읽기 전용으로 참조한다.
    FROZEN_DIR.mkdir(parents=True, exist_ok=True)
    if FROZEN_VALID_CSV.exists():
        say(f"      [경고] 이미 동결본이 있습니다: data/frozen/{FROZEN_VALID_CSV.name}")
        say(f"             기존 해시 {sha256_of(FROZEN_VALID_CSV)[:16]}… 를 덮어씁니다.")
        say("             (심사 대응용으로 이전 동결본을 보관해야 한다면 먼저 백업하세요)")
        try:
            FROZEN_VALID_CSV.chmod(0o644)
        except Exception:  # noqa: BLE001
            pass
    test_df.to_csv(FROZEN_VALID_CSV, index=False, encoding="utf-8-sig")
    digest = sha256_of(FROZEN_VALID_CSV)
    FROZEN_VALID_SHA.write_text(
        f"{digest}  {FROZEN_VALID_CSV.name}\n"
        f"n_rows={len(test_df)}\n"
        f"n_positive={int((test_df['label'] == 1).sum())}\n"
        f"n_negative={int((test_df['label'] == 0).sum())}\n"
        f"split_mode={split_mode}\n",
        encoding="utf-8",
    )
    try:  # 실수로 덮어쓰지 못하게 읽기 전용으로 잠근다
        FROZEN_VALID_CSV.chmod(0o444)
    except Exception:  # noqa: BLE001
        pass
    say(f"      검증셋 동결: data/frozen/{FROZEN_VALID_CSV.name}   ← 논문 external validation set")
    say(f"      SHA-256    : {digest}")
    say(f"      건수       : {len(test_df):,}건 "
        f"(양성 {int((test_df['label'] == 1).sum()):,} / 음성 {int((test_df['label'] == 0).sum()):,})")
    say("      ※ 이후 s4/s5는 이 파일을 읽기 전용으로 사용합니다. 지표를 본 뒤")
    say("        검증셋을 바꾸는 일이 없도록 해시로 고정합니다.")
    say()

    # --- 상식 점검 -------------------------------------------------
    say(BAR)
    say("  상식 점검 — 알려진 변이원성 물질")
    say(BAR)
    dp._model_ready = False
    dp._model = None
    dp._model_error = None
    checks = [
        ("CN(C)N=O", "NDMA 니트로사민", "양성"),
        ("CCN(CC)N=O", "NDEA 니트로사민", "양성"),
        ("Nc1ccccc1", "아닐린", "양성"),
        ("C1CO1", "에틸렌옥사이드", "양성"),
        ("CC(=O)O", "아세트산", "음성"),
        ("CC(=O)Oc1ccccc1C(=O)O", "아스피린", "음성"),
    ]
    res = dp.predict_advanced_properties([c[0] for c in checks])
    hit = 0
    for smi, name, expect in checks:
        p = res.get(smi, {}).get("Ames_Mutagenicity_Prob")
        if p is None:
            say(f"   {name:22s}  예측 실패")
            continue
        got = "양성" if p >= 0.5 else "음성"
        mark = "O" if got == expect else "X"
        if got == expect:
            hit += 1
        say(f"   [{mark}] {name:22s}  P = {p:.3f}   예상 {expect} / 예측 {got}")
    say()
    say(f"   {hit}/{len(checks)} 일치")
    say()

    if hit >= len(checks) - 1:
        say("Layer C가 이제 실제 Ames 변이원성을 예측합니다.")
    else:
        say("일부 상식 점검에 실패했습니다. 데이터 파일과 위 지표를 확인하세요.")

    say()
    say("PRIME-DS를 실행하면 새 모델이 자동으로 사용됩니다.")
    say(f"되돌리려면 models/{MODEL_PKL.name} 파일을 지우면 됩니다.")

    REPORT_TXT.write_text("\n".join(_log_lines), encoding="utf-8")
    print(f"\n기록 저장: data/{REPORT_TXT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
