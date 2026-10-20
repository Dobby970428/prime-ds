"""
Layer C 적용범위(applicability domain) 판정.

ICH M7(R2)은 통계 모델의 적용범위를 벗어나는 화합물을 전문가 검토로 회부할 것을
요구한다. 이 모듈은 질의 구조가 Layer C 학습셋에서 얼마나 떨어져 있는지를
d_max 로 반환한다.

    d_max(x) = max over t in TRAIN  Tanimoto( ECFP4(x), ECFP4(t) )
    d_max < 0.30  →  out-of-domain

정의·파라미터는 s11_stats_ad_analysis.py 와 **완전히 동일**하다. 논문에 보고된
수치와 같은 척도를 유지하기 위함이며, 둘 중 하나를 고치면 반드시 다른 하나도
함께 고쳐야 한다. 표준화는 p74_structure_standardize 가 단일 기준을 갖는다.

    학습셋   : 논문자료/01_데이터셋/ames_train.csv
               2026-09-14 스캐폴드 분할 재실행 후 5,309건
    표준화   : p74 (염 제거 → 전하 중화 → 호변이성체 정규화) — s1 과 동일
    지문     : Morgan radius 2 · 2048 bits (ECFP4), 키랄성 미부호화
    임계값   : 0.30

설계 원칙
    · 이 모듈은 판정 결과를 바꾸지 않는다. 등급 부여에 관여하지 않고
      「이 예측을 얼마나 믿을 수 있는가」만 별도로 보고한다.
    · 어떤 이유로든 산출이 불가능하면 예외를 올리지 않고 available=False 로
      반환한다. 호출 측의 기존 동작이 영향을 받아서는 안 된다.
    · 학습셋 지문은 최초 1회만 만들고 pickle 로 캐시한다. 이후 호출은
      BulkTanimotoSimilarity 한 번으로 끝난다.

단독 실행
    python p73_applicability_domain.py
"""

from __future__ import annotations

import pickle
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent

TRAIN_CSV = ROOT / "논문자료" / "01_데이터셋" / "ames_train.csv"
CACHE_DIR = ROOT / "cache"
CACHE_FILE = CACHE_DIR / "ad_train_fps.pkl"

# s11_stats_ad_analysis.py 와 동일해야 하는 값들
AD_THRESHOLD = 0.30
ECFP_RADIUS = 2
ECFP_BITS = 2048

SMILES_COLUMN_HINTS = ("smiles", "canonical_smiles", "structure")

_lock = threading.RLock()
_train_fps: list | None = None
_train_error: str | None = None
_generator = None


# ---------------------------------------------------------------------------
# 결과 표현
# ---------------------------------------------------------------------------
@dataclass
class ApplicabilityDomainResult:
    """d_max 와 그 해석. available=False 이면 나머지 값은 의미가 없다."""

    d_max: float | None = None
    in_domain: bool | None = None
    threshold: float = AD_THRESHOLD
    n_train: int = 0
    available: bool = False
    reason: str = ""

    @property
    def label(self) -> str:
        if not self.available or self.d_max is None:
            return "AD not evaluated"
        return "In domain" if self.in_domain else "Out of domain"

    @property
    def summary(self) -> str:
        """보고서에 한 줄로 넣을 문자열."""
        if not self.available or self.d_max is None:
            return f"AD not evaluated ({self.reason or 'unavailable'})"
        return (
            f"d_max = {self.d_max:.3f} vs training set (n = {self.n_train:,}); "
            f"{self.label.lower()} at threshold {self.threshold:.2f}"
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "Applicability Domain d_max": round(self.d_max, 4) if self.d_max is not None else None,
            "Applicability Domain Status": self.label,
            "Applicability Domain Threshold": self.threshold,
            "Applicability Domain Training Size": self.n_train,
        }


# ---------------------------------------------------------------------------
# 지문
# ---------------------------------------------------------------------------
def _get_generator():
    global _generator
    if _generator is None:
        from rdkit.Chem import rdFingerprintGenerator

        _generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=ECFP_RADIUS, fpSize=ECFP_BITS
        )
    return _generator


def fingerprint(smiles: str):
    """
    SMILES → ECFP4 비트벡터. 실패 시 None. (s11 의 fingerprint 와 동일)

    2026-09-14 수정 — p74 로 구조를 표준화한 뒤 지문을 만든다. s1 이
    학습셋을 표준화해 저장하므로, 질의 구조도 같은 처리를 거쳐야 같은
    척도가 된다. 표준화 없이 재면 생성 구조의 d_max 가 낮게 나온다.
    """
    from p74_structure_standardize import standardized_mol_from_smiles

    mol = standardized_mol_from_smiles(smiles)
    if mol is None:
        return None
    return _get_generator().GetFingerprint(mol)


# ---------------------------------------------------------------------------
# 학습셋 로드 및 캐시
# ---------------------------------------------------------------------------
def _cache_signature() -> tuple:
    """
    학습 파일이 바뀌거나 지문 파라미터가 바뀌면 캐시를 버리기 위한 서명.

    2026-09-14 — 표준화 버전을 서명에 넣었다. 표준화 정의가 바뀌면 지문도
    달라지는데, 파일 mtime 만 보면 낡은 캐시를 그대로 쓰게 된다.
    """
    from p74_structure_standardize import STANDARDIZER_VERSION

    try:
        stat = TRAIN_CSV.stat()
        return (TRAIN_CSV.name, stat.st_size, int(stat.st_mtime),
                ECFP_RADIUS, ECFP_BITS, STANDARDIZER_VERSION)
    except OSError:
        return ("missing", 0, 0, ECFP_RADIUS, ECFP_BITS, STANDARDIZER_VERSION)


def _pick_smiles_column(columns: list[str]) -> str | None:
    lowered = {c.lower(): c for c in columns}
    for hint in SMILES_COLUMN_HINTS:
        if hint in lowered:
            return lowered[hint]
    for col in columns:
        if "smi" in col.lower():
            return col
    return None


def _load_from_cache() -> list | None:
    if not CACHE_FILE.exists():
        return None
    try:
        with open(CACHE_FILE, "rb") as handle:
            payload = pickle.load(handle)
        if not isinstance(payload, dict):
            return None
        if tuple(payload.get("signature", ())) != _cache_signature():
            return None
        fps = payload.get("fps")
        return fps if fps else None
    except Exception:
        return None


def _save_to_cache(fps: list) -> None:
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with open(CACHE_FILE, "wb") as handle:
            pickle.dump({"signature": _cache_signature(), "fps": fps}, handle)
    except Exception:
        # 캐시 저장 실패는 기능에 영향을 주지 않는다.
        pass


def _build_train_fps() -> list:
    import pandas as pd
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")

    if not TRAIN_CSV.exists():
        raise FileNotFoundError(f"학습셋을 찾을 수 없습니다: {TRAIN_CSV}")

    frame = pd.read_csv(TRAIN_CSV)
    column = _pick_smiles_column(list(frame.columns))
    if column is None:
        raise ValueError(f"구조 열을 찾지 못했습니다. 열 목록: {list(frame.columns)}")

    fps = []
    for value in frame[column].tolist():
        fp = fingerprint(value)
        if fp is not None:
            fps.append(fp)

    if not fps:
        raise ValueError("학습셋에서 유효한 구조를 하나도 읽지 못했습니다.")
    return fps


def _ensure_train_fps() -> list | None:
    """학습셋 지문을 한 번만 준비한다. 스레드 안전. 실패는 기록만 하고 None."""
    global _train_fps, _train_error

    with _lock:
        if _train_fps is not None:
            return _train_fps
        if _train_error is not None:
            return None

        try:
            cached = _load_from_cache()
            if cached is not None:
                _train_fps = cached
                return _train_fps

            fps = _build_train_fps()
            _save_to_cache(fps)
            _train_fps = fps
            return _train_fps
        except Exception as exc:  # noqa: BLE001
            _train_error = f"{type(exc).__name__}: {exc}"
            return None


# ---------------------------------------------------------------------------
# 공개 API
# ---------------------------------------------------------------------------
def applicability_domain(smiles: str) -> ApplicabilityDomainResult:
    """
    질의 구조의 d_max 를 산출한다.

    어떤 이유로도 예외를 올리지 않는다. 산출 불가 시 available=False.
    """
    train_fps = _ensure_train_fps()
    if train_fps is None:
        return ApplicabilityDomainResult(
            available=False,
            reason=_train_error or "training fingerprints unavailable",
        )

    try:
        from rdkit.DataStructs import BulkTanimotoSimilarity

        query_fp = fingerprint(smiles)
        if query_fp is None:
            return ApplicabilityDomainResult(
                n_train=len(train_fps),
                available=False,
                reason="invalid SMILES",
            )

        value = float(max(BulkTanimotoSimilarity(query_fp, train_fps)))
        return ApplicabilityDomainResult(
            d_max=value,
            in_domain=value >= AD_THRESHOLD,
            threshold=AD_THRESHOLD,
            n_train=len(train_fps),
            available=True,
        )
    except Exception as exc:  # noqa: BLE001
        return ApplicabilityDomainResult(
            n_train=len(train_fps),
            available=False,
            reason=f"{type(exc).__name__}: {exc}"[:200],
        )


def d_max(smiles: str) -> float | None:
    """편의 함수. 산출 불가 시 None."""
    return applicability_domain(smiles).d_max


def is_out_of_domain(smiles: str) -> bool:
    """
    적용범위를 벗어났으면 True.

    ⚠️ 산출이 불가능한 경우에도 False 를 반환한다. 「모른다」와 「범위 안이다」는
    다르므로, 그 구분이 필요하면 applicability_domain() 의 available 을 볼 것.
    """
    result = applicability_domain(smiles)
    return bool(result.available and result.in_domain is False)


def status() -> dict[str, Any]:
    """진단용. 엔진 정보 화면에 그대로 실을 수 있다."""
    train_fps = _ensure_train_fps()
    return {
        "ready": train_fps is not None,
        "n_train": len(train_fps) if train_fps else 0,
        "train_csv": str(TRAIN_CSV),
        "train_csv_exists": TRAIN_CSV.exists(),
        "cache_file": str(CACHE_FILE),
        "cache_exists": CACHE_FILE.exists(),
        "threshold": AD_THRESHOLD,
        "fingerprint": f"Morgan radius {ECFP_RADIUS}, {ECFP_BITS} bits (ECFP4)",
        "error": _train_error,
    }


# ---------------------------------------------------------------------------
# 단독 실행 — 자가 점검
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import json

    print("=" * 70)
    print("  p73_applicability_domain  자가 점검")
    print("=" * 70)
    print()

    info = status()
    print(json.dumps(info, ensure_ascii=False, indent=2))
    print()

    if not info["ready"]:
        print("[중단] 학습셋 지문을 준비하지 못했습니다. 위 error 를 확인하세요.")
        raise SystemExit(1)

    samples = [
        ("CN(C)N=O", "NDMA — 니트로사민, 학습셋에 있을 가능성이 높음"),
        ("CC(=O)O", "아세트산 — 단순 구조"),
        ("c1ccccc1[N+](=O)[O-]", "니트로벤젠"),
        ("CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)N=O)o1", "라니티딘 NDSRI — 생성 구조"),
        ("not-a-smiles", "잘못된 입력 — available=False 여야 정상"),
    ]

    print(f"  {'d_max':>8}  {'판정':<16} 구조")
    print("  " + "-" * 66)
    for smi, note in samples:
        res = applicability_domain(smi)
        shown = f"{res.d_max:.3f}" if res.d_max is not None else "   n/a"
        print(f"  {shown:>8}  {res.label:<16} {note}")
    print()
    print("  ※ s11_stats_ad_analysis.py 산출과의 정합 기준")
    print("       검증셋 1,440건      중앙값 0.646  (p25 0.538 · p75 0.769)")
    print("       니트로사민 4건       최소 0.286 · 중앙값 0.742 · 최대 1.000")
    print("       부형제 부가물 44건   최소 0.250 · 중앙값 0.391")
    print("       이원 배합 산물 53건  최소 0.250 · 중앙값 0.388")
    print("     학습셋에 있는 화합물(NDMA 등)이 1.000 으로 나오는 것은 정상입니다.")
