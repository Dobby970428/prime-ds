"""
Layer C 통계 모델 — Ames 변이원성 확률 예측 (RDKit + scikit-learn).

이 모듈은 원래 DeepChem에 의존했으나, DeepChem이 최신 Python(3.13+)을
지원하지 않아 동일한 계산을 RDKit과 scikit-learn만으로 수행하도록 교체되었습니다.

계산 내용은 동일합니다:
  - 분자 지문: ECFP / Morgan fingerprint (radius=2, 2048 bits)
    (DeepChem의 CircularFingerprint(size=2048)는 내부적으로 RDKit의 동일 함수를 호출)
  - 분류기: RandomForestClassifier(n_estimators=100, random_state=0)
  - 학습 데이터: Tox21 SR-p53 endpoint (models/tox21.csv.gz)

원본(DeepChem 판)은 deepchem_predictor.py.backup 에 보존되어 있었습니다.
  ※ 2026-09-03 확인: 이 저장소 스냅샷에 해당 백업 파일은 존재하지 않습니다
     (`*.backup` 0건). 별도 보관 위치를 확인하거나 이 안내를 삭제하세요.

모델 재학습:  python s6_retrain_layer_c_model.py
"""
# 2026-09-03 수정: 위 docstring 의 재학습 안내가 구식 파일명(retrain_layer_c_model.py)
# 이었습니다. 실제 파일은 s6_retrain_layer_c_model.py 입니다.

from __future__ import annotations

import gzip
import os
import pickle
import threading
from pathlib import Path
from typing import Any
from urllib.request import urlretrieve

import numpy as np
import pandas as pd

_BASE_DIR = Path(__file__).resolve().parent
_MODEL_DIR = _BASE_DIR / "models"

# 모델 우선순위: Ames 학습 모델이 있으면 그것을, 없으면 Tox21 대체 모델을 쓴다.
_AMES_MODEL_PATH = _MODEL_DIR / "ames_ecfp_rf.pkl"
_TOX21_MODEL_PATH = _MODEL_DIR / "tox21_ecfp_rf_sklearn.pkl"
_LEGACY_MODEL_PATH = _MODEL_DIR / "tox21_ecfp_rf.pkl"  # DeepChem 형식 (참고용, 사용 안 함)

_MODEL_PATH = _TOX21_MODEL_PATH  # 하위 호환 (재학습 스크립트가 참조)


def _forced_model_choice() -> str:
    """PRIME_LAYER_C_MODEL 환경변수로 강제 지정된 모델 ('' = 미지정)."""
    forced = os.environ.get("PRIME_LAYER_C_MODEL", "").strip().lower()
    if forced == "ames":
        return "ames"
    if forced in ("tox21", "tox"):
        return "tox21"
    return ""


def _resolve_model_path() -> Path:
    """
    사용할 모델 파일을 결정한다.

    환경변수 PRIME_LAYER_C_MODEL 로 강제 지정할 수 있다 (비교 실험용):
        ames   → models/ames_ecfp_rf.pkl        (Ames 실측 학습)
        tox21  → models/tox21_ecfp_rf_sklearn.pkl (Tox21 SR-p53 학습)
    지정하지 않으면 Ames 모델이 있으면 그것을, 없으면 Tox21을 쓴다.
    """
    forced = _forced_model_choice()
    if forced == "ames":
        return _AMES_MODEL_PATH
    if forced == "tox21":
        return _TOX21_MODEL_PATH
    if _AMES_MODEL_PATH.exists():
        return _AMES_MODEL_PATH
    return _TOX21_MODEL_PATH

_TOX21_CSV = _MODEL_DIR / "tox21.csv.gz"
_TOX21_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz"

_TOX21_TASKS = [
    "NR-AR",
    "NR-AR-LBD",
    "NR-AhR",
    "NR-Aromatase",
    "NR-ER",
    "NR-ER-LBD",
    "NR-PPAR-gamma",
    "SR-ARE",
    "SR-ATAD5",
    "SR-HSE",
    "SR-MMP",
    "SR-p53",
]
# 변이원성/Ames와 가장 관련 높은 Tox21 endpoint (stress response p53 pathway)
_MUTAGEN_TASK = "SR-p53"

MODEL_LABEL = "ECFP-2048 + RandomForest (Tox21 SR-p53)"
AMES_MODEL_LABEL = "ECFP-2048 + RandomForest (Ames bacterial mutagenicity)"
MODEL_UNAVAILABLE_LABEL = "unavailable"

# ---------------------------------------------------------------------------
# 엔드포인트 표기
#
# ⚠️ Tox21 대체 모델의 학습 엔드포인트는 SR-p53(진핵세포 p53 스트레스 응답)이며
#    세균 변이원성(Ames)이 아닙니다. 자체 검증 산출물
#    `논문자료/02_검증결과/report_tox21.txt`가 이 대체의 Layer C 단독
#    Sensitivity를 5.1%로 보고합니다. 확률값을 «Ames 변이원성 확률»로
#    부르면 안 됩니다.
# ---------------------------------------------------------------------------
AMES_ENDPOINT = "Ames bacterial reverse mutation"
SURROGATE_ENDPOINT = "Tox21 SR-p53 (eukaryotic p53 stress response)"
SURROGATE_ENDPOINT_NOTE = "surrogate endpoint — not bacterial mutagenicity"


def _endpoint_is_ames(endpoint: str) -> bool:
    return "ames" in str(endpoint or "").lower()


def _endpoint_for_path(path: Path) -> str:
    return AMES_ENDPOINT if path == _AMES_MODEL_PATH else SURROGATE_ENDPOINT


def active_model_label() -> str:
    """현재 사용 중인 모델의 표시용 이름 (리포트에 그대로 기재된다).

    2026-09-03 수정: `_model_meta`가 비면 **파일 존재 여부만** 보고 라벨을
    돌려주어, 로드에 실패했거나 파일이 아예 없어도 그럴듯한 모델명이
    보고서에 찍혔습니다. 실제 로드 성공 여부를 반영하도록 변경.
    """
    model = _ensure_model()
    if model is None:
        return MODEL_UNAVAILABLE_LABEL
    label = (_model_meta or {}).get("model_label")
    if label:
        return str(label)
    path = _resolve_model_path()
    return AMES_MODEL_LABEL if path == _AMES_MODEL_PATH else MODEL_LABEL


def active_endpoint() -> str:
    """실제로 로드된 모델의 학습 엔드포인트 ('' = 사용 불가)."""
    model = _ensure_model()
    if model is None:
        return ""
    endpoint = (_model_meta or {}).get("endpoint")
    if endpoint:
        return str(endpoint)
    return _endpoint_for_path(_resolve_model_path())

FP_RADIUS = 2
FP_BITS = 2048

_model_lock = threading.RLock()
_model: Any = None
_model_meta: dict[str, Any] = {}
_model_ready = False
_model_error: str | None = None

_fp_generator = None


# --------------------------------------------------------------------------
# 분자 지문 (RDKit Morgan / ECFP)
# --------------------------------------------------------------------------
def _get_fp_generator():
    """Morgan fingerprint 생성기 (radius=2, 2048 bits)."""
    global _fp_generator
    if _fp_generator is None:
        from rdkit.Chem import rdFingerprintGenerator

        _fp_generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=FP_RADIUS, fpSize=FP_BITS
        )
    return _fp_generator


def featurize_smiles_one(smiles: str) -> np.ndarray | None:
    """SMILES 하나를 2048비트 지문 배열로 변환. 실패 시 None."""
    from rdkit import Chem

    if not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None

    gen = _get_fp_generator()
    try:
        arr = gen.GetFingerprintAsNumPy(mol)
    except AttributeError:  # 구버전 RDKit 대비
        from rdkit import DataStructs

        bv = gen.GetFingerprint(mol)
        arr = np.zeros((FP_BITS,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(bv, arr)

    arr = np.asarray(arr, dtype=float).ravel()
    if arr.size != FP_BITS:
        return None
    return arr


def _featurize_smiles(smiles_list: list[str]) -> tuple[list[np.ndarray], list[str]]:
    features: list[np.ndarray] = []
    valid_smiles: list[str] = []
    for smi in smiles_list:
        feat = featurize_smiles_one(smi)
        if feat is not None:
            features.append(feat)
            valid_smiles.append(smi)
    return features, valid_smiles


# --------------------------------------------------------------------------
# 학습
# --------------------------------------------------------------------------
def _download_tox21_csv(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        urlretrieve(_TOX21_URL, dest)
    return dest


def train_and_save_model(verbose: bool = False) -> Any:
    """Tox21 SR-p53로 RandomForest를 학습하고 저장한다.

    ⚠️ CLI 전용입니다 (`python s6_retrain_layer_c_model.py`). 웹 요청 처리 중에는
    호출하지 마세요 — Tox21 데이터셋을 내려받느라 수 분간 블로킹됩니다.
    """
    from sklearn.ensemble import RandomForestClassifier

    _download_tox21_csv(_TOX21_CSV)

    with gzip.open(_TOX21_CSV, "rt") as handle:
        df = pd.read_csv(handle)

    task = _MUTAGEN_TASK
    df = df.dropna(subset=["smiles", task])
    df[task] = pd.to_numeric(df[task], errors="coerce")
    df = df.dropna(subset=[task])
    df = df[df[task].isin([0.0, 1.0])]

    if verbose:
        print(f"[1/3] 학습 대상 {len(df):,}건 (SR-p53 라벨 보유)")

    X_rows: list[np.ndarray] = []
    y_rows: list[int] = []
    skipped = 0
    for _, row in df.iterrows():
        feat = featurize_smiles_one(str(row["smiles"]))
        if feat is None:
            skipped += 1
            continue
        X_rows.append(feat)
        y_rows.append(int(row[task]))

    if len(X_rows) < 50:
        raise RuntimeError("Tox21 학습 데이터가 충분하지 않습니다.")

    X = np.vstack(X_rows)
    y = np.asarray(y_rows, dtype=int)

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if verbose:
        print(f"[2/3] 지문 변환 완료 — 사용 {len(X_rows):,}건 / 제외 {skipped:,}건")
        print(f"       양성 {n_pos:,} · 음성 {n_neg:,} (양성 비율 {n_pos / len(y) * 100:.1f}%)")
        print("[3/3] RandomForest 학습 중...")

    model = RandomForestClassifier(n_estimators=100, random_state=0)
    model.fit(X, y)

    meta = {
        "model_label": MODEL_LABEL,
        # 2026-09-03 추가: 확률값의 실제 엔드포인트를 메타에 남긴다 (Ames 아님).
        "endpoint": SURROGATE_ENDPOINT,
        "endpoint_note": SURROGATE_ENDPOINT_NOTE,
        "task": task,
        "fp_radius": FP_RADIUS,
        "fp_bits": FP_BITS,
        "n_train": int(len(y)),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "n_skipped": skipped,
        "n_estimators": 100,
        "random_state": 0,
        "sklearn_native": True,
    }

    _MODEL_DIR.mkdir(parents=True, exist_ok=True)
    with open(_MODEL_PATH, "wb") as handle:
        pickle.dump({"model": model, "meta": meta}, handle)

    if verbose:
        size_mb = _MODEL_PATH.stat().st_size / (1024 * 1024)
        print(f"저장 완료: {_MODEL_PATH.name} ({size_mb:.1f} MB)")

    return model


# 하위 호환 (구 이름)
_train_and_save_model = train_and_save_model


def _ensure_model():
    """저장된 모델을 한 번만 로드한다. 스레드 안전.

    2026-09-03 수정: 파일이 없을 때 `train_and_save_model()`로 넘어가던 경로를
    제거했습니다. ① `PRIME_LAYER_C_MODEL=ames`로 강제 지정했는데 파일이 없으면
    엉뚱하게 Tox21 SR-p53을 학습해 **다른 경로**에 저장했고, ② 웹 요청 처리 중
    `urlretrieve`로 외부 데이터셋을 받아 수 분간 블로킹됐습니다.
    자동 학습은 CLI(`s6_retrain_layer_c_model.py` → `train_and_save_model`)에서만.
    """
    global _model, _model_meta, _model_ready, _model_error

    with _model_lock:
        if _model_ready:
            return _model
        if _model_error:
            return None

        model_path = _resolve_model_path()
        forced = _forced_model_choice()
        try:
            if not model_path.exists():
                if forced:
                    _model_error = (
                        f"PRIME_LAYER_C_MODEL={forced} was requested but "
                        f"{model_path.name} does not exist — refusing to train a "
                        "different endpoint during inference "
                        "(run s1_setup_layer_c_ames.py / s6_retrain_layer_c_model.py)"
                    )
                else:
                    _model_error = (
                        f"no Layer C model file ({model_path.name}) — statistical "
                        "(Q)SAR not executed; run s6_retrain_layer_c_model.py "
                        "(training is CLI-only, never during a request)"
                    )
                return None

            with open(model_path, "rb") as handle:
                payload = pickle.load(handle)
            if isinstance(payload, dict) and "model" in payload:
                _model = payload["model"]
                _model_meta = dict(payload.get("meta", {}) or {})
            else:  # 순수 estimator만 저장된 경우
                _model = payload
                _model_meta = {}
            # 메타에 엔드포인트가 없으면 파일 경로로 추정해 채운다.
            if not _model_meta.get("endpoint"):
                _model_meta["endpoint"] = _endpoint_for_path(model_path)
            _model_ready = True
            return _model
        except Exception as exc:  # noqa: BLE001
            _model_error = f"{type(exc).__name__}: {exc}"
            return None


def model_status() -> dict[str, Any]:
    """모델 상태 조회 (진단용)."""
    model = _ensure_model()
    path = _resolve_model_path()
    endpoint = active_endpoint()
    return {
        "ready": model is not None,
        "model_label": active_model_label(),
        "model_path": str(path),
        "model_exists": path.exists(),
        "trained_on_ames": model is not None and _endpoint_is_ames(endpoint),
        # 2026-09-03 추가: 학습 엔드포인트를 명시 (SR-p53은 Ames 대체 지표).
        "endpoint": endpoint,
        "endpoint_is_ames": _endpoint_is_ames(endpoint),
        "endpoint_note": (
            "" if (not endpoint or _endpoint_is_ames(endpoint)) else SURROGATE_ENDPOINT_NOTE
        ),
        "qsar_models_executed": 1 if model is not None else 0,
        "error": _model_error,
        "meta": dict(_model_meta),
    }


# --------------------------------------------------------------------------
# 예측
# --------------------------------------------------------------------------
def _positive_column_index(model: Any) -> int:
    """predict_proba 출력에서 '양성(1)' 확률이 있는 열 번호."""
    classes = getattr(model, "classes_", None)
    if classes is None:
        return 1
    classes = list(classes)
    if 1 in classes:
        return classes.index(1)
    if 1.0 in classes:
        return classes.index(1.0)
    return len(classes) - 1


def _endpoint_context() -> dict[str, Any]:
    """결과 dict에 병기할 엔드포인트 메타데이터."""
    endpoint = active_endpoint()
    is_ames = _endpoint_is_ames(endpoint)
    return {
        "Endpoint": endpoint or "none",
        "Endpoint_Is_Ames": is_ames,
        "Endpoint_Note": (
            "" if (not endpoint or is_ames) else SURROGATE_ENDPOINT_NOTE
        ),
        "Model": active_model_label(),
        "(Q)SAR Models Executed": 1 if endpoint else 0,
    }


def _empty_row(reason: str) -> dict[str, Any]:
    ctx = _endpoint_context()
    row: dict[str, Any] = {
        "Ames_Mutagenicity_Prob": None,
        "SRp53_StressResponse_Prob": None,
        "Prediction_Class": reason,
    }
    row.update(ctx)
    row["(Q)SAR Models Executed"] = 0
    return row


def _prob_row(prob: float) -> dict[str, Any]:
    """확률 하나를 엔드포인트 표기와 함께 결과 dict으로 만든다.

    2026-09-03 수정: 학습 엔드포인트가 SR-p53(진핵세포 p53 스트레스 응답)인데
    반환 키가 `Ames_Mutagenicity_Prob`(세균 변이원성) 하나뿐이라 대체 지표가
    Ames 예측으로 둔갑했습니다. 엔드포인트별 키를 분리하되, 기존 키를 읽는
    상위 모듈(p41·p50·p64)이 깨지지 않도록 값은 유지하고 엔드포인트 정보를
    병기합니다.
    """
    ctx = _endpoint_context()
    is_ames = bool(ctx["Endpoint_Is_Ames"])
    label = "Positive (위험)" if prob >= 0.5 else "Negative (안전)"
    if not is_ames:
        label = f"{label} — {SURROGATE_ENDPOINT_NOTE}"
    row: dict[str, Any] = {
        # 하위 호환 키 — 값은 실제 예측 확률, 다만 엔드포인트를 반드시 함께 읽을 것
        "Ames_Mutagenicity_Prob": prob,
        "Ames_Mutagenicity_Prob_Is_Surrogate": not is_ames,
        "SRp53_StressResponse_Prob": None if is_ames else prob,
        "Prediction_Class": label,
    }
    row.update(ctx)
    return row


def predict_advanced_properties(smiles_list: list[str]) -> dict[str, dict[str, Any]]:
    """
    SMILES 리스트를 입력받아 (Q)SAR 확률을 예측한다.

    ECFP-2048 + RandomForest 분류기를 쓰며, 학습 엔드포인트는 로드된 모델에
    따라 Ames(실측) 또는 Tox21 SR-p53(대체 지표)입니다. 반환 dict은 기존
    키(`Ames_Mutagenicity_Prob`, `Prediction_Class`)를 유지하고
    `Endpoint`·`Endpoint_Note`·`SRp53_StressResponse_Prob`를 추가합니다.
    """
    model = _ensure_model()
    if model is None:
        err = _model_error or "model load failed"
        return {smi: _empty_row(f"Unavailable ({err})") for smi in smiles_list}

    features, valid_smiles = _featurize_smiles(smiles_list)
    if not valid_smiles:
        return {smi: _empty_row("Invalid SMILES") for smi in smiles_list}

    X = np.vstack(features)
    proba = model.predict_proba(X)
    pos_idx = _positive_column_index(model)

    prob_by_smiles: dict[str, float] = {}
    for smi, row in zip(valid_smiles, proba):
        row = np.asarray(row, dtype=float).ravel()
        prob_by_smiles[smi] = float(row[pos_idx]) if row.size > pos_idx else float(row[-1])

    results: dict[str, dict[str, Any]] = {}
    for smi in smiles_list:
        prob = prob_by_smiles.get(smi)
        if prob is None:
            results[smi] = _empty_row("Invalid SMILES")
            continue
        results[smi] = _prob_row(prob)
    return results


def predict_single_smiles(smiles: str) -> dict[str, Any] | None:
    """단일 SMILES 예측 (Streamlit 연동용)."""
    row = predict_advanced_properties([smiles]).get(smiles)
    if not row or row.get("Ames_Mutagenicity_Prob") is None:
        return None
    return row


if __name__ == "__main__":
    print("모델 상태:", model_status())
    test_smiles = ["CN(C)N=O", "CC(=O)O"]
    for smi, row in predict_advanced_properties(test_smiles).items():
        print(f"  {smi:20s} -> {row}")
