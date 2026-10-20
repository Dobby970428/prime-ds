"""
Layer C — AttentiveFP / graph deep-learning mutagenicity inference (DeepChem).

Thread-safe singleton with automatic fallback to ECFP+RandomForest when
AttentiveFP checkpoint is absent or PyTorch/TensorFlow memory fails.

Set PRIME_GNN_BACKEND=rf|attentivefp|auto (default auto).
Build AttentiveFP checkpoint offline: python p97_build_attentivefp_model.py
"""
# 2026-09-03 수정: 위 docstring 및 아래 fallback 메시지의 구식 파일명
# (build_attentivefp_model.py)을 실제 파일명 p97_build_attentivefp_model.py 로 정정.

from __future__ import annotations

import logging
import os
import pickle
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger("prime_ds.gnn")

_MODEL_DIR = Path(__file__).resolve().parent / "models"
_ATTENTIVEFP_DIR = _MODEL_DIR / "attentivefp_ames"
_ATTENTIVEFP_META = _MODEL_DIR / "attentivefp_ames_meta.pkl"
_MUTAGEN_TASK = "SR-p53"

_model_lock = threading.RLock()
_attentivefp_model: Any = None
_attentivefp_ready = False
_attentivefp_load_attempted = False
_attentivefp_error: str | None = None
_backend_name = "none"


MODEL_UNAVAILABLE_LABEL = "unavailable"


@dataclass(frozen=True)
class GnnMutagenicityResult:
    """P(Mutagenic) from graph deep-learning layer (or RF fallback).

    2026-09-03 수정: 체크포인트도 RF pkl도 없어 `available=False, p=None`인데도
    `model_name`에 "ECFP-2048 + RandomForest (Tox21 SR-p53)" 같은 모델명이
    그대로 실려 보고서에 «GNN 포함 앙상블»로 기재됐습니다. 통계 (Q)SAR이
    0개 실행된 사실이 표기 없이 Class 확정으로 이어지지 않도록,
    사용 불가일 때는 모델명을 강제로 `unavailable`로 만들고 실행된 (Q)SAR
    개수·엔드포인트·fallback 로그를 결과 dict에 필수로 담습니다.
    """

    p_mutagenic: float | None
    model_name: str
    available: bool
    fallback_used: bool
    fallback_reason: str | None = None
    inference_ms: float = 0.0
    endpoint: str = ""
    endpoint_note: str = ""
    models_executed: int = 0

    def __post_init__(self) -> None:
        if not self.available and self.model_name != MODEL_UNAVAILABLE_LABEL:
            object.__setattr__(self, "model_name", MODEL_UNAVAILABLE_LABEL)
        if not self.available and self.models_executed:
            object.__setattr__(self, "models_executed", 0)

    def to_dict(self) -> dict[str, Any]:
        return {
            "P(Mutagenic)": self.p_mutagenic,
            "GNN Model": self.model_name,
            "GNN Available": self.available,
            "GNN Fallback": self.fallback_used,
            "GNN Fallback Reason": self.fallback_reason or "",
            # 2026-09-03 추가 — ICH M7은 상보적 (Q)SAR 2개를 요구합니다.
            # 통계 모델이 몇 개 실행됐는지, 어떤 엔드포인트인지 명시합니다.
            "(Q)SAR Models Executed": self.models_executed,
            "Layer C Endpoint": self.endpoint or "none",
            "Layer C Endpoint Note": self.endpoint_note,
            "Layer C Fallback Log": self.fallback_reason or "",
        }


def _gnn_backend_mode() -> str:
    return os.environ.get("PRIME_GNN_BACKEND", "auto").strip().lower() or "auto"


def _checkpoint_available() -> bool:
    return _ATTENTIVEFP_META.is_file() and _ATTENTIVEFP_DIR.is_dir()


def _load_attentivefp() -> tuple[Any, dict] | None:
    """Load pre-trained AttentiveFP checkpoint only — never trains during inference."""
    global _attentivefp_model, _attentivefp_ready, _attentivefp_error
    global _attentivefp_load_attempted, _backend_name

    with _model_lock:
        if _attentivefp_ready and _attentivefp_model is not None:
            with open(_ATTENTIVEFP_META, "rb") as handle:
                return _attentivefp_model, pickle.load(handle)

        if _attentivefp_load_attempted:
            return None

        _attentivefp_load_attempted = True

        if not _checkpoint_available():
            # 2026-09-03 수정: RF pkl조차 없을 수 있으므로 근거 없는 "validated"
            # 표현을 제거했습니다 (검증된 대체 모델이 있다는 뜻으로 읽힘).
            _attentivefp_error = (
                "AttentiveFP checkpoint not found — falling back to ECFP+RF "
                "statistical model if one is available "
                # 2026-09-03 수정: 런타임에 그대로 출력되는 안내 문구의 파일명 정정.
                "(run p97_build_attentivefp_model.py to enable graph model)"
            )
            return None

        try:
            import deepchem as dc

            with open(_ATTENTIVEFP_META, "rb") as handle:
                meta = pickle.load(handle)
            model = dc.models.AttentiveFPModel(
                mode="classification",
                n_tasks=int(meta["n_tasks"]),
                model_dir=str(_ATTENTIVEFP_DIR),
            )
            model.reload()
            _attentivefp_model = model
            _attentivefp_ready = True
            _backend_name = "AttentiveFP (DeepChem GraphConv / Tox21 SR-p53)"
            return model, meta
        except Exception as exc:
            _attentivefp_error = str(exc)[:240]
            logger.warning("AttentiveFP unavailable: %s", _attentivefp_error)
            return None


def _predict_attentivefp(smiles: str) -> GnnMutagenicityResult | None:
    import time

    t0 = time.perf_counter()
    loaded = _load_attentivefp()
    if loaded is None:
        return None

    model, meta = loaded
    import deepchem as dc

    featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
    feat = featurizer.featurize([smiles])[0]
    if feat is None or (hasattr(feat, "num_atoms") and feat.num_atoms == 0):
        return None

    dataset = dc.data.NumpyDataset(X=np.array([feat], dtype=object))
    with _model_lock:
        preds = model.predict(dataset)
    arr = np.asarray(preds, dtype=float)
    task_idx = int(meta["task_idx"])
    if arr.ndim == 3 and arr.shape[-1] >= 2:
        prob = float(arr[0, task_idx, 1])
    elif arr.ndim == 2 and arr.shape[-1] >= 2:
        prob = float(arr[task_idx, 1]) if arr.shape[0] > task_idx else float(arr[0, 1])
    elif arr.size >= 2:
        prob = float(arr.flat[1])
    else:
        prob = float(arr.flat[0])

    return GnnMutagenicityResult(
        p_mutagenic=max(0.0, min(1.0, prob)),
        model_name=_backend_name,
        available=True,
        fallback_used=False,
        inference_ms=(time.perf_counter() - t0) * 1000,
        endpoint=str(meta.get("endpoint") or _MUTAGEN_TASK),
        endpoint_note=(
            "" if "ames" in str(meta.get("endpoint") or "").lower()
            else "surrogate endpoint — not bacterial mutagenicity"
        ),
        models_executed=1,
    )


def _predict_rf_fallback(smiles: str, *, reason: str) -> GnnMutagenicityResult:
    import time

    from p63_deepchem_predictor import (
        active_endpoint,
        active_model_label,
        predict_advanced_properties,
    )

    t0 = time.perf_counter()
    try:
        row = predict_advanced_properties([smiles]).get(smiles, {})
        prob = row.get("Ames_Mutagenicity_Prob")
        label = active_model_label()
        endpoint = str(row.get("Endpoint") or active_endpoint() or "")
        endpoint_note = str(row.get("Endpoint_Note") or "")
    except Exception as exc:
        logger.warning("RF fallback failed: %s", exc)
        return GnnMutagenicityResult(
            p_mutagenic=None,
            # 2026-09-03 수정: 실행되지 않은 모델의 이름을 남기지 않습니다.
            model_name=MODEL_UNAVAILABLE_LABEL,
            available=False,
            fallback_used=True,
            fallback_reason=f"{reason}; RF error: {exc}"[:240],
            inference_ms=(time.perf_counter() - t0) * 1000,
            models_executed=0,
        )

    available = prob is not None
    if not available:
        # 2026-09-03 수정: 모델이 하나도 실행되지 않았음을 fallback 로그에 남깁니다.
        detail = str(row.get("Prediction_Class") or "statistical model unavailable")
        reason = f"{reason}; Layer C not executed: {detail}"[:400]
    return GnnMutagenicityResult(
        p_mutagenic=float(prob) if available else None,
        model_name=label if available else MODEL_UNAVAILABLE_LABEL,
        available=available,
        fallback_used=True,
        fallback_reason=reason,
        inference_ms=(time.perf_counter() - t0) * 1000,
        endpoint=endpoint if available else "",
        endpoint_note=endpoint_note if available else "",
        models_executed=1 if available else 0,
    )


def _last_attentivefp_error() -> str | None:
    """Module-level error string (safe read — avoids UnboundLocalError in callers)."""
    return _attentivefp_error


def predict_mutagenicity_probability(smiles: str) -> GnnMutagenicityResult:
    """
    Layer C inference — AttentiveFP when checkpoint exists, RF fallback otherwise.
    Thread-safe for multi-component process lines.
    """
    smi = (smiles or "").strip()
    if not smi:
        return GnnMutagenicityResult(
            p_mutagenic=None,
            model_name=MODEL_UNAVAILABLE_LABEL,
            available=False,
            fallback_used=True,
            fallback_reason="Empty SMILES",
            models_executed=0,
        )

    mode = _gnn_backend_mode()
    if mode == "rf":
        return _predict_rf_fallback(smi, reason="PRIME_GNN_BACKEND=rf")

    inference_error: str | None = None
    if mode in ("auto", "attentivefp"):
        try:
            result = _predict_attentivefp(smi)
            if result is not None and result.available:
                return result
        except Exception as exc:
            logger.warning("AttentiveFP inference failed: %s", exc)
            inference_error = str(exc)[:240]

    # 2026-09-03 수정: RF 대체 모델의 존재·검증 여부를 확인하지 않은 채
    # "validated"라고 적던 문구를 제거했습니다.
    reason = (
        inference_error
        or _last_attentivefp_error()
        or "AttentiveFP unavailable — attempting ECFP+RF statistical model"
    )
    return _predict_rf_fallback(smi, reason=reason)
