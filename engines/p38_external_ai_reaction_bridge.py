"""
PRIME-DS — External AI Reaction Bridge (ASKCOS · AiZynthFinder · Unified Gateway).

로컬 Docker ASKCOS/AiZynthFinder 또는 PRIME 통합 게이트웨이로
출발물질 + 시약 환경을 전송하여 generative 부반응을 예측한다.
"""

# 경로 B
# ⚠️ ASKCOS 연결

from __future__ import annotations

import os
import re
from typing import Any
from urllib.parse import urljoin

import requests
from rdkit import Chem

from engines.p03_branding import log_process_warning
from engines.p33_process_impurity_engine import mol_from_smiles

SCREENING_ENGINE_TAG = "External_AI_Reaction_Bridge"

_DEFAULT_UNIFIED_URL = "http://localhost:5000/predict/forward"
_DEFAULT_ASKCOS_URL = "http://localhost:9100"
_DEFAULT_AIZYNTH_URL = "http://localhost:5001"
_DEFAULT_TIMEOUT = 45.0

_MECHANISM_HEURISTICS: tuple[tuple[str, str], ...] = (
    (r"homo|dimer|coupling", "Homo-coupling Byproduct (AI)"),
    (r"dehal|des.?halo|des.?bromo|des.?chloro", "Des-halogenated Impurity (AI)"),
    (r"deboro|protodeboro|des.?boro", "Des-boronated Impurity (AI)"),
    (r"solvent|esterif", "Solvent-Interacted Impurity (AI)"),
    (r"alkyl|substitut", "Reagent-Alkylating Impurity (AI)"),
)


def _canonical_smiles(smiles: str) -> str:
    mol = mol_from_smiles(str(smiles or "").strip(), log_bypass=False)
    if mol is None:
        return str(smiles or "").strip()
    return Chem.MolToSmiles(mol, canonical=True)


def _env_str(name: str, default: str) -> str:
    return os.environ.get(name, default).strip()


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        return default


def _env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in ("1", "true", "yes", "on")


class ExternalAIReactionBridge:
    """
    ASKCOS / AiZynthFinder / PRIME Unified Gateway 순방향 부반응 브릿지.

    환경 변수
    ---------
    PRIME_AI_REACTION_ENABLED : auto|on|off (default auto)
    PRIME_AI_REACTION_BACKEND : auto|askcos|aizynth|unified (default auto)
    PRIME_AI_FORWARD_URL      : unified gateway (default http://localhost:5000/predict/forward)
    PRIME_ASKCOS_BASE_URL     : default http://localhost:9100
    PRIME_AIZYNTH_BASE_URL    : default http://localhost:5001
    PRIME_AI_TIMEOUT_SEC      : default 45
    """

    def __init__(
        self,
        *,
        unified_url: str | None = None,
        askcos_base: str | None = None,
        aizynth_base: str | None = None,
        backend: str | None = None,
        timeout_sec: float | None = None,
        enabled: bool | None = None,
    ) -> None:
        self.unified_url = unified_url or _env_str("PRIME_AI_FORWARD_URL", _DEFAULT_UNIFIED_URL)
        self.askcos_base = (askcos_base or _env_str("PRIME_ASKCOS_BASE_URL", _DEFAULT_ASKCOS_URL)).rstrip("/")
        self.aizynth_base = (aizynth_base or _env_str("PRIME_AIZYNTH_BASE_URL", _DEFAULT_AIZYNTH_URL)).rstrip("/")
        self.backend = (backend or _env_str("PRIME_AI_REACTION_BACKEND", "auto")).casefold()
        self.timeout_sec = timeout_sec if timeout_sec is not None else _env_float("PRIME_AI_TIMEOUT_SEC", _DEFAULT_TIMEOUT)
        enabled_mode = _env_str("PRIME_AI_REACTION_ENABLED", "auto").casefold()
        if enabled is not None:
            self.enabled = enabled
        else:
            self.enabled = enabled_mode != "off"

    def build_request_payload(
        self,
        reactants_smiles: list[str],
        target_smiles: str,
        selected_environment: dict[str, list[str]],
    ) -> dict[str, Any]:
        """PRIME 통합 요청 스키마 — 게이트웨이·ASKCOS·AiZynth 공통."""
        reactants = [_canonical_smiles(s) for s in reactants_smiles if str(s or "").strip()]
        env = {
            k: [
                str(s).strip()
                for s in (selected_environment.get(k) or [])
                if str(s or "").strip()
            ]
            for k in ("Solvents", "Acids", "Bases", "Catalysts")
        }
        reagents_flat: list[str] = []
        for cat in ("Solvents", "Acids", "Bases", "Catalysts"):
            for smi in selected_environment.get(cat) or []:
                cs = _canonical_smiles(str(smi))
                if cs:
                    reagents_flat.append(cs)
        dot_smiles = ".".join(reactants + [r for r in reagents_flat if r not in reactants])
        return {
            "reactants": reactants,
            "reagents": reagents_flat,
            "target_smiles": _canonical_smiles(target_smiles) if target_smiles else "",
            "environment": env,
            "smiles": dot_smiles,
            "include_side_products": True,
            "top_k": 25,
            "backend": self.backend,
            "model": _env_str("PRIME_ASKCOS_FORWARD_MODEL", "wldn5"),
        }

    def predict_side_products(
        self,
        reactants_smiles: list[str],
        *,
        target_smiles: str = "",
        selected_environment: dict[str, list[str]] | None = None,
    ) -> dict[str, Any]:
        """
        Generative 부반응 예측 — predicted_impurities 리스트 반환.

        Returns
        -------
        dict: predicted_impurities, status, backend_used, error (optional)
        """
        empty: dict[str, Any] = {
            "predicted_impurities": [],
            "status": "disabled",
            "backend_used": None,
            "scenario_count": 0,
        }
        if not self.enabled:
            return empty

        cleaned = [str(s).strip() for s in reactants_smiles if str(s or "").strip()]
        if not cleaned:
            return {**empty, "status": "no_reactants"}

        env = selected_environment or {}
        payload = self.build_request_payload(cleaned, target_smiles, env)
        backends = self._resolve_backend_order()

        last_error: str | None = None
        for backend_name in backends:
            try:
                raw = self._dispatch(backend_name, payload)
                preds = self._parse_response(
                    raw,
                    backend_name,
                    reactants_canon={_canonical_smiles(s) for s in cleaned},
                    target_canon=_canonical_smiles(target_smiles) if target_smiles else "",
                )
                if preds:
                    return {
                        "predicted_impurities": preds,
                        "status": "ok",
                        "backend_used": backend_name,
                        "scenario_count": len(preds),
                        "raw_response_keys": list(raw.keys()) if isinstance(raw, dict) else [],
                    }
            except Exception as exc:
                last_error = f"{backend_name}: {exc}"
                log_process_warning(f"External AI bridge skipped ({last_error})")

        return {
            **empty,
            "status": "unavailable",
            "error": last_error or "no AI backend reachable",
        }

    def is_available(self) -> bool:
        """헬스체크 — unified → askcos → aizynth 순."""
        if not self.enabled:
            return False
        for backend_name in self._resolve_backend_order():
            try:
                if backend_name == "unified":
                    r = requests.get(self.unified_url.rsplit("/", 1)[0] + "/health", timeout=3)
                    if r.status_code < 500:
                        return True
                elif backend_name == "askcos":
                    r = requests.get(f"{self.askcos_base}/docs", timeout=3)
                    if r.status_code < 500:
                        return True
                elif backend_name == "aizynth":
                    r = requests.get(f"{self.aizynth_base}/health", timeout=3)
                    if r.status_code < 500:
                        return True
            except Exception:
                continue
        return False

    def _resolve_backend_order(self) -> list[str]:
        if self.backend == "askcos":
            return ["askcos", "unified"]
        if self.backend == "aizynth":
            return ["aizynth", "unified"]
        if self.backend == "unified":
            return ["unified"]
        return ["unified", "askcos", "aizynth"]

    def _dispatch(self, backend_name: str, payload: dict[str, Any]) -> Any:
        if backend_name == "unified":
            return self._post_json(self.unified_url, payload)
        if backend_name == "askcos":
            return self._post_askcos(payload)
        if backend_name == "aizynth":
            return self._post_aizynth(payload)
        raise ValueError(f"Unknown backend: {backend_name}")

    def _post_askcos(self, payload: dict[str, Any]) -> Any:
        """ASKCOS v2 forward — call-sync-without-token 우선."""
        endpoints = (
            "/api/forward/call-sync-without-token",
            "/api/forward/call-sync",
            "/api/legacy/forward/",
        )
        body = {
            "smiles": payload.get("smiles") or ".".join(payload.get("reactants") or []),
            "reagents": payload.get("reagents") or [],
            "model": payload.get("model", "wldn5"),
            "top_k": payload.get("top_k", 25),
            "include_side_products": True,
        }
        last_exc: Exception | None = None
        for path in endpoints:
            url = urljoin(self.askcos_base + "/", path.lstrip("/"))
            try:
                return self._post_json(url, body)
            except Exception as exc:
                last_exc = exc
        raise RuntimeError(f"ASKCOS forward failed: {last_exc}")

    def _post_aizynth(self, payload: dict[str, Any]) -> Any:
        """AiZynthFinder HTTP 확장 — 표준 게이트웨이 또는 expansion API."""
        paths = (
            "/api/forward/side-products",
            "/predict/forward",
            "/api/expansion",
        )
        body = {
            "smiles": payload.get("target_smiles") or payload.get("smiles"),
            "reactants": payload.get("reactants"),
            "environment": payload.get("environment"),
            "return_n": payload.get("top_k", 25),
            "include_side_products": True,
        }
        last_exc: Exception | None = None
        for path in paths:
            url = urljoin(self.aizynth_base + "/", path.lstrip("/"))
            try:
                return self._post_json(url, body)
            except Exception as exc:
                last_exc = exc
        raise RuntimeError(f"AiZynth forward failed: {last_exc}")

    def _post_json(self, url: str, body: dict[str, Any]) -> Any:
        resp = requests.post(
            url,
            json=body,
            headers={"Accept": "application/json", "Content-Type": "application/json"},
            timeout=self.timeout_sec,
        )
        resp.raise_for_status()
        return resp.json()

    def _parse_response(
        self,
        data: Any,
        backend_name: str,
        *,
        reactants_canon: set[str],
        target_canon: str,
    ) -> list[dict[str, Any]]:
        """다양한 ASKCOS/AiZynth/Unified JSON → predicted_impurities."""
        if isinstance(data, dict) and isinstance(data.get("predicted_impurities"), list):
            return self._normalize_prediction_list(
                data["predicted_impurities"],
                backend_name,
                reactants_canon=reactants_canon,
                target_canon=target_canon,
            )
        if isinstance(data, dict) and isinstance(data.get("predictions"), list):
            return self._normalize_prediction_list(
                data["predictions"],
                backend_name,
                reactants_canon=reactants_canon,
                target_canon=target_canon,
            )

        extracted = self._extract_smiles_scores_recursive(data)
        return self._normalize_prediction_list(
            extracted,
            backend_name,
            reactants_canon=reactants_canon,
            target_canon=target_canon,
        )

    def _extract_smiles_scores_recursive(self, obj: Any, depth: int = 0) -> list[dict[str, Any]]:
        if depth > 12:
            return []
        found: list[dict[str, Any]] = []
        if isinstance(obj, dict):
            smi = (
                obj.get("smiles")
                or obj.get("product")
                or obj.get("product_smiles")
                or obj.get("outcome")
            )
            if isinstance(smi, str) and smi.strip():
                canon = _canonical_smiles(smi)
                if mol_from_smiles(canon, log_bypass=False) is not None:
                    score = self._pick_confidence(obj)
                    mech = (
                        obj.get("mechanism")
                        or obj.get("reaction_name")
                        or obj.get("reaction_class")
                        or obj.get("template")
                        or ""
                    )
                    found.append(
                        {
                            "smiles": canon,
                            "ai_confidence": score,
                            "mechanism": str(mech),
                            "reaction_name": obj.get("reaction_name"),
                        }
                    )
            for val in obj.values():
                found.extend(self._extract_smiles_scores_recursive(val, depth + 1))
        elif isinstance(obj, list):
            for item in obj:
                found.extend(self._extract_smiles_scores_recursive(item, depth + 1))
        elif isinstance(obj, str) and "." in obj and len(obj) < 500:
            for part in obj.split("."):
                part = part.strip()
                if part and mol_from_smiles(part, log_bypass=False):
                    found.append({"smiles": _canonical_smiles(part), "ai_confidence": None})
        return found

    @staticmethod
    def _pick_confidence(obj: dict[str, Any]) -> float | None:
        for key in (
            "confidence",
            "score",
            "probability",
            "normalized_score",
            "cum_prob",
            "ai_confidence",
        ):
            val = obj.get(key)
            if val is None:
                continue
            try:
                f = float(val)
                if f > 1.0 and f <= 100.0:
                    f /= 100.0
                return max(0.0, min(1.0, f))
            except (TypeError, ValueError):
                continue
        return None

    def _normalize_prediction_list(
        self,
        items: list[dict[str, Any]],
        backend_name: str,
        *,
        reactants_canon: set[str],
        target_canon: str,
    ) -> list[dict[str, Any]]:
        seen: set[tuple[str, str]] = set()
        normalized: list[dict[str, Any]] = []
        for item in items:
            smi = _canonical_smiles(str(item.get("smiles") or ""))
            if not smi or smi in reactants_canon:
                continue
            rxn_name = self._resolve_reaction_name(item, smi, target_canon)
            key = (smi, rxn_name)
            if key in seen:
                continue
            seen.add(key)
            ai_conf = item.get("ai_confidence")
            if ai_conf is None:
                ai_conf = self._pick_confidence(item)
            mech_text = str(item.get("mechanism") or item.get("mechanism_ui_text") or "")
            is_main = bool(target_canon and smi == target_canon and "main" in rxn_name.casefold())
            normalized.append(
                {
                    "smiles": smi,
                    "reaction_name": rxn_name,
                    "origin_reaction": rxn_name,
                    "mechanism": mech_text or f"{backend_name.upper()} generative forward prediction",
                    "ai_confidence": ai_conf,
                    "ai_backend": backend_name,
                    "screening_engine": SCREENING_ENGINE_TAG,
                    "impurity_role": "Main" if is_main else "Side_Impurity",
                    "impurity_classification": (
                        "Main Product (Target)" if is_main else "Side Impurity"
                    ),
                    "mechanism_ui_text": mech_text or self._mechanism_ui_for_name(rxn_name),
                    "impurity_warn": not is_main,
                    "is_main_product": is_main,
                }
            )
        return normalized

    def _resolve_reaction_name(
        self,
        item: dict[str, Any],
        smiles: str,
        target_canon: str,
    ) -> str:
        raw = str(
            item.get("reaction_name")
            or item.get("origin_reaction")
            or item.get("reaction_class")
            or ""
        ).strip()
        if raw:
            return raw if "(AI)" in raw else f"{raw} (AI)"
        if target_canon and smiles == target_canon:
            return "Cross-Coupling Product (AI)"
        mech = str(item.get("mechanism") or "").casefold()
        for pattern, label in _MECHANISM_HEURISTICS:
            if re.search(pattern, mech, re.I):
                return label
        return "AI-Generated Side Product"

    @staticmethod
    def _mechanism_ui_for_name(reaction_name: str) -> str:
        name = reaction_name.casefold()
        if "homo" in name:
            return "Homo-coupling 불순물 — AI 특허/반응 DB 기반 예측"
        if "des-halo" in name or "dehal" in name:
            return "Des-halogenated 불순물 — AI 생성 부반응"
        if "des-boro" in name or "deboro" in name:
            return "Des-boronated 불순물 — AI 생성 부반응"
        if "solvent" in name:
            return "Solvent 부반응 산물 — AI 생성"
        return "AI 생성 불순물 — ASKCOS/AiZynth generative forward"


_BRIDGE: ExternalAIReactionBridge | None = None


def get_external_ai_reaction_bridge() -> ExternalAIReactionBridge:
    global _BRIDGE
    if _BRIDGE is None:
        _BRIDGE = ExternalAIReactionBridge()
    return _BRIDGE
