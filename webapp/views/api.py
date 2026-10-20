"""AJAX · 이미지 엔드포인트."""

from __future__ import annotations

import hashlib
import io
import threading

from flask import Blueprint, Response, jsonify, request, send_file

from .. import state

bp = Blueprint("api", __name__)

# 구조 이미지 캐시 — 같은 SMILES·크기 조합은 다시 그리지 않습니다.
_IMG_LOCK = threading.Lock()
_IMG_CACHE: dict[str, bytes] = {}
_IMG_CACHE_MAX = 600

_PLACEHOLDER = (
    b'<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d">'
    b'<rect width="100%%" height="100%%" fill="#f8fafc" stroke="#cbd5e1" '
    b'stroke-dasharray="4 3"/><text x="50%%" y="50%%" font-family="sans-serif" '
    b'font-size="12" fill="#94a3b8" text-anchor="middle" dominant-baseline="middle">'
    b"2D \xea\xb5\xac\xec\xa1\xb0 \xec\x97\x86\xec\x9d\x8c</text></svg>"
)


def _parse_size(raw: str) -> tuple[int, int]:
    try:
        w, h = raw.lower().split("x", 1)
        return max(40, min(1200, int(w))), max(40, min(1200, int(h)))
    except Exception:
        return 300, 220


@bp.get("/structure.png")
def structure_png() -> Response:
    """SMILES를 2D 구조 PNG로. 없거나 실패하면 SVG 자리표시자."""
    smiles = (request.args.get("smiles") or "").strip()
    size = (request.args.get("size") or "300x220").strip()
    highlight = (request.args.get("hl") or "").strip()
    width, height = _parse_size(size)

    if not smiles:
        svg = _PLACEHOLDER % (width, height)
        return Response(svg, mimetype="image/svg+xml")

    key = hashlib.md5(f"{smiles}|{width}x{height}|{highlight}".encode()).hexdigest()
    with _IMG_LOCK:
        hit = _IMG_CACHE.get(key)
    if hit:
        return send_file(io.BytesIO(hit), mimetype="image/png")

    data: bytes | None = None
    try:
        if highlight:
            from p84_structure_images import rdkit_alert_highlight_png

            atoms = [int(a) for a in highlight.split(",") if a.strip().isdigit()]
            data = rdkit_alert_highlight_png(smiles, atoms, f"{width}x{height}")
        else:
            from p84_structure_images import rdkit_clean_structure_png

            data = rdkit_clean_structure_png(smiles, f"{width}x{height}")
    except Exception:
        data = None

    if not data:
        try:
            from p84_structure_images import rdkit_png

            data = rdkit_png(smiles, f"{width}x{height}")
        except Exception:
            data = None

    if not data:
        svg = _PLACEHOLDER % (width, height)
        return Response(svg, mimetype="image/svg+xml")

    with _IMG_LOCK:
        if len(_IMG_CACHE) > _IMG_CACHE_MAX:
            _IMG_CACHE.clear()
        _IMG_CACHE[key] = data

    resp = send_file(io.BytesIO(data), mimetype="image/png")
    resp.headers["Cache-Control"] = "public, max-age=86400"
    return resp


@bp.post("/canvas")
def save_canvas() -> Response:
    """Ketcher 캔버스 SMILES를 서버 세션에 저장합니다."""
    payload = request.get_json(silent=True) or {}
    name = (payload.get("name") or "").strip()
    smiles = (payload.get("smiles") or "").strip()
    if not name:
        return jsonify({"ok": False, "error": "name 누락"}), 400
    state.set_canvas(name, smiles)
    return jsonify({"ok": True, "name": name, "smiles": smiles})


@bp.get("/canvas")
def read_canvas() -> Response:
    return jsonify({"ok": True, "canvas": dict(state.canvas())})


@bp.post("/validate")
def validate_smiles() -> Response:
    """SMILES 유효성 + 분자식·분자량을 즉시 돌려줍니다."""
    payload = request.get_json(silent=True) or {}
    smiles = (payload.get("smiles") or "").strip()
    if not smiles:
        return jsonify({"ok": False, "error": "빈 구조"})
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"ok": False, "error": "RDKit이 해석하지 못한 구조입니다"})
        return jsonify(
            {
                "ok": True,
                "canonical": Chem.MolToSmiles(mol),
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "mw": round(Descriptors.ExactMolWt(mol), 4),
                "atoms": mol.GetNumAtoms(),
            }
        )
    except Exception as exc:  # pragma: no cover
        return jsonify({"ok": False, "error": str(exc)})


@bp.post("/registry/search")
def registry_search() -> Response:
    """공정 시약 동적 등록 — PubChem 조회 후 카테고리 자동 분류.

    Streamlit 판에서는 렌더 함수 안에서 네트워크 조회를 하고 st.rerun()을
    불렀습니다. 여기서는 별도 엔드포인트로 분리했습니다.
    """
    payload = request.get_json(silent=True) or {}
    query = (payload.get("query") or "").strip()
    prefix = (payload.get("prefix") or "process").strip()
    if not query:
        return jsonify({"ok": False, "message": "검색어를 입력하세요"})

    try:
        from engines.p14_dynamic_registry import get_dynamic_registry

        result = get_dynamic_registry().get_unlimited_entry(query) or {}
    except Exception as exc:
        return jsonify({"ok": False, "message": f"조회 실패: {exc}"})

    if not result.get("success"):
        return jsonify({"ok": False, "message": result.get("message") or "찾지 못했습니다"})

    entry = {
        "id": result.get("id"),
        "name": result.get("name"),
        "smiles": result.get("smiles"),
        "category": result.get("category"),
    }
    catalog = state.get(f"{prefix}_dynamic_catalog", [])
    if not any(e.get("id") == entry["id"] for e in catalog):
        catalog = list(catalog) + [entry]
        state.put(f"{prefix}_dynamic_catalog", catalog)

    selected = list(state.get(f"{prefix}_registry_{_cat_slug(entry['category'])}", []))
    if entry["id"] and entry["id"] not in selected:
        selected.append(entry["id"])
        state.put(f"{prefix}_registry_{_cat_slug(entry['category'])}", selected)

    return jsonify({"ok": True, "entry": entry, "catalog": catalog})


def _cat_slug(category: str | None) -> str:
    return (category or "").strip().lower() or "solvents"


@bp.post("/registry/remove")
def registry_remove() -> Response:
    payload = request.get_json(silent=True) or {}
    entry_id = (payload.get("id") or "").strip()
    prefix = (payload.get("prefix") or "process").strip()
    catalog = [e for e in state.get(f"{prefix}_dynamic_catalog", []) if e.get("id") != entry_id]
    state.put(f"{prefix}_dynamic_catalog", catalog)
    for cat in ("solvents", "acids", "bases", "catalysts"):
        key = f"{prefix}_registry_{cat}"
        sel = [i for i in state.get(key, []) if i != entry_id]
        state.put(key, sel)
    return jsonify({"ok": True, "catalog": catalog})
