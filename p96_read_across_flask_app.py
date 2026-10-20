"""
PRIME-DS Read-Across 검증 Flask API + HTML 리포트 대시보드.

Streamlit 앱(app.py)과 별도 포트(기본 5000)에서 실행합니다.
  python read_across_flask_app.py
"""

from __future__ import annotations

import os

from flask import Flask, jsonify, render_template, request
from rdkit import Chem

from p51_prime_backend_core import READ_ACROSS_MIN_SIMILARITY, ReadAcrossEngine
from p84_structure_images import smiles_to_base64_img
from p62_toxicity_alerts import scan_ashby_tennant_alerts

try:
    from p66_ich_m7 import classify_ich_m7
except ImportError:
    classify_ich_m7 = None  # type: ignore

DEFAULT_TARGET_SMILES = "Nc1ccc(Cl)cc1"
DEFAULT_TOP_N = 3
IMAGE_SIZE = (220, 220)

app = Flask(__name__, template_folder="templates")


def _structure_base64(smiles: str) -> str:
    """SMILES → PNG Base64 (2D 구조식, 220×220)."""
    if not smiles:
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return ""
        from rdkit.Chem import AllChem

        AllChem.Compute2DCoords(mol)
        return smiles_to_base64_img(Chem.MolToSmiles(mol), size=IMAGE_SIZE)
    except Exception:
        return smiles_to_base64_img(smiles, size=IMAGE_SIZE)


def _target_prediction_label(smiles: str) -> str:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return "SMILES 파싱 실패"
    try:
        alerts, _ = scan_ashby_tennant_alerts(mol)
        if classify_ich_m7 is not None:
            m7 = classify_ich_m7(mol, alerts)
            tier = m7.class_name_ko or f"Class {m7.class_number}"
            if m7.class_number in (1, 2, 3):
                return f"ICH M7 {tier} — Mutagenic concern"
            if m7.class_number == 5:
                return f"ICH M7 {tier}"
            return f"ICH M7 {tier}"
    except Exception:
        pass
    return "QSAR 예측: 검토 필요"


def _build_predict_payload(
    target_smiles: str,
    *,
    top_n: int = DEFAULT_TOP_N,
    min_similarity: float | None = None,
) -> dict:
    smi = (target_smiles or "").strip() or DEFAULT_TARGET_SMILES
    if Chem.MolFromSmiles(smi) is None:
        return {"error": f"유효하지 않은 SMILES: {smi}"}

    min_sim = (
        float(min_similarity)
        if min_similarity is not None
        else float(READ_ACROSS_MIN_SIMILARITY)
    )
    hits = ReadAcrossEngine.find_analogues(
        smi,
        top_n=max(1, min(top_n, 20)),
        min_similarity=max(0.0, min(1.0, min_sim)),
    )

    analogues = [
        {
            "name": h["name"],
            "smiles": h["smiles"],
            "ames_result": h["ames_result"],
            "similarity": h["similarity"],
            "source": h.get("source", "Hansen"),
            "img_base64": _structure_base64(h["smiles"]),
        }
        for h in hits
    ]

    return {
        "target": {
            "smiles": smi,
            "prediction": _target_prediction_label(smi),
            "img_base64": _structure_base64(smi),
        },
        "analogues": analogues,
        "min_similarity": min_sim,
    }


@app.route("/api/predict", methods=["GET"])
def predict_impurity():
    """예측 대상 + Read-Across analogues JSON."""
    target_smiles = request.args.get("smiles", DEFAULT_TARGET_SMILES)
    try:
        top_n = int(request.args.get("top_n", DEFAULT_TOP_N))
    except (TypeError, ValueError):
        top_n = DEFAULT_TOP_N

    min_raw = request.args.get("min_similarity")
    min_similarity = None
    if min_raw is not None:
        try:
            min_similarity = float(min_raw)
        except ValueError:
            return jsonify({"error": "min_similarity는 0~1 숫자여야 합니다."}), 400

    payload = _build_predict_payload(
        target_smiles,
        top_n=top_n,
        min_similarity=min_similarity,
    )
    if "error" in payload:
        return jsonify(payload), 400
    return jsonify(payload)


@app.route("/")
def index():
    return render_template("read_across_report.html")


if __name__ == "__main__":
    host = os.environ.get("PRIME_FLASK_HOST", "127.0.0.1")
    port = int(os.environ.get("PRIME_FLASK_PORT", "5000"))
    debug = os.environ.get("PRIME_FLASK_DEBUG", "0") == "1"
    print(f"PRIME-DS Read-Across API → http://{host}:{port}/")
    app.run(host=host, port=port, debug=debug)
