"""구조식 이미지 — PubChem / RDKit, Ashby-Tennant 하이라이트 지원."""

# 구조식 그림
# structure_png: 기본 그림
# rdkit_clean_structure_png: 깔끔한 버전
# rdkit_alert_highlight_png: 경보 부위를 빨갛게 표시

from __future__ import annotations

import base64
import urllib.parse
from io import BytesIO

import requests

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    from rdkit.Chem.Draw import rdMolDraw2D

    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False

from p62_toxicity_alerts import get_alert_highlight_indices, scan_report_structural_alerts

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "PRIME-DS/1.0"})

HIGHLIGHT_ALERT_RGB = (1.0, 0.15, 0.15)  # vivid red for alert zones
DEFAULT_REPORT_IMAGE_SIZE = "480x360"

__all__ = [
    "structure_png",
    "rdkit_png",
    "rdkit_clean_structure_png",
    "rdkit_alert_highlight_png",
    "load_image_bytes",
    "smiles_to_png_bytes",
    "smiles_to_base64_img",
    "smiles_to_data_url",
]


def structure_png(cid: int | str | None, smiles: str | None, size: str = "300x300") -> bytes | None:
    if cid is not None and str(cid).strip().isdigit():
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{int(cid)}/PNG"
        params = {"image_size": size}
        try:
            r = SESSION.get(url, params=params, timeout=15)
            if r.status_code == 200 and r.content[:8] == b"\x89PNG\r\n\x1a\n":
                return r.content
        except requests.RequestException:
            pass

    if smiles and isinstance(smiles, str) and smiles.strip():
        encoded = urllib.parse.quote(smiles.strip(), safe="")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/PNG"
        try:
            r = SESSION.get(url, params={"image_size": size}, timeout=20)
            if r.status_code == 200 and r.content[:8] == b"\x89PNG\r\n\x1a\n":
                return r.content
        except requests.RequestException:
            pass
    return None


def _parse_size(size: str) -> tuple[int, int]:
    try:
        w, h = size.lower().split("x", 1)
        return int(w), int(h)
    except (ValueError, AttributeError):
        return 300, 300


def _prepare_2d_mol(mol):
    """2D 좌표 생성·정리 — 결합이 끊겨 보이는 문제 방지."""
    rw = Chem.RWMol(mol)
    Chem.SanitizeMol(rw)
    mol2 = rw.GetMol()
    AllChem.Compute2DCoords(mol2)
    try:
        rdMolDraw2D.PrepareMolForDrawing(mol2, addChiralHs=False)
    except Exception:
        pass
    return mol2


def rdkit_png(
    smiles: str | None,
    size: str = "300x300",
    highlight_atoms: list[int] | None = None,
    highlight_bonds: list[int] | None = None,
    highlight_color: tuple[float, float, float] = HIGHLIGHT_ALERT_RGB,
) -> bytes | None:
    if not _HAS_RDKIT or not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    w, h = _parse_size(size)
    try:
        mol = _prepare_2d_mol(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
        opts = drawer.drawOptions()
        opts.bondLineWidth = 2.0
        opts.padding = 0.12
        opts.additionalAtomLabelPadding = 0.06
        opts.fixedBondLength = 30
        opts.fixedScale = 0.85
        opts.atomLabelDeuteriumTritium = False
        opts.highlightBondWidthMultiplier = 2.2

        atoms = highlight_atoms or []
        bonds = highlight_bonds or []
        if atoms or bonds:
            atom_colors = {a: highlight_color for a in atoms}
            bond_colors = {b: highlight_color for b in bonds}
            drawer.DrawMolecule(
                mol,
                highlightAtoms=atoms,
                highlightBonds=bonds,
                highlightAtomColors=atom_colors,
                highlightBondColors=bond_colors,
            )
        else:
            drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        try:
            kwargs: dict = {"size": (w, h)}
            if highlight_atoms:
                colors = {a: highlight_color for a in highlight_atoms}
                kwargs["highlightAtoms"] = highlight_atoms
                kwargs["highlightAtomColors"] = colors
            img = Draw.MolToImage(_prepare_2d_mol(mol), **kwargs)
            buf = BytesIO()
            img.save(buf, format="PNG")
            return buf.getvalue()
        except Exception:
            return None


def rdkit_clean_structure_png(
    smiles: str | None,
    size: str = DEFAULT_REPORT_IMAGE_SIZE,
) -> bytes | None:
    """평가 대상 구조 — 하이라이트 없이 깨끗한 2D 렌더링."""
    return rdkit_png(smiles, size=size)


def rdkit_alert_highlight_png(
    smiles: str | None,
    alert_names: list[str] | None,
    size: str = DEFAULT_REPORT_IMAGE_SIZE,
) -> bytes | None:
    """동일 구조에 alert SMARTS 매칭 영역을 빨간색으로 표시."""
    if not _HAS_RDKIT or not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    mol = _prepare_2d_mol(mol)
    highlight_atoms, highlight_bonds = get_alert_highlight_indices(mol, alert_names)
    if not highlight_atoms:
        scanned_names, highlight_atoms, highlight_bonds = scan_report_structural_alerts(mol)
        if not highlight_atoms:
            return rdkit_png(smiles, size=size)
        alert_names = scanned_names
    return rdkit_png(
        smiles,
        size=size,
        highlight_atoms=highlight_atoms,
        highlight_bonds=highlight_bonds,
        highlight_color=HIGHLIGHT_ALERT_RGB,
    )


def load_image_bytes(
    cid,
    smiles,
    size: str = "300x300",
    highlight_atoms: list[int] | None = None,
):
    # RDKit 2D 배치 우선 — 결합·원자 위치가 일관되게 표시됨
    data = rdkit_png(smiles, size=size, highlight_atoms=highlight_atoms)
    if data:
        return BytesIO(data)
    data = structure_png(cid, smiles, size=size)
    return BytesIO(data) if data else None


def smiles_to_png_bytes(
    smiles: str | None,
    *,
    size: str | tuple[int, int] = (200, 200),
) -> bytes | None:
    """SMILES → PNG bytes (2D coords 우선, MolToImage fallback)."""
    if isinstance(size, tuple):
        size_str = f"{size[0]}x{size[1]}"
    else:
        size_str = size or "200x200"

    png = rdkit_clean_structure_png(smiles, size=size_str)
    if png:
        return png
    if not _HAS_RDKIT or not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    try:
        if isinstance(size, tuple):
            wh = size
        else:
            wh = _parse_size(size_str)
        img = Draw.MolToImage(_prepare_2d_mol(mol), size=wh)
        buf = BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    except Exception:
        return None


def smiles_to_base64_img(
    smiles: str | None,
    size: str | tuple[int, int] = (200, 200),
) -> str:
    """SMILES → Base64 PNG (`data:image/png;base64,...` 바인딩용)."""
    if not _HAS_RDKIT or not smiles or not str(smiles).strip():
        return ""
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return ""
    if isinstance(size, tuple):
        wh = size
    else:
        wh = _parse_size(size or "200x200")
    try:
        # 2026-09-03 수정: 이 함수만 `_prepare_2d_mol`(= Compute2DCoords)을 거치지
        #   않고 `Draw.MolToImage` 를 직접 호출해, 같은 파일 77행 주석이 스스로
        #   "결합이 끊겨 보이는 문제 방지"라고 적은 좌표 계산을 빠뜨리고 있었습니다.
        #   다른 렌더 경로(rdkit_png · smiles_to_png_bytes)와 동일하게 맞춥니다.
        img = Draw.MolToImage(_prepare_2d_mol(mol), size=wh)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode("utf-8")
    except Exception:
        return ""


def smiles_to_data_url(
    smiles: str | None,
    size: str | tuple[int, int] = (200, 200),
) -> str:
    b64 = smiles_to_base64_img(smiles, size=size)
    return f"data:image/png;base64,{b64}" if b64 else ""
