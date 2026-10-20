"""
논문 Phase 5 — Scheme 2 · 3 · 4 작도.

원칙: **코드에 실제로 들어 있는 SMIRKS 문자열을 그대로 싣고, 그 문자열을
RDKit 으로 실행해 얻은 산물을 그린다.** 손으로 그린 반응식이 코드와 어긋날
여지를 없앤다.

출처
----
    p41_degradation_smirks_rules.DEGRADATION_RULES   (가수분해 · 산화)
    p30_reaction_engine.SECONDARY_AMINE_NITROSATION_SMIRKS
    p30_reaction_engine.TERTIARY_AMINE_NITROSATION_SMIRKS

출력 (600 dpi PNG)
    figures/Scheme2_hydrolysis.png
    figures/Scheme3_nitrosation.png
    figures/Scheme4_oxidation.png
"""

from __future__ import annotations

import io
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({
    "figure.dpi": 600, "savefig.dpi": 600,
    "font.family": "DejaVu Sans", "font.size": 8,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})

from p41_degradation_smirks_rules import DEGRADATION_RULES
from p30_reaction_engine import (
    SECONDARY_AMINE_NITROSATION_SMIRKS,
    TERTIARY_AMINE_NITROSATION_SMIRKS,
)

NITRITE = "[N+](=O)[O-]"


def draw_mol(smiles: str, w: int = 300, h: int = 210) -> Image.Image:
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    o = d.drawOptions()
    o.padding = 0.10
    o.bondLineWidth = 2
    o.fixedFontSize = 22
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))


def run_rule(smirks: str, reactants: list[str]) -> list[str]:
    rxn = AllChem.ReactionFromSmarts(smirks.replace(" ", ""))
    mols = tuple(Chem.MolFromSmiles(s) for s in reactants)
    out: list[str] = []
    for pset in rxn.RunReactants(mols):
        got = []
        for p in pset:
            try:
                Chem.SanitizeMol(p)
                got.append(Chem.MolToSmiles(p))
            except Exception:
                got = []
                break
        if got:
            return got
    return out


def panel(ax, tag: str, title: str, smirks: str,
          reactants: list[str], products: list[str], caption: str) -> None:
    """한 패널 = 제목 · 반응식 그림 · SMIRKS 문자열."""
    ax.axis("off")
    species = [(s, False) for s in reactants] + [(s, True) for s in products]
    n = len(species)
    # 좌표 배치 — 분자 사이에 글리프 자리를 둔다
    slot_w = 1.0 / n
    for i, (smi, is_prod) in enumerate(species):
        x0 = i * slot_w
        sub = ax.inset_axes([x0 + 0.010, 0.34, slot_w - 0.020, 0.56])
        sub.imshow(draw_mol(smi))
        sub.axis("off")
        if i < n - 1:
            glyph = "→" if i == len(reactants) - 1 else "+"
            ax.text(x0 + slot_w - 0.005, 0.62, glyph, ha="center", va="center",
                    fontsize=13, color=INK2, transform=ax.transAxes)
    ax.text(0.0, 1.0, f"({tag})  {title}", ha="left", va="top", fontsize=8.6,
            weight="bold", color=INK, transform=ax.transAxes)
    lines = textwrap.wrap(smirks.replace(" ", ""), width=104,
                          break_long_words=True, break_on_hyphens=False)
    ax.text(0.0, 0.265, "SMIRKS", ha="left", va="top", fontsize=6.6,
            color=INK2, weight="bold", transform=ax.transAxes)
    for k, ln in enumerate(lines):
        ax.text(0.0, 0.205 - 0.052 * k, ln, ha="left", va="top", fontsize=6.3,
                family="DejaVu Sans Mono", color=INK, transform=ax.transAxes)
    ax.text(0.0, 0.205 - 0.052 * len(lines) - 0.02, caption, ha="left", va="top",
            fontsize=7, color=INK2, transform=ax.transAxes, wrap=True)


def build(fname: str, panels: list[dict], height_per: float = 2.25) -> None:
    fig, axes = plt.subplots(len(panels), 1,
                             figsize=(7.4, height_per * len(panels)))
    if len(panels) == 1:
        axes = [axes]
    for ax, spec in zip(axes, panels):
        prods = run_rule(spec["smirks"], spec["reactants"])
        if not prods:
            raise RuntimeError(f"산물 없음: {spec['title']}")
        panel(ax, spec["tag"], spec["title"], spec["smirks"],
              spec["reactants"], prods, spec["caption"])
        spec["products"] = prods
    fig.tight_layout(pad=0.5, h_pad=1.2)
    fig.savefig(OUT / fname, bbox_inches="tight")
    plt.close(fig)
    print(f"저장: {fname}")
    for spec in panels:
        print(f"    ({spec['tag']}) {spec['title']}")
        print(f"        {' + '.join(spec['reactants'])}  →  {' + '.join(spec['products'])}")


# ── Scheme 2 — 가수분해 ─────────────────────────────────────────────────
build("Scheme2_hydrolysis.png", [
    dict(tag="a", title="Ester hydrolysis",
         smirks=DEGRADATION_RULES["Hydrolysis_Ester"]["smirks"],
         reactants=["CC(=O)Oc1ccccc1C(=O)O"],
         caption="Worked example: acetylsalicylic acid → salicylic acid + acetic acid."),
    dict(tag="b", title="Amide hydrolysis",
         smirks=DEGRADATION_RULES["Hydrolysis_Amide"]["smirks"],
         reactants=["CC(=O)Nc1ccc(O)cc1"],
         caption="Worked example: paracetamol → acetic acid + 4-aminophenol."),
    dict(tag="c", title="Carbamate hydrolysis",
         smirks=DEGRADATION_RULES["Hydrolysis_Carbamate"]["smirks"],
         reactants=["CNC(=O)Oc1ccccc1"],
         caption="Worked example: phenyl N-methylcarbamate → methylamine + "
                 "carbon dioxide + phenol, i.e. the carbamic acid intermediate is "
                 "taken to decarboxylate immediately."),
])

# ── Scheme 3 — 니트로소화 ───────────────────────────────────────────────
build("Scheme3_nitrosation.png", [
    dict(tag="a", title="Secondary amine nitrosation",
         smirks=SECONDARY_AMINE_NITROSATION_SMIRKS,
         reactants=["CNC", NITRITE],
         caption="Worked example: dimethylamine + nitrite → "
                 "N-nitrosodimethylamine (NDMA)."),
    dict(tag="b", title="Tertiary amine nitrosative dealkylation",
         smirks=TERTIARY_AMINE_NITROSATION_SMIRKS,
         reactants=["CN(C)Cc1ccccc1", NITRITE],
         caption="Worked example: N,N-dimethylbenzylamine + nitrite → NDMA. "
                 "The departing alkyl group is not mapped into the product "
                 "template and is therefore not emitted as a separate species."),
])

# ── Scheme 4 — 산화 ─────────────────────────────────────────────────────
build("Scheme4_oxidation.png", [
    dict(tag="a", title="Thioether oxidation to sulfoxide",
         smirks=DEGRADATION_RULES["Oxidation_Thioether"]["smirks"],
         reactants=["CSCc1ccccc1"],
         caption="Worked example: benzyl methyl sulfide → the corresponding sulfoxide."),
    dict(tag="b", title="Tertiary amine oxidation to N-oxide",
         smirks=DEGRADATION_RULES["Oxidation_Tertiary_Amine"]["smirks"],
         reactants=["CN(C)Cc1ccccc1"],
         caption="Worked example: N,N-dimethylbenzylamine → the corresponding N-oxide."),
])
