"""
논문 Phase 6 — 보충자료 워크북(Table S1–S4) 생성.

출력
    supplementary/PRIME-DS_Supplementary_Tables.xlsx
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

BASE = Path(__file__).resolve().parent
SRC_V = Path("/mnt/user-data/uploads/PRIME_DS(New)/논문자료/02_검증결과")
SRC_C = Path("/mnt/user-data/uploads/PRIME_DS(New)/논문자료/03_사례연구")
S1_MD = Path("/tmp/phase1/표S1_경보예외_추가안.md")
OUT = BASE / "supplementary"
OUT.mkdir(exist_ok=True)
XLSX = OUT / "PRIME-DS_Supplementary_Tables.xlsx"

FONT = "Arial"
HEAD_FILL = PatternFill("solid", fgColor="1B6B3A")
HEAD_FONT = Font(name=FONT, size=10, bold=True, color="FFFFFF")
BODY_FONT = Font(name=FONT, size=10)
THIN = Side(style="thin", color="D8D7D2")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

# ── Table S1 — 구조 경보 ────────────────────────────────────────────────
rows = []
for line in S1_MD.read_text(encoding="utf-8").splitlines():
    m = re.match(r"^\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*([\d,]+)\s*\|\s*([\d.]+)%\s*\|\s*(.+?)\s*\|$",
                 line)
    if m and not line.startswith("| Alert"):
        name, source, fires, prec, exc = m.groups()
        rows.append({
            "Alert": name.replace("**", ""),
            "Source": source.replace("**", ""),
            "Fires (training split)": int(fires.replace(",", "")),
            "Precision (%)": float(prec),
            "Exception": exc.replace("**", "").split(" (아래")[0],
        })
s1 = pd.DataFrame(rows)

# ── Table S2 — 반응 템플릿 ──────────────────────────────────────────────
s2 = pd.read_csv(SRC_V / "rule_inventory.csv").rename(columns={
    "module": "Module", "group": "Group", "container": "Container",
    "rule": "Rule", "status": "Status",
    "n_reactant": "Reactant templates", "n_product": "Product templates",
    "smirks": "SMIRKS",
})
s2["Status"] = s2["Status"].replace({"정상": "active"})

# ── Table S3 — 구조별 검증 결과 ─────────────────────────────────────────
pv = pd.read_csv(SRC_V / "pipeline_validation_results.csv")
ad = pd.read_csv(SRC_V / "ad_tanimoto_validation.csv")[
    ["smiles", "max_tanimoto", "in_domain"]]
s3 = pv.merge(ad, on="smiles", how="left").rename(columns={
    "smiles": "SMILES", "actual_ames": "Ames (observed)",
    "layer_a_db_hit": "Layer A record hit", "layer_b_alert": "Layer B alert",
    "layer_b_alerts": "Layer B alert names", "layer_c_prob": "Layer C p(mutagenic)",
    "rule_class": "Rule-layer class", "final_class": "Ensemble class",
    "final_class_label": "Ensemble class label", "ai_inferred": "AI-inferred",
    "pred_ensemble": "Ensemble call", "pred_layer_b": "Layer B call",
    "pred_layer_c": "Layer C call", "max_tanimoto": "d_max",
    "in_domain": "Inside applicability domain",
})

# ── Table S4 — 사례연구 생성 구조 ───────────────────────────────────────
ds = pd.read_csv(BASE / "data" / "ds_process_case_study_results.csv")
s4a = pd.DataFrame({
    "Case study": "1 — DS process (losartan route)",
    "Context": ds["rule"],
    "Role": ds["role"],
    "SMILES": ds["smiles"],
    "MW": ds["mw"],
    "ICH M7 class": ds["final_class"],
    "Layer C p(mutagenic)": ds["layer_c_prob"],
    "Layer B alerts": ds["layer_b_alerts"].fillna(""),
    "d_max": ds["ad_dmax"],
    "Documented occurrence": ds["reported_source"].fillna(""),
})

cs = pd.read_csv(SRC_C / "case_study_results.csv")
cs = cs[cs["product_smiles"].notna()]
s4b = pd.DataFrame({
    "Case study": "2 — NDSRI formation",
    "Context": cs["api_name"],
    "Role": "Generated nitrosamine",
    "SMILES": cs["product_smiles"],
    "MW": cs["product_mw"],
    "ICH M7 class": cs["final_class"].astype(str).str.replace("Class ", "", regex=False),
    "Layer C p(mutagenic)": cs["layer_c_prob"],
    "Layer B alerts": cs["layer_b_alerts"].fillna(""),
    "d_max": "",
    "Documented occurrence": "",
})

ex = pd.read_csv(SRC_C / "excipient_case_study_results.csv")
s4c = pd.DataFrame({
    "Case study": "3 — Excipient matrix",
    "Context": ex["api_name"] + " × " + ex["excipient"].astype(str),
    "Role": ex["rule"],
    "SMILES": ex["smiles"],
    "MW": ex["mw"],
    "ICH M7 class": ex["ich_m7"].astype(str).str.replace("Class ", "", regex=False),
    "Layer C p(mutagenic)": ex["prob"],
    "Layer B alerts": ex["alerts"].fillna(""),
    "d_max": "",
    "Documented occurrence": "",
})

bi = pd.read_csv(SRC_C / "binary_compatibility_pairs.csv")
bi = bi[bi["Impurity_SMILES"].notna()]
s4d = pd.DataFrame({
    "Case study": "4 — Binary 1:1 pairs",
    "Context": bi["API"] + " × " + bi["Excipient"].astype(str),
    "Role": bi["Rule"].fillna(""),
    "SMILES": bi["Impurity_SMILES"],
    "MW": "",
    "ICH M7 class": bi["ICH_M7_hint"].fillna(""),
    "Layer C p(mutagenic)": "",
    "Layer B alerts": "",
    "d_max": "",
    "Documented occurrence": "",
})
s4 = pd.concat([s4a, s4b, s4c, s4d], ignore_index=True)

# ── README ──────────────────────────────────────────────────────────────
readme = pd.DataFrame([
    ("Table S1", "Structural alerts of the rule layer",
     f"{len(s1)} rows — alert name, source rule base, firing count and precision "
     "measured on the training split (5,309 compounds), and curated exceptions. "
     "Precision is the basis on which alerts were selected and is not a performance "
     "figure."),
    ("Table S2", "Reaction template inventory",
     f"{len(s2)} rows — every SMIRKS template reachable by the generation engine, "
     "with its module, container and reactant/product template counts. Produced by "
     "s19_rule_inventory.py."),
    ("Table S3", "Per-compound validation results",
     f"{len(s3)} rows — the scaffold-separated external validation cohort with the "
     "observed Ames outcome, each layer's call, the ensemble class and the "
     "applicability-domain distance d_max."),
    ("Table S4", "Structures generated in the four case studies",
     f"{len(s4)} rows — every structure proposed by the generation engine, with its "
     "ICH M7 assignment where one was made. The binary compatibility study reports a "
     "risk band rather than an ICH M7 class."),
    ("", "", ""),
    ("Provenance", "",
     "All sheets are generated from the committed analysis outputs of the PRIME-DS "
     "repository; no value was entered by hand. Split basis: Murcko scaffold "
     "separation with InChIKey14 connectivity de-duplication. Applicability-domain "
     "threshold: d_max = 0.30 on ECFP4-2048 Tanimoto similarity."),
    ("Generated", "", "2026-09-30"),
], columns=["Sheet", "Title", "Description"])


def write_sheet(writer, name: str, frame: pd.DataFrame, widths: dict | None = None):
    frame.to_excel(writer, sheet_name=name, index=False)
    ws = writer.sheets[name]
    ws.freeze_panes = "A2"
    for c in range(1, frame.shape[1] + 1):
        cell = ws.cell(row=1, column=c)
        cell.fill = HEAD_FILL
        cell.font = HEAD_FONT
        cell.alignment = Alignment(vertical="center", wrap_text=True)
        cell.border = BORDER
    for r in range(2, frame.shape[0] + 2):
        for c in range(1, frame.shape[1] + 1):
            cl = ws.cell(row=r, column=c)
            cl.font = BODY_FONT
            cl.border = BORDER
            cl.alignment = Alignment(vertical="top", wrap_text=False)
    for i, col in enumerate(frame.columns, start=1):
        w = (widths or {}).get(col)
        if w is None:
            longest = max([len(str(col))] +
                          [len(str(v)) for v in frame[col].head(400)])
            w = min(max(longest + 2, 10), 60)
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.row_dimensions[1].height = 30


with pd.ExcelWriter(XLSX, engine="openpyxl") as xw:
    write_sheet(xw, "README", readme, {"Description": 90, "Title": 42, "Sheet": 14})
    write_sheet(xw, "Table S1", s1, {"Exception": 46, "Alert": 34})
    write_sheet(xw, "Table S2", s2, {"SMIRKS": 70, "Module": 32, "Rule": 34})
    write_sheet(xw, "Table S3", s3, {"SMILES": 52, "Layer B alert names": 52})
    write_sheet(xw, "Table S4", s4, {"SMILES": 52, "Context": 34,
                                     "Documented occurrence": 52, "Role": 34})

print(f"저장: {XLSX}")
for nm, f in (("S1", s1), ("S2", s2), ("S3", s3), ("S4", s4)):
    print(f"    Table {nm}: {len(f):,}행 × {f.shape[1]}열")
