"""
논문 Phase 5 — Figure 1 (플랫폼 아키텍처) · Scheme 1 (ICH M7 결정 트리) 작도.

Scheme 1 의 분기는 p66_ich_m7.route_m7_consensus_class_detailed() 의 실제
분기 순서를 그대로 옮긴 것이다.

출력 (600 dpi PNG)
    figures/Figure1_architecture.png
    figures/Scheme1_ich_m7_decision_tree.png
"""

from __future__ import annotations

from pathlib import Path

import cairosvg

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

INK = "#0b0b0b"
INK2 = "#52514e"
LINE = "#8a8984"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
FILL_N = "#f4f4f2"
FONT = "DejaVu Sans, Helvetica, Arial, sans-serif"
MONO = "DejaVu Sans Mono, monospace"


def _wrap(text: str, width: int) -> list[str]:
    words, lines, cur = text.split(), [], ""
    for w_ in words:
        trial = f"{cur} {w_}".strip()
        if len(trial) <= width:
            cur = trial
        else:
            lines.append(cur); cur = w_
    if cur:
        lines.append(cur)
    return lines


def esc(t: str) -> str:
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def box(x, y, w, h, title, lines=(), accent=INK2, fill="#ffffff", r=5):
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" '
           f'fill="{fill}" stroke="{accent}" stroke-width="1.4"/>']
    out.append(f'<rect x="{x}" y="{y}" width="4.5" height="{h}" rx="2" fill="{accent}"/>')
    ty = y + 17
    out.append(f'<text x="{x + 14}" y="{ty}" font-family="{FONT}" font-size="11.5" '
               f'font-weight="700" fill="{INK}">{esc(title)}</text>')
    for i, ln in enumerate(lines):
        out.append(f'<text x="{x + 14}" y="{ty + 15 + i * 13}" font-family="{FONT}" '
                   f'font-size="9.6" fill="{INK2}">{esc(ln)}</text>')
    return "".join(out)


def arrow(x1, y1, x2, y2, label="", dash=False, colour=LINE):
    d = ' stroke-dasharray="5 3"' if dash else ""
    s = (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{colour}" '
         f'stroke-width="1.4" marker-end="url(#ar)"{d}/>')
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        s += (f'<text x="{mx + 7}" y="{my - 3}" font-family="{FONT}" font-size="9" '
              f'fill="{INK2}">{esc(label)}</text>')
    return s


def svg(w, h, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}"><defs>'
            f'<marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
            f'markerHeight="6" orient="auto-start-reverse">'
            f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{LINE}"/></marker></defs>'
            f'<rect width="{w}" height="{h}" fill="white"/>{body}</svg>')


def write(name: str, markup: str, width_px: int):
    p = OUT / name
    cairosvg.svg2png(bytestring=markup.encode(), write_to=str(p),
                     output_width=width_px, background_color="white")
    print(f"저장: {name}")


# ===========================================================================
#  Figure 1 — 플랫폼 아키텍처
# ===========================================================================
W, H = 980, 610
b = []
COL = 300
X0, XW = 40, 900

b.append(f'<text x="{X0}" y="26" font-family="{FONT}" font-size="13" '
         f'font-weight="700" fill="{INK}">PRIME-DS — platform architecture</text>')

# 1 입력
b.append(box(X0, 44, XW, 54, "1 · Input (any structure in the synthetic route)",
             ["Drug substance · starting material · intermediate · reagent · "
              "formulation context (excipient selection, nitrite source)"],
             accent=INK2, fill=FILL_N))
b.append(arrow(X0 + XW / 2, 98, X0 + XW / 2, 118))

# 2 생성
b.append(box(X0, 120, 430, 96, "2 · Generation engine",
             ["112 SMIRKS reaction templates",
              "Hydrolysis · oxidation · nitrosation cascade",
              "Named reactions and a mechanism-specific",
              "cross-coupling rule set (Schemes 2–4)"], accent=BLUE))
b.append(box(X0 + 470, 120, 430, 96, "2′ · Formulation channel",
             ["143-entry excipient library",
              "14 interaction rules (Maillard, formylation,",
              "peroxide oxidation, ester/amide condensation)",
              "Whole-matrix and binary 1:1 modes"], accent=BLUE))
b.append(arrow(X0 + 215, 216, X0 + XW / 2, 240))
b.append(arrow(X0 + 685, 216, X0 + XW / 2, 240))

# 3 물성 필터
b.append(box(X0, 242, XW, 56, "3 · Physicochemical screening and plausibility filters",
             ["Molecular-weight bound · suppression of radical species and of Phase I "
              "metabolic products (pKa, LogD and boiling point are shown as "
              "indicative values only)"], accent=AQUA))
b.append(arrow(X0 + XW / 2, 298, X0 + XW / 2, 318))

# 4 4계층
b.append(f'<rect x="{X0}" y="320" width="{XW}" height="150" rx="6" fill="#fbfbfa" '
         f'stroke="{INK2}" stroke-width="1.4"/>')
b.append(f'<text x="{X0 + 14}" y="340" font-family="{FONT}" font-size="11.5" '
         f'font-weight="700" fill="{INK}">4 · Four-layer ICH M7 assessment</text>')
lay = [
    ("Layer A", "Experimental record index", "Ames / carcinogenicity — Class 1, 2"),
    ("Layer B", "23 alerts + 4 fragments", "21 SMARTS + 2 algorithmic PAH + 4 metabolic; exceptions"),
    ("Layer C", "ECFP4-2048 random forest", "n = 5,309 training compounds; p(mutagenic)"),
    ("Layer D", "Read-across arbitration", "invoked when B and C disagree"),
]
for i, (tag, name, sub) in enumerate(lay):
    x = X0 + 14 + i * 222
    b.append(f'<rect x="{x}" y="350" width="206" height="106" rx="4" fill="white" '
             f'stroke="{LINE}" stroke-width="1"/>')
    b.append(f'<rect x="{x}" y="350" width="206" height="4" fill="{ORANGE}"/>')
    b.append(f'<text x="{x + 11}" y="374" font-family="{FONT}" font-size="10.5" '
             f'font-weight="700" fill="{INK}">{tag}</text>')
    for j, ln in enumerate(_wrap(name, 26)):
        b.append(f'<text x="{x + 11}" y="{390 + j * 12}" font-family="{FONT}" '
                 f'font-size="9.4" fill="{INK}">{esc(ln)}</text>')
    for j, ln in enumerate(_wrap(sub, 30)):
        b.append(f'<text x="{x + 11}" y="{418 + j * 11.5}" font-family="{FONT}" '
                 f'font-size="8.6" fill="{INK2}">{esc(ln)}</text>')
b.append(arrow(X0 + XW / 2, 470, X0 + XW / 2, 492))

# 5 합의
b.append(box(X0, 494, 560, 56, "5 · Upgrade-only consensus (ICH M7)",
             ["A statistical or read-across call may raise a class, never lower one "
              "assigned by an alert"], accent=ORANGE))
b.append(box(X0 + 600, 494, 300, 56, "6 · Output",
             ["ICH M7 Class 1–5 with d_max applicability-domain",
              "similarity on every prediction; TXT / DOCX / PDF dossier"],
             accent=INK2, fill=FILL_N))
b.append(arrow(X0 + 560, 522, X0 + 598, 522))

b.append(f'<text x="{X0}" y="578" font-family="{FONT}" font-size="9" fill="{INK2}">'
         f'Numbers shown are the counts verified by the rule inventory '
         f'(s19_rule_inventory.py) and the Layer C model metadata.</text>')

# ===========================================================================
#  Scheme 1 — ICH M7 결정 트리
# ===========================================================================
SW, SH = 980, 636
c = []
c.append(f'<text x="40" y="26" font-family="{FONT}" font-size="13" font-weight="700" '
         f'fill="{INK}">ICH M7 class assignment — decision order as implemented</text>')
c.append(f'<text x="40" y="43" font-family="{FONT}" font-size="9" fill="{INK2}">'
         f'Branches follow route_m7_consensus_class_detailed() in p66_ich_m7.py, '
         f'evaluated top to bottom; the first matching branch returns.</text>')

STEPS = [
    ("Empirical record: carcinogenicity positive AND Ames positive",
     "Class 1", "verified mutagenic carcinogen", ORANGE),
    ("Empirical record: carcinogenicity positive AND Ames negative",
     "Class 5 + expert review", "non-DNA-reactive mechanism, outside M7 scope", AQUA),
    ("Empirical record: carcinogenicity positive, Ames data absent",
     "Class 3 + expert review", "mutagenic mechanism unconfirmed; Ames testing required", ORANGE),
    ("Empirical record: Ames positive",
     "Class 2", "known mutagen, carcinogenic potential unknown", ORANGE),
    ("Structural alert shared with the parent API, parent Ames negative",
     "Class 4", "alerting but same alert as the drug substance; Q3A/Q3B limits", BLUE),
    ("Structural alert present, compound in record index with Ames negative",
     "Class 4", "alerting structure with negative experimental data", BLUE),
    ("Structural alert present, no experimental record",
     "Class 3 (preliminary)", "alerting structure of unknown potential", ORANGE),
    ("No structural alert",
     "Class 5 (preliminary)",
     "no alert identified; final unless a negative record exists", AQUA),
]

y = 62
for i, (cond, cls, note, accent) in enumerate(STEPS):
    c.append(f'<rect x="40" y="{y}" width="560" height="46" rx="4" fill="{FILL_N}" '
             f'stroke="{LINE}" stroke-width="1"/>')
    c.append(f'<text x="54" y="{y + 19}" font-family="{FONT}" font-size="9.2" '
             f'fill="{INK2}" font-weight="700">{i + 1}</text>')
    for j, ln in enumerate(_wrap(cond, 72)):
        c.append(f'<text x="70" y="{y + 19 + j * 12}" font-family="{FONT}" '
                 f'font-size="9.6" fill="{INK}">{esc(ln)}</text>')
    c.append(f'<line x1="600" y1="{y + 23}" x2="632" y2="{y + 23}" stroke="{LINE}" '
             f'stroke-width="1.4" marker-end="url(#ar)"/>')
    c.append(f'<rect x="634" y="{y}" width="306" height="46" rx="4" fill="white" '
             f'stroke="{accent}" stroke-width="1.4"/>')
    c.append(f'<rect x="634" y="{y}" width="4.5" height="46" rx="2" fill="{accent}"/>')
    for j, ln in enumerate(_wrap(cls, 40)):
        c.append(f'<text x="648" y="{y + 17 + j * 12}" font-family="{FONT}" '
                 f'font-size="10" font-weight="700" fill="{INK}">{esc(ln)}</text>')
    c.append(f'<text x="648" y="{y + 38}" font-family="{FONT}" font-size="8.6" '
             f'fill="{INK2}">{esc(note[:64])}</text>')
    if i < len(STEPS) - 1:
        c.append(f'<line x1="52" y1="{y + 46}" x2="52" y2="{y + 60}" stroke="{LINE}" '
                 f'stroke-width="1.2" stroke-dasharray="3 3"/>')
        c.append(f'<text x="60" y="{y + 57}" font-family="{FONT}" font-size="8" '
                 f'fill="{INK2}">no</text>')
    y += 60

c.append(f'<rect x="40" y="{y + 6}" width="900" height="52" rx="5" fill="white" '
         f'stroke="{ORANGE}" stroke-width="1.4"/>')
c.append(f'<rect x="40" y="{y + 6}" width="4.5" height="52" rx="2" fill="{ORANGE}"/>')
c.append(f'<text x="54" y="{y + 25}" font-family="{FONT}" font-size="10.5" '
         f'font-weight="700" fill="{INK}">Consensus overlay — upgrade only</text>')
c.append(f'<text x="54" y="{y + 42}" font-family="{FONT}" font-size="9.4" fill="{INK2}">'
         f'A Layer C probability above 0.70, or a Layer D read-across call, may raise '
         f'the class assigned above; neither may lower it.</text>')


write("Figure1_architecture.png", svg(W, H, "".join(b)), 4200)
write("Scheme1_ich_m7_decision_tree.png", svg(SW, SH, "".join(c)), 4200)
