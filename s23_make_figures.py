"""
논문 Phase 5 — Figure 2 · Figure 3 작도.

자료
----
    논문자료/02_검증결과/ad_tanimoto_validation.csv   (구조별 d_max · 라벨 · Layer C 확률)
    논문자료/02_검증결과/statistics_and_ad.json        (임계값 스윕 · 적용범위 요약)

출력 (600 dpi PNG · MDPI 기준 폭 충족)
    figures/Figure2_layerC_discrimination.png
    figures/Figure3_applicability_domain.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter
from sklearn.metrics import roc_auc_score, roc_curve

SRC = Path("/mnt/user-data/uploads/PRIME_DS(New)/논문자료/02_검증결과")
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

# 검증된 팔레트 (dataviz 기준 3슬롯 · all-pairs 통과)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d8d7d2"

plt.rcParams.update({
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.edgecolor": INK2,
    "axes.linewidth": 0.7,
    "axes.labelcolor": INK,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.frameon": False,
    "legend.fontsize": 7.5,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})


def tidy(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, color=GRID, linewidth=0.5, alpha=0.9)
    ax.set_axisbelow(True)


ad = pd.read_csv(SRC / "ad_tanimoto_validation.csv")
stats = json.load((SRC / "statistics_and_ad.json").open(encoding="utf-8"))

# ===========================================================================
#  Figure 2 — Layer C 의 판별력
# ===========================================================================
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.5))

# (a) ROC — 적용범위 안 / 밖
ins = ad[ad["in_domain"] == 1]
out = ad[ad["in_domain"] == 0]
allc = ad

for frame, colour, style, label in [
    (allc, INK2, (0, (1, 1)), "All compounds"),
    (ins, BLUE, "-", "Inside domain"),
    (out, ORANGE, (0, (4, 2)), "Outside domain"),
]:
    fpr, tpr, _ = roc_curve(frame["actual_ames"], frame["layer_c_prob"])
    auc = roc_auc_score(frame["actual_ames"], frame["layer_c_prob"])
    axA.plot(fpr, tpr, color=colour, linestyle=style, linewidth=1.6,
             label=f"{label}  (n = {len(frame):,}, AUC {auc:.3f})")

axA.plot([0, 1], [0, 1], color=GRID, linewidth=0.8)
axA.set_xlim(0, 1); axA.set_ylim(0, 1)
axA.set_xlabel("False positive rate")
axA.set_ylabel("True positive rate")
axA.set_title("(a)  Discrimination by applicability domain", fontsize=8.5, pad=6)
axA.legend(loc="lower right", handlelength=2.4, borderpad=0.4)
tidy(axA)

# (b) 임계값 스윕
sweep = pd.DataFrame(stats["threshold_sweep"]).dropna(subset=["mcc"])
LABEL_AT = {"sensitivity": 0.25, "specificity": 0.80, "mcc": 0.40}
for col, colour, style, label in [
    ("sensitivity", BLUE, "-", "Sensitivity"),
    ("specificity", ORANGE, (0, (4, 2)), "Specificity"),
    ("mcc", AQUA, (0, (1, 1.4)), "MCC"),
]:
    axB.plot(sweep["threshold"], sweep[col], color=colour, linestyle=style,
             linewidth=1.6, label=label)
    # 곡선 위 직접 라벨 — 대비가 낮은 슬롯의 구제 규칙
    xt = LABEL_AT[col]
    yt = float(sweep.loc[np.isclose(sweep["threshold"], xt), col].iloc[0])
    axB.annotate(label, (xt, yt), color=colour, fontsize=7.5, weight="bold",
                 ha="center", va="center",
                 bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none"))

esc = 0.70
axB.axvline(esc, color=INK2, linewidth=0.8, linestyle=(0, (2, 2)))
axB.annotate("escalation threshold 0.70", (esc - 0.025, 0.055), ha="right",
             va="bottom", fontsize=7, color=INK2)
best = sweep.loc[sweep["mcc"].idxmax()]
axB.plot([best["threshold"]], [best["mcc"]], "o", color=AQUA, markersize=5,
         markeredgecolor="white", markeredgewidth=1.0)
axB.annotate(f"max MCC {best['mcc']:.3f} at {best['threshold']:.2f}",
             xy=(best["threshold"], best["mcc"]), xytext=(0.60, 0.26),
             fontsize=7, color=INK2, ha="left",
             arrowprops=dict(arrowstyle="-", color=INK2, linewidth=0.6,
                             shrinkA=0, shrinkB=3))
axB.set_xlim(0, 1); axB.set_ylim(0, 1.02)
axB.set_xticks(np.arange(0, 1.01, 0.2))
axB.set_xlabel("Layer C probability threshold")
axB.set_ylabel("Value")
axB.set_title("(b)  Operating characteristics vs. threshold", fontsize=8.5, pad=6)
axB.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=3,
           handlelength=2.2, columnspacing=1.6)
tidy(axB)

fig.tight_layout(pad=0.7)
fig.savefig(OUT / "Figure2_layerC_discrimination.png", bbox_inches="tight")
plt.close(fig)
print("저장: Figure2_layerC_discrimination.png")

# ===========================================================================
#  Figure 3 — 적용범위
# ===========================================================================
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.4),
                               gridspec_kw={"width_ratios": [1.3, 1]})

# (a) d_max 분포 + 임계값
thr = stats["meta"]["ad_threshold"]
dist = stats["ad_validation"]["distribution"]
bins = np.linspace(0.1, 1.0, 37)
axA.hist(ad["max_tanimoto"], bins=bins, color=BLUE, alpha=0.85,
         edgecolor="white", linewidth=0.4)
axA.axvline(thr, color=ORANGE, linewidth=1.4, linestyle=(0, (4, 2)))
n_out = int(stats["ad_validation"]["n_out_of_domain"])
n_all = int(stats["ad_validation"]["validation_n"])
axA.annotate(f"domain threshold {thr:.2f}\n{n_out} of {n_all:,} "
             f"({100 * n_out / n_all:.1f}%) below",
             (thr + 0.015, axA.get_ylim()[1] * 0.93), fontsize=7,
             color=ORANGE, va="top", linespacing=1.3, weight="bold")
axA.axvline(dist["median"], color=INK2, linewidth=0.9)
axA.annotate(f"median {dist['median']:.3f}", (dist["median"] + 0.014,
             axA.get_ylim()[1] * 0.55), fontsize=7, color=INK2,
             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
axA.set_xlabel("$d_{max}$ — maximum Tanimoto similarity to the training set")
axA.set_ylabel("Validation compounds")
axA.set_title("(a)  Distance to training set, validation cohort", fontsize=8.5, pad=6)
tidy(axA)

# (b) 적용범위 밖 비중 — 검증셋 대 생성 구조
cases = stats["ad_case_study"]


def share(key_part: str):
    for k, v in cases.items():
        if key_part in k:
            return 100.0 * v["fraction_out_of_domain"], v["n_structures"]
    raise KeyError(key_part)


labels, values, ns = [], [], []
labels.append("Validation\ncohort")
values.append(100 * n_out / n_all)
ns.append(n_all)
for key, lab in [("니트로사민", "Nitrosamines\n(NDSRIs)"),
                 ("부형제", "Excipient\nadducts"),
                 ("이원", "Binary blend\nproducts")]:
    v, n = share(key)
    labels.append(lab); values.append(v); ns.append(n)

colours = [INK2, AQUA, BLUE, ORANGE]
bars = axB.bar(range(len(values)), values, color=colours, width=0.62,
               edgecolor="white", linewidth=1.2)
ymax = max(values) * 1.28
for b, v in zip(bars, values):
    axB.annotate(f"{v:.1f}%", (b.get_x() + b.get_width() / 2, v + ymax * 0.025),
                 ha="center", fontsize=7.8, weight="bold", color=INK)
axB.axhline(values[0], color=INK2, linewidth=0.8, linestyle=(0, (2, 2)))
axB.annotate("validation baseline", (len(values) - 0.45, values[0] + ymax * 0.012),
             ha="right", fontsize=6.8, color=INK2)
axB.set_xticks(range(len(labels)))
axB.set_xticklabels([f"{l}\nn = {n:,}" for l, n in zip(labels, ns)], fontsize=7,
                    linespacing=1.35)
axB.set_ylim(0, ymax)
axB.yaxis.set_major_formatter(PercentFormatter())
axB.set_ylabel("Share outside the applicability domain")
axB.set_title("(b)  Generated structures vs. validation cohort", fontsize=8.5, pad=6)
tidy(axB)
axB.grid(axis="x", visible=False)

fig.tight_layout(pad=0.6)
fig.savefig(OUT / "Figure3_applicability_domain.png", bbox_inches="tight")
plt.close(fig)
print("저장: Figure3_applicability_domain.png")

# 수치 확인용 출력
print()
print(f"  전체 AUC {roc_auc_score(allc['actual_ames'], allc['layer_c_prob']):.4f}"
      f"  (JSON {stats['threshold_free']['roc_auc']:.4f})")
print(f"  적용범위 안 AUC {roc_auc_score(ins['actual_ames'], ins['layer_c_prob']):.4f}"
      f"  (JSON {stats['ad_stratified_auc']['in_domain']['roc_auc']:.4f})")
print(f"  적용범위 밖 AUC {roc_auc_score(out['actual_ames'], out['layer_c_prob']):.4f}"
      f"  (JSON {stats['ad_stratified_auc']['out_of_domain']['roc_auc']:.4f})")
print(f"  최대 MCC {best['mcc']:.4f} @ {best['threshold']:.2f}"
      f"  (JSON {stats['threshold_selection']['max_mcc']['mcc']:.4f} @ "
      f"{stats['threshold_selection']['max_mcc']['threshold']:.2f})")
