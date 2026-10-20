"""Figure 2 - sex bias in the magnitude of differential expression."""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from adjustText import adjust_text

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *  # noqa: E402,F401
from figdata import *  # noqa: E402,F401

apply_style()
plt.rcParams["mathtext.default"] = "regular"
M = lambda v, f="+.3f": format(v, f).replace("-", "−")
os.makedirs(SRC, exist_ok=True)

FW, FH = W2, 9.0
fig = plt.figure(figsize=(FW, FH))


def ax_in(x, y, w, h, **kw):
    return fig.add_axes([x / FW, y / FH, w / FW, h / FH], **kw)


def plabel(letter, x, y):
    fig.text(x / FW, y / FH, letter, fontsize=11, fontweight="bold", va="top", ha="left", color=INK)


def fmt_p(p):
    if p >= 0.01:
        return f"{p:.3f}" if p < 0.1 else f"{p:.2f}"
    m, e = f"{p:.2e}".split("e")
    return f"{m} × 10$^{{{int(e)}}}$"


def ref_lines(ax, lo, hi, zero=True):
    if zero:
        ax.axhline(0, color=GRID, lw=0.6, zorder=0)
        ax.axvline(0, color=GRID, lw=0.6, zorder=0)
    ax.plot([lo, hi], [lo, hi], color=INK2, lw=0.6, zorder=1)


hpc, pfc, pfcu = load_hpc(), load_pfc(), load_pfc_unfiltered()
core = core109()
stats_rows = []

# =============================================================================
# A  male vs female log2FC of sex-shared genes, per region
# =============================================================================
plabel("A", 0.0, FH - 0.02)
a_src = []
for k, (reg, d, lim) in enumerate((("HPC", hpc, (-3.5, 11.2)), ("PFC", pfc, (-4.0, 7.5)))):
    s = d[d.category == "sex-shared"]
    ax = ax_in(0.42, FH - 0.32 - 1.45 - k * 1.95, 1.45, 1.45)
    ref_lines(ax, *lim)
    ax.scatter(s.lfc_m, s.lfc_f, s=5, color=SHARED, alpha=0.45, linewidths=0, zorder=2)
    rho, prho = stats.spearmanr(s.lfc_m, s.lfc_f)
    r, pr = stats.pearsonr(s.lfc_m, s.lfc_f)
    stats_rows.append(dict(panel="2A", region=reg, n=len(s), spearman_rho=rho, spearman_p=prho, pearson_r=r, pearson_p=pr))
    ax.text(0.03, 0.97, f"{reg}, n = {len(s):,}\nSpearman ρ = {rho:.3f}\nPearson r = {r:.3f}",
            transform=ax.transAxes, fontsize=6, va="top", ha="left", color=INK, linespacing=1.25)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("log₂FC, men", fontsize=7, labelpad=1)
    ax.set_ylabel("log₂FC, women", fontsize=7, labelpad=1)
    ax.tick_params(labelsize=6)
    a_src.append(s.assign(region=reg)[["region", "gene", "lfc_m", "lfc_f", "padj_m", "padj_f"]])
pd.concat(a_src).to_csv(f"{SRC}/Figure2A.csv", index=False)

# =============================================================================
# B  |log2FC| women vs men for the 109 core genes, per region (log-log)
# =============================================================================
plabel("B", 2.08, FH - 0.02)
mk = core.state_marker != ""
fb = core.feedback_regulator
core["marker_class"] = np.select([mk & fb, mk, fb], ["state marker + feedback regulator", "state marker",
                                                     "feedback regulator"], "other")
SHAPES = {"other": ("o", 9), "state marker": ("^", 18), "feedback regulator": ("D", 13),
          "state marker + feedback regulator": ("*", 34)}
b_src = []
BLIM = (0.07, 14)
BXLIM = (0.07, 40)
for k, reg in enumerate(("HPC", "PFC")):
    ax = ax_in(2.36 + k * 2.45, FH - 0.45 - 3.0, 2.35, 3.0)
    xm, xf = core[f"{reg}_men"].abs(), core[f"{reg}_women"].abs()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(BXLIM, BXLIM, color=INK2, lw=0.6, zorder=1)
    hi = core.stratum == "high"
    for cls, (m, sz) in SHAPES.items():
        sel = core.marker_class == cls
        for strat, col in (("low", "#bdbbb6"), ("high", INK2 if cls == "other" else INK)):
            ss = sel & (core.stratum == strat)
            if ss.any():
                ax.scatter(xm[ss], xf[ss], marker=m, s=sz, color=col, linewidths=0.4, edgecolors="white", zorder=3)
    ax.set_xlim(BXLIM)
    ax.set_ylim(BLIM)
    # deterministic callout columns: F>=M genes left of the cloud, M>F genes right of it
    ax_h_in = 3.0
    pitch = 0.092 * np.log10(BLIM[1] / BLIM[0]) / ax_h_in      # label pitch in decades
    for side, sel, xa, ylo, yhi in (("L", hi & (xf >= xm), 1.08, 0.72, 13.0),
                                    ("R", hi & (xf < xm), 8.2, 0.085, 7.5)):
        idx = core.index[sel]
        order = sorted(idx, key=lambda i: -xf[i])
        want = [np.log10(xf[i]) for i in order]
        top, bot = np.log10(yhi), np.log10(ylo)
        pos = []
        for wv in want:
            v = min(wv, top) if not pos else min(wv, pos[-1] - pitch)
            pos.append(v)
        for j in range(len(pos) - 1, -1, -1):
            lim = bot if j == len(pos) - 1 else pos[j + 1] + pitch
            pos[j] = max(pos[j], lim)
        assert pos[0] <= top + 1e-9, (reg, side)
        for i, py in zip(order, pos):
            ly = 10 ** py
            ax.text(xa, ly, core.gene[i], fontsize=6, style="italic", color=INK, zorder=4,
                    ha="right" if side == "L" else "left", va="center")
            x0 = xa * (1.06 if side == "L" else 1 / 1.06)
            ax.plot([x0, xm[i]], [ly, xf[i]], color=MUTED, lw=0.4, zorder=1.5)
    d = (xf - xm)[hi]
    w = stats.wilcoxon(xf[hi], xm[hi])
    stats_rows.append(dict(panel="2B", region=reg, n=int(hi.sum()), median_diff_F_minus_M=np.median(d),
                           n_F_gt_M=int((d > 0).sum()), n_M_gt_F=int((d < 0).sum()), wilcoxon_p=w.pvalue))
    ax.text(0.03, 0.97, f"{reg}", transform=ax.transAxes, fontsize=7, fontweight="bold", va="top", color=INK)
    ax.set_title(f"High-effect genes (n = 38): women > men {int((d > 0).sum())}, men > women {int((d < 0).sum())}\n"
                 f"median Δ = {M(np.median(d))}, Wilcoxon P = {fmt_p(w.pvalue)}",
                 fontsize=6, fontweight="normal", color=INK2, loc="center", pad=3, linespacing=1.25)
    ticks = [0.1, 0.5, 1.5, 5, 10]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks], fontsize=6)
    ax.set_yticklabels([f"{t:g}" for t in ticks] if k == 0 else [], fontsize=6)
    ax.minorticks_off()
    ax.set_xlabel("|log₂FC|, men", fontsize=7, labelpad=1)
    ax.set_ylabel("|log₂FC|, women" if k == 0 else "", fontsize=7, labelpad=1)
    b_src.append(pd.DataFrame({"region": reg, "gene": core.gene, "abs_lfc_men": xm, "abs_lfc_women": xf,
                               "stratum": core.stratum, "marker_class": core.marker_class,
                               "labeled": hi}))
pd.concat(b_src).to_csv(f"{SRC}/Figure2B.csv", index=False)
handles = [Line2D([], [], marker=SHAPES[c][0], ls="none", color=INK, ms=np.sqrt(SHAPES[c][1]) * 1.05,
                  label=lab) for c, lab in [("state marker", "State marker"), ("feedback regulator", "Feedback regulator"),
                                            ("state marker + feedback regulator", "Both")]]
handles += [Line2D([], [], marker="o", ls="none", color=INK2, ms=3, label="Other, high effect"),
            Line2D([], [], marker="o", ls="none", color="#bdbbb6", ms=3, label="Other, low effect")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(2.4 / FW, (FH - 3.82) / FH), ncol=5, fontsize=6,
           handletextpad=0.2, columnspacing=1.0, borderaxespad=0, labelspacing=0.4)

# =============================================================================
# C  sex bias of the 38 high-effect core genes, HPC vs PFC (dumbbell)
# =============================================================================
hc = core[core.stratum == "high"].copy()
hc["bias_HPC"] = hc.HPC_women.abs() - hc.HPC_men.abs()
hc["bias_PFC"] = hc.PFC_women.abs() - hc.PFC_men.abs()
hc["switch_F_to_M"] = (hc.bias_HPC > 0) & (hc.bias_PFC < 0)
hc = hc.sort_values("bias_HPC", ascending=False).reset_index(drop=True)
wc = stats.wilcoxon(hc.bias_HPC, hc.bias_PFC)
stats_rows.append(dict(panel="2C", region="HPC vs PFC", n=len(hc), median_bias_HPC=hc.bias_HPC.median(),
                       median_bias_PFC=hc.bias_PFC.median(), n_switch=int(hc.switch_F_to_M.sum()),
                       wilcoxon_p=wc.pvalue))
hc[["gene", "HPC_men", "HPC_women", "PFC_men", "PFC_women", "bias_HPC", "bias_PFC", "switch_F_to_M"]].to_csv(
    f"{SRC}/Figure2C.csv", index=False)

plabel("C", 0.0, 5.0)
axc = ax_in(0.47, 3.15, 6.68, 1.62)
x = np.arange(len(hc))
axc.axhline(0, color=INK2, lw=0.6, zorder=1)
for i, r in hc.iterrows():
    sw = r.switch_F_to_M
    axc.plot([i, i], [r.bias_HPC, r.bias_PFC], color=MEN if sw else "#c9c7c2", lw=1.6 if sw else 1.0,
             zorder=2, solid_capstyle="butt")
axc.scatter(x, hc.bias_HPC, s=16, facecolors="white", edgecolors=INK, linewidths=0.8, zorder=3, label="HPC")
axc.scatter(x, hc.bias_PFC, s=16, color=INK, linewidths=0, zorder=3, label="PFC")
axc.set_xlim(-0.7, len(hc) - 0.3)
axc.set_ylim(-2.6, 3.8)
axc.set_xticks(x)
axc.set_xticklabels(hc.gene, rotation=90, fontsize=6, style="italic")
for lab, sw in zip(axc.get_xticklabels(), hc.switch_F_to_M):
    if sw:
        lab.set_fontweight("bold")
        lab.set_color(INK)
    else:
        lab.set_color(INK2)
axc.tick_params(axis="x", length=0, pad=2)
axc.set_ylabel("|log₂FC_F| − |log₂FC_M|", fontsize=7)
axc.tick_params(axis="y", labelsize=6)
axc.text(len(hc) - 0.5, 3.75, "↑ larger in women", fontsize=6, color=INK2, ha="right", va="top")
axc.text(len(hc) - 0.5, -2.55, "↓ larger in men", fontsize=6, color=INK2, ha="right", va="bottom")
axc.text(len(hc) - 0.5, 2.95,
         f"Paired Wilcoxon P = {fmt_p(wc.pvalue)}\nmedian HPC {M(hc.bias_HPC.median())}, "
         f"PFC {M(hc.bias_PFC.median())}",
         fontsize=6, color=INK, ha="right", va="top", linespacing=1.25)
hl = [Line2D([], [], marker="o", ls="none", mfc="white", mec=INK, mew=0.8, ms=4, label="HPC"),
      Line2D([], [], marker="o", ls="none", color=INK, ms=4, label="PFC"),
      Line2D([], [], color=MEN, lw=1.6,
             label=f"Female-biased in HPC → male-biased in PFC (n = {int(hc.switch_F_to_M.sum())}; bold)")]
axc.legend(handles=hl, loc="upper left", bbox_to_anchor=(0.30, 1.0), ncol=3, fontsize=6, handletextpad=0.3,
           columnspacing=1.2, borderaxespad=0.1, handlelength=1.4)

# =============================================================================
# D  58 non-shared HPC genes, coloured by category, shape by PFC status
# =============================================================================
plabel("D", 0.0, 2.42)
ns = hpc[hpc.category != "sex-shared"].copy()
pf = pfc.set_index("gene").category
pu = pfcu.set_index("gene").category


def pfc_status(g):
    if g in pf.index:
        return "PFC sex-shared (filtered list)" if pf[g] == "sex-shared" else "PFC other category"
    if g in pu.index:
        return "PFC sex-shared, unfiltered list only" if pu[g] == "sex-shared" else "PFC other category (unfiltered only)"
    return "absent from PFC lists"


ns["pfc_status"] = ns.gene.map(pfc_status)
STATUS = {"PFC sex-shared (filtered list)": "o", "PFC sex-shared, unfiltered list only": "s",
          "PFC other category": "^", "absent from PFC lists": "X"}
assert set(ns.pfc_status) <= set(STATUS), set(ns.pfc_status)
ns.to_csv(f"{SRC}/Figure2D.csv", index=False)
CATS3 = ["female-specific", "male-specific", "sex-dimorphic"]
axd = ax_in(0.47, 0.36, 2.25, 1.95)
ref_lines(axd, -2.3, 2.5)
for c in CATS3:
    for st, m in STATUS.items():
        ss = ns[(ns.category == c) & (ns.pfc_status == st)]
        if len(ss):
            axd.scatter(ss.lfc_m, ss.lfc_f, marker=m, s=22 if m != "X" else 26, color=CAT_COLORS[c],
                        edgecolors="white", linewidths=0.5, zorder=3)
LBL_D = ["PRDM11", "CXCL12", "DHRS9", "MTRNR2L8", "OTUD1", "GABBR1", "BCL10", "TICAM1", "SLC27A6", "BTBD11",
         "INO80D", "SLC26A3"]
nsl = ns.set_index("gene")
texts = [axd.text(nsl.loc[g, "lfc_m"], nsl.loc[g, "lfc_f"], g, fontsize=6, style="italic", color=INK, zorder=4)
         for g in LBL_D]
adjust_text(texts, x=ns.lfc_m.values, y=ns.lfc_f.values, ax=axd, expand=(1.2, 1.4), force_text=(0.3, 0.5),
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.4, shrinkA=0, shrinkB=2))
axd.set_xlim(-2.3, 2.5)
axd.set_ylim(-0.95, 2.2)
axd.set_xlabel("log₂FC, men (HPC)", fontsize=7, labelpad=1)
axd.set_ylabel("log₂FC, women (HPC)", fontsize=7, labelpad=1)
axd.tick_params(labelsize=6)
# legends: colour = HPC category; shape = PFC status
lh = [Line2D([], [], marker="o", ls="none", color=CAT_COLORS[c], ms=4.5,
             label=f"{c} ({int((ns.category == c).sum())})") for c in CATS3]
fig.legend(handles=lh, loc="upper left", bbox_to_anchor=(2.8 / FW, 2.32 / FH), fontsize=6, handletextpad=0.3,
           borderaxespad=0, labelspacing=0.45, title="HPC category (color)", title_fontsize=6, alignment="left")
lh2 = [Line2D([], [], marker=m, ls="none", color=INK2, ms=4.5 if m != "X" else 5,
              label=f"{st} ({int((ns.pfc_status == st).sum())})") for st, m in STATUS.items()]
fig.legend(handles=lh2, loc="upper left", bbox_to_anchor=(2.8 / FW, 1.55 / FH), fontsize=6, handletextpad=0.3,
           borderaxespad=0, labelspacing=0.45, title="Status in PFC (shape)", title_fontsize=6, alignment="left")

# =============================================================================
# E  PFC non-shared genes after the filter (8) and removed by it (24)
# =============================================================================
plabel("E", 4.52, 2.42)
kept = pfc[pfc.category != "sex-shared"].assign(list="retained after filter")
removed = pfcu[(pfcu.category != "sex-shared") & ~pfcu.gene.isin(pfc.gene)].assign(list="removed by filter")
assert len(kept) == 8 and len(removed) == 24, (len(kept), len(removed))
e = pd.concat([kept, removed])
e.to_csv(f"{SRC}/Figure2E.csv", index=False)
axe = ax_in(4.98, 0.36, 2.17, 1.95)
axe.set_xscale("symlog", linthresh=1, linscale=1.2)
axe.set_yscale("symlog", linthresh=1, linscale=1.2)
ref_lines(axe, -4, 18)
for c in CATS3:
    k_ = kept[kept.category == c]
    r_ = removed[removed.category == c]
    axe.scatter(k_.lfc_m, k_.lfc_f, marker="o", s=22, color=CAT_COLORS[c], edgecolors="white", linewidths=0.5,
                zorder=3)
    axe.scatter(r_.lfc_m, r_.lfc_f, marker="s", s=16, facecolors="white", edgecolors=CAT_COLORS[c],
                linewidths=0.9, zorder=3)
LBL_E = ["ZBTB16", "VWC2", "EDIL3", "HLA-DQB1", "CCL4", "ZBTB20", "SMAP2"]
el = e.set_index("gene")
OFF_E = {"ZBTB20": (-34, 10, "right"), "SMAP2": (-34, -9, "right"), "VWC2": (8, 4, "left"),
         "EDIL3": (8, 0, "left"), "ZBTB16": (8, -6, "left"), "HLA-DQB1": (8, -6, "left"),
         "CCL4": (-7, 0, "right")}
for g in LBL_E:
    dx, dy, ha = OFF_E[g]
    axe.annotate(g, (el.loc[g, "lfc_m"], el.loc[g, "lfc_f"]), xytext=(dx, dy), textcoords="offset points",
                 fontsize=6, style="italic", color=INK, ha=ha, va="center", zorder=4,
                 arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.4, shrinkA=0, shrinkB=2.5)
                 if abs(dx) > 10 else None)
tk = [-3, -1, 0, 1, 3, 10]
axe.set_xlim(-4, 18)
axe.set_ylim(-4, 18)
axe.set_xticks(tk)
axe.set_yticks(tk)
axe.set_xticklabels([f"{t:g}".replace("-", "−") for t in tk], fontsize=6)
axe.set_yticklabels([f"{t:g}".replace("-", "−") for t in tk], fontsize=6)
axe.minorticks_off()
axe.set_xlabel("log₂FC, men (PFC)", fontsize=7, labelpad=1)
axe.set_ylabel("log₂FC, women (PFC)", fontsize=7, labelpad=1)
le = [Line2D([], [], marker="o", ls="none", color=INK2, ms=4.5, label="Retained (8)"),
      Line2D([], [], marker="s", ls="none", mfc="white", mec=INK2, mew=0.9, ms=4, label="Removed by filter (24)")]
axe.legend(handles=le, loc="upper left", fontsize=6, handletextpad=0.2, borderaxespad=0.2, labelspacing=0.3,
           title="Colors as in D", title_fontsize=6, alignment="left")
axe.text(0.98, 0.02, "symmetric-log axes\n(linear within ±1)", transform=axe.transAxes, fontsize=6,
         color=MUTED, ha="right", va="bottom")

pd.DataFrame(stats_rows).to_csv(f"{SRC}/Figure2_statistics.csv", index=False)
plt.rcParams["savefig.pad_inches"] = 0.03
save(fig, f"{OUT}/Figure2")
print(pd.DataFrame(stats_rows).to_string())
print(hc.loc[hc.switch_F_to_M, "gene"].tolist())
print(ns.pfc_status.value_counts())
