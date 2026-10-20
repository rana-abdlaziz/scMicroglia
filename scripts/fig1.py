"""Figure 1 - study design, DEG categories, cross-region overlap and the 109-gene core."""
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *  # noqa: E402,F401
from figdata import *  # noqa: E402,F401

apply_style()
os.makedirs(SRC, exist_ok=True)

FW, FH = W2, 8.7
fig = plt.figure(figsize=(FW, FH))


def ax_in(x, y, w, h, **kw):
    """Axes placed in inches from the bottom-left corner of the figure."""
    return fig.add_axes([x / FW, y / FH, w / FW, h / FH], **kw)


def plabel(ax, letter, x, y):
    """Panel letter at a figure position in inches (top-left)."""
    fig.text(x / FW, y / FH, letter, fontsize=11, fontweight="bold", va="top", ha="left", color=INK)


# =============================================================================
# A  study design schematic
# =============================================================================
axA = ax_in(0, 6.62, FW, 2.0)
axA.set_xlim(0, FW)
axA.set_ylim(0, 2.0)
axA.axis("off")
plabel(axA, "A", 0.0, FH - 0.02)

BW, BH = 1.62, 0.74
xs = [0.28, 2.04, 3.80, 5.56]
y1, y2 = 1.13, 0.08
FILL = "#f7f6f3"

boxes = [
    ("ROSMAP microglia atlas", ["Sun et al. 2023", "snRNA-seq, microglial nuclei", "AD and non-AD donors"]),
    ("Per-region objects", ["HPC: 86 samples (66 identified donors)", "PFC: 469 samples (423 donors)", "Women and men separately"]),
    ("State annotation", ["12 microglial states", "Module scores per nucleus", "Hungarian cluster→state assignment"]),
    ("Within-sex contrast", ["Inflammatory MG2+MG8+MG10 (AD)", "vs homeostatic MG0 (non-AD)", "Poisson GLM; age as latent variable"]),
    ("Categorization", None),
    ("Effect-size filter", ["Shared: |log₂FC| >1.5 or <0.5", "  in both sexes; specific: ≥2×", "  other sex; PFC 1,967 → 1,136"]),
    ("Cross-region core", ["HPC sex-shared ∩ PFC sex-shared", "109 genes", "(101 increased, 8 decreased)"]),
    ("Enrichment & networks", ["ORA (GO, KEGG), gseGO", "MetaCore TF network", "decoupleR TF activity"]),
]
pos = [(x, y1) for x in xs] + [(x, y2) for x in xs]
for i, ((title, body), (x, y)) in enumerate(zip(boxes, pos)):
    axA.add_patch(FancyBboxPatch((x, y), BW, BH, boxstyle="round,pad=0,rounding_size=0.06",
                                 fc=FILL, ec=INK2, lw=0.6))
    axA.add_patch(Circle((x + 0.11, y + BH - 0.11), 0.072, fc=INK, ec="none"))
    axA.text(x + 0.11, y + BH - 0.112, str(i + 1), color="white", fontsize=6, fontweight="bold",
             ha="center", va="center")
    axA.text(x + 0.22, y + BH - 0.11, title, fontsize=7, fontweight="bold", va="center", ha="left", color=INK)
    if body:
        for j, line in enumerate(body):
            axA.text(x + 0.08, y + BH - 0.27 - j * 0.125, line, fontsize=6, va="center", ha="left", color=INK2)
    else:  # categorisation swatches
        cats = [("sex-shared", SHARED), ("female-specific", WOMEN), ("male-specific", MEN), ("sex-dimorphic", DIMORPHIC)]
        for k, (lab, col) in enumerate(cats):
            cx = x + 0.1 + (k % 2) * 0.78
            cy = y + BH - 0.30 - (k // 2) * 0.16
            axA.add_patch(Rectangle((cx, cy - 0.04), 0.08, 0.08, fc=col, ec="none"))
            axA.text(cx + 0.12, cy, lab, fontsize=6, va="center", ha="left", color=INK2)
        axA.text(x + 0.08, y + 0.08, "P_adj, |log₂FC| and sign per sex", fontsize=6, va="center",
                 ha="left", color=MUTED)

arrow_kw = dict(arrowstyle="-|>", color=INK2, lw=0.8, mutation_scale=7, shrinkA=0, shrinkB=0)
for row_y in (y1, y2):
    for a, b in zip(xs[:-1], xs[1:]):
        axA.annotate("", xy=(b - 0.01, row_y + BH / 2), xytext=(a + BW + 0.01, row_y + BH / 2), arrowprops=arrow_kw)
# return connector 4 -> 5
mid = (y2 + BH + y1) / 2
axA.plot([xs[3] + BW / 2, xs[3] + BW / 2, xs[0] + BW / 2], [y1, mid, mid], color=INK2, lw=0.8,
         solid_capstyle="butt")
axA.annotate("", xy=(xs[0] + BW / 2, y2 + BH + 0.005), xytext=(xs[0] + BW / 2, mid + 0.004), arrowprops=arrow_kw)

pd.DataFrame({"step": range(1, 9), "title": [b[0] for b in boxes],
              "details": ["; ".join(b[1]) if b[1] else "sex-shared; female-specific; male-specific; sex-dimorphic"
                          for b in boxes]}).to_csv(f"{SRC}/Figure1A.csv", index=False)

# =============================================================================
# B  DEG counts per category, region and list (log-scale dot plot)
# =============================================================================
hpc, pfc, pfcu = load_hpc(), load_pfc(), load_pfc_unfiltered()
lists = [("HPC", "filtered", hpc), ("PFC", "before filter", pfcu), ("PFC", "after filter", pfc)]
rows = []
for reg, lab, d in lists:
    d = d.assign(dir=direction(d))
    for c in CAT_ORDER:
        sub = d[d.category == c]
        rows.append(dict(region=reg, list=lab, category=c, increased=int((sub.dir > 0).sum()),
                         decreased=int((sub.dir < 0).sum()), total=len(sub)))
cnt = pd.DataFrame(rows)
cnt.to_csv(f"{SRC}/Figure1B.csv", index=False)

bx0, bw, bgap, by0, bh = 0.98, 0.86, 0.09, 4.62, 1.30
plabel(None, "B", 0.0, 6.48)
OFF = 0.2
for k, (reg, lab, d) in enumerate(lists):
    ax = ax_in(bx0 + k * (bw + bgap), by0, bw, bh)
    sub = cnt[(cnt.region == reg) & (cnt.list == lab)].set_index("category").loc[CAT_ORDER]
    ax.set_xscale("log")
    ax.set_xlim(0.45, 6000)
    ax.set_ylim(len(CAT_ORDER) - 0.5, -0.55)
    for yy in range(len(CAT_ORDER)):
        ax.axhline(yy, color=GRID, lw=0.5, zorder=0)
    for i, c in enumerate(CAT_ORDER):
        for key, off, col, mk in (("increased", -OFF, UP, "^"), ("decreased", OFF, DOWN, "v")):
            n = sub.loc[c, key]
            yy = i + off
            if n == 0:
                ax.text(0.62, yy, "0", fontsize=6, color=MUTED, ha="center", va="center")
                continue
            ax.plot(n, yy, marker=mk, ms=4.6, color=col, mec="white", mew=0.5, ls="none", zorder=3)
            if n >= 100:
                ax.text(n / 1.7, yy, f"{n:,}", fontsize=6, color=INK, ha="right", va="center")
            else:
                ax.text(n * 1.7, yy, f"{n:,}", fontsize=6, color=INK, ha="left", va="center")
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(["1", "10", "100", "1,000"], fontsize=6)
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_yticks(range(len(CAT_ORDER)))
    ax.set_yticklabels(CAT_ORDER if k == 0 else [], fontsize=7)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_title(f"{reg}, {lab}\nn = {int(sub.total.sum()):,}", fontsize=7, fontweight="normal", loc="center",
                 color=INK, pad=3)
    if k == 1:
        ax.set_xlabel("Number of DEGs (log scale)", fontsize=7, labelpad=2)
# legend + direction note
handles = [Line2D([], [], marker="^", color=UP, ls="none", ms=4.6, label="Increased"),
           Line2D([], [], marker="v", color=DOWN, ls="none", ms=4.6, label="Decreased")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.98 / FW, 6.47 / FH), ncol=2, fontsize=6.5,
           handletextpad=0.2, columnspacing=1.0, borderaxespad=0)
fig.text(0.25 / FW, 4.08 / FH, "Direction refers to: sex-shared, both sexes; female-specific, women;\nmale-specific and sex-dimorphic, men. "
         "Filter applied to PFC only (unfiltered\nHPC list not retained).",
         fontsize=6, color=INK2, va="top", ha="left", linespacing=1.3)

# =============================================================================
# C  UpSet plot of exclusive intersections (computed here)
# =============================================================================
sets = {}
for reg, d in (("HPC", hpc), ("PFC", pfc)):
    for c in CAT_ORDER:
        sets[f"{reg} {c}"] = set(d.gene[d.category == c])
set_names = list(sets)
allg = set().union(*sets.values())
memb = Counter(tuple(n for n in set_names if g in sets[n]) for g in allg)
inter = sorted(memb.items(), key=lambda kv: (-kv[1], len(kv[0])))
up_rows = []
for key, n in inter:
    genes = sorted(g for g in allg if tuple(s for s in set_names if g in sets[s]) == key)
    up_rows.append(dict(intersection=" & ".join(key), size=n, genes=";".join(genes)))
pd.DataFrame(up_rows).to_csv(f"{SRC}/Figure1C.csv", index=False)

plabel(None, "C", 0.0, 3.72)
mx0, mw = 1.42, 2.25
ncol = len(inter)
colw = mw / ncol
axbar = ax_in(mx0, 1.86, mw, 1.62)
axmat = ax_in(mx0, 0.22, mw, 1.56)
xpos = np.arange(ncol)
sizes = [n for _, n in inter]
is_core = [k == ("HPC sex-shared", "PFC sex-shared") for k, _ in inter]
axbar.bar(xpos, sizes, width=0.62, color=[INK if c else INK2 for c in is_core], ec="none")
for x, n in zip(xpos, sizes):
    axbar.text(x, n + 18, f"{n:,}", ha="center", va="bottom", fontsize=6, color=INK)
axbar.set_xlim(-0.75, ncol - 0.4)
axbar.set_ylim(0, 1150)
axbar.set_xticks([])
axbar.set_yticks([0, 500, 1000])
axbar.set_yticklabels(["0", "500", "1,000"], fontsize=6)
axbar.set_ylabel("Genes in exclusive\nintersection", fontsize=7)
axbar.spines["bottom"].set_visible(False)
ci = is_core.index(True)
axbar.annotate("cross-region core", xy=(ci + 0.15, sizes[ci] + 95), xytext=(ci + 1.4, 560), fontsize=6,
               color=INK, ha="left", va="center",
               arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5, shrinkA=1, shrinkB=0))

nset = len(set_names)
axmat.set_xlim(-0.75, ncol - 0.4)
axmat.set_ylim(nset - 0.5, -0.5)
axmat.axis("off")
for r in range(nset):
    if set_names[r].startswith("PFC"):
        axmat.add_patch(Rectangle((-0.75, r - 0.5), ncol + 0.35, 1, fc="#f4f3f0", ec="none", zorder=0))
for x, (key, _) in enumerate(inter):
    ys = [set_names.index(s) for s in key]
    axmat.scatter([x] * nset, range(nset), s=16, color=GRID, zorder=1, linewidths=0)
    axmat.scatter([x] * len(ys), ys, s=18, color=INK, zorder=3, linewidths=0)
    if len(ys) > 1:
        axmat.plot([x, x], [min(ys), max(ys)], color=INK, lw=1.2, zorder=2)
# set labels with category swatch + set size
axlab = ax_in(0.05, 0.22, mx0 - 0.08, 1.56)
axlab.set_xlim(0, 1)
axlab.set_ylim(nset - 0.5, -0.5)
axlab.axis("off")
for r, s in enumerate(set_names):
    reg, cat = s.split(" ", 1)
    if reg == "PFC":
        axlab.add_patch(Rectangle((0, r - 0.5), 1.05, 1, fc="#f4f3f0", ec="none", zorder=0, clip_on=False))
    axlab.text(0.93, r, f"{reg} {cat} ({len(sets[s]):,})", ha="right", va="center", fontsize=6.5, color=INK)
    axlab.add_patch(Rectangle((0.955, r - 0.17), 0.04, 0.34, fc=CAT_COLORS[cat], ec="none"))

# =============================================================================
# D  heatmap of the 109 core genes x 4 contrasts
# =============================================================================
core = core109()
GROUP_ORDER = [
    ("Chemokines/cytokines & receptors", "Chemokines, cytokines & receptors"),
    ("Alarmins/complement/immune receptors & signalling", "Alarmins & immune signalling"),
    ("AP-1/IEG transcription factors", "AP-1 immediate-early TFs"),
    ("EGR family", "EGR family"),
    ("NR4A family", "NR4A family"),
    ("Other transcriptional/chromatin regulators", "Other transcriptional regulators"),
    ("Apoptosis (BCL2 family) & autophagy", "Apoptosis (BCL2) & autophagy"),
    ("Lipid / proteoglycan / glycan metabolism", "Lipid, proteoglycan & glycan"),
    ("Cytoskeleton / small-GTPase regulation / motility", "Cytoskeleton, GTPases & motility"),
    ("Vesicle trafficking / endolysosome / nuclear transport", "Vesicle trafficking & endolysosome"),
    ("Solute/amino-acid/ion transport", "Solute & ion transport"),
    ("Proteostasis / translation / DNA repair / other", "Proteostasis, translation & other"),
]
gorder = {g: i for i, (g, _) in enumerate(GROUP_ORDER)}
core["group_rank"] = core.functional_group.map(gorder)
core = core.sort_values(["group_rank", "sumFC"], ascending=[True, False]).reset_index(drop=True)
CONTR = ["HPC_men", "HPC_women", "PFC_men", "PFC_women"]
CLIP = 4.0
core_out = core[["gene", "functional_group", "stratum", "state_marker", "feedback_regulator"] + CONTR + ["sumFC"]]
core_out.to_csv(f"{SRC}/Figure1D.csv", index=False)

blocks = [GROUP_ORDER[:8], GROUP_ORDER[8:]]
plabel(None, "D", 3.80, 6.48)
DX0 = 3.86            # left edge of right column
LBLW = 0.50           # gene-label width
CW, CGAP = 0.165, 0.05
GW = 0.11
blk_w = LBLW + 4 * CW + CGAP + 0.06 + 3 * GW
top, bottom = 5.93, 0.18
nlines = max(sum(len(core[core.functional_group == g]) + 1 for g, _ in b) for b in blocks)
pitch = (top - bottom) / nlines
for b_i, blk in enumerate(blocks):
    x0 = DX0 + b_i * (blk_w + 0.12)
    ax = ax_in(x0, bottom, blk_w, top - bottom)
    ax.set_xlim(0, blk_w)
    ax.set_ylim(top - bottom, 0)
    ax.axis("off")
    xcell = [LBLW + i * CW + (CGAP if i >= 2 else 0) for i in range(4)]
    xg = [LBLW + 4 * CW + CGAP + 0.06 + j * GW + GW / 2 for j in range(3)]
    yy = 0.0
    for g, gname in blk:
        ax.text(0.0, yy + pitch * 0.55, gname, fontsize=6, fontweight="bold", color=INK, va="center", ha="left")
        yy += pitch
        sub = core[core.functional_group == g]
        for _, r in sub.iterrows():
            ax.text(LBLW - 0.03, yy + pitch / 2, r.gene, fontsize=6, style="italic", color=INK,
                    va="center", ha="right")
            for i, c in enumerate(CONTR):
                v = r[c]
                col = DIVERGING((np.clip(v, -CLIP, CLIP) + CLIP) / (2 * CLIP))
                ax.add_patch(Rectangle((xcell[i], yy + 0.004), CW - 0.012, pitch - 0.008, fc=col, ec="none"))
                if abs(v) > CLIP:
                    ax.text(xcell[i] + (CW - 0.012) / 2, yy + pitch / 2, f"{v:.1f}", fontsize=6,
                            color="white", ha="center", va="center")
            # annotation glyphs: high effect, state marker, feedback regulator
            if r.stratum == "high":
                ax.plot(xg[0], yy + pitch / 2, marker="s", ms=3.0, color=INK, ls="none")
            if r.state_marker:
                ax.plot(xg[1], yy + pitch / 2, marker="o", ms=3.2, color=INK, ls="none")
            if r.feedback_regulator:
                ax.plot(xg[2], yy + pitch / 2, marker="D", ms=2.8, color=INK, ls="none")
            yy += pitch
    # column headers
    hy = -0.04
    for i, c in enumerate(CONTR):
        ax.text(xcell[i] + (CW - 0.012) / 2, hy, "M" if c.endswith("men") and not c.endswith("women") else "F",
                fontsize=6.5, ha="center", va="bottom", color=INK, clip_on=False)
    for i0, reg in ((0, "HPC"), (2, "PFC")):
        xa, xb = xcell[i0], xcell[i0 + 1] + CW - 0.012
        ax.plot([xa, xb], [hy - 0.15, hy - 0.15], color=INK2, lw=0.6, clip_on=False)
        ax.text((xa + xb) / 2, hy - 0.17, reg, fontsize=7, fontweight="bold", ha="center", va="bottom",
                color=INK, clip_on=False)
    for j, lab in enumerate(["High effect", "State marker", "Feedback reg."]):
        ax.text(xg[j], hy, lab, rotation=90, fontsize=6, ha="center", va="bottom", color=INK2, clip_on=False)

# legend for D (below the shorter second block)
lines2 = sum(len(core[core.functional_group == g]) + 1 for g, _ in blocks[1])
leg_top = top - lines2 * pitch - 0.12
lx = DX0 + blk_w + 0.12
cax = ax_in(lx + 0.08, leg_top - 0.36, 1.35, 0.08)
grad = np.linspace(-CLIP, CLIP, 256)[None, :]
cax.imshow(grad, aspect="auto", cmap=DIVERGING, vmin=-CLIP, vmax=CLIP, extent=[-CLIP, CLIP, 0, 1])
cax.set_yticks([])
cax.set_xticks([-4, -2, 0, 2, 4])
cax.set_xticklabels(["≤−4", "−2", "0", "2", "≥4"], fontsize=6)
cax.tick_params(length=2, pad=1)
for s in cax.spines.values():
    s.set_visible(False)
fig.text((lx + 0.08) / FW, (leg_top - 0.02) / FH, "log₂FC, inflammatory-AD vs\nhomeostatic non-AD",
         fontsize=6, ha="left", va="top", color=INK)
ky = leg_top - 0.62
for j, (mk, ms, lab) in enumerate([("s", 3.0, "High effect (|log₂FC| > 1.5)"),
                                   ("o", 3.2, "State marker (MG8/MG10 list)"),
                                   ("D", 2.8, "Negative-feedback regulator")]):
    fig.add_artist(Line2D([(lx + 0.12) / FW], [(ky - j * 0.12) / FH], marker=mk, ms=ms, color=INK, ls="none",
                          transform=fig.transFigure))
    fig.text((lx + 0.22) / FW, (ky - j * 0.12) / FH, lab, fontsize=6, va="center", ha="left", color=INK2)
fig.text((lx + 0.08) / FW, (ky - 3 * 0.12 + 0.03) / FH,
         "M, men; F, women. Colors clipped at ±4;\nvalues beyond ±4 are printed.",
         fontsize=6, va="top", ha="left", color=INK2, linespacing=1.3)

plt.rcParams["savefig.pad_inches"] = 0.03
save(fig, f"{OUT}/Figure1")
print("state markers in core:", core.loc[core.state_marker != "", ["gene", "state_marker"]].values.tolist())
print(cnt)
print(inter)
