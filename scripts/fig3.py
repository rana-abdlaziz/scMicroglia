"""Figure 3 - pathway enrichment of sex-shared, sex-specific and sex-dimorphic DEGs.

Panels
  A  GO-BP / GO-MF dot plots of HPC and PFC sex-shared genes (top terms, FDR < 0.05) with the
     number of increased / decreased member genes per term (F_updown_recomputed.csv).
  B  KEGG pathways of the 109 core genes (existing clusterProfiler table).
  C  GO-BP of the 109 core genes (top 20) - RECOMPUTED (fig3_enrichGO_recompute.R).
  D  KEGG of HPC male- and female-specific genes (existing tables).
  E  GO of HPC / PFC sex-dimorphic genes - RECOMPUTED (fig3_enrichGO_recompute.R).
Run fig3_enrichGO_recompute.R first.
"""
import sys, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D

sys.path.insert(0, "/home/claude/scm/figures/scripts")
from figstyle import *  # noqa

apply_style()
T = "/home/claude/scm/agents_out/tables/"
E = T + "enrichment/"
OUT = "/home/claude/scm/figures/out/"
SD = OUT + "source_data/"
FS = 6.3  # term-label font size (pt)

# --------------------------------------------------------------------------- helpers
ABBR = [("RNA polymerase II-specific", "RNA Pol II-specific"), ("DNA-binding transcription factor", "DNA-binding TF"),
        ("transcription factor", "TF"), ("positive regulation of", "pos. reg. of"),
        ("negative regulation of", "neg. reg. of"), ("regulation of", "reg. of"),
        ("involved in", "in"), ("establishment of", "estab. of"), ("nucleoside-triphosphatase", "NTPase"),
        (" signaling pathway", " signaling"), ("AGE-RAGE signaling in diabetic complications", "AGE-RAGE signaling (diabetic compl.)"),
        ("structural constituent", "struct. constituent"), ("transmembrane transporter", "transm. transporter"),
        ("- Homo sapiens (human)", "")]


def short(s, n=46):
    for a, b in ABBR:
        s = s.replace(a, b)
    s = s.strip()
    return s if len(s) <= n else s[: n - 1].rstrip() + "…"


def nlog(p):
    return -np.log10(p)


def ratio(gr):
    a, b = gr.split("/")
    return int(a) / int(b)


def dedup(df):
    """drop terms whose member-gene set is identical to a better-ranked term (e.g. taxis = chemotaxis)."""
    seen, keep = set(), []
    for g in df.geneID:
        key = frozenset(str(g).split("/"))
        keep.append(key not in seen)
        seen.add(key)
    return df[keep]


# --------------------------------------------------------------------------- data A
ud = pd.read_csv(T + "F_updown_recomputed.csv")
ud["Count"] = ud.n_up + ud.n_down
ud["ratio"] = ud.GeneRatio.map(ratio)
selA = []
NTOP = {"BP": 8, "MF": 8}
for reg in ["HPC", "PFC"]:
    for ont in ["BP", "MF"]:
        d = ud[(ud.set == f"{reg}_{ont}_GenderShared") & (ud["p.adjust"] < 0.05)].sort_values("p.adjust")
        d = dedup(d)
        top = d.head(NTOP[ont]).copy()
        top["group"] = "top"
        if ont == "BP":
            # also show the best terms dominated by decreased genes (only HPC has any)
            dn = d[(d.upDown_correct == "down") & ~d.ID.isin(top.ID)].head(3).copy()
            dn["group"] = "down-dominated"
            if len(dn) == 0:  # keep facet heights equal: next-best terms
                dn = d[~d.ID.isin(top.ID)].head(3).copy()
                dn["group"] = "top"
            top = pd.concat([top, dn])
        top["region"], top["ont"] = reg, ont
        selA.append(top)
A = pd.concat(selA)
A[["region", "ont", "group", "ID", "Description", "GeneRatio", "Count", "p.adjust", "n_up", "n_down", "upDown_correct", "geneID"]] \
    .to_csv(SD + "Fig3A_GO_shared_HPC_PFC_dotplot.csv", index=False)

# --------------------------------------------------------------------------- data B
B = pd.read_csv(T + "F_shared_pfc_hipo__shared_DEGs_enrichment_KEGG_FDRadj_symbols.csv")
B["Count"] = B.GeneRatio.str.split("/").str[0].astype(int)
B.to_csv(SD + "Fig3B_KEGG_core109.csv", index=False)

# --------------------------------------------------------------------------- data C (recomputed)
C = pd.read_csv(SD + "Fig3_recomputed_enrichGO_core109_BP.csv").sort_values("p.adjust").head(20)
C.to_csv(SD + "Fig3C_GOBP_core109_top20_recomputed.csv", index=False)

# --------------------------------------------------------------------------- data D
km = pd.read_csv(E + "HPC_Pathway_Enrichment_Analysis_KEGG_Male.csv")
km = km[km["p.adjust"] < 0.05].copy()
km["genes"] = "MED12L, PFKFB2"  # Entrez 116931 / 5208
assert set(km.geneID) == {"116931/5208"}
kf = pd.read_csv(T + "F_hipo_analysis__DEGs_enrichment_KEGG_FDRadj_female_symbols.csv")
kf = kf[kf["p.adjust"] < 0.05].copy()
km["sex"], kf["sex"] = "HPC male-specific", "HPC female-specific"
D = pd.concat([km[["sex", "ID", "Description", "GeneRatio", "p.adjust", "genes"]],
               kf[["sex", "ID", "Description", "GeneRatio", "p.adjust", "genes"]]])
D["genes"] = D.genes.str.replace(",", ", ").str.replace(",  ", ", ")
D.to_csv(SD + "Fig3D_KEGG_HPC_male_female_specific.csv", index=False)

# --------------------------------------------------------------------------- data E (recomputed)
Elist = []
for reg, ont, n in [("HPC", "BP", 4), ("HPC", "MF", 3), ("PFC", "BP", 3), ("PFC", "MF", 2)]:
    d = pd.read_csv(SD + f"Fig3_recomputed_enrichGO_{reg}_dimorphic_{ont}.csv").sort_values(["p.adjust", "pvalue"])
    d = d.head(n).copy()
    d["region"], d["ont"] = reg, ont
    Elist.append(d)
Ed = pd.concat(Elist)
Ed.to_csv(SD + "Fig3E_GO_dimorphic_recomputed.csv", index=False)

# =========================================================================== figure
FW, FH = W2, 9.1
fig = plt.figure(figsize=(FW, FH))


def axin(x, y, w, h):
    """axes from inch coordinates measured from the top-left corner."""
    return fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH])


def plabel(letter, x, y):
    fig.text(x / FW, 1 - y / FH, letter, fontsize=11, fontweight="bold", va="top", ha="left", color=INK)


def ptitle(text, x, y, **kw):
    fig.text(x / FW, 1 - y / FH, text, fontsize=8, fontweight="bold", va="top", ha="left", color=INK, **kw)


RH = 0.112  # row height (in)

# ---------------- A ----------------------------------------------------------------
norm = Normalize(vmin=1.3, vmax=12, clip=True)
cmap = SEQUENTIAL
SIZE = lambda c: 2.0 + np.asarray(c, float) * 1.1  # marker area (pt^2) from gene count
plabel("A", 0.0, 0.02)
nrow = {"BP": 11, "MF": 8}
ytop = {"BP": 0.30, "MF": 0.30 + 11.4 * RH + 0.40}
xcol = {"HPC": 1.72, "PFC": 5.32}  # left edge of dot axes
DW, UW = 0.92, 0.30
for reg in ["HPC", "PFC"]:
    ptitle(f"{reg} sex-shared genes", xcol[reg] - 1.50, 0.04)
    for ont in ["BP", "MF"]:
        d = A[(A.region == reg) & (A.ont == ont)].reset_index(drop=True)
        n = len(d)
        y = np.arange(n)[::-1].astype(float)
        ndn = int((d.group == "down-dominated").sum())
        span = n + (1.0 if ndn else 0)
        ax = axin(xcol[reg], ytop[ont], DW, span * RH)
        axu = axin(xcol[reg] + DW + 0.06, ytop[ont], UW, span * RH)
        if ndn:
            y[-ndn:] -= 1.0
        ax.scatter(d.ratio, y, s=SIZE(d.Count), c=nlog(d["p.adjust"]), cmap=cmap, norm=norm,
                   edgecolor=INK2, linewidth=0.3, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels([short(s, 40) for s in d.Description], fontsize=FS)
        ax.set_ylim(y.min() - 0.55, y.max() + 0.55)
        xmax = d.ratio.max()
        ax.set_xlim(0, xmax * 1.2)
        ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        ax.grid(axis="x", color=GRID, lw=0.4)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0, pad=2)
        ax.tick_params(axis="x", labelsize=6, pad=1.5)
        ax.set_xlabel("GeneRatio", fontsize=6, labelpad=1)
        ax.text(1.0, 1.0, f"GO-{ont}", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=6.3, color=INK2, fontweight="bold")
        if ndn:
            yy = y[-ndn:].max() + 1.0
            ax.axhline(yy - 0.15, color=MUTED, lw=0.4)
            ax.text(-0.02, yy, "best terms dominated by decreased genes", transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=6, color=INK2, style="italic")
        fu = d.n_up / d.Count
        axu.barh(y, fu, height=0.64, color=UP, lw=0)
        axu.barh(y, 1 - fu, left=fu, height=0.64, color=DOWN, lw=0)
        for yi, nu, nd in zip(y, d.n_up, d.n_down):
            axu.text(1.05, yi, f"{nu}/{nd}", va="center", ha="left", fontsize=6, color=INK2)
        axu.set_xlim(0, 1)
        axu.set_ylim(ax.get_ylim())
        axu.axis("off")
        if ont == "BP":
            axu.text(0.0, 1.0 + 0.1 / (span * RH), "up/down", transform=axu.transAxes, fontsize=6, color=INK2, va="bottom")

# legend strip for A
ly = ytop["MF"] + 8 * RH + 0.40
cax = axin(1.05, ly + 0.02, 1.1, 0.07)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", extend="max")
cb.ax.tick_params(labelsize=6, width=0.4, length=2, pad=1)
cb.outline.set_linewidth(0.4)
fig.text(1.0 / FW, 1 - (ly + 0.055) / FH, "−log$_{10}$ FDR", fontsize=6, ha="right", va="center", color=INK2)
lk = axin(2.65, ly - 0.06, 1.75, 0.22); lk.axis("off"); lk.set_xlim(0, 1); lk.set_ylim(0, 1)
lk.text(0.0, 0.5, "genes", fontsize=6, color=INK2, va="center")
for i, c in enumerate([5, 25, 50, 75]):
    xx = 0.22 + i * 0.2
    lk.scatter(xx, 0.5, s=SIZE(c), color="#d8d6d0", edgecolor=INK2, lw=0.3, clip_on=False)
    lk.text(xx + 0.055, 0.5, str(c), va="center", fontsize=6, color=INK2)
uk = axin(4.40, ly - 0.06, 2.5, 0.22); uk.axis("off"); uk.set_xlim(0, 1); uk.set_ylim(0, 1)
uk.add_patch(plt.Rectangle((0.0, 0.3), 0.05, 0.4, color=UP, lw=0))
uk.text(0.065, 0.5, "increased", fontsize=6, va="center", color=INK2)
uk.add_patch(plt.Rectangle((0.25, 0.3), 0.05, 0.4, color=DOWN, lw=0))
uk.text(0.315, 0.5, "decreased member genes (n up/down)", fontsize=6, va="center", color=INK2)

# ---------------- B ----------------------------------------------------------------
yB = ly + 0.30
plabel("B", 0.0, yB)
ptitle("KEGG pathways of the 109 core genes (77 KEGG-annotated)", 0.22, yB + 0.02)
d = B.sort_values("p.adjust").reset_index(drop=True)
axB = axin(1.72, yB + 0.30, 1.25, len(d) * RH)
y = np.arange(len(d))[::-1]
axB.barh(y, nlog(d["p.adjust"]), height=0.66, color=SHARED, lw=0)
axB.set_yticks(y)
axB.set_yticklabels([short(s, 40) for s in d.Description], fontsize=FS)
axB.tick_params(axis="y", length=0, pad=2)
axB.set_xlim(0, 3.0)
axB.set_ylim(-0.55, y.max() + 0.55)
axB.set_xlabel("−log$_{10}$ FDR", fontsize=6, labelpad=1)
axB.tick_params(axis="x", labelsize=6, pad=1.5)
axB.axvline(nlog(0.05), color=MUTED, lw=0.5, ls=(0, (2, 1.5)))
for yi, (_, r) in zip(y, d.iterrows()):
    axB.text(3.08, yi, r.GeneRatio, va="center", ha="left", fontsize=6, color=INK, clip_on=False)
    axB.text(3.62, yi, r.genes.replace(",", ", "), va="center", ha="left", fontsize=6, color=INK2,
             style="italic", clip_on=False)
axB.text(3.08, y.max() + 0.7, "n/N", fontsize=6, color=INK2, va="bottom", fontweight="bold")
axB.text(3.62, y.max() + 0.7, "member genes", fontsize=6, color=INK2, va="bottom", fontweight="bold")
axB.set_xlabel("\u2212log$_{10}$ FDR (dashed line: FDR = 0.05)", fontsize=6, labelpad=1)

# ---------------- C ----------------------------------------------------------------
yC = yB + 0.30 + len(d) * RH + 0.36
plabel("C", 0.0, yC)
ptitle("GO-BP of the 109 core genes (top 20, recomputed)", 0.22, yC + 0.02)
d = C.reset_index(drop=True)
axC = axin(1.72, yC + 0.30, 1.30, len(d) * RH)
y = np.arange(len(d))[::-1]
axC.barh(y, nlog(d["p.adjust"]), height=0.66, color=SHARED, lw=0)
for yi, (_, r) in zip(y, d.iterrows()):
    axC.text(nlog(r["p.adjust"]) + 0.12, yi, r.GeneRatio, va="center", fontsize=6, color=INK2)
axC.set_yticks(y)
axC.set_yticklabels([short(s, 40) for s in d.Description], fontsize=FS)
axC.tick_params(axis="y", length=0, pad=2)
axC.set_xlim(0, 9.0)
axC.set_ylim(-0.55, y.max() + 0.55)
axC.set_xlabel("−log$_{10}$ FDR", fontsize=6, labelpad=1)
axC.tick_params(axis="x", labelsize=6, pad=1.5)
fig.text(0.22 / FW, 1 - (yC + 0.30 + len(d) * RH + 0.30) / FH,
         "org.Hs.eg.db/GO.db 3.18; 105 of 107 mapped genes GO-BP-annotated;\nno GO-MF term at FDR < 0.05",
         fontsize=6, color=INK2, va="top")

# ---------------- D ----------------------------------------------------------------
xL = 4.95  # bar-axes left edge for D/E
plabel("D", 3.55, yC)
ptitle("KEGG, HPC sex-specific genes", 3.77, yC + 0.02)
d = D.reset_index(drop=True)
yv = np.array([3.6, 1.0, 0.0])
axD = axin(xL, yC + 0.30, 0.62, 4.6 * RH)
cols = [MEN if s.startswith("HPC male") else WOMEN for s in d.sex]
axD.barh(yv, nlog(d["p.adjust"]), height=0.66, color=cols, lw=0)
axD.set_yticks(yv)
axD.set_yticklabels([short(s, 34) for s in d.Description], fontsize=FS)
axD.tick_params(axis="y", length=0, pad=2)
axD.set_xlim(0, 2.0)
axD.set_ylim(-0.55, 4.2)
axD.axvline(nlog(0.05), color=MUTED, lw=0.5, ls=(0, (2, 1.5)))
for yi, (_, r) in zip(yv, d.iterrows()):
    axD.text(2.08, yi, f"{r.GeneRatio}  {r.genes}", va="center", fontsize=6, color=INK2, style="italic", clip_on=False)
tr = axD.get_yaxis_transform()
axD.text(-1.95, 4.45, "male-specific (9 genes)", transform=tr, fontsize=6, color=INK, fontweight="bold", va="center")
axD.text(-1.95, 1.85, "female-specific (45 genes)", transform=tr, fontsize=6, color=INK, fontweight="bold", va="center")
axD.set_xlabel("−log$_{10}$ FDR", fontsize=6, labelpad=1)
axD.tick_params(axis="x", labelsize=6, pad=1.5)
# small colour keys next to the group headers
for yy, cc in [(4.45, MEN), (1.85, WOMEN)]:
    axD.scatter([-2.05], [yy], transform=tr, marker="s", s=10, color=cc, clip_on=False)
axD.set_ylim(-0.55, 4.2)

# ---------------- E ----------------------------------------------------------------
yE = yC + 0.30 + 4.6 * RH + 0.42
plabel("E", 3.55, yE)
ptitle("GO, sex-dimorphic genes (recomputed)", 3.77, yE + 0.02)
fig.text(3.77 / FW, 1 - (yE + 0.20) / FH, "HPC: SLC27A6, SLC26A3, INO80D (BTBD11 unmapped); PFC: ZBTB16.\n1\u20132 genes per term: restates single-gene annotations",
         fontsize=6, color=INK, va="top", style="italic")
d = Ed.reset_index(drop=True)
ys, yv, heads, prev = [], 0.0, [], None
for _, r in d.iloc[::-1].iterrows():
    key = (r.region, r.ont)
    if prev is not None and key != prev:
        heads.append((yv - 0.1, prev)); yv += 1.25
    ys.append(yv); yv += 1.0; prev = key
heads.append((yv - 0.1, prev))
ys = np.array(ys)[::-1]
axE = axin(xL, yE + 0.42, 0.62, (yv + 0.2) * RH)
axE.barh(ys, nlog(d["p.adjust"]), height=0.66, color=DIMORPHIC, lw=0)
axE.set_yticks(ys)
axE.set_yticklabels([short(s, 34) for s in d.Description], fontsize=FS)
axE.tick_params(axis="y", length=0, pad=2)
axE.set_xlim(0, 3.2)
axE.set_ylim(-0.55, yv + 0.3)
axE.axvline(nlog(0.05), color=MUTED, lw=0.5, ls=(0, (2, 1.5)))
for yi, (_, r) in zip(ys, d.iterrows()):
    axE.text(3.3, yi, f"{r.GeneRatio}  {r.geneID.replace('/', ', ')}", va="center", fontsize=6,
             color=INK2, style="italic", clip_on=False)
tr = axE.get_yaxis_transform()
names = {"HPC": "HPC", "PFC": "PFC"}
for yh, (reg, ont) in heads:
    axE.text(-1.95, yh + 0.05, f"{names[reg]} · GO-{ont}", transform=tr, fontsize=6, color=INK,
             fontweight="bold", va="center")
axE.set_xlabel("−log$_{10}$ FDR", fontsize=6, labelpad=1)
axE.tick_params(axis="x", labelsize=6, pad=1.5)

plt.rcParams["savefig.bbox"] = None  # layout is set explicitly in inches; keep exact 7.2-in width
save(fig, OUT + "Figure3")
print("A rows", len(A), "| C FDR range", C["p.adjust"].min(), C["p.adjust"].max())
