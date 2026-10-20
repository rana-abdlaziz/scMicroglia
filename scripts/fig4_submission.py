"""Figure 4 - prior-knowledge regulatory networks.

  A  MetaCore pathway map (existing image, metacore/pathway_white_background.png).
  B  Curated MetaCore transcription-regulation network of the 109 core genes (+ JunB, JunD).
  C  TF-to-TF subnetwork among 14 AP-1/EGR/NR4A/other regulators (72 interaction records).
  D  MetaCore network expansions of HPC male- and female-specific genes (existing images, cropped).
  E  Genetic-algorithm Boolean networks (adjacency.txt) for HPC male- and female-specific genes.
"""
import sys
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.lines import Line2D
from PIL import Image

sys.path.insert(0, "<BUILD_DIR>/figures/scripts")
from figstyle import *  # noqa

apply_style()
T = "<BUILD_DIR>/agents_out/tables/"
U = "<PROJECT_ROOT>/"
OUT = "<BUILD_DIR>/figures/out/"
SD = OUT + "source_data/"

ACT, INH, UNS, MIX = UP, DOWN, "#c9c7c1", DIMORPHIC
ECOL = {"Activation": ACT, "Inhibition": INH, "Unspecified": UNS, "Both": MIX}
NODE_TF, NODE_T = "#8f8d88", "#e6e5e1"
FEEDBACK = {"NR4A1", "NR4A3", "BCL6", "BCL3", "SOCS3", "DUSP1", "ZFP36L1"}
NONCORE = {"JUNB", "JUND", "KDM3A"}
ALIAS = {"SMARCA4": "SMARCA4\n(BRG1)", "FOSL2": "FOSL2\n(Fra-2)", "KDM3A": "KDM3A\n(JMJD1A)"}

# =========================================================================== data
core_df = pd.read_csv(T + "C_109_shared_both_regions.csv")
core = set(core_df.Gene_Symbols)
E = pd.read_csv(T + "G_shared_network_edges_curated.csv")
E = E[E.effect != "Technical"]
rows = []
for _, r in E.iterrows():
    for s in str(r.src_sym).split("/"):
        rows.append(dict(n=r.n, src=s, tgt=r.tgt_sym, effect=r.effect, mech=r.mech))
X = pd.DataFrame(rows)

# verify regulator numbers quoted in the text
st = pd.read_csv(T + "G_shared_network_regulator_stats.csv").set_index("src")
quoted = {"JUND": 80, "JUNB": 68, "FOSL2": 58, "JUN": 40, "SMARCA4": 32, "BCL6": 15, "EGR1": 14, "MKL1": 8}
recount = {k: X[(X.src == k) & X.tgt.isin(core)].tgt.nunique() for k in quoted}
ver = pd.DataFrame({"quoted": quoted, "regulator_stats": st.loc[list(quoted), "n_unique_targets109"], "recount": recount})
print(ver)
ver.to_csv(SD + "Fig4B_regulator_target_counts_check.csv")


def collapse(df):
    g = df.groupby(["src", "tgt"]).agg(effects=("effect", lambda v: set(v)), n_records=("n", "nunique")).reset_index()

    def eff(s):
        if "Activation" in s and "Inhibition" in s:
            return "Both"
        if "Activation" in s:
            return "Activation"
        if "Inhibition" in s:
            return "Inhibition"
        return "Unspecified"
    g["effect"] = g.effects.map(eff)
    g["effects"] = g.effects.map(lambda s: ";".join(sorted(s)))
    return g


# ---- B: core network
inside = X[X.src.isin(core) & X.tgt.isin(core)]
nodes_core = set(inside.src) | set(inside.tgt)
jj = X[X.src.isin(["JUNB", "JUND"]) & X.tgt.isin(core)]
nodes_core = nodes_core | set(jj.tgt)  # core genes reached only by JunB/JunD (e.g. SRGN) are kept
EB = collapse(pd.concat([inside, jj]))
EB.to_csv(SD + "Fig4B_core_network_edges.csv", index=False)
nB_records_inside = inside.n.nunique()
selfB = EB[EB.src == EB.tgt]
EBd = EB[EB.src != EB.tgt]
GB = nx.DiGraph()
GB.add_nodes_from(sorted(nodes_core | {"JUNB", "JUND"}))
for _, r in EBd.iterrows():
    GB.add_edge(r.src, r.tgt, effect=r.effect)
print(f"B: {GB.number_of_nodes()} nodes, {GB.number_of_edges()} drawn pairs "
      f"({nB_records_inside} core-core records; JunB/JunD pairs {len(jj[['src','tgt']].drop_duplicates())}; "
      f"self-loops omitted {list(selfB.src)}); effects {EBd.effect.value_counts().to_dict()}")
pd.DataFrame({"node": list(GB.nodes), "out_degree": [GB.out_degree(n) for n in GB.nodes],
              "in_degree": [GB.in_degree(n) for n in GB.nodes],
              "in_core": [n in core for n in GB.nodes], "feedback": [n in FEEDBACK for n in GB.nodes]}) \
    .to_csv(SD + "Fig4B_core_network_nodes.csv", index=False)

# ---- C: TF-TF subnetwork (14 regulators; records counted as in G_shared_network_analyze.py)
R14 = ["JUN", "JUNB", "JUND", "FOSL2", "EGR1", "EGR3", "NR4A1", "NR4A3", "BCL6", "BCL3", "SMARCA4", "MKL1",
       "ZFP36L1", "KDM3A"]
rec = E[E.src_sym.map(lambda x: any(p in R14 for p in str(x).split("/"))) &
        E.tgt_sym.map(lambda x: any(p in R14 for p in str(x).split("/")))]
n_tftf_records = len(rec)
XC = X[X.src.isin(R14) & X.tgt.isin(R14)]
EC = collapse(XC)
EC.to_csv(SD + "Fig4C_TF_TF_edges.csv", index=False)
print(f"C: {n_tftf_records} TF->TF interaction records; {len(EC)} unique directed pairs after splitting complexes "
      f"(incl. {int((EC.src == EC.tgt).sum())} self-loops); effects {EC.effect.value_counts().to_dict()}")

# ---- E: GA Boolean networks


def read_ga(sex):
    d = U + f"metacore/{sex}/"
    adj = pd.read_csv(d + "adjacency.txt", sep="\t", index_col=0)
    expr = pd.read_csv(d + "Expression.txt", sep="\t", header=None, names=["gene", "sign"])
    expr["gene"] = expr.gene.str.upper()  # C5orf15 vs C5ORF15
    expr = expr.set_index("gene").sign
    diff = pd.read_csv(d + "DifferentialNetworkGenerator_Output.txt", sep="\t", header=None,
                       names=["src", "type", "tgt", "code"])
    filt = pd.read_csv(d + "FilteredMappedInteractionsNoDuplicates.txt", sep="\t", header=None, names=["src", "type", "tgt"])
    code2type = {1: "-->", -1: "--|", 2: "--?"}
    ed = [(s, t, int(adj.loc[s, t])) for s in adj.index for t in adj.columns if adj.loc[s, t] != 0]
    ed = pd.DataFrame(ed, columns=["src", "tgt", "code"])
    # verify adjacency codes against the filtered MetaCore interactions
    m = ed.merge(filt, on=["src", "tgt"], how="outer", indicator=True)
    ok = (m._merge == "both").all() and (m.code.map(code2type) == m.type).all()
    print(f"E {sex}: {adj.shape[0]} nodes, {len(ed)} edges ({(ed.src == ed.tgt).sum()} self-loops); "
          f"adjacency codes match FilteredMappedInteractions: {ok}; differential: {diff[['src','type','tgt']].values.tolist()}")
    ed["differential"] = [((diff.src == s) & (diff.tgt == t)).any() for s, t in zip(ed.src, ed.tgt)]
    ed["diff_type"] = [diff.type[(diff.src == s) & (diff.tgt == t)].iloc[0] if f else "" for s, t, f in
                       zip(ed.src, ed.tgt, ed.differential)]
    ed["expr_src"] = ed.src.map(expr)
    ed["expr_tgt"] = ed.tgt.map(expr)
    ed.to_csv(SD + f"Fig4E_GA_{sex}_edges.csv", index=False)
    return list(adj.index), ed, expr


ga = {s: read_ga(s) for s in ["male", "female"]}

# =========================================================================== drawing helpers


def draw_edges(ax, pos, edges, node_r, lw=0.5, alpha=1.0, rad=0.08, zorder=1, shrink_pt=None):
    """edges: iterable of (u, v, effect, highlight). node_r: dict node -> radius (points) for shrinking."""
    for u, v, eff, hl in edges:
        col = ECOL[eff]
        style = "-|>" if eff in ("Activation", "Both") else ("-[" if eff == "Inhibition" else "-")
        if eff == "Inhibition":
            style = "-[,widthB=2.2,lengthB=0"
        ms = 4 if eff in ("Activation", "Both") else 3
        if u == v:
            continue
        sa = node_r.get(u, 3) + 0.5
        sb = node_r.get(v, 3) + 1.0
        if hl:
            ax.add_patch(FancyArrowPatch(pos[u], pos[v], arrowstyle="-", connectionstyle=f"arc3,rad={rad}",
                                         lw=lw + 3.2, color=HIGHLIGHT, alpha=0.55, shrinkA=sa, shrinkB=sb,
                                         zorder=zorder - 0.5))
        ax.add_patch(FancyArrowPatch(pos[u], pos[v], arrowstyle=style, mutation_scale=ms,
                                     connectionstyle=f"arc3,rad={rad}", lw=lw, color=col, alpha=alpha,
                                     shrinkA=sa, shrinkB=sb, zorder=zorder))


def self_loop(ax, xy, r_pt, eff, ang=90, size=0.07, hl=False, lw=0.6):
    """small self-loop drawn as an arc above the node (data units)."""
    x, y = xy
    a = np.deg2rad(ang)
    cx, cy = x + np.cos(a) * size, y + np.sin(a) * size
    if hl:
        ax.add_patch(Circle((cx, cy), size * 0.55, fill=False, lw=lw + 3, color=HIGHLIGHT, alpha=0.55, zorder=0.5))
    ax.add_patch(Circle((cx, cy), size * 0.55, fill=False, lw=lw, color=ECOL[eff], zorder=0.8))


# =========================================================================== figure
FW, FH = W2, 9.0
fig = plt.figure(figsize=(FW, FH))


def axin(x, y, w, h):
    return fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH])


def plabel(letter, x, y):
    fig.text(x / FW, 1 - y / FH, letter, fontsize=11, fontweight="bold", va="top", ha="left", color=INK)


def ptitle(text, x, y, size=8, **kw):
    fig.text(x / FW, 1 - y / FH, text, fontsize=size, fontweight="bold", va="top", ha="left", color=INK, **kw)


# ---------------- A: simplified MetaCore map (redrawn; DEG nodes coloured by measured log2FC) -------------
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.colors import Normalize
plabel("A", 0.0, 0.02)
ptitle("MetaCore map: substance P-mediated inflammation", 0.22, 0.04)
fig.text(0.22 / FW, 1 - 0.22 / FH, "simplified redraw; DEG nodes show measured log₂FC", fontsize=6, color=INK2, va="top")
S2 = pd.read_csv(r"<PROJECT_ROOT>/Supplementary_Tables/csv/S2.csv")
S2f = S2[S2.list == "filtered"]
def lfc(gene):
    out = {}
    for reg in ["HPC", "PFC"]:
        r = S2f[(S2f.region == reg) & (S2f.gene == gene)]
        out[reg] = (r.log2FC_men.iloc[0], r.log2FC_women.iloc[0]) if len(r) else (np.nan, np.nan)
    return out
AW, AH = 3.22, 3.20
axA = axin(0.02, 0.30, AW, AH)
axA.set_xlim(0, AW); axA.set_ylim(0, AH); axA.axis("off")
cnorm = Normalize(-4, 4)
def heat(x, y, gene, cw=0.215, ch=0.115):
    v = lfc(gene)
    for i, reg in enumerate(["HPC", "PFC"]):
        for j in range(2):
            val = v[reg][j]
            xx, yy = x + j * cw, y - (i + 1) * ch
            if np.isnan(val):
                axA.add_patch(Rectangle((xx, yy), cw, ch, fc="#f4f3f0", ec="white", lw=0.5, zorder=4))
                axA.text(xx + cw / 2, yy + ch / 2, "–", fontsize=5.5, ha="center", va="center", color=MUTED, zorder=5)
            else:
                c = DIVERGING(cnorm(np.clip(val, -4, 4)))
                axA.add_patch(Rectangle((xx, yy), cw, ch, fc=c, ec="white", lw=0.5, zorder=4))
                axA.text(xx + cw / 2, yy + ch / 2 - 0.004, f"{val:.1f}", fontsize=5.5, ha="center", va="center",
                         color="white" if abs(val) > 2.2 else INK, zorder=5)
def node(x, y, label, gene=None, w=0.56, h=0.44):
    measured = gene is not None and not all(np.isnan(t) for p in lfc(gene).values() for t in p)
    axA.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.04",
                                 fc="white", ec=INK if measured else MUTED, lw=0.8 if measured else 0.5,
                                 ls="-" if measured else (0, (2, 1.2)), zorder=3))
    if measured:
        axA.text(x, y + h / 2 - 0.035, label, fontsize=6, ha="center", va="top", color=INK, fontweight="bold", zorder=5)
        heat(x - 0.215, y + h / 2 - 0.165, gene)
    else:
        axA.text(x, y, label, fontsize=6, ha="center", va="center", color=INK2, zorder=5, linespacing=1.25)
# upstream signalling (map context)
axA.add_patch(FancyBboxPatch((0.18, 2.72), 2.86, 0.38, boxstyle="round,pad=0,rounding_size=0.05",
                             fc="#f7f6f3", ec=MUTED, lw=0.5, zorder=2))
axA.text(1.61, 2.91, "Substance P → NK1 receptor → Gαq / PLC-β → PKC, Ca²⁺\n→ ERK1/2, JNK, p38 MAPK and IKK (map context)",
         fontsize=6, ha="center", va="center", color=INK2, linespacing=1.15)
# TFs
xj, xn, yt = 0.95, 2.27, 2.25
for xx in (xj, xn):
    axA.annotate("", xy=(xx, yt + 0.23), xytext=(xx, 2.72), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=INK2, mutation_scale=6))
node(xj, yt, "c-Jun (JUN)", "JUN", w=0.70)
node(xn, yt, "NF-κB\np50/p65", None, w=0.70)
# transcriptional-regulation bus
ybus = 1.86
for xx in (xj, xn):
    axA.plot([xx, xx], [yt - 0.22, ybus], color=ACT, lw=0.8, zorder=1)
axA.plot([0.32, 3.17], [ybus, ybus], color=ACT, lw=0.8, zorder=1)
axA.text(1.61, ybus + 0.03, "transcriptional regulation", fontsize=5.5, ha="center", va="bottom", color=ACT,
         bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none"), zorder=2)
row1 = [("CCL3L1", "CCL3L1"), ("MIP-1β/CCL4", "CCL4"), ("CCL13", None), ("COX-2\n(PTGS2)", None), ("TNF-α", None)]
row2 = [("CCL2", "CCL2"), ("IL-8\n(CXCL8)", None), ("IL-6", None), ("GRO-2\n(CXCL2)", None), ("GRO-3\n(CXCL3)", None)]
xs1 = [0.32 + i * 0.645 for i in range(5)]
y1, y2 = 1.40, 0.80
for (lab, g), xx in zip(row1, xs1):
    axA.annotate("", xy=(xx, y1 + 0.225), xytext=(xx, ybus), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=ACT, mutation_scale=5))
    node(xx, y1, lab, g)
ybus2 = y2 + 0.33
axA.plot([3.17, 3.17], [ybus, ybus2], color=ACT, lw=0.8, zorder=1)
axA.plot([xs1[0], 3.17], [ybus2, ybus2], color=ACT, lw=0.8, zorder=1)
axA.plot([xs1[-1] + 0.28, 3.17], [ybus, ybus], color=ACT, lw=0.8, zorder=1)
for (lab, g), xx in zip(row2, xs1):
    axA.annotate("", xy=(xx, y2 + 0.225), xytext=(xx, ybus2), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=ACT, mutation_scale=5))
    node(xx, y2, lab, g)
# output
axA.plot([0.32, 2.90], [0.50, 0.50], color=MUTED, lw=0.6)
for xx in xs1:
    axA.plot([xx, xx], [y2 - 0.22, 0.50], color=MUTED, lw=0.6)
axA.annotate("", xy=(1.61, 0.36), xytext=(1.61, 0.50), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=MUTED, mutation_scale=6))
axA.add_patch(FancyBboxPatch((0.95, 0.14), 1.32, 0.22, boxstyle="round,pad=0,rounding_size=0.04", fc="#f7f6f3", ec=MUTED, lw=0.5))
axA.text(1.61, 0.25, "Inflammatory response", fontsize=6, ha="center", va="center", color=INK2)
# key for the mini-heatmaps
kx, ky = 2.45, 0.30
axA.text(kx - 0.02, ky + 0.08, "M   F", fontsize=5.5, ha="left", va="bottom", color=INK2)
for i, reg in enumerate(["HPC", "PFC"]):
    for j in range(2):
        axA.add_patch(Rectangle((kx + j * 0.12, ky - (i + 1) * 0.085 + 0.08), 0.12, 0.085, fc="#f4f3f0", ec=INK2, lw=0.4))
    axA.text(kx + 0.27, ky - (i + 1) * 0.085 + 0.08 + 0.042, reg, fontsize=5.5, va="center", color=INK2)
yA_end = 0.30 + AH

# ---------------- B: core network -------------------------------------------------------
plabel("B", 3.25, 0.02)
ptitle("Curated MetaCore network of the 109 core genes", 3.47, 0.04)
outdeg = dict(GB.out_degree())
regs = [n for n in GB.nodes if outdeg[n] > 0]
inner = ["JUND", "JUNB", "FOSL2", "JUN", "SMARCA4", "BCL6", "EGR1", "MKL1"]
TFset = set(R14) | {"KDM2A", "ZHX2", "POU2F2"}
middle = [n for n in regs if n not in inner and n in TFset]
outer_src = [n for n in regs if n not in inner and n not in TFset]
targets = [n for n in GB.nodes if outdeg[n] == 0]


def dominant(n):
    preds = [p for p in GB.predecessors(n)]
    for p in inner:
        if p in preds:
            return inner.index(p)
    return len(inner)


pos = {}
for i, n in enumerate(inner):
    a = 2 * np.pi * i / len(inner) + 0.2
    pos[n] = (0.30 * np.cos(a), 0.30 * np.sin(a))
middle = sorted(middle, key=lambda n: (dominant(n), n))
for i, n in enumerate(middle):
    a = 2 * np.pi * i / len(middle) + 0.15
    pos[n] = (0.62 * np.cos(a), 0.62 * np.sin(a))
outer = sorted(targets + outer_src, key=lambda n: (min([inner.index(p) for p in GB.predecessors(n) if p in inner] +
                                                       [len(inner)]), sorted(GB.predecessors(n))[:1], n))
for i, n in enumerate(outer):
    a = 2 * np.pi * i / len(outer)
    pos[n] = (np.cos(a), np.sin(a))

axB = axin(3.42, 0.40, 3.45, 3.12)
axB.set_xlim(-1.42, 1.42)
axB.set_ylim(-1.48, 1.48)
axB.set_aspect("equal")
axB.axis("off")
size_pt = {n: 2.2 + 1.35 * np.sqrt(outdeg[n]) for n in GB.nodes}  # marker radius (pt)
order = {"Unspecified": 0, "Activation": 1, "Inhibition": 2, "Both": 2}
edges = sorted([(u, v, d["effect"], False) for u, v, d in GB.edges(data=True)], key=lambda e: order[e[2]])
for e_ in edges:
    draw_edges(axB, pos, [e_], size_pt, lw=0.35 if e_[2] == "Unspecified" else 0.5,
               alpha=0.55 if e_[2] == "Unspecified" else 0.85, rad=0.06, zorder=1 + order[e_[2]])
LAB = set(inner) | set(middle) | FEEDBACK | {"CCL2", "CCL3L1", "IL1B", "SOCS3", "BCL2", "CXCR4", "CD83", "SRGN", "ACSL1"}
for n in GB.nodes:
    x, y = pos[n]
    r = size_pt[n]
    fc = NODE_TF if n in TFset else NODE_T
    ec, lw, ls = INK2, 0.4, "-"
    if n in FEEDBACK:
        ec, lw = HIGHLIGHT, 1.4
    if n in NONCORE:
        fc, ec, lw, ls = "white", INK, 0.8, (0, (2, 1))
    axB.scatter([x], [y], s=(2 * r) ** 2, color=fc, edgecolor=ec, linewidth=lw, linestyle=ls, zorder=5)
    if n in LAB:
        if n in inner:
            axB.text(x, y, ALIAS.get(n, n), fontsize=6, ha="center", va="center", zorder=6, color=INK,
                     fontweight="bold", linespacing=0.9,
                     bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
        else:
            rr = np.hypot(x, y)
            ux, uy = x / rr, y / rr
            off = 0.075 if n in outer else 0.075
            axB.text(x + ux * off, y + uy * off, ALIAS.get(n, n).replace("\n", " "), fontsize=6,
                     ha="left" if ux > 0.25 else ("right" if ux < -0.25 else "center"),
                     va="bottom" if uy > 0.25 else ("top" if uy < -0.25 else "center"), zorder=6, color=INK,
                     bbox=dict(boxstyle="round,pad=0.05", fc="white", ec="none", alpha=0.7))
# legend B
lg = [Line2D([0], [0], color=ACT, lw=1, label="activation"), Line2D([0], [0], color=INH, lw=1, label="inhibition"),
      Line2D([0], [0], color=UNS, lw=1, label="unspecified"), Line2D([0], [0], color=MIX, lw=1, label="both signs"),
      Line2D([0], [0], marker="o", ls="", ms=5, mfc=NODE_TF, mec=INK2, mew=0.4, label="TF / regulator"),
      Line2D([0], [0], marker="o", ls="", ms=5, mfc="white", mec=HIGHLIGHT, mew=1.4, label="feedback regulator"),
      Line2D([0], [0], marker="o", ls="", ms=5, mfc="white", mec=INK, mew=0.8, label="not a core gene")]
axB.legend(handles=lg, loc="upper left", bbox_to_anchor=(-0.04, 1.03), fontsize=6, handlelength=1.4,
           handletextpad=0.4, borderaxespad=0, labelspacing=0.25)
axB.text(1.04, 0.07, f"{GB.number_of_nodes()} nodes, {GB.number_of_edges()} directed pairs\n(1 self-loop not drawn); node size ∝ out-degree",
         transform=axB.transAxes, ha="right", va="bottom", fontsize=6, color=INK2)

# ---------------- C: TF-TF ---------------------------------------------------------------
yC = max(yA_end, 0.40 + 3.12) + 0.10
plabel("C", 0.0, yC)
ptitle(f"TF–TF edges among 14 regulators", 0.22, yC + 0.02)
fig.text(0.22 / FW, 1 - (yC + 0.20) / FH, f"{n_tftf_records} interaction records, {len(EC)} directed pairs",
         fontsize=6, color=INK2, va="top")
famorder = ["JUN", "JUNB", "JUND", "FOSL2", "BCL3", "BCL6", "ZFP36L1", "NR4A1", "NR4A3", "EGR1", "EGR3",
            "MKL1", "SMARCA4", "KDM3A"]
GC = nx.DiGraph()
GC.add_nodes_from(famorder)
for _, r in EC.iterrows():
    GC.add_edge(r.src, r.tgt, effect=r.effect)
posC = {n: (np.cos(np.pi / 2 + 2 * np.pi * i / 14), np.sin(np.pi / 2 + 2 * np.pi * i / 14)) for i, n in enumerate(famorder)}
axC = axin(0.40, yC + 0.42, 2.15, 2.15)
axC.set_xlim(-1.38, 1.38); axC.set_ylim(-1.38, 1.38); axC.set_aspect("equal"); axC.axis("off")
odC = dict(GC.out_degree())
rC = {n: 3.2 + 1.0 * np.sqrt(odC[n]) for n in GC.nodes}
edC = sorted([(u, v, d["effect"], False) for u, v, d in GC.edges(data=True) if u != v], key=lambda e: order[e[2]])
for e_ in edC:
    draw_edges(axC, posC, [e_], rC, lw=0.45 if e_[2] == "Unspecified" else 0.6,
               alpha=0.7 if e_[2] == "Unspecified" else 0.9, rad=0.12, zorder=1 + order[e_[2]])
for u, v, d in GC.edges(data=True):
    if u == v:
        x, y = posC[u]
        self_loop(axC, (x, y), rC[u], d["effect"], ang=np.degrees(np.arctan2(y, x)), size=0.11)
for n in GC.nodes:
    x, y = posC[n]
    fc, ec, lw, ls = NODE_TF, INK2, 0.4, "-"
    if n in FEEDBACK:
        ec, lw = HIGHLIGHT, 1.4
    if n in NONCORE:
        fc, ec, lw, ls = "white", INK, 0.8, (0, (2, 1))
    axC.scatter([x], [y], s=(2 * rC[n]) ** 2, color=fc, edgecolor=ec, linewidth=lw, linestyle=ls, zorder=5)
    ux, uy = x, y
    axC.text(x * 1.2, y * 1.2, ALIAS.get(n, n).replace("\n", " ") if n != "SMARCA4" else "SMARCA4\n(BRG1)",
             fontsize=6, ha="left" if ux > 0.3 else ("right" if ux < -0.3 else "center"),
             va="bottom" if uy > 0.3 else ("top" if uy < -0.3 else "center"), color=INK, zorder=6, linespacing=0.9)

# ---------------- E: GA networks ----------------------------------------------------------
xE = 0.0
yE = yC + 2.56
plabel("E", xE, yE)
ptitle("Genetic-algorithm Boolean networks, HPC sex-specific genes", xE + 0.22, yE + 0.02)

CODE = {1: "Activation", -1: "Inhibition", 2: "Unspecified"}
FILL = {1: "#f2a3a2", -1: "#86b6ef"}


def ga_node(ax, xy, n, expr, r, fs=6, label_xy=None, ha="center", va="center", rot=0):
    s = expr.get(n, np.nan)
    if np.isnan(s):
        fc, ec, ls, lw = "white", INK, (0, (2, 1)), 0.8
    else:
        fc, ec, ls, lw = FILL[int(s)], INK2, "-", 0.4
    ax.scatter([xy[0]], [xy[1]], s=(2 * r) ** 2, color=fc, edgecolor=ec, linewidth=lw, linestyle=ls, zorder=5)
    lx, ly_ = label_xy if label_xy is not None else xy
    ax.text(lx, ly_, {"C5ORF15": "C5orf15"}.get(n, n), fontsize=fs, ha=ha, va=va, zorder=6, rotation=rot, rotation_mode="anchor", color=INK)


# male: hand layout (8 nodes)
nodes_m, ed_m, ex_m = ga["male"]
pm = {"RUNX1": (-0.25, 0.2), "FOXP1": (0.65, 0.2), "ESRRB": (0.2, -0.45), "PRDM11": (0.2, 0.85),
      "CXCL12": (1.05, -0.5), "TATDN3": (1.05, 0.85), "DHRS9": (-0.95, 0.85), "PFKFB2": (-0.95, -0.45)}
axEm = axin(xE + 0.10, yE + 0.30, 1.85, 2.0)
axEm.set_xlim(-1.3, 1.4); axEm.set_ylim(-1.45, 1.48); axEm.set_aspect("equal"); axEm.axis("off")
axEm.text(-1.3, 1.48, "Male (8 nodes, 15 edges; best fitness 6.0)", fontsize=6.5, fontweight="bold", va="top", color=INK)
rm = {n: 13 for n in nodes_m}
for _, r in ed_m.iterrows():
    eff = CODE[r.code]
    if r.differential and r.diff_type == "--|":
        eff = "Inhibition"  # GA-resolved sign of the unspecified database edge
    if r.src == r.tgt:
        x, y = pm[r.src]
        self_loop(axEm, (x, y), 10, eff, ang={"RUNX1": 200, "FOXP1": 20, "ESRRB": 225}[r.src], size=0.3)
        continue
    rad = 0.18 if ((ed_m.src == r.tgt) & (ed_m.tgt == r.src)).any() else 0.05
    draw_edges(axEm, pm, [(r.src, r.tgt, eff, r.differential)], rm, lw=0.7, rad=rad, zorder=2)
for n in nodes_m:
    ga_node(axEm, pm[n], n, ex_m, 13, fs=5.6)
axEm.annotate("differential:\nESRRB —| CXCL12\n(unspecified in DB;\nGA-resolved)", xy=(0.62, -0.48), xytext=(-1.3, -1.45),
              fontsize=6, color=INK2, va="bottom", ha="left",
              arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED, shrinkB=2))

# female: RUNX1 star
nodes_f, ed_f, ex_f = ga["female"]
axEf = axin(4.70, yE + 0.22, 2.30, 2.30)
axEf.set_xlim(-1.62, 1.62); axEf.set_ylim(-1.62, 1.62); axEf.set_aspect("equal"); axEf.axis("off")
axEf.text(-1.62, 1.62, "Female (36 nodes, 39 edges; best fitness 4.0)", fontsize=6.5, fontweight="bold", va="top", color=INK)
others = [n for n in nodes_f if n != "RUNX1"]
# place SMC1A next to OTUD1 and FZR1 near the top so the extra edges stay short
pref = ["OTUD1", "SMC1A"]
others = pref + [n for n in others if n not in pref]
pf = {"RUNX1": (0.0, 0.0)}
for i, n in enumerate(others):
    a = np.pi / 2 - 2 * np.pi * (i + 0.5) / len(others)
    pf[n] = (0.8 * np.cos(a), 0.8 * np.sin(a))
rf = {n: 3.4 for n in nodes_f}
rf["RUNX1"] = 8
for _, r in ed_f.sort_values("code", key=lambda c: c.map({2: 0, 1: 1, -1: 2})).iterrows():
    eff = CODE[r.code]
    if r.src == r.tgt:
        x, y = pf[r.src]
        ang = np.degrees(np.arctan2(y, x)) if r.src != "RUNX1" else 200
        self_loop(axEf, (x, y), rf[r.src], eff, ang=ang, size=0.1 if r.src != "RUNX1" else 0.14, hl=r.differential)
        continue
    rad = 0.25 if r.src == "FZR1" else (0.35 if r.src == "SMC1A" else 0.0)
    draw_edges(axEf, pf, [(r.src, r.tgt, eff, r.differential)], rf, lw=0.45 if eff == "Unspecified" else 0.7,
               rad=rad, zorder=2 + (eff != "Unspecified"))
for n in nodes_f:
    x, y = pf[n]
    if n == "RUNX1":
        ga_node(axEf, (x, y), n, ex_f, 8, fs=6, label_xy=(0, -0.2), va="top")
        continue
    a = np.degrees(np.arctan2(y, x))
    rot = a if -90 <= a <= 90 else a - 180
    ga_node(axEf, (x, y), n, ex_f, 3.4, fs=6, label_xy=(x * 1.11, y * 1.11),
            ha="left" if -90 <= a <= 90 else "right", va="center", rot=rot)
fig.text(4.78 / FW, 1 - (yE + 2.50) / FH,
         "differential: SMC1A → OTUD1; BCL10 self-loop", fontsize=6, color=INK2, va="top")

# shared legend for E
lgE = axin(2.30, yE + 1.00, 1.8, 0.45); lgE.axis("off")
hE = [Line2D([0], [0], marker="o", ls="", ms=5, mfc=FILL[1], mec=INK2, mew=0.4, label="input +1 (up)"),
      Line2D([0], [0], marker="o", ls="", ms=5, mfc=FILL[-1], mec=INK2, mew=0.4, label="input −1 (down)"),
      Line2D([0], [0], marker="o", ls="", ms=5, mfc="white", mec=INK, mew=0.8, label="not in the expression input"),
      Line2D([0], [0], color=HIGHLIGHT, lw=3, alpha=0.6, label="differential edge")]
lgE.legend(handles=hE, loc="upper left", fontsize=6, ncol=1, handlelength=1.4, handletextpad=0.4,
           labelspacing=0.2, borderaxespad=0)

# ---------------- D: MetaCore expansions ---------------------------------------------------
yD = yC
xD = 2.72
plabel("D", xD, yD)
ptitle("MetaCore network expansions, HPC sex-specific genes", xD + 0.22, yD + 0.02)
fig.text((xD + 0.22) / FW, 1 - (yD + 0.20 + 1.80) / FH, "node fill: HPC male log₂FC (pink up, blue down);\ndashed: added by MetaCore, not measured", fontsize=5.6, color=INK2, va="top")
# male expansion redrawn as a vector graph (edges transcribed from the MetaCore export)
fig.text((xD + 0.22) / FW, 1 - (yD + 0.22) / FH, "Male-specific (Oct-3/4–SMAD2/3 hub)", fontsize=6.5, color=INK, va="top", fontweight="bold")
PM = {"NODAL": (385, 70), "FGF4": (805, 105), "Esrrb": (1040, 85), "CDX2": (1175, 215), "MED12L": (80, 305),
      "SDF-1\n(CXCL12)": (760, 300), "DBL\n(MCF2)": (1330, 455), "SMAD3": (510, 485), "ALK-4": (225, 555),
      "Oct-3/4": (955, 540), "KCNH5": (1340, 800), "SMAD2": (590, 785), "PRDM11": (880, 965), "DHRS9": (1140, 905)}
PM = {k: (x / 1500 * 1.85, (1060 - y) / 1060 * 1.32) for k, (x, y) in PM.items()}
MALE_DIR = {"MED12L": -1, "SDF-1\n(CXCL12)": -1, "Esrrb": -1, "DBL\n(MCF2)": 1, "KCNH5": 1, "PRDM11": -1, "DHRS9": -1}
A_, I_, U_ = "Activation", "Inhibition", "Unspecified"
EDM = [("NODAL", "FGF4", A_), ("FGF4", "Esrrb", A_), ("FGF4", "CDX2", A_), ("NODAL", "ALK-4", A_),
       ("ALK-4", "SMAD3", A_), ("ALK-4", "SMAD2", A_), ("SMAD3", "Oct-3/4", A_), ("SMAD2", "Oct-3/4", A_),
       ("SMAD3", "MED12L", U_), ("SMAD3", "SDF-1\n(CXCL12)", U_), ("Oct-3/4", "SDF-1\n(CXCL12)", A_),
       ("Oct-3/4", "FGF4", A_), ("Oct-3/4", "Esrrb", A_), ("Oct-3/4", "CDX2", I_), ("CDX2", "Oct-3/4", I_),
       ("CDX2", "Esrrb", A_), ("Esrrb", "SDF-1\n(CXCL12)", U_), ("Oct-3/4", "DBL\n(MCF2)", U_),
       ("Oct-3/4", "KCNH5", I_), ("Oct-3/4", "PRDM11", U_), ("Oct-3/4", "DHRS9", U_)]
axDm = axin(xD + 0.08, yD + 0.40, 1.95, 1.56)
axDm.set_xlim(-0.06, 1.89); axDm.set_ylim(-0.14, 1.42); axDm.set_aspect("equal"); axDm.axis("off")
rDm = {n: 5.0 for n in PM}
for u, v, e in sorted(EDM, key=lambda t: order[t[2]]):
    bidir = any((a == v and b == u) for a, b, _ in EDM)
    draw_edges(axDm, PM, [(u, v, e, False)], rDm, lw=0.45 if e == U_ else 0.6, alpha=0.75 if e == U_ else 0.95,
               rad=0.18 if bidir else 0.04, zorder=2 + order[e])
self_loop(axDm, PM["Esrrb"], 5, U_, ang=60, size=0.07)
for n, (x, y) in PM.items():
    d = MALE_DIR.get(n)
    if d is None:
        fc, ec, ls, lw = "white", INK, (0, (2, 1)), 0.7
    else:
        fc, ec, ls, lw = FILL[d], INK2, "-", 0.4
    axDm.scatter([x], [y], s=(2 * 5.0) ** 2, color=fc, edgecolor=ec, linewidth=lw, linestyle=ls, zorder=6)
    LOFF = {"SDF-1\n(CXCL12)": (-0.07, 0.0, "right", "center"), "Esrrb": (0.0, 0.075, "center", "bottom"),
            "KCNH5": (0.0, -0.075, "center", "top"),
            "CDX2": (0.07, 0.0, "left", "center")}
    dx, dy, ha_, va_ = LOFF.get(n, (0.0, -0.075, "center", "top"))
    axDm.text(x + dx, y + dy, n, fontsize=5.6, ha=ha_, va=va_, color=INK, zorder=7, linespacing=0.9,
              bbox=dict(boxstyle="round,pad=0.04", fc="white", ec="none", alpha=0.8))
# female expansion: MetaCore raster with legible hub call-outs
fig.text((xD + 2.42) / FW, 1 - (yD + 0.22) / FH, "Female-specific (E2F1/E2F4–p53 hub)", fontsize=6.5, color=INK, va="top", fontweight="bold")
im = Image.open(U + "female_network.png")
crop = im.crop((1359, 0, 4453, 2140))
bg = Image.new("RGB", crop.size, "white"); bg.paste(crop, mask=crop.split()[3])
cw, ch = bg.size
Wf = 2.40; Hf = Wf * ch / cw
axDf = axin(xD + 2.06, yD + 0.40, Wf, Hf)
axDf.imshow(np.asarray(bg), interpolation="lanczos"); axDf.axis("off")
scf = 1500 / cw  # call-out coordinates were read on a 1500-px-wide rendering of this crop
CALL = [("E2F1/E2F4–p53\nhub", (595, 560), (1130, 440)), ("CDK1", (400, 322), (300, 230)),
        ("CDC25A", (483, 285), (560, 170)), ("p21", (535, 252), (720, 210)), ("CDK2", (395, 92), (230, 60)),
        ("Rb", (880, 752), (1010, 840)), ("HDAC1", (848, 990), (1000, 1000)), ("PPAR-γ", (1045, 690), (1150, 600)),
        ("E2F5", (660, 890), (480, 960))]
for lab, (px, py), (tx, ty) in CALL:
    axDf.annotate(lab, xy=(px / scf, py / scf), xytext=(tx / scf, ty / scf), fontsize=5.6, color=INK, ha="center",
                  va="center", fontweight="bold", linespacing=0.9, zorder=8,
                  bbox=dict(boxstyle="round,pad=0.1", fc="white", ec=INK2, lw=0.3, alpha=0.92),
                  arrowprops=dict(arrowstyle="-", lw=0.4, color=INK, shrinkA=0, shrinkB=1))
axDf.set_xlim(0, cw); axDf.set_ylim(ch, 0)
fig.text((xD + 2.42) / FW, 1 - (yD + 0.42 + Hf + 0.02) / FH, "MetaCore-added objects were not measured",
         fontsize=5.6, color=INK2, va="top")

r = fig.canvas.get_renderer()
bb = fig.get_tightbbox(r)
print("tight bbox (in):", bb.bounds)
for a in fig.axes:
    eb = a.get_tightbbox(r)
    if eb.x1 > 7.1 * fig.dpi or eb.x0 < 0: print("overflow axes", a.get_position().bounds, eb.x0 / fig.dpi, eb.x1 / fig.dpi)
for t in fig.texts:
    eb = t.get_window_extent(r)
    if eb.x1 > 7.2 * fig.dpi: print("overflow text", t.get_text()[:30])
plt.rcParams["savefig.bbox"] = None  # layout is set explicitly in inches; keep exact 7.2-in width
save(fig, OUT + "Figure4")
print("FH", FH, "yC", yC, "yE", yE)
