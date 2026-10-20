"""Supplementary Figures S1-S5 built from Supplementary_Tables/csv (S2-S8) of the scMicroglia project.
Run inside the analysis kernel after `S` (dict of CSV DataFrames) is loaded; uses figure-style helpers."""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
import os
apply_figure_style(frame="open", sizes=(8, 7, 6))
plt.rcParams["mathtext.default"] = "regular"
W = 7.2
RET, REM = "#4d4d4d", "#d95f02"          # retained vs removed by filter
WOM, MEN = "#c2410c", "#2563eb"           # women / men (not red-green)
ACT, INH, UNS = "#b2182b", "#2166ac", "#bdbdbd"
os.makedirs("out", exist_ok=True)

def neg(s):
    return s.replace("-", "\u2212")

# ---------------------------------------------------------------- S1: sensitivity to the PFC filter (163 genes)
s4 = S["S4"].copy()
s4["HPC_mean"] = s4[["log2FC_HPC_men", "log2FC_HPC_women"]].mean(axis=1)
s4["PFC_mean"] = s4[["log2FC_PFC_men", "log2FC_PFC_women"]].mean(axis=1)
ret = s4[s4.in_filtered_core_109 == "yes"]; rem = s4[s4.in_filtered_core_109 != "yes"]
assert len(ret) == 109 and len(rem) == 54, (len(ret), len(rem))
fig = plt.figure(figsize=(W, 3.4))
axA = fig.add_axes([0.08, 0.16, 0.40, 0.72]); axB = fig.add_axes([0.60, 0.16, 0.37, 0.72])
axA.axhline(0, color="#cccccc", lw=0.6, zorder=0); axA.axvline(0, color="#cccccc", lw=0.6, zorder=0)
axA.plot([-3, 10.5], [-3, 10.5], color="#999999", lw=0.6, ls="--", zorder=0)
axA.scatter(ret.HPC_mean, ret.PFC_mean, s=10, color=RET, alpha=0.8, lw=0, label=f"retained by PFC filter (filtered core, n = {len(ret)})")
axA.scatter(rem.HPC_mean, rem.PFC_mean, s=14, facecolor="none", edgecolor=REM, lw=0.8, label=f"removed by PFC filter (intermediate PFC effect, n = {len(rem)})")
for g, tx, ty in [("C22orf34", 2.2, -1.4), ("RUNX1", 5.2, 0.4), ("FOXP1", 5.2, 1.3), ("CDKN1A", 6.4, 2.6), ("RELB", 6.4, 3.6)]:
    r = s4[s4.gene == g].iloc[0]
    axA.annotate(g, (r.HPC_mean, r.PFC_mean), xytext=(tx, ty), fontsize=6, style="italic",
                 arrowprops=dict(arrowstyle="-", lw=0.5, color="#555555"), ha="left", va="center")
axA.set_xlabel("mean log$_2$FC in HPC (men, women)"); axA.set_ylabel("mean log$_2$FC in PFC (men, women)")
axA.set_title("Direction concordant for 162 of 163 genes\nwhen the PFC filter is relaxed", loc="left")
axA.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=6, handletextpad=0.4)
axA.margins(0.05); panel_letter(axA, "A", case="upper", dx=-0.16)
# B: PFC |log2FC| of removed genes per sex vs retained
vals = [np.abs(rem.log2FC_PFC_men), np.abs(rem.log2FC_PFC_women), np.abs(ret.log2FC_PFC_men), np.abs(ret.log2FC_PFC_women)]
labels = ["removed\nmen", "removed\nwomen", "retained\nmen", "retained\nwomen"]
rng = np.random.default_rng(1)
for i, (v, c) in enumerate(zip(vals, [REM, REM, RET, RET])):
    x = i + rng.uniform(-0.15, 0.15, len(v))
    axB.scatter(x, v, s=6, color=c, alpha=0.7, lw=0)
    axB.hlines(np.median(v), i - 0.25, i + 0.25, color="black", lw=1.4)
axB.axhspan(0.5, 1.5, color="#f0f0f0", zorder=0)
axB.text(3.45, 1.0, "intermediate\nstratum\n(0.5–1.5)", fontsize=6, va="center", ha="left", color="#555555")
axB.set_yscale("log"); axB.set_yticks([0.1, 0.3, 1, 3, 10]); axB.set_yticklabels(["0.1", "0.3", "1", "3", "10"])
axB.set_xticks(range(4)); axB.set_xticklabels(labels); axB.set_ylabel("|log$_2$FC| in PFC")
axB.set_title("Removed genes have intermediate\nPFC effects (0.5–1.5) by construction", loc="left"); axB.set_xlim(-0.5, 4.3)
panel_letter(axB, "B", case="upper", dx=-0.2)
fig.savefig("out/Supplementary_Figure_S1.png", dpi=600)
fig.savefig("out/Supplementary_Figure_S1.pdf")
plt.close(fig)

# ---------------------------------------------------------------- S2: overlap of DEGs with state-marker lists (circularity)
s2 = S["S2"].copy(); s1 = S["S1"]
fig = plt.figure(figsize=(W, 3.6))
axA = fig.add_axes([0.08, 0.20, 0.36, 0.66]); axB = fig.add_axes([0.58, 0.20, 0.39, 0.66])
sets = [("HPC", "filtered", "HPC\nsex-shared\n(253)"), ("PFC", "filtered", "PFC filtered\nsex-shared\n(1,128)"), ("PFC", "unfiltered", "PFC unfiltered\nsex-shared\n(1,935)")]
infl, homeo, other = [], [], []
for reg, lst, _ in sets:
    d = s2[(s2.region == reg) & (s2.list == lst) & (s2.category == "sex-shared")]
    infl.append((d.marker_of_contrast_group == "Inflammatory").sum()); homeo.append((d.marker_of_contrast_group == "Homeostatic").sum())
    other.append(d.state_marker_lists.notna().sum() - infl[-1] - homeo[-1])
x = np.arange(3)
axA.bar(x, infl, color="#b2182b", width=0.6, label="marker of MG2/MG8/MG10 (inflammatory, contrast group)")
axA.bar(x, homeo, bottom=infl, color="#2166ac", width=0.6, label="marker of MG0 (homeostatic, reference group)")
axA.bar(x, other, bottom=np.array(infl) + np.array(homeo), color="#bdbdbd", width=0.6, label="marker of another state list")
for i in range(3):
    axA.text(i, infl[i] + homeo[i] + other[i] + 1.5, f"{infl[i]} / {homeo[i]} / {other[i]}", ha="center", fontsize=6)
axA.set_xticks(x); axA.set_xticklabels([s[2] for s in sets]); axA.set_ylabel("sex-shared DEGs that are curated state markers (n)")
axA.set_title("Few shared DEGs are annotation markers", loc="left"); axA.set_ylim(0, max(np.array(infl) + np.array(homeo) + np.array(other)) * 1.35)
axA.legend(loc="upper left", fontsize=6, handlelength=1.0)
panel_letter(axA, "A", case="upper", dx=-0.2)
# B: the 9 female-biased HPC genes (>1 log2 unit) + 2 male-biased, with marker list
g = S["S5_G_sexdiff_gt1"]; g = g[g.region == "HPC"].sort_values("sexbias_absF_minus_absM", ascending=False)
mk = s2[(s2.region == "HPC")].set_index("gene").state_marker_lists
y = np.arange(len(g))[::-1]
cols = [WOM if v > 0 else MEN for v in g.sexbias_absF_minus_absM]
axB.hlines(y, 0, g.sexbias_absF_minus_absM, color=cols, lw=1.2)
axB.scatter(g.sexbias_absF_minus_absM, y, color=cols, s=22, zorder=3)
axB.axvline(0, color="#999999", lw=0.6)
axB.set_yticks(y); axB.set_yticklabels([f"{gn}" for gn in g.gene], style="italic")
for yy, gn in zip(y, g.gene):
    m = mk.get(gn, np.nan)
    axB.text(3.95, yy, "–" if pd.isna(m) else f"{m} list", fontsize=6, va="center", ha="left", color="#333333")
axB.text(3.95, len(g) - 0.45, "marker list", fontsize=6, ha="left", va="bottom", fontweight="bold")
axB.set_xlim(-2.3, 5.2); axB.set_xlabel(neg("|log$_2$FC$_F$| − |log$_2$FC$_M$| in HPC")); axB.set_ylim(-0.7, len(g) - 0.3)
axB.set_title("7 of 9 female-biased HPC genes are MG8/MG10 markers", loc="left")
axB.text(-2.2, -0.55, "larger in men ←", fontsize=6, color=MEN, ha="left"); axB.text(3.8, -0.55, "→ larger in women", fontsize=6, color=WOM, ha="right")
panel_letter(axB, "B", case="upper", dx=-0.24)
fig.savefig("out/Supplementary_Figure_S2.png", dpi=600); fig.savefig("out/Supplementary_Figure_S2.pdf"); plt.close(fig)

# ---------------------------------------------------------------- S3: a priori gene sets and per-gene sex bias
f = S["S5_F_program_genes"]; s5 = S["S5"]
fig = plt.figure(figsize=(W, 3.4))
axA = fig.add_axes([0.08, 0.20, 0.40, 0.64]); axB = fig.add_axes([0.60, 0.20, 0.37, 0.64])
groups = [("chemokine/cytokine", "HPC"), ("IEG/AP-1/EGR/NR4A", "HPC"), ("chemokine/cytokine", "PFC"), ("IEG/AP-1/EGR/NR4A", "PFC")]
rng = np.random.default_rng(2)
for i, (prog, reg) in enumerate(groups):
    d = f[(f.program == prog)][f"sexbias_{reg}_absF_minus_absM"].dropna()
    # restrict to the a priori sets named in Methods (10 chemokine genes; 6 IEG genes) when available in that region
    x = i + rng.uniform(-0.15, 0.15, len(d))
    c = WOM if np.median(d) > 0 else MEN
    axA.scatter(x, d, s=14, color="#4d4d4d", alpha=0.8, lw=0)
    axA.hlines(np.median(d), i - 0.28, i + 0.28, color="black", lw=1.6)
    key = ("chemokine_cytokine" if prog.startswith("chem") else "IEG_AP1_EGR_NR4A") + f"_{reg}"
    row = s5[s5.set_id.str.startswith(key)].iloc[0]
    axA.text(i, 4.55, f"n = {int(row.n)}\nmedian {neg(f'{row.median_absF_minus_absM:+.2f}')}\nP = {row.Wilcoxon_P:.2g}", ha="center", fontsize=6, va="top")
axA.axhline(0, color="#999999", lw=0.6)
axA.set_xticks(range(4)); axA.set_xticklabels(["chemokine/\ncytokine\nHPC", "IEG/AP-1/\nEGR/NR4A\nHPC", "chemokine/\ncytokine\nPFC", "IEG/AP-1/\nEGR/NR4A\nPFC"])
axA.set_ylabel(neg("|log$_2$FC$_F$| − |log$_2$FC$_M$|")); axA.set_ylim(-2.6, 4.6)
axA.set_title("A priori gene sets: only the PFC\nchemokine/cytokine set is sex-biased", loc="left")
axA.text(3.45, -2.5, "larger in men", fontsize=6, color=MEN, va="bottom", ha="right"); axA.text(3.45, 2.3, "larger in women", fontsize=6, color=WOM, va="top", ha="right")
panel_letter(axA, "A", case="upper", dx=-0.16)
# B: genes with sex difference > 1 log2 unit in both regions
g = S["S5_G_sexdiff_gt1"].copy(); g["lab"] = g.gene + " (" + g.region + ")"
g = g.sort_values(["region", "sexbias_absF_minus_absM"], ascending=[True, False])
y = np.arange(len(g))[::-1]; cols = [WOM if v > 0 else MEN for v in g.sexbias_absF_minus_absM]
axB.hlines(y, 0, g.sexbias_absF_minus_absM, color=cols, lw=1.2); axB.scatter(g.sexbias_absF_minus_absM, y, color=cols, s=20, zorder=3)
axB.axvline(0, color="#999999", lw=0.6)
axB.set_yticks(y); axB.set_yticklabels(g.lab, style="italic")
for yy, (fm, ff) in zip(y, zip(g.log2FC_men, g.log2FC_women)):
    axB.text(3.75, yy, f"{fm:.2f} / {ff:.2f}", fontsize=6, va="center", ha="left", color="#333333")
axB.text(3.75, len(g) - 0.45, "log$_2$FC men / women", fontsize=6, ha="left", va="bottom", fontweight="bold")
axB.set_xlim(-2.6, 5.6); axB.set_xlabel(neg("|log$_2$FC$_F$| − |log$_2$FC$_M$|")); axB.set_ylim(-0.7, len(g) - 0.3)
axB.set_title("Genes with a sex difference > 1 log$_2$ unit\n(high-effect shared genes)", loc="left", pad=10)
panel_letter(axB, "B", case="upper", dx=-0.3)
fig.savefig("out/Supplementary_Figure_S3.png", dpi=600); fig.savefig("out/Supplementary_Figure_S3.pdf"); plt.close(fig)

# ---------------------------------------------------------------- S4: gseGO of sex-shared genes (nominal P)
s6 = S["S6"]; gse = s6[(s6.analysis == "gseGO") & (s6.category_list == "sex-shared")].copy()
fig = plt.figure(figsize=(W, 5.0))
pos = {("HPC", "GO-BP"): [0.30, 0.57, 0.17, 0.34], ("HPC", "GO-MF"): [0.30, 0.09, 0.17, 0.34], ("PFC", "GO-BP"): [0.80, 0.57, 0.17, 0.34], ("PFC", "GO-MF"): [0.80, 0.09, 0.17, 0.34]}
letters = {("HPC", "GO-BP"): "A", ("HPC", "GO-MF"): "B", ("PFC", "GO-BP"): "C", ("PFC", "GO-MF"): "D"}
def short(s, n=46):
    return s if len(s) <= n else s[: n - 1] + "…"
for (reg, db), p in pos.items():
    d = gse[(gse.region == reg) & (gse.database == db)].sort_values("P_value").head(12).iloc[::-1]
    ax = fig.add_axes(p)
    cols = [WOM if v > 0 else MEN for v in d.NES]  # positive NES = larger female log2FC ranking (HPC/PFC shared ranked by female log2FC)
    ax.barh(np.arange(len(d)), d.NES, color=["#4d4d4d" if v > 0 else "#9ecae1" for v in d.NES], height=0.65)
    ax.set_yticks(np.arange(len(d))); ax.set_yticklabels([short(x) for x in d.Description])
    for i, (nes, pv) in enumerate(zip(d.NES, d.P_value)):
        ax.text(nes + (0.04 if nes > 0 else -0.04), i, f"P = {pv:.1e}".replace("e-0", "e−").replace("e-", "e−"), fontsize=5.5, va="center", ha="left" if nes > 0 else "right")
    ax.axvline(0, color="#999999", lw=0.6)
    n_all = len(gse[(gse.region == reg) & (gse.database == db)]); n_pos = int((gse[(gse.region == reg) & (gse.database == db)].NES > 0).sum())
    ax.set_title(f"{reg} sex-shared genes, {db}\n{n_all} terms, nominal P < 0.05 ({n_pos} with NES > 0)", loc="left", fontsize=7)
    ax.set_xlabel("normalized enrichment score (NES)"); ax.set_xlim(min(-0.3, d.NES.min() - 0.9), d.NES.max() + 1.3)
    fig.text(p[0] - 0.285, p[1] + p[3] + 0.03, letters[(reg, db)], fontweight="bold", fontsize=10, ha="left", va="bottom")
fig.savefig("out/Supplementary_Figure_S4.png", dpi=600); fig.savefig("out/Supplementary_Figure_S4.pdf"); plt.close(fig)

# ---------------------------------------------------------------- S5: evidence behind the MetaCore network and decoupleR edges
s7 = S["S7"]; rs = S["S7_B_regulator_stats"].sort_values("n_unique_core_targets", ascending=False).head(12)
s8 = S["S8"]
fig = plt.figure(figsize=(W, 4.6))
axA = fig.add_axes([0.13, 0.60, 0.36, 0.33]); axB = fig.add_axes([0.62, 0.60, 0.35, 0.33]); axC = fig.add_axes([0.13, 0.08, 0.84, 0.36])
y = np.arange(len(rs))[::-1]
nt = s7[(s7.effect != "Technical") & (s7.in_core_target == "yes")]
eff = nt.groupby(["source_HGNC", "effect"]).size().unstack(fill_value=0)
rows = []
for r in rs.regulator_HGNC:
    sub = nt[nt.source_HGNC.str.split("/").apply(lambda xs: r in xs)]
    c = sub.effect.value_counts()
    rows.append((r, c.get("Activation", 0), c.get("Inhibition", 0), c.get("Unspecified", 0)))
ev = pd.DataFrame(rows, columns=["reg", "act", "inh", "uns"])
axA.barh(y, ev.act, color=ACT, height=0.7, label="activation"); axA.barh(y, ev.inh, left=ev.act, color=INH, height=0.7, label="inhibition")
axA.barh(y, ev.uns, left=ev.act + ev.inh, color=UNS, height=0.7, label="unspecified")
axA.set_yticks(y); axA.set_yticklabels([f"{r} ({'core' if ic == 'yes' else 'not core'})" for r, ic in zip(rs.regulator_HGNC, rs.regulator_in_core)], style="italic")
axA.set_xlabel("non-technical edges to core genes (n; composite sources counted for each member)"); axA.set_title("Most JunD/JunB edges to the core\nare of unspecified sign", loc="left")
axA.legend(loc="lower right", fontsize=6, handlelength=1.0); panel_letter(axA, "A", case="upper", dx=-0.42)
# B: number of supporting references per edge, JunD / JunB vs all other regulators
def refs(sub):
    return sub.n_references.fillna(0).astype(int).values
allnt = s7[s7.effect != "Technical"]
jd = allnt[allnt.source_HGNC == "JUND"]; jb = allnt[allnt.source_HGNC == "JUNB"]
oth = allnt[~allnt.source_HGNC.str.split("/").apply(lambda xs: ("JUND" in xs) or ("JUNB" in xs))]
bins = [0.5, 1.5, 2.5, 4.5, 9.5, 1000]; blab = ["1", "2", "3–4", "5–9", "≥10"]
for i, (name, sub, c) in enumerate([("JunD edges", jd, "#4d4d4d"), ("JunB edges", jb, "#969696"), ("all other regulators", oth, "#cccccc")]):
    h, _ = np.histogram(refs(sub), bins=bins); h = h / h.sum() * 100
    axB.bar(np.arange(5) + (i - 1) * 0.27, h, width=0.26, color=c, label=f"{name} (n = {len(sub)})")
axB.set_xticks(range(5)); axB.set_xticklabels(blab); axB.set_xlabel("supporting references per edge (n)"); axB.set_ylabel("edges (%)")
single = (jd.n_references.fillna(0).astype(int) == 1).sum()
pm = jd.references_PMID.fillna("").astype(str).str.contains("40455848").sum()
axB.set_title(f"JunD edges rest on one study:\n{pm} of {len(jd)} single-source JunD edges cite PMID 40455848", loc="left")
axB.legend(loc="upper right", fontsize=6, handlelength=1.0); panel_letter(axB, "B", case="upper", dx=-0.2)
# C: decoupleR/CollecTRI edges (S8): TF (y) x target (x), marker by sex, color by sign; targets not in final list hatched
tfs = ["SP1", "STAT3", "JUN", "NFKB1", "RELA", "HIF1A", "NFKB", "CREB1", "AP1"]
tfs = [t for t in tfs if t in set(s8.TF_source)] + [t for t in sorted(set(s8.TF_source)) if t not in tfs]
tg_m = list(dict.fromkeys(s8[s8.sex == "men"].target)); tg_w = list(dict.fromkeys(s8[s8.sex == "women"].target))
targets = tg_m + [t for t in tg_w if t not in tg_m]
xi = {t: i for i, t in enumerate(targets)}; yi = {t: i for i, t in enumerate(tfs[::-1])}
for _, r in s8.iterrows():
    col = ACT if r.sign == "activating" else INH
    mk = "o" if r.sex == "men" else "s"
    face = col if r.target_in_final_HPC_list == "yes" else "white"
    axC.scatter(xi[r.target], yi[r.TF_source], marker=mk, s=42, facecolor=face, edgecolor=col, lw=1.0, zorder=3)
axC.set_xticks(range(len(targets))); axC.set_xticklabels(targets, rotation=60, ha="right", style="italic")
axC.set_yticks(range(len(tfs))); axC.set_yticklabels(tfs[::-1]); axC.set_xlim(-0.6, len(targets) + 4.5); axC.set_ylim(-0.6, len(tfs) - 0.4)
axC.grid(True, color="#eeeeee", lw=0.5); axC.set_axisbelow(True)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", color=ACT, ls="", mfc=ACT, label="men, activating"), Line2D([], [], marker="o", color=INH, ls="", mfc=INH, label="men, repressing"),
     Line2D([], [], marker="s", color=ACT, ls="", mfc=ACT, label="women, activating"), Line2D([], [], marker="o", color="#4d4d4d", ls="", mfc="white", label="target not in final HPC list")]
axC.legend(handles=h, loc="center right", fontsize=6, ncol=1, handletextpad=0.3)
nm = (s8.sex == "men").sum(); nw = (s8.sex == "women").sum()
axC.set_title(f"CollecTRI edges from the 10 top-ranked decoupleR TFs to HPC sex-specific targets (men {nm} edges, women {nw}; run on an earlier, unfiltered HPC list)", loc="left", fontsize=7)
axC.set_ylabel("transcription factor"); axC.set_xlabel("sex-specific target gene")
panel_letter(axC, "C", case="upper", dx=-0.09)
fig.savefig("out/Supplementary_Figure_S5.png", dpi=600); fig.savefig("out/Supplementary_Figure_S5.pdf"); plt.close(fig)
print("done", ev.to_string(), "\nJunD single-ref", single, "of", len(jd), "; PMID 40455848:", pm, "; JunB n", len(jb))