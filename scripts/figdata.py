"""Shared data loaders for Figures 1 and 2 (read-only inputs; 'gender' relabelled 'sex')."""
import re
import numpy as np
import pandas as pd

ROOT = "<BUILD_DIR>"
DATA = f"{ROOT}/data/microglia-subtypes"
TAB = f"{ROOT}/agents_out/tables"
OUT = f"{ROOT}/figures/out"
SRC = f"{OUT}/source_data"

COLS = ["category", "gene", "lfc_m", "lfc_f", "padj_m", "padj_f"]
CAT_ORDER = ["sex-shared", "female-specific", "male-specific", "sex-dimorphic"]
FEEDBACK = ["NR4A1", "NR4A3", "BCL6", "BCL3", "SOCS3", "DUSP1", "ZFP36L1"]


def _relabel(df):
    df = df.copy()
    df.columns = COLS
    df["category"] = df["category"].str.replace("gender", "sex", regex=False)
    return df


def load_hpc():
    return _relabel(pd.read_csv(f"{DATA}/hipo_analysis/dge_filterd.csv"))


def load_pfc():
    return _relabel(pd.read_csv(f"{DATA}/pfc_result/dge_filterd.csv"))


def load_pfc_unfiltered():
    return _relabel(pd.read_csv(f"{DATA}/pfc_result/pfc_degs.csv", index_col=0))


def direction(df):
    """+1/-1: female-specific -> sign in women; all others -> sign in men
    (sex-shared genes have the same sign in both sexes)."""
    s = np.where(df["category"] == "female-specific", np.sign(df["lfc_f"]), np.sign(df["lfc_m"]))
    return s.astype(int)


def state_markers():
    """Parse MG0/MG2/MG8/MG10 marker lists from script.R."""
    txt = open(f"{ROOT}/data/script.R").read()
    out = {}
    for m in re.finditer(r"(MG\d+)_\w+\s*=\s*c\(([^)]*)\)", txt, flags=re.S):
        out[m.group(1)] = re.findall(r'"([^"]+)"', m.group(2))
    return {k: out[k] for k in ["MG0", "MG2", "MG8", "MG10"]}


def core109():
    """109 core genes with full-precision log2FC from the filtered lists."""
    h = load_hpc().set_index("gene")
    p = load_pfc().set_index("gene")
    genes = sorted(set(h.index[h.category == "sex-shared"]) & set(p.index[p.category == "sex-shared"]))
    df = pd.DataFrame({
        "gene": genes,
        "HPC_men": h.loc[genes, "lfc_m"].values, "HPC_women": h.loc[genes, "lfc_f"].values,
        "PFC_men": p.loc[genes, "lfc_m"].values, "PFC_women": p.loc[genes, "lfc_f"].values,
    })
    fg = pd.read_csv(f"{TAB}/C_109_functional_groups.csv")
    fg = fg[fg.group != "group"]
    df = df.merge(fg.rename(columns={"group": "functional_group"}), on="gene", how="left")
    assert df.functional_group.notna().all() and len(df) == 109
    hi = lambda a, b: (a.abs() > 1.5) & (b.abs() > 1.5)
    df["stratum"] = np.where(hi(df.HPC_men, df.HPC_women), "high", "low")
    assert (df["stratum"] == np.where(hi(df.PFC_men, df.PFC_women), "high", "low")).all()
    mk = state_markers()
    df["state_marker"] = df.gene.map(lambda g: "/".join(k for k in ["MG0", "MG2", "MG8", "MG10"] if g in mk[k]))
    df["feedback_regulator"] = df.gene.isin(FEEDBACK)
    df["sumFC"] = df[["HPC_men", "HPC_women", "PFC_men", "PFC_women"]].sum(axis=1)
    return df
