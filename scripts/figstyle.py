"""Shared figure style for the scMicroglia manuscript (print figures).

Palette = validated reference palette (dataviz skill, light mode):
  men #2a78d6 / women #eb6834 / dimorphic #4a3aa7 pass all-pairs CVD + normal-vision checks.
  HPC #1baf7a (aqua) is < 3:1 on white -> always direct-label it.
Usage:  from figstyle import *; apply_style()
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ---- semantic colours -------------------------------------------------------
MEN = "#2a78d6"        # categorical slot 1 (blue)
WOMEN = "#eb6834"      # categorical slot 2 (orange)
DIMORPHIC = "#4a3aa7"  # slot 7 (violet)
SHARED = "#52514e"     # neutral (text-secondary ink) for sex-shared
HPC = "#1baf7a"        # slot 3 (aqua) - needs direct labels
PFC = "#4a3aa7"        # violet (only used where dimorphic is not shown)
UP = "#e34948"         # diverging warm pole (increased)
DOWN = "#2a78d6"       # diverging cool pole (decreased)
NEUTRAL = "#f0efec"    # diverging midpoint
GRID = "#e6e5e1"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#9a9893"
HIGHLIGHT = "#eda100"  # slot 4 (yellow) - outline/accent only, never alone

CAT_COLORS = {"sex-shared": SHARED, "female-specific": WOMEN,
              "male-specific": MEN, "sex-dimorphic": DIMORPHIC}

# diverging blue -> grey -> red for log2FC heatmaps
DIVERGING = LinearSegmentedColormap.from_list(
    "bluegreyred", ["#104281", "#2a78d6", "#86b6ef", NEUTRAL, "#f2a3a2", "#e34948", "#a32a29"])
# sequential blue for magnitude (-log10 FDR etc.)
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "blueseq", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])

def apply_style():
    mpl.rcParams.update({
        "font.family": ["Arial", "Segoe UI", "DejaVu Sans"],
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.edgecolor": INK2, "axes.linewidth": 0.6, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.5,
        "legend.frameon": False, "figure.dpi": 150, "savefig.dpi": 600,
        "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.titleweight": "bold", "axes.titlelocation": "left",
    })

def panel_label(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left", color=INK)

def save(fig, path_noext):
    """Save PNG (600 dpi) + PDF (vector)."""
    fig.savefig(path_noext + ".png")
    fig.savefig(path_noext + ".pdf")
    plt.close(fig)

# journal widths (inches): single column 3.5, 1.5 col 5.5, double column 7.2
W1, W15, W2 = 3.5, 5.5, 7.2
