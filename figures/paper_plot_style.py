"""Shared figure style (Scientific Reports). Sans-serif (Helvetica-metric Nimbus Sans) per Nature-family figure
guidance -- a deliberate override of the serif default; vector PDF; no in-figure titles; Okabe-Ito palette."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
SRC = os.path.join(ROOT, "figures", "source_data")
FIG_DIR = os.path.join(ROOT, "figures")
FS = 8  # Sci Rep: 5-8 pt at final size; figures drawn at final width (89 mm single / 183 mm double)
MM = 1 / 25.4
W1, W2 = 89 * MM, 183 * MM

matplotlib.rcParams.update({
    "font.size": FS, "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "axes.labelsize": FS, "xtick.labelsize": FS - 1, "ytick.labelsize": FS - 1, "legend.fontsize": FS - 1,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 1.0, "lines.markersize": 4,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "axes.grid": False, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "mathtext.fontset": "dejavusans", "legend.frameon": False,
})

# Okabe-Ito (colour-blind safe)
OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "verm": "#D55E00", "purple": "#CC79A7", "grey": "#7F7F7F"}
# construction -> (label, colour); u = open marker, c = filled
CONS = ["mean_u", "mean_c", "adm_u", "adm_c", "pereye_u", "pereye_c"]
FAM_COL = {"mean": OI["verm"], "adm": OI["blue"], "pereye": OI["green"]}
FAM_LAB = {"mean": "mean-speed", "adm": "admitted", "pereye": "per-eye"}


def cons_label(c):
    f, k = c.split("_")
    return f"{FAM_LAB[f]}, {'centred' if k == 'c' else 'unfiltered'}"


def cons_style(c):
    f, k = c.split("_")
    col = FAM_COL[f]
    return dict(color=col, mfc=col if k == "c" else "white", mec=col, marker="o" if k == "c" else "s")


def res(name):
    return pd.read_csv(os.path.join(RES, name))


def src(name):
    return pd.read_csv(os.path.join(SRC, name))


def panel(ax, letter, x=-0.16, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=FS + 1, fontweight="bold", va="bottom", ha="left")


def f3(x):
    """3-decimal, leading zero dropped (journal style for bounded statistics)."""
    s = f"{x:.3f}"
    if s in ("0.000", "-0.000") and x != 0:   # avoid a bare signed zero: one more decimal
        s = f"{x:.4f}"
        if s in ("0.0000", "-0.0000"):        # numerically zero (float noise)
            s = "0.000"
    return s.replace("0.", ".", 1) if s.startswith(("0.", "-0.")) else s


def save_fig(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"))
    plt.close(fig)
    print("Saved:", name)


def fm(x):
    """f3 for figure text: typographic minus (U+2212)."""
    return f3(x).replace("-", "−")
