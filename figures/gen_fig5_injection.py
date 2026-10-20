"""Fig 5: synthetic one-eye fault injection (B3a) -- abstention recall among admitted events, and excess rejection under
mild coherent asymmetry (N1). Source: results/b3a_injection.csv (frozen criterion: recall >= .80 with lower bound > .70
for judged conditions)."""
from paper_plot_style import *
import matplotlib.ticker
import numpy as np

B = res("b3a_injection.csv").set_index("cond")
GROUPS = [("F1\none-eye step", [("F1_a2", "2°"), ("F1_a4", "4°"), ("F1_a8", "8°")]),
          ("F2a\none-eye slowing", [("F2a_f1.5", "×1.5"), ("F2a_f2", "×2"), ("F2a_f3", "×3")]),
          ("F2b\npeak glitch", [("F2b_b3", "3°"), ("F2b_b6", "6°")]),
          ("F3\ndropout", [("F3", "")]),
          ("contamination\n(slowing)", [("C_p5", "5%"), ("C_p10", "10%"), ("C_p20", "20%"), ("C_p30", "30%")])]
CRIT = 0.80

fig, axs = plt.subplots(1, 2, figsize=(W2, 2.5), gridspec_kw=dict(width_ratios=[4.2, 1], wspace=0.35))
ax = axs[0]
x = 0; ticks, labs, gx = [], [], []
for g, conds in GROUPS:
    x0 = x
    for c, lab in conds:
        r = B.loc[c]
        judged, inc = bool(r.judged), bool(r.inconclusive)
        col = OI["grey"] if c.startswith("C_") else OI["blue"]
        if inc:
            ax.plot(x, r.recall, "x", color=OI["verm"], ms=5, mew=1.0)
            ax.annotate(f"not interpreted:\n{int(r.n_admitted)} of {int(r.injections)}\nreached centring", (x, r.recall),
                        xytext=(x, 0.33), ha="center", fontsize=FS - 2, color=OI["verm"],
                        arrowprops=dict(arrowstyle="-", lw=0.4, color=OI["verm"]))
        else:
            ax.errorbar(x, r.recall, yerr=[[r.recall - r.recall_lo], [r.recall_hi - r.recall]], fmt="o", color=col,
                        mfc=col if judged else "white", mec=col, ms=4, capsize=1.5, elinewidth=0.8)
        ax.plot(x, r.admitted_rate, "_", color=OI["black"], ms=7, mew=1.0)
        ticks.append(x); labs.append(lab); x += 1
    gx.append(((x0 + x - 1) / 2, g)); x += 0.8
ax.axhline(CRIT, color=OI["verm"], lw=0.7, ls="--")
ax.set_xticks(ticks); ax.set_xticklabels(labs); ax.set_ylim(-0.03, 1.05); ax.set_xlim(-0.6, ticks[-1] + 0.6)
ax.set_ylabel("proportion")
for xc, g in gx:
    ax.text(xc, -0.11, g, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=FS - 1)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", color=OI["blue"], ls="none", ms=4, label="recall (of reaching centring), judged"),
     Line2D([], [], marker="o", color=OI["blue"], mfc="white", ls="none", ms=4, label="recall, unjudged"),
     Line2D([], [], marker="o", color=OI["grey"], mfc="white", ls="none", ms=4, label="recall, descriptive"),
     Line2D([], [], marker="_", color=OI["black"], ls="none", ms=7, mew=1.0, label="reached centring (of all injected)"),
     Line2D([], [], marker="x", color=OI["verm"], ls="none", ms=5, mew=1.0, label="recall not interpreted (< 30% reached centring)"),
     Line2D([], [], color=OI["verm"], lw=0.7, ls="--", label="criterion: recall ≥ .80 and lower bound > .70")]
ax.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=FS - 1.5, columnspacing=1.0,
          handletextpad=0.3)
panel(ax, "a", x=-0.07, y=1.17)

# (b) N1: excess rejection of untouched-but-scaled events (false abstention)
ax = axs[1]
for i, (c, lab) in enumerate([("N1_k0.85", "κ .85"), ("N1_k1.15", "κ 1.15")]):
    r = B.loc[c]
    ax.errorbar(i, r.excess, yerr=[[r.excess - r.excess_lo], [r.excess_hi - r.excess]], fmt="o", color=OI["black"],
                mfc="white", ms=4, capsize=1.5, elinewidth=0.8)
ax.set_xticks([0, 1]); ax.set_xticklabels(["κ .85", "κ 1.15"]); ax.set_xlim(-0.6, 1.6)
ax.set_ylim(0, 0.02); ax.set_yticks([0, .005, .01, .015, .02]); ax.set_ylabel("excess rejection")
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: fm(v)))
ax.text(0.5, -0.11, "N1 mild one-eye\nscaling (control)", transform=ax.transAxes, ha="center", va="top", fontsize=FS - 1)
panel(ax, "b", x=-0.6, y=1.17)
save_fig(fig, "fig5_injection")
