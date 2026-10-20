"""Fig 4: between-construction agreement of recording-level A (Bland-Altman), all eligible recordings, floor 5.
Agreement between constructions, not between-session repeatability. LoA/bias/r from results/pw_p2_agreement.csv
(participant-clustered CIs); empirical coverage and 2.5th/97.5th percentiles from results/rv_v2_ba_coverage.csv;
points from figures/source_data/recording_A.csv."""
from paper_plot_style import *
from matplotlib.lines import Line2D
import numpy as np

R = src("recording_A.csv").query("floor == 5")
AG = res("pw_p2_agreement.csv").query("floor == 5")
CV = res("rv_v2_ba_coverage.csv").query("floor == 5")
PAIRS = [("mean_c", "pereye_c", "a"), ("mean_u", "pereye_u", "b")]
PCT_LS = (0, (1, 1.5))

fig, axs = plt.subplots(1, 2, figsize=(W2 * 0.78, 2.9), gridspec_kw=dict(wspace=0.42))
for ax, (x, y, lab) in zip(axs, PAIRS):
    a = AG.query("x == @x and y == @y").iloc[0]; cv = CV.query("x == @x and y == @y").iloc[0]
    w = R.query("construction in [@x, @y]").pivot_table(index=["pid", "round", "session"], columns="construction",
                                                        values="A").dropna()
    assert len(w) == int(a.n_rec) == int(cv.n_rec) and w.index.get_level_values("pid").nunique() == int(a.n_pid)
    d = w[y] - w[x]; m = (w[x] + w[y]) / 2
    assert abs(d.mean() - a.bias) < 1e-9 and abs(d.mean() - cv.bias) < 1e-9 and abs(cv.loa_lo - a.loa_lo) < 1e-9

    def draw(ax_, small=False):
        ax_.scatter(m, d, s=3 if small else 4, facecolor="none", edgecolor=OI["blue"], lw=0.4, alpha=0.6)
        for v, ls in ((a.bias, "-"), (a.loa_lo, "--"), (a.loa_hi, "--")):
            ax_.axhline(v, color=OI["black"], lw=0.7, ls=ls)
        for v in (cv.p025, cv.p975):
            ax_.axhline(v, color=OI["verm"], lw=0.8, ls=PCT_LS)
    draw(ax)
    ax.axhline(0, color=OI["grey"], lw=0.5, ls=":")
    ax.set_xlabel("mean of the two A values")
    ax.set_ylabel(f"A({cons_label(y)}) −\nA({cons_label(x)})")
    rng = d.max() - d.min(); ax.set_ylim(d.min() - (0.42 if lab == "a" else 0.30) * rng, d.max() + 0.05 * rng)
    ax.text(0.03, 0.03, f"{int(a.n_rec)} recordings, {int(a.n_pid)} participants\n"
            f"bias {('.000' if abs(a.bias) < 0.0005 else fm(a.bias))}; LoA [{fm(a.loa_lo)}, {fm(a.loa_hi)}]\n"
            f"inside LoA: {int(cv.n_inside)} ({100 * cv.coverage:.1f}%)\n"
            f"2.5th–97.5th pct [{fm(cv.p025)}, {fm(cv.p975)}]\nr = {fm(a.pearson)} [{fm(a.pearson_lo)}, {fm(a.pearson_hi)}]",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=FS - 1.5)
    panel(ax, lab, x=-0.36)
    if lab == "b":
        n_out = int((d.abs() > 0.1).sum())
        ins = ax.inset_axes([0.47, 0.36, 0.50, 0.34])
        draw(ins, small=True)
        ins.set_ylim(-0.1, 0.1); ins.set_xlim(m.quantile(0.005) - 0.02, m.quantile(0.995) + 0.02)
        ins.tick_params(labelsize=FS - 2.5, length=1.5, pad=1)
        ins.set_title(f"central range (|difference| ≤ .1);\n{n_out} recordings outside", fontsize=FS - 2.5, pad=2)
fig.legend(handles=[Line2D([], [], color=OI["black"], lw=0.7, ls="-", label="mean difference (bias)"),
                    Line2D([], [], color=OI["black"], lw=0.7, ls="--", label="bias ± 1.96 SD (limits of agreement)"),
                    Line2D([], [], color=OI["verm"], lw=0.8, ls=PCT_LS, label="empirical 2.5th / 97.5th percentiles")],
           loc="lower center", bbox_to_anchor=(0.5, 0.98), ncol=3, fontsize=FS - 1.5, handlelength=2.2, columnspacing=1.2)
save_fig(fig, "fig4_bland_altman")
