"""Fig 6: amplitude-matched left/right peak-velocity difference D (B5) and its association with A.
Sources: figures/source_data/participant_D.csv (from results/passes/ran/pairs_*_primary.parquet via b5.participant_D),
results/b5_matched.csv, results/b5_transfer.csv, results/pw_p8_cross_transfer.csv"""
from paper_plot_style import *
import numpy as np
from scipy.stats import spearmanr

P = src("participant_D.csv")
M = res("b5_matched.csv").query("variant == 'primary'").set_index("task")
T = res("b5_transfer.csv").query("variant == 'primary'").set_index("task")
X = res("pw_p8_cross_transfer.csv").query("task == 'TEX' and variant == 'primary'")
COL = {"TEX": OI["blue"], "RAN": OI["orange"]}
LAB = {"TEX": "reading (TEX)", "RAN": "random saccades (RAN)"}

fig = plt.figure(figsize=(W2, 2.5))
axs = [fig.add_axes(r) for r in ([0.065, 0.17, 0.235, 0.78], [0.36, 0.17, 0.215, 0.78], [0.815, 0.17, 0.175, 0.78])]
# (a) participant D distributions
ax = axs[0]
LO, HI = -0.2, 0.25  # values beyond are collected in the edge bins (counted on the panel)
bins = np.linspace(LO, HI, 46)
for t in ("TEX", "RAN"):
    d = P.query("task == @t").D
    assert len(d) == int(M.loc[t].n_pid) and abs(d.mean() - M.loc[t].mean_D) < 1e-9
    k = int(((d < LO) | (d > HI)).sum())
    edge = f" ({k} in edge bins)" if k else ""
    ax.hist(d.clip(LO + 1e-9, HI - 1e-9), bins=bins, density=True, histtype="step", color=COL[t], lw=1.0,
            label=f"{t}: n = {len(d)}{edge}\nmean {fm(M.loc[t].mean_D)} [{fm(M.loc[t].mean_D_lo)}, {fm(M.loc[t].mean_D_hi)}]")
ax.axvline(0, color=OI["grey"], lw=0.5, ls=":")
ax.set_xlabel("participant D"); ax.set_ylabel("density"); ax.legend(loc="upper left", fontsize=FS - 2, handlelength=1.2)
ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
panel(ax, "a", x=-0.25)

# (b) D vs A (A: R1 two-session mean of centred mean-speed A, reading task); y-range set to the central spread,
# values below it are drawn as triangles at the lower edge and counted in the panel
ax = axs[1]
YLO, YHI = -0.3, 0.6
n_below = 0
for t in ("TEX", "RAN"):
    q = P.query("task == @t").dropna(subset=["A_tex_mean_c"])
    rho = spearmanr(q.D, q.A_tex_mean_c).statistic
    assert len(q) == int(T.loc[t].n) and abs(rho - T.loc[t].rho) < 1e-6
    assert (q.D <= YHI).all()
    inn, low = q[q.D >= YLO], q[q.D < YLO]; n_below += len(low)
    ax.scatter(inn.A_tex_mean_c, inn.D, s=4, facecolor="none" if t == "TEX" else COL[t], edgecolor=COL[t], lw=0.4,
               alpha=0.7, label=f"{t} D: ρ = {fm(rho)}, n = {len(q)}")
    ax.scatter(low.A_tex_mean_c, np.full(len(low), YLO + 0.02), s=10, marker="v",
               facecolor="none" if t == "TEX" else COL[t], edgecolor=COL[t], lw=0.6)
ax.set_ylim(YLO, YHI)
ax.text(0.98, 0.075, f"{n_below} values < −.3\nshown as ▽ at edge", transform=ax.transAxes, ha="right", va="bottom", fontsize=FS - 2)
ax.axhline(0, color=OI["grey"], lw=0.5, ls=":"); ax.axvline(0, color=OI["grey"], lw=0.5, ls=":")
ax.set_xlabel("A (reading, centred mean-speed)"); ax.set_ylabel("participant D")
ax.legend(loc="upper left", fontsize=FS - 2, handletextpad=0.2, markerscale=1.5)
panel(ax, "b", x=-0.27)

# (c) Spearman rho with 95% CI: same-session, cross-session and pooled associations
ax = axs[2]
rows = []
for _, r in X.iterrows():
    if r.kind in ("same", "cross"):
        rows.append((f"A S{int(r.A_session)} – D S{int(r.D_session)} (n = {int(r.n)})", r.rho, r.rho_lo, r.rho_hi, COL["TEX"]))
for t in ("TEX", "RAN"):
    r = T.loc[t]
    rows.append((f"{t} D – A, pooled (n = {int(r.n)})", r.rho, r.rho_lo, r.rho_hi, COL[t]))
nss = sum(X.kind.isin(["same", "cross"]))
ax.axhline(nss - 0.5, color=OI["grey"], lw=0.4, ls="-")  # single-session rows | pooled rows
for i, (lab, v, lo, hi, c) in enumerate(rows):
    ax.errorbar(v, i, xerr=[[v - lo], [hi - v]], fmt="o", color=c, ms=3.5, capsize=1.5, elinewidth=0.8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=FS - 0.5)
ax.set_ylim(len(rows) - 0.5, -1.9)
ax.axvline(0, color=OI["grey"], lw=0.5, ls=":"); ax.set_xlim(-0.05, 0.75); ax.set_xlabel("Spearman ρ (95% CI)")
m = M.loc["TEX"]
ax.text(0.08, 0.985, f"TEX D repeatability:\nICC = {fm(m.icc_a1)} [{fm(m.icc_a1_lo)}, {fm(m.icc_a1_hi)}]\nn = {int(m.icc_n)}",
        transform=ax.transAxes, ha="left", va="top", fontsize=FS - 1)
panel(ax, "c", x=-1.12)
save_fig(fig, "fig6_matched")
