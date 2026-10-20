"""Fig 3: participant influence (sensitivity analysis, not a cleaning rule).
Top row: own-cohort influence for mean-speed events; bottom row: the admission contrast in the joint cohort and the
diagnostics of its most influential participant. No new statistics: every plotted value is read from
figures/source_data/recording_A.csv (from results/passes via make_source_data.py), results/pw_p4_loo.csv,
results/pw_p4_loo_top.csv, results/rv_v4_loo_summary.csv and results/rv_v4_loo_admission.csv; the joint-cohort
ICCs are recomputed only to assert that the plotted cohort reproduces rv_v4_loo_summary.csv."""
import os, sys
from paper_plot_style import *
import numpy as np

sys.path.insert(0, os.path.join(ROOT, "experiments"))
from common import icc  # noqa: E402

R = src("recording_A.csv").query("floor == 5 and round == 1")
LOO = res("pw_p4_loo.csv").query("floor == 5 and comparison == 'R1_S1S2'")
TOP = res("pw_p4_loo_top.csv").query("floor == 5 and comparison == 'R1_S1S2' and method == 'mean_u' and cohort == 'own_y'")
SUM = res("rv_v4_loo_summary.csv").query("floor == 5").iloc[0]
LA = res("rv_v4_loo_admission.csv")


def s1s2(c, var="A"):
    w = R.query("construction == @c").pivot(index="pid", columns="session", values=var).dropna()
    return w[1], w[2]


lu = LOO.query("method == 'mean_u' and cohort == 'own_y'").iloc[0]
lc = LOO.query("method == 'mean_c' and cohort == 'own_x'").iloc[0]
x_u, y_u = s1s2("mean_u")
assert len(x_u) == int(lu.n) == 379
sq = ((x_u - y_u) ** 2).sort_values(ascending=False); share = sq / sq.sum()
assert abs(share.iloc[:3].sum() - lu.top3_sqdiff_share) < 1e-9
# participants ranked by |leave-one-out shift| (as in the deletion analysis); for mean_u at floor 5 the top two coincide
top2 = list(TOP.reindex(TOP["shift"].abs().sort_values(ascending=False).index).pid.iloc[:2])
assert top2 == list(share.index[:2])

# joint cohort: both R1 sessions eligible under all six constructions
joint = sorted(set.intersection(*[set(s1s2(c)[0].index) for c in CONS]))
assert len(joint) == int(SUM.n) == 376
J = {c: [s.loc[joint] for s in s1s2(c)] for c in ("mean_u", "adm_u")}
for c, ref in (("mean_u", SUM.icc_mean_u), ("adm_u", SUM.icc_adm_u)):
    assert abs(icc(*J[c])[0] - ref) < 1e-9
P = int(SUM.top_pid)
assert P == int(LA.loc[LA.sq_share.idxmax(), "pid"]) == top2[1]

fig, axs = plt.subplots(2, 3, figsize=(W2, 4.7), gridspec_kw=dict(wspace=0.55, hspace=0.62))
axs = axs.ravel()
SMALL = FS - 1.5

# (a) share of summed squared session differences, own cohort, top 10
ax = axs[0]
k = 10
ax.bar(range(1, k + 1), share.iloc[:k].values * 100, color=[OI["black"]] * 3 + [OI["grey"]] * (k - 3), width=0.7)
ax.set_xticks(range(1, k + 1)); ax.set_xlabel("participant rank"); ax.set_ylabel("share of Σ(S1 − S2)² (%)")
ax.text(0.97, 0.95, f"top 3 (black): {100 * lu.top3_sqdiff_share:.0f}%\n{cons_label('mean_u')}\nown cohort, n = {int(lu.n)}",
        transform=ax.transAxes, ha="right", va="top", fontsize=SMALL)
panel(ax, "a", x=-0.36)

# (b) ICC after deleting the k most influential participants (own cohort of each construction)
ax = axs[1]
ks = [0, 1, 2, 3, 5]
for r, c in ((lu, "mean_u"), (lc, "mean_c")):
    v = [r.icc, r.icc_drop_top1, r.icc_drop_top2, r.icc_drop_top3, r.icc_drop_top5]
    st = cons_style(c)
    ax.plot(ks, v, marker=st["marker"], color=st["color"], mfc=st["mfc"], mec=st["mec"], ms=3.5,
            label=f"{cons_label(c).split(', ')[1]}, n = {int(r.n)}")
ax.set_xticks(ks); ax.set_xlabel("most influential removed (k)"); ax.set_ylabel("ICC(A,1), R1 S1–S2")
ax.set_ylim(0, 0.85)
ax.legend(loc="lower right", bbox_to_anchor=(1, 0.0), fontsize=SMALL, handlelength=1.5, title="mean-speed, own cohort",
          title_fontsize=SMALL)
panel(ax, "b", x=-0.36)

# (c) leave-one-out admission contrast in the joint cohort
ax = axs[2]
oth = LA[LA.pid != P]; me = LA[LA.pid == P].iloc[0]
ax.scatter(oth.sq_share * 100, oth.delta, s=5, facecolor="white", edgecolor=OI["grey"], lw=0.4)
ax.errorbar([me.sq_share * 100], [me.delta], yerr=[[me.delta - SUM.delta_wo_top_lo], [SUM.delta_wo_top_hi - me.delta]],
            fmt="*", ms=7, color=OI["black"], elinewidth=0.6, capsize=2, zorder=3)
ax.axhline(SUM.delta_full, color=OI["blue"], lw=0.6, ls="--")
ax.axhline(0, color=OI["grey"], lw=0.5, ls=":")
ax.text(0.97, SUM.delta_full + 0.02, f"all {int(SUM.n)}: {fm(SUM.delta_full)}", transform=ax.get_yaxis_transform(),
        ha="right", va="bottom", fontsize=SMALL, color=OI["blue"])
ax.annotate(f"participant omitted\n{fm(me.delta)} [{fm(SUM.delta_wo_top_lo)}, {fm(SUM.delta_wo_top_hi)}]",
            (me.sq_share * 100, me.delta), xytext=(-6, 4), textcoords="offset points", ha="right", va="bottom",
            fontsize=SMALL)
ax.set_xlim(-3, 75); ax.set_ylim(-0.05, 0.45)
ax.set_xlabel("omitted participant's\nshare (%)")
ax.set_ylabel("admission contrast (ΔICC)")
panel(ax, "c", x=-0.36)

# (d, e) joint-cohort S1 vs S2 under mean_u and adm_u, shared axes
allv = np.concatenate([np.r_[J[c][0].values, J[c][1].values] for c in J])
lo_, hi_ = allv.min(), allv.max(); pad = 0.05 * (hi_ - lo_)
for ax, c, lab, ic in ((axs[3], "mean_u", "d", SUM.icc_mean_u), (axs[4], "adm_u", "e", SUM.icc_adm_u)):
    x, y = J[c]; st = cons_style(c); hi = x.index == P
    ax.scatter(x[~hi], y[~hi], s=5, facecolor=st["mfc"], edgecolor=st["color"], lw=0.4, alpha=0.7)
    ax.scatter(x[hi], y[hi], s=30, marker="*", color=OI["black"], lw=0, zorder=3)
    ax.plot([lo_ - pad, hi_ + pad], [lo_ - pad, hi_ + pad], color=OI["grey"], lw=0.5, ls=":")
    ax.set_xlim(lo_ - pad, hi_ + pad); ax.set_ylim(lo_ - pad, hi_ + pad); ax.set_aspect("equal")
    ax.set_xlabel("A, session 1"); ax.set_ylabel("A, session 2")
    ax.text(0.03, 0.97, f"{cons_label(c)}\njoint cohort\nn = {len(x)}; ICC = {fm(ic)}", transform=ax.transAxes,
            ha="left", va="top", fontsize=SMALL)
    panel(ax, lab, x=-0.36)

# (f) diagnostics of the starred participant: class medians of l and RS counts per session
ax = axs[5]
pr = R.query("pid == @P and construction in ['mean_u', 'adm_u']").set_index(["construction", "session"])
xs = {(1, "rs"): 0, (1, "prog"): 1, (2, "rs"): 2.4, (2, "prog"): 3.4}
for c, dx in (("mean_u", -0.12), ("adm_u", 0.12)):
    st = cons_style(c)
    for (s, cl), xp in xs.items():
        v = pr.loc[(c, s), cl]
        ax.plot(xp + dx, v, marker=st["marker"], ms=4, color=st["color"], mfc=st["mfc"], mec=st["mec"], ls="none",
                label=FAM_LAB[c.split("_")[0]] + ", unfiltered" if (s, cl) == (1, "rs") else None)
        if cl == "rs":
            up = v > 1; side = "center" if up else ("right" if dx < 0 else "left")
            ax.text(xp + dx + (0 if up else 1.2 * dx), v + (0.12 if up else -0.12), f"{int(pr.loc[(c, s), 'n_rs'])}", ha=side,
                    va="bottom" if v > 1 else "top", fontsize=SMALL - 0.5, color=st["color"])
ax.axhline(0, color=OI["grey"], lw=0.5, ls=":")
ax.set_xticks(list(xs.values())); ax.set_xticklabels(["RS", "prog.", "RS", "prog."])
ax.set_xlim(-0.5, 3.9); ax.set_ylim(-0.6, 2.6)
for xc, s in ((0.5, "session 1"), (2.9, "session 2")):
    ax.text(xc, -0.17, s, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=SMALL)
ax.set_ylabel("class median of $l$")
ax.legend(loc="center right", bbox_to_anchor=(1, 0.62), fontsize=SMALL, handletextpad=0.2, borderaxespad=0.1)
ax.text(0.97, 0.40, "starred participant\nnumbers: RS count", transform=ax.transAxes, ha="right", va="top",
        fontsize=SMALL)
panel(ax, "f", x=-0.36)
save_fig(fig, "fig3_influence")
