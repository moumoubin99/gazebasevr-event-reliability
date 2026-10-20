"""Fig 2: short-term (R1 S1-S2) repeatability by event construction, path-specific contrasts with and without the
most influential participant, eligibility, and detector-setting sensitivity.
Sources: results/pw_p1_joint_icc.csv, results/pw_p1_joint_delta.csv, results/b1_reliability.csv, results/b2_selection.csv,
results/rv_v7_ceiling_icc.csv, results/rv_v7_paths.csv, results/rv_v9_sens_summary.csv (experiments/r3_revision.py)"""
from paper_plot_style import *
import numpy as np

J = res("pw_p1_joint_icc.csv").query("floor == 5 and comparison == 'R1_S1S2' and cohort == 'joint6'").set_index("method")
D = res("pw_p1_joint_delta.csv").query("floor == 5 and comparison == 'R1_S1S2' and cohort == 'joint6'").set_index("contrast")
B1 = res("b1_reliability.csv").query("cond == 'base' and comparison == 'R1_S1S2'")
B2 = res("b2_selection.csv").set_index("floor")
N = int(J.n.iloc[0]); assert N == 376
CI = res("rv_v7_ceiling_icc.csv").query("tag == 'floor5_ceil1000'").set_index("method")
PT = res("rv_v7_paths.csv")
SV = res("rv_v9_sens_summary.csv").set_index("variant")
assert int(CI.n.iloc[0]) == N and all(abs(CI.loc[c].icc_a1 - J.loc[c].icc_a1) < 1e-9 for c in CONS)
CEIL = dict(color=OI["purple"], mfc="white", mec=OI["purple"], marker="D")

fig = plt.figure(figsize=(W2, 4.9))
axs = [fig.add_axes([0.17, 0.60, 0.22, 0.38]), fig.add_axes([0.615, 0.60, 0.12, 0.38]),
       fig.add_axes([0.17, 0.10, 0.22, 0.38]), fig.add_axes([0.79, 0.10, 0.15, 0.38])]

# (a) ICC(A,1) per construction, joint cohort
ax = axs[0]
for i, c in enumerate(CONS):
    r = J.loc[c]; st = cons_style(c)
    ax.errorbar(r.icc_a1, i, xerr=[[r.icc_a1 - r.icc_a1_lo], [r.icc_a1_hi - r.icc_a1]], fmt=st["marker"], color=st["color"],
                mfc=st["mfc"], mec=st["mec"], ms=4, capsize=1.5, elinewidth=0.8)
    ax.text(1.08, i, fm(r.icc_a1), va="center", ha="right", fontsize=FS - 1)
r = CI.loc["ceil1000_u"]
ax.errorbar(r.icc_a1, 6, xerr=[[r.icc_a1 - r.icc_a1_lo], [r.icc_a1_hi - r.icc_a1]], fmt=CEIL["marker"], color=CEIL["color"],
            mfc=CEIL["mfc"], mec=CEIL["mec"], ms=3.5, capsize=1.5, elinewidth=0.8)
ax.text(1.08, 6, fm(r.icc_a1), va="center", ha="right", fontsize=FS - 1)
ax.axhline(5.5, color=OI["grey"], lw=0.4, ls=":")
ax.set_yticks(range(7)); ax.set_yticklabels([cons_label(c) for c in CONS] + ["ceiling 1,000°/s, unfiltered‡"]); ax.invert_yaxis()
ax.set_xlim(0, 1.1); ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xlabel(f"ICC(A,1), R1 S1–S2 (joint cohort, n = {N})")
panel(ax, "a", x=-0.55)

# (b) path-specific DeltaICC (start -> end): all participants (filled) and omitting the most influential one (open)
PATHS = [("adm_u-mean_u", "mean-speed u → admitted u"), ("ceil1000_u-mean_u", "mean-speed u → ceiling u‡"),
         ("adm_u-ceil1000_u", "ceiling u‡ → admitted u"), ("pereye_u-adm_u", "admitted u → per-eye u"),
         ("adm_c-adm_u", "admitted u → admitted c"), ("mean_c-mean_u", "mean-speed u → mean-speed c"),
         ("pereye_c-pereye_u", "per-eye u → per-eye c"), ("pereye_c-mean_c", "mean-speed c → per-eye c")]
ax = axs[1]
P_all = PT.query("subset == 'all'").set_index("contrast"); P_wo = PT.query("subset == 'without_top'").set_index("contrast")
for k in ("adm_u-mean_u", "pereye_u-adm_u", "adm_c-adm_u", "mean_c-mean_u", "pereye_c-pereye_u", "pereye_c-mean_c"):
    assert abs(P_all.loc[k].delta - D.loc[k].delta) < 1e-9 and abs(P_all.loc[k].lo - D.loc[k].lo) < 1e-9
NW = int(P_wo.n.iloc[0])
for i, (k, lab) in enumerate(PATHS):
    r, w = P_all.loc[k], P_wo.loc[k]; assert int(r.n) == N
    col = OI["black"] if i < 3 else OI["grey"]
    ax.errorbar(r.delta, i - 0.17, xerr=[[r.delta - r.lo], [r.hi - r.delta]], fmt="D", color=col, ms=3.2, capsize=1.3,
                elinewidth=0.8)
    ax.errorbar(w.delta, i + 0.17, xerr=[[w.delta - w.lo], [w.hi - w.delta]], fmt="o", color=col, mfc="white", ms=3.2,
                capsize=1.3, elinewidth=0.8)
    for xt, (v, lo, hi) in ((1.05, (r.delta, r.lo, r.hi)), (2.15, (w.delta, w.lo, w.hi))):
        ax.text(xt, i, f"{fm(v)} [{fm(lo)}, {fm(hi)}]", va="center", ha="left", fontsize=FS - 1.5,
                transform=ax.get_yaxis_transform())
ax.text(1.05, -0.75, f"all (n = {N})", fontsize=FS - 1.5, fontweight="bold", transform=ax.get_yaxis_transform(), ha="left")
ax.text(2.15, -0.75, f"one omitted (n = {NW})", fontsize=FS - 1.5, fontweight="bold", transform=ax.get_yaxis_transform(),
        ha="left")
ax.axvline(0, color=OI["grey"], lw=0.6, ls=":")
ax.set_yticks(range(len(PATHS))); ax.set_yticklabels([p[1] for p in PATHS]); ax.set_ylim(len(PATHS) - 0.5, -0.5)
ax.set_xlim(-0.15, 0.6); ax.set_xlabel("ΔICC (end − start)")
panel(ax, "b", x=-1.75)

# (c) eligibility and random-exclusion comparator, floors 5 and 10 (mean-speed detector)
ax = axs[2]
rows = []
for k, fl in enumerate((5, 10)):
    own = B1.query("floor == @fl and filt == 'unfiltered' and sample == 'own'").iloc[0]
    inter = B1.query("floor == @fl and filt == 'unfiltered' and sample == 'intersection'").iloc[0]
    cen = B1.query("floor == @fl and filt == 'centred' and sample == 'intersection'").iloc[0]
    b = B2.loc[fl]; assert abs(cen.icc_a1 - b.centred) < 1e-9
    assert int(inter.n) == int(b.n)
    y0 = k * 5
    items = [(y0, own.icc_a1, own.icc_a1_lo, own.icc_a1_hi, f"unfiltered, own (n = {int(own.n)})", "s", "white"),
             (y0 + 1, inter.icc_a1, inter.icc_a1_lo, inter.icc_a1_hi, f"unfiltered, common (n = {int(inter.n)})", "s", "white"),
             (y0 + 2, b.random_mean, b.random_rep_lo, b.random_rep_hi, f"random excl.†, common (n = {int(b.n)})", "^", OI["grey"]),
             (y0 + 3, cen.icc_a1, cen.icc_a1_lo, cen.icc_a1_hi, f"centred, common (n = {int(b.n)})", "o", OI["verm"])]
    for y, v, lo, hi, lab, m, fc in items:
        col = OI["grey"] if m == "^" else OI["verm"]
        if m == "^":  # repetition range of the random comparator: dashed, no caps (not a CI)
            ax.hlines(y, lo, hi, color=col, lw=0.9, ls=(0, (2, 1.2)))
            ax.plot(v, y, m, color=col, mfc=fc, mec=col, ms=4)
        elif np.isfinite(lo):
            ax.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt=m, color=col, mfc=fc, mec=col, ms=4, capsize=1.5,
                        elinewidth=0.8, ls="none")
        else:
            ax.plot(v, y, m, color=col, mfc=fc, mec=col, ms=4)
        rows.append((y, lab))
    ax.text(1.02, y0 + 1.5, f"floor {fl}", fontsize=FS - 1, fontweight="bold", transform=ax.get_yaxis_transform(),
            ha="left", va="center")
    ax.text(1.02, y0 + 2.4, f"centred − random: {fm(b.d_c_r)} [{fm(b.d_c_r_lo)}, {fm(b.d_c_r_hi)}]\n"
            f"centred − unfiltered: {fm(b.d_c_u)} [{fm(b.d_c_u_lo)}, {fm(b.d_c_u_hi)}]",
            fontsize=FS - 1.5, transform=ax.get_yaxis_transform(), ha="left", va="center")
    ax.axhline(y0 + 0.5, color=OI["grey"], lw=0.4, ls=":")  # own cohort | common-cohort comparisons
    if k == 0:
        ax.axhline(3.9 + 0.1, color=OI["grey"], lw=0.4, ls="-", alpha=0.5)
ax.set_yticks([r[0] for r in rows]); ax.set_yticklabels([r[1] for r in rows]); ax.invert_yaxis()
ax.set_xlim(-0.05, 1.0); ax.set_xlabel("ICC(A,1), R1 S1–S2")
panel(ax, "c", x=-0.62)

# (d) detector-setting sensitivity (all R1 recordings re-extracted; one setting changed at a time)
ax = axs[3]
VAR = [("ref", "frozen settings"), ("vth30", "threshold 30°/s"), ("vth50", "threshold 50°/s"),
       ("sg5", "SG window 5"), ("sg9", "SG window 9"), ("gap1", "pairing gap 1"), ("gap5", "pairing gap 5")]
r0 = SV.loc["ref"]
assert int(r0.n_joint) == N and abs(r0.icc_mean_u - J.loc["mean_u"].icc_a1) < 1e-9 and abs(r0.icc_adm_u - J.loc["adm_u"].icc_a1) < 1e-9
for i, (v, lab) in enumerate(VAR):
    r = SV.loc[v]
    for x, st, dy in ((r.icc_mean_u, cons_style("mean_u"), -0.2), (r.icc_adm_u, cons_style("adm_u"), 0.0), (r.icc_ceil_u, CEIL, 0.2)):
        ax.plot(x, i + dy, st["marker"], color=st["color"], mfc=st["mfc"], mec=st["mec"], ms=3.5)
    ax.text(1.03, i, f"n = {int(r.n_joint)}", va="center", ha="left", fontsize=FS - 1.5, transform=ax.get_yaxis_transform())
ax.axhline(0.5, color=OI["grey"], lw=0.4, ls=":")
ax.set_yticks(range(len(VAR))); ax.set_yticklabels([v[1] for v in VAR]); ax.set_ylim(len(VAR) - 0.5, -0.5)
ax.set_xlim(0, 1.0); ax.set_xlabel("ICC(A,1), R1 S1–S2 (joint cohort)")
h = [plt.Line2D([], [], ls="none", **{k: st[k] for k in ("marker", "color", "mfc", "mec")}, ms=3.5)
     for st in (cons_style("mean_u"), cons_style("adm_u"), CEIL)]
ax.legend(h, ["mean-speed u", "admitted u", "ceiling u‡"], loc="upper center", bbox_to_anchor=(0.35, -0.2), ncol=3,
          fontsize=FS - 1.5, handletextpad=0.1, columnspacing=0.6, borderaxespad=0.1)
panel(ax, "d", x=-0.75)
save_fig(fig, "fig2_icc")
