"""Figure environments + captions for the main text (Figs 1-6, Table 1) -> figures/latex_includes.tex.
Every number in a caption is read from the named results file (contract C2); nothing is hand-typed."""
from paper_plot_style import res, f3
from pathlib import Path
import pandas as pd


def t3(x): return f3(x).replace("-", "$-$")
def t4(x): return f"{x:.4f}".replace("0.", ".", 1).replace("-", "$-$")
def b3(x): return "$-$.000" if -0.0005 < x < 0 else (".000" if 0 <= x < 0.0005 else t3(x))
def ci(v, lo, hi): return f"{t3(v)} [{t3(lo)}, {t3(hi)}]"
def one(df, **kw):
    m = pd.Series(True, index=df.index)
    for k, v in kw.items(): m &= df[k] == v
    assert m.sum() == 1, (kw, m.sum()); return df[m].iloc[0]


J = res("pw_p1_joint_icc.csv"); D = res("pw_p1_joint_delta.csv"); B2 = res("b2_selection.csv")
L = res("pw_p4_loo.csv"); AG = res("pw_p2_agreement.csv"); CV = res("rv_v2_ba_coverage.csv"); INJ = res("b3a_injection.csv").set_index("cond")
LS = res("rv_v4_loo_summary.csv").set_index("floor").loc[5]
M = res("b5_matched.csv"); T = res("b5_transfer.csv"); X = res("pw_p8_cross_transfer.csv")

j = lambda m: one(J, floor=5, comparison="R1_S1S2", cohort="joint6", method=m)
d = lambda c: one(D, floor=5, comparison="R1_S1S2", cohort="joint6", contrast=c)
N = int(j("mean_u").n)
mu, au, mc = j("mean_u"), j("adm_u"), j("mean_c")
ac, pc = j("adm_c"), j("pereye_c")
dadm, dcen = d("adm_u-mean_u"), d("mean_c-mean_u")
b5, b10 = B2.set_index("floor").loc[5], B2.set_index("floor").loc[10]
lu = one(L, floor=5, comparison="R1_S1S2", method="mean_u", cohort="own_y")
lc = one(L, floor=5, comparison="R1_S1S2", method="mean_c", cohort="own_x")
a1 = one(AG, floor=5, x="mean_c", y="pereye_c"); a2 = one(AG, floor=5, x="mean_u", y="pereye_u")
c1 = one(CV, floor=5, x="mean_c", y="pereye_c"); c2 = one(CV, floor=5, x="mean_u", y="pereye_u")
def cov(c): return (f"{int(c.n_inside):,} of {int(c.n_rec):,} recordings ({100 * c.coverage:.1f}\\%, 95\\% CI {100 * c.coverage_lo:.1f}--{100 * c.coverage_hi:.1f}\\%) inside the limits; "
                    f"2.5th and 97.5th percentiles of the differences {t3(c.p025)} [{t3(c.p025_lo)}, {t3(c.p025_hi)}] and {t3(c.p975)} [{t3(c.p975_lo)}, {t3(c.p975_hi)}]")
mt, mr = one(M, task="TEX", variant="primary"), one(M, task="RAN", variant="primary")
tt, tr = one(T, task="TEX", variant="primary"), one(T, task="RAN", variant="primary")
xc = X[(X.task == "TEX") & (X.variant == "primary") & X.kind.isin(["same", "cross"])]
n1a, n1b = INJ.loc["N1_k0.85"], INJ.loc["N1_k1.15"]
f3d = INJ.loc["F3"]
CE = res("rv_v7_ceiling_icc.csv").query("tag == 'floor5_ceil1000'").set_index("method")
PT = res("rv_v7_paths.csv"); pa = PT.query("subset == 'all'").set_index("contrast"); pw = PT.query("subset == 'without_top'").set_index("contrast")
SV = res("rv_v9_sens_summary.csv").set_index("variant")
ce = CE.loc["ceil1000_u"]; assert int(ce.n) == N and abs(pa.loc["adm_u-mean_u"].delta - dadm.delta) < 1e-9

S = []
S.append(rf"""% === Fig 1: schematic (FigEdit; source figures/fig1_schematic.svg) ===
\begin{{figure}}[t]
\centering
\includegraphics[width=\textwidth]{{figures/fig1_schematic.pdf}}
\caption{{\textbf{{Event construction for the interocular reading-event contrast $A$.}}
Left- and right-eye gaze traces are converted to per-eye angular velocity. Saccades are detected either on the mean of the two per-eye speed
magnitudes (mean-speed detector) or separately in each eye (per-eye detector). The admission gate keeps a mean-speed event only when it overlaps a
matched left/right per-eye event pair, whether or not that pair passes its own exclusion rules; detection windows, peak samples and saccade classes of admitted events are unchanged. Centring is an optional
branch that keeps events whose log left/right peak-velocity ratio $l$ lies within $\log 2$ of the recording median. $A$ is half the difference between
the median $l$ of return sweeps and of progressive saccades in a recording. A recording enters only when it meets the event floor; at floor 5 it
must contain at least 5 return sweeps and 20 progressive saccades. Inset: R1 S1--S2 ICC(A,1) of unfiltered mean-speed events, of admitted events and of mean-speed events with a post hoc 1,000\textdegree/s peak-velocity ceiling,
in the joint cohort, the participants whose two R1 sessions were eligible under all six constructions ($n$ = {N}, floor 5), and the first two again with the
single most influential participant omitted ($n$ = {N - 1}; Results).
Detector thresholds and the six-construction lookup are given in Methods.}}
\label{{fig:schematic}}
\end{{figure}}
% inset values: unfiltered mean-speed {t3(mu.icc_a1)}, admitted {t3(au.icc_a1)} (pw_p1_joint_icc.csv), ceiling {t3(ce.icc_a1)} (rv_v7_ceiling_icc.csv), omitted {t3(LS.icc_mean_u_wo_top)} vs {t3(LS.icc_adm_u_wo_top)} (rv_v4_loo_summary.csv)
""")

S.append(rf"""% === Fig 2: repeatability by construction (gen_fig2_icc.py) ===
\begin{{figure}}[t]
\centering
\includegraphics[width=\textwidth]{{figures/fig2_icc.pdf}}
\caption{{\textbf{{Short-term repeatability of $A$ under each event construction, with one participant omitted, under recording eligibility and under detector settings.}}
\textbf{{a}}, ICC(A,1) between R1 sessions 1 and 2 with 95\% participant-bootstrap CI (5,000 draws) for the six constructions, in the joint
cohort of participants eligible under all six ($n$ = {N}, floor 5); the unit is the recording-level $A$ of each session. Open squares, unfiltered; filled circles, centred;
orange, mean-speed events; blue, admitted mean-speed events; green, per-eye events. $\ddagger$, post hoc comparator not in the frozen plan: mean-speed events after removing those with
either per-eye peak velocity above 1,000\textdegree/s (same cohort; {t3(ce.icc_a1)} [{t3(ce.icc_a1_lo)}, {t3(ce.icc_a1_hi)}]).
\textbf{{b}}, Path-specific contrasts $\Delta$ICC (end $-$ start; u, unfiltered; c, centred) with 95\% CI from paired participant-bootstrap draws. Filled diamonds, all participants ($n$ = {N});
open circles, the same contrasts with the participant who contributed most to the session differences under unfiltered mean-speed events omitted ($n$ = {N - 1}). The contrasts share endpoints and are not
an additive partition of reliability. With windows, peaks and classes unchanged, admission raised ICC from {t3(mu.icc_a1)} to {t3(au.icc_a1)},
$\Delta$ICC = {ci(dadm.delta, dadm.lo, dadm.hi)}, or {ci(pw.loc["adm_u-mean_u"].delta, pw.loc["adm_u-mean_u"].lo, pw.loc["adm_u-mean_u"].hi)} with that participant omitted;
the ceiling comparator reproduced the increase, and admitted minus ceiling was {ci(pa.loc["adm_u-ceil1000_u"].delta, pa.loc["adm_u-ceil1000_u"].lo, pa.loc["adm_u-ceil1000_u"].hi)}. \textbf{{c}}, Mean-speed detector under different eligibility, floors 5 and 10:
each construction's own cohort, the common cohort eligible both unfiltered and centred, and a count-matched random-exclusion comparator ($\dagger$)
that keeps, in each participant, session and class, the number of events retained by centring (mean of 200 repetitions; its dashed bars without caps show the
2.5--97.5 percentile range of repetitions, not a CI). Dotted lines separate the own-cohort row from the common-cohort comparisons. Right-hand text: contrasts with 95\% CI in the common cohort ($n$ = {int(b5.n)} at floor 5;
$n$ = {int(b10.n)} at floor 10). \textbf{{d}}, Post hoc detector-setting sensitivity: all R1 recordings re-extracted with one detector setting changed at a time (velocity threshold, Savitzky--Golay window, maximum
left/right pairing gap in samples), each in its own joint cohort (right-hand $n$; floor 5; point estimates). Unfiltered mean-speed ICC ranged from {t3(SV.icc_mean_u.min())} to {t3(SV.icc_mean_u.max())},
admitted from {t3(SV.icc_adm_u.min())} to {t3(SV.icc_adm_u.max())} and ceiling from {t3(SV.icc_ceil_u.min())} to {t3(SV.icc_ceil_u.max())}; re-extraction with the frozen settings reproduced the analysed events exactly.}}
\label{{fig:icc}}
\end{{figure}}
""")

S.append(rf"""% === Fig 3: participant influence (gen_fig3_influence.py) ===
\begin{{figure}}[t]
\centering
\includegraphics[width=\textwidth]{{figures/fig3_influence.pdf}}
\caption{{\textbf{{The low ICC of unfiltered mean-speed events and the admission contrast were carried by a few participants.}}
\textbf{{a}}, Share of the summed squared session differences of $A$, $\Sigma(S1 - S2)^2$ (R1, floor 5), for the ten largest contributors under
unfiltered mean-speed events in their own cohort ($n$ = {int(lu.n)}); the top three (black) account for {100 * lu.top3_sqdiff_share:.0f}\%.
\textbf{{b}}, ICC(A,1) after removing the $k$ most influential participants (largest absolute leave-one-out shift in ICC), ranked separately within
each mean-speed construction's own cohort ($n$ before removal: unfiltered {int(lu.n)}, centred {int(lc.n)}); the unfiltered ICC changed from
{t3(lu.icc)} to {t3(lu.icc_drop_top2)} after two removals, and the centred ICC stayed between {t3(lc.icc)} (no removal) and
{t3(max(lc.icc, lc.icc_drop_top1, lc.icc_drop_top2, lc.icc_drop_top3, lc.icc_drop_top5))}. Because cohorts and rankings differ, a and b do not isolate an effect of centring.
\textbf{{c}}, Admission contrast (ICC of admitted minus unfiltered mean-speed events) in the joint cohort ($n$ = {int(LS.n)}) after each single-participant omission,
against the omitted participant's share of the summed squared deviations of the session differences from their mean under unfiltered mean-speed events.
Dashed line, full-cohort contrast {ci(LS.delta_full, LS.delta_full_lo, LS.delta_full_hi)}; star, the omission that gave {ci(LS.delta_wo_top, LS.delta_wo_top_lo, LS.delta_wo_top_hi)}
($n$ = {int(LS.n) - 1}); brackets, 95\% participant-cluster bootstrap CIs. \textbf{{d}}, \textbf{{e}}, Session 2 against session 1 $A$ in the joint cohort for
unfiltered mean-speed (d, ICC = {t3(LS.icc_mean_u)}) and admitted unfiltered events (e, ICC = {t3(LS.icc_adm_u)}) on common axes; star, the participant omitted in c;
dotted line, identity. \textbf{{f}}, For that participant, the median $l = \ln(\mathrm{{PV}}_L/\mathrm{{PV}}_R)$ of return sweeps (RS) and of progressive saccades
in each R1 session under the two constructions; numbers, RS counts. Intervals are descriptive. Participant deletion is a sensitivity analysis here, not a cleaning rule.}}
\label{{fig:influence}}
\end{{figure}}
""")

S.append(rf"""% === Fig 4: between-construction agreement (gen_fig4_bland_altman.py) ===
\begin{{figure}}[t]
\centering
\includegraphics[width=\textwidth]{{figures/fig4_bland_altman.pdf}}
\caption{{\textbf{{Agreement of recording-level $A$ between constructions.}} Bland--Altman plots of all recordings eligible under both
constructions (floor 5, all rounds): vertical axis, per-eye minus mean-speed $A$; horizontal axis, mean of the two. \textbf{{a}}, Centred per-eye against centred mean-speed construction ({int(a1.n_rec):,} recordings,
{int(a1.n_pid)} participants): bias {t4(a1.bias)} [{t4(c1.bias_lo)}, {t4(c1.bias_hi)}], nominal 95\% limits of agreement {t3(a1.loa_lo)} [{t3(c1.loa_lo_lo)}, {t3(c1.loa_lo_hi)}] and {t3(a1.loa_hi)} [{t3(c1.loa_hi_lo)}, {t3(c1.loa_hi_hi)}];
{cov(c1)}.
\textbf{{b}}, Unfiltered per-eye against unfiltered mean-speed construction ({int(a2.n_rec):,} recordings, {int(a2.n_pid)} participants):
bias {t3(a2.bias)} [{t3(c2.bias_lo)}, {t3(c2.bias_hi)}], limits {t3(a2.loa_lo)} [{t3(c2.loa_lo_lo)}, {t3(c2.loa_lo_hi)}] and {t3(a2.loa_hi)} [{t3(c2.loa_hi_lo)}, {t3(c2.loa_hi_hi)}]; {cov(c2)}.
The few large differences widen the limits in b beyond the central spread; the inset shows the recordings with $|$difference$|$ $\leq$ .1. This panel compares detectors at the same filtering.
Solid line, bias; dashed lines, limits of agreement (bias $\pm$ 1.96 SD of the differences); dotted lines, empirical 2.5th and 97.5th percentiles.
All summaries are pooled over recordings, so participants with more eligible recordings carry more weight; brackets are 95\% participant-cluster bootstrap CIs (5,000 draws; Supplementary Table S13).
The recording sets differ between panels, so the panel widths do not isolate the effect of filtering. Pearson $r$ (in-panel text) carries a 95\% participant-cluster bootstrap CI. These are agreement between constructions, not
between-session repeatability.}}
\label{{fig:agreement}}
\end{{figure}}
""")

S.append(rf"""% === Fig 5: synthetic one-eye faults (gen_fig5_injection.py) ===
\begin{{figure}}[t]
\centering
\includegraphics[width=\textwidth]{{figures/fig5_injection.pdf}}
\caption{{\textbf{{Abstention by the centring rule under synthetic one-eye faults.}}
\textbf{{a}}, Recall: fraction of injections reaching centring (overlapped by a detected event that passed the class and exclusion rules) that the centring rule abstained from, with
95\% participant-cluster bootstrap CI; denominator, injections reaching centring. Dashes show the fraction reaching centring; denominator, all injections. The frozen criterion for judged conditions (filled) was
a point estimate of at least .80 (dashed line) with lower bound $>$ .70; open blue markers are unjudged magnitudes, grey markers are descriptive contamination
conditions (one-eye slowing by a factor of 3 applied to 5--30\% of eligible events). Of the judged one-eye conditions, only the tested 4$^\circ$ and 8$^\circ$ steps met the criterion; one-eye slowing (factor 3) and a
6$^\circ$ peak glitch did not. For one-eye dropout, only {int(f3d.n_admitted)} of {int(f3d.injections):,} injections reached centring (most were removed at the missingness layer),
so its recall was not interpreted. \textbf{{b}}, Specificity control: the within-saccade displacement of one randomly chosen eye was scaled by
$\kappa$, a negative-control perturbation ($|\Delta l| \approx .15 < \log 2$). Excess rejection is the rejection rate of the same
target events reaching centring after scaling minus that in the uninjected baseline pass ($\kappa$ = 0.85: {ci(n1a.excess, n1a.excess_lo, n1a.excess_hi)};
$\kappa$ = 1.15: {ci(n1b.excess, n1b.excess_lo, n1b.excess_hi)}; participant-cluster bootstrap).
Injection counts, the fraction reaching centring, recall or excess rejection, untouched-event changes and outcome layers per condition are in Supplementary Table S9.}}
\label{{fig:injection}}
\end{{figure}}
""")

S.append(rf"""% === Fig 6: amplitude-matched peak-velocity difference (gen_fig6_matched.py) ===
\begin{{figure}}[t]
\centering
\includegraphics[width=\textwidth]{{figures/fig6_matched.pdf}}
\caption{{\textbf{{Amplitude-matched left/right peak-velocity difference $D$ and its association with $A$.}}
$D = \tfrac12(l_{{\mathrm{{leftward}}}} - l_{{\mathrm{{rightward}}}})$ for a leftward and a rightward saccade matched on amplitude and position within a
recording (dimensionless, log units); positive $D$ means a relatively faster left eye in leftward than in rightward saccades. Participant $D$ requires
both R1 sessions and at least 5 matched pairs in total; its repeatability requires at least 5 pairs in each session, hence the smaller $n$ below.
\textbf{{a}}, Participant $D$ (mean of the two R1 session medians) in the reading (TEX; $n$ = {int(mt.n_pid)}, mean {ci(mt.mean_D, mt.mean_D_lo, mt.mean_D_hi)})
and random-saccade tasks (RAN; $n$ = {int(mr.n_pid)}, mean {ci(mr.mean_D, mr.mean_D_lo, mr.mean_D_hi)}); histograms are normalized to density; values beyond the axis range are placed in
the edge bins and counted in the legend. Different TEX and RAN sample sizes do not by themselves establish a task difference. \textbf{{b}}, Participant $D$ against TEX $A$ (R1 two-session mean, centred mean-speed construction):
TEX $D$, Spearman $\rho$ = {ci(tt.rho, tt.rho_lo, tt.rho_hi)} ($n$ = {int(tt.n)}); RAN $D$, $\rho$ = {ci(tr.rho, tr.rho_lo, tr.rho_hi)}
($n$ = {int(tr.n)}); $D$ values below $-$.3 are drawn as triangles at the lower edge, and their number is given in the panel. \textbf{{c}}, Spearman $\rho$ with 95\% participant-bootstrap CI between single-session TEX $A$ and single-session $D$ (same- and cross-session pairings,
$n$ = {int(xc.n.min())}--{int(xc.n.max())}) and, below the line, for the pooled associations of panel b; cohorts and averaging differ between the two groups, so their magnitudes are
compared descriptively. TEX $D$ repeatability (R1 S1--S2):
ICC(A,1) = {ci(mt.icc_a1, mt.icc_a1_lo, mt.icc_a1_hi)}, $n$ = {int(mt.icc_n)}. RAN repeatability was not estimable, because no participant had at least 5 pairs in each session, and the RAN session-specific associations
rest on small samples (Supplementary Table S8).}}
\label{{fig:matched}}
\end{{figure}}
""")

S.append(r"""% === Table 1: longitudinal repeatability (gen_table1_longitudinal.py) ===
\input{figures/TABLE_1.tex}
""")

Path(__file__).with_name("latex_includes.tex").write_text("\n".join(S))
print("Saved: figures/latex_includes.tex")
