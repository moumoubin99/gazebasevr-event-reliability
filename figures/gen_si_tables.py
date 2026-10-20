"""Supplementary Tables S1-S13, generated from results/*.csv (no hand-typed values).
Writes figures/si/TABLE_S<k>.tex and figures/SI_TABLES.tex (\\input list)."""
from paper_plot_style import res, f3, CONS, cons_label
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path(__file__).with_name("si"); OUT.mkdir(exist_ok=True)


def t3(x):
    return "--" if pd.isna(x) else f3(x).replace("-", "$-$")


def ci(v, lo, hi):
    if pd.isna(v): return "--"
    return f"{t3(v)} [{t3(lo)}, {t3(hi)}]" if not (pd.isna(lo) and pd.isna(hi)) else f"{t3(v)} [--]"


def ni(x):
    return "--" if pd.isna(x) else f"{int(x):,}".replace(",", "{,}")


def pct(x, d=1):
    return "--" if pd.isna(x) else f"{100 * x:.{d}f}"


def table(k, caption, blocks, note=None, size=r"\footnotesize"):
    """blocks: list of (colspec, header_row, body_rows[, subtitle])."""
    s = [r"\begin{table}[!htbp]", r"\centering", rf"\caption{{{caption}}}", rf"\label{{tab:S{k}}}", size]
    for b in blocks:
        spec, head, body = b[:3]
        if len(b) > 3: s += [rf"\par\smallskip\textit{{{b[3]}}}\par\smallskip"]
        s += [r"\begin{adjustbox}{max width=\linewidth}", rf"\begin{{tabular}}{{{spec}}}", r"\toprule", head + r" \\", r"\midrule"] + [r + r" \\" for r in body] + [r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}"]
    if note: s += [r"\par\smallskip", rf"\begin{{minipage}}{{0.95\linewidth}}\scriptsize {note}\end{{minipage}}"]
    s += [r"\end{table}"]
    (OUT / f"TABLE_S{k}.tex").write_text("\n".join(s) + "\n")
    return k


def row(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items(): m &= (df[k] == v).values
    assert m.sum() == 1, kw
    return df[m].iloc[0]


made = []
PATHS = ["adm_u-mean_u", "pereye_u-adm_u", "adm_c-adm_u", "mean_c-mean_u", "pereye_c-pereye_u", "pereye_c-mean_c"]
OTHER = ["pereye_u-mean_u", "adm_c-mean_c", "pereye_u-mean_c"]
CL = {c: cons_label(c) for c in CONS}
def clab(contrast):
    a, b = contrast.split("-"); return f"{CL[a]} $-$ {CL[b]}"

# ---------------- S1 random-exclusion comparator (b2_selection)
B2 = res("b2_selection.csv")
made.append(table(1, "Count-matched random-exclusion comparator (mean-speed detector, R1 S1--S2, common cohort per floor). "
    "In each of 200 repetitions, random exclusion keeps, within each participant $\\times$ session $\\times$ class (RS, progressive), a random subset of "
    "unfiltered events of the same size as the set retained by the centred rule; `random' is the mean ICC over repetitions, with the 2.5th--97.5th "
    "percentile range of repetitions in brackets (a repetition range, not a confidence interval). The cohort is the participants eligible both unfiltered and "
    "centred at that floor (common cohort). Contrasts are ICC(A,1) differences with 95\\% participant-bootstrap intervals; the random "
    "comparator enters the centred $-$ random contrast as its repetition mean.",
    [("rcccccc", r"Floor & $n$ & Unfiltered & Centred & Random [rep.\ range] & Centred $-$ unfiltered & Centred $-$ random",
      [f"{int(r.floor)} & {int(r.n)} & {t3(r.unfiltered)} & {t3(r.centred)} & {t3(r.random_mean)} [{t3(r.random_rep_lo)}, {t3(r.random_rep_hi)}] & "
       f"{ci(r.d_c_u, r.d_c_u_lo, r.d_c_u_hi)} & {ci(r.d_c_r, r.d_c_r_lo, r.d_c_r_hi)}" for r in B2.itertuples()])]))

# ---------------- S2 floor sensitivity of joint-cohort ICCs and path contrasts (pw_p1_joint_*)
JI = res("pw_p1_joint_icc.csv").query("comparison == 'R1_S1S2' and cohort == 'joint6'")
JD = res("pw_p1_joint_delta.csv").query("comparison == 'R1_S1S2' and cohort == 'joint6'")
FLS = sorted(JI.floor.unique())
nf = {fl: int(JI[JI.floor == fl].n.iloc[0]) for fl in FLS}
head = "Construction / contrast & " + " & ".join(f"floor {fl} ($n$ = {nf[fl]})" for fl in FLS)
body = [f"{CL[m]} & " + " & ".join(ci(*row(JI, floor=fl, method=m)[["icc_a1", "icc_a1_lo", "icc_a1_hi"]]) for fl in FLS) for m in CONS]
body += [r"\midrule \multicolumn{%d}{l}{\textit{Path-specific contrasts ($\Delta$ICC)}}" % (len(FLS) + 1)]
for c in PATHS + OTHER:
    body.append(f"{clab(c)} & " + " & ".join(ci(*row(JD, floor=fl, contrast=c)[["delta", "lo", "hi"]]) for fl in FLS))
made.append(table(2, "Eligibility-floor sensitivity of short-term repeatability (R1 S1--S2, joint cohort eligible under all six constructions at each floor). "
    "Entries are ICC(A,1) [95\\% CI]; contrasts are differences on the same participants and bootstrap draws. Contrasts are path-specific and do not "
    "form an additive partition. Floor 5 is the primary analysis.", [("l" + "c" * len(FLS), head, body)]))

# ---------------- S3 mean-speed vs per-eye detection: denominators (b1_denominators*)
D0, D1 = res("b1_denominators.csv").iloc[0], res("b1_denominators_pereye.csv").iloc[0]
rows = [("Recordings", "recordings"), ("Events detected", "events_detected"), ("Events in RS/progressive pools", "events_pool"),
        ("RS pool", "rs_pool"), ("Progressive pool", "prog_pool"), ("Events kept by centring", "events_kept")]
body = [f"{lab} & {ni(D0[c])} & {ni(D1[c])}" for lab, c in rows] + [f"Retention after centring & {t3(D0.retention)} & {t3(D1.retention)}"]
for filt in ("unfilt", "centred"):
    for fl in (1, 5, 10):
        body.append(f"{'Unfiltered' if filt == 'unfilt' else 'Centred'}, floor {fl}: eligible recordings / R1 pairs & "
                    f"{ni(D0[f'eligible_rec_{filt}_f{fl}'])} / {ni(D0[f'pairs_R1_{filt}_f{fl}'])} & {ni(D1[f'eligible_rec_{filt}_f{fl}'])} / {ni(D1[f'pairs_R1_{filt}_f{fl}'])}")
made.append(table(3, "Event and eligibility denominators for the mean-speed detector and per-eye detection (all rounds, TEX). "
    "Event counts are pooled over all recordings. Floor $k$: a recording is eligible when it has at least $k$ return sweeps and at least 20 progressive "
    "saccades (after centring, for the centred rows). R1 pairs: participants with both R1 sessions eligible (each detector's own cohort).",
    [("lcc", "Quantity & Mean-speed & Per-eye", body)]))

# ---------------- S4 event composition and centring by pool (pw_p3_admission_events)
P3 = res("pw_p3_admission_events.csv")
PL = {"mean_pool": "mean-speed", "adm_pool": "admitted", "pereye_pool": "per-eye",
      "mean_pool|pe_support=True": "mean-speed, with per-eye support", "mean_pool|pe_support=False": "mean-speed, without per-eye support"}
made.append(table(4, "Event composition by construction pool (all rounds, TEX). $|l|$: absolute value of $l = \\log(\\text{left}/\\text{right peak velocity})$, that is, the distance of the interocular log ratio from equality, "
    "summarised over all events of the pool without recording centring; kept: fraction retained by the centred rule ($|l - \\mathrm{median}(l)| < \\log 2$, "
    "with the recording median computed within each pool; the admitted pool's median is recomputed from admitted events, so its kept fractions can differ "
    "slightly from those of the mean-speed events it contains). Per-eye events are matched left/right pairs of "
    "single-eye detections. The admitted pool is the mean-speed pool restricted to events whose window overlaps the window of at least one matched left/right "
    "per-eye event in the same recording (per-eye support); its windows, peaks and classes are those of the mean-speed events.",
    [("llrccc", r"Pool & Class & Events & Kept & Median $|l|$ & P95 $|l|$",
      [f"{PL[r.pool]} & {'RS' if r.cls == 'rs' else 'progressive'} & {ni(r.n)} & {t3(r.kept_frac)} & {t3(r.med_absl)} & {t3(r.p95_absl)}" for r in P3.itertuples()])]))

# ---------------- S5 retention and eligibility transitions (pw_p5_*)
R5, E5 = res("pw_p5_retention.csv"), res("pw_p5_eligibility.csv")
b1 = [f"{'mean-speed' if r.detector == 'mean' else 'per-eye'} & {r.cls if r.cls != 'prog' else 'progressive'} & {ni(r.n_rec)} & {t3(r['mean'])} & "
      f"{t3(r.p05)} & {t3(r["median"])} & {t3(r.p95)} & {t3(r.frac_rec_below_80)} & {t3(r.frac_rec_below_50)}" for _, r in R5.iterrows()]
b2 = [f"{'mean-speed' if r.detector == 'mean' else 'per-eye'} & {int(r.floor)} & R{int(r['round'])} & {'recording' if r.unit == 'recording' else 'R1 pair'} & "
      f"{ni(r.n)} & {ni(r.both)} & {ni(r.unfilt_only)} & {ni(r.centred_only)} & {ni(r.neither)}" for _, r in E5.iterrows()]
made.append(table(5, "Per-recording event retention under centring, and eligibility transitions between the unfiltered and centred markers "
    "(TEX; mean-speed detector and per-eye detection, not the admitted construction). (a) Fraction of a recording's RS or progressive events retained by the "
    "centred rule, over all recordings (all rounds); columns $<$.80 and $<$.50 give the fraction of recordings below that retention. (b) Counts of units "
    "eligible under both markers, only one, or neither, at floor $k$ ($\\ge k$ RS and $\\ge$ 20 progressive saccades); unit `recording' counts recordings of the "
    "given round, unit `R1 pair' counts participants with both R1 sessions eligible.",
    [("llccccccc", r"Detector & Class & Recordings & Mean & P5 & Median & P95 & $<$.80 & $<$.50", b1, "(a) Retention fraction per recording"),
     ("lclcccccc", r"Detector & Floor & Round & Unit & $n$ & Both & Unfiltered only & Centred only & Neither", b2, "(b) Eligibility under unfiltered and centred markers")]))

# ---------------- S6 mixed synthetic contamination (REC) detail (pw_p6_*, b3a_recovery)
P6 = res("pw_p6_rec_injected.csv"); RL = res("pw_p6_rec_layers.csv"); RV = res("b3a_recovery.csv")
e = P6[P6.injected_recs.notna()]
b1 = [f"{int(r.floor)} & {r.filt} & {ni(r.injected_recs)} & {ni(r.eligible_base)} & {ni(r.eligible_rec)} & {ni(r.lost)} & {ni(r.gained)} & "
      f"{t3(r.absdA_med)} & {t3(r.absdA_max)}" for r in e.itertuples()]
x = P6[P6.filt.str.startswith("ICC_injected") | P6.filt.str.startswith("dICC_injected")]
XL = {"ICC_injected_pids base_c": "ICC, uncontaminated, centred", "ICC_injected_pids base_u": "ICC, uncontaminated, unfiltered",
      "ICC_injected_pids rec_c": "ICC, contaminated, centred", "ICC_injected_pids rec_u": "ICC, contaminated, unfiltered",
      "dICC_injected_pids rec-base c": r"$\Delta$ICC contaminated $-$ uncontaminated, centred",
      "dICC_injected_pids rec-base u": r"$\Delta$ICC contaminated $-$ uncontaminated, unfiltered"}
b2 = [f"{XL[r.filt]} & {int(r.n_shared)} & {ci(r.icc_a1, r.icc_a1_lo, r.icc_a1_hi)}" for r in x.itertuples()]
rv = RV[RV.run.isin(["base", "rec"]) & RV.n.notna()]
b3 = [f"{'uncontaminated' if r.run == 'base' else 'contaminated'} & {r.filt} & {int(r.n)} & {ci(r.icc_a1, r.icc_a1_lo, r.icc_a1_hi)}" for r in rv.itertuples()]
b3 += [f"{'contaminated $-$ uncontaminated' if r.run == 'rec-base' else 'contaminated'} & {'centred' if r.filt == 'centred' else 'centred $-$ unfiltered'} & {int(r.d_n)} & "
       f"{ci(r.d_delta, r.d_lo, r.d_hi)}" for r in RV[RV.d_delta.notna()].itertuples()]
lay_cols = ["retained", "abstained", "not_detected", "out_of_class", "duration", "missingness"]
b4 = [f"{r.kind} & {int(r.n)} & " + " & ".join(t3(r[c]) for c in lay_cols) for _, r in RL.iterrows()]
made.append(table(6, "Mixed synthetic contamination (REC): one R1 session of each of 20 participants received one-eye slowing (factor 3) "
    "on 30\\% of eligible events and 4$^\\circ$ one-eye steps in 10\\% of fixations. Four runs: uncontaminated and contaminated data, each with the centred and the unfiltered marker (mean-speed detector, floor 5 unless stated). "
    "(a) Eligibility and marker displacement $|\\Delta A|$ (contaminated minus uncontaminated $A$ of the same recording) in injected recordings; "
    "(b) R1 S1--S2 ICC restricted to injected participants with S1--S2 pairs in all four runs (exposed-participant protection was not demonstrated); "
    "(c) R1 S1--S2 ICCs in the four-run common cohort, i.e. all participants with S1--S2 pairs in all four runs (this cohort differs in membership from the "
    "joint cohort of Fig.~2 even where $n$ coincides, so its baseline ICCs differ from Fig.~2a); "
    "(d) outcome layer of each injected fault, as fractions of all injections of that kind (rows sum to 1).",
    [("clccccccc", r"Floor & Marker & Injected & Eligible (base) & Eligible (REC) & Lost & Gained & Median $|\Delta A|$ & Max $|\Delta A|$", b1, "(a)"),
     ("lcc", r"Quantity & $n$ & Estimate [95\% CI]", b2, "(b)"),
     ("llcc", r"Run & Marker & $n$ & ICC(A,1) or $\Delta$ [95\% CI]", b3, "(c)"),
     ("lccccccc", r"Fault & $n$ & Retained & Abstained & Not detected & Out of class & Duration & Missingness", b4, "(d)")]))

# ---------------- S7 overlap of injected faults (pw_p7_overlap)
P7 = res("pw_p7_overlap.csv")
a = P7[P7.kind.isin(["F1", "F2a"])]
b1 = [f"{r.cond.replace('C_p', 'contamination ') + ('\\%' if r.cond.startswith('C_p') else '')} & {r.kind} & {r.subset.replace('_', '-')} & {ni(r.injections)} & {t3(r.admitted_rate)} & "
      f"{ci(r.recall, r.recall_lo, r.recall_hi)} & {ni(r.n_admitted)}" for r in a.itertuples()]
d = P7[P7.kind.str.startswith("absdA")]
b2 = [f"{r.cond.replace('C_p', 'contamination ')}\\% & {'centred' if 'centred' in r.kind else 'unfiltered'} & {r.subset.replace('rec_', '').replace('_', ' ')} & {ni(r.injections)} & "
      f"{t3(r.absdA_med)} & {t3(r.absdA_p95)}" for r in d.itertuples()]
made.append(table(7, "Sensitivity of abstention recall and marker displacement to overlapping injections. Contamination $p$\\%: one-eye slowing (F2a, factor 3) applied to $p$\\% of a recording's eligible events; REC: mixed synthetic contamination (Supplementary Table S6). Overlapping: the injected fault's "
    "modified samples or event window overlap another injection. Reached: fraction of all injections that reached the centring step; recall: fraction of "
    "injections reaching centring that the centred rule abstained from, with 95\\% participant-cluster bootstrap CI. (b) $|\\Delta A|$: absolute difference between "
    "the contaminated and the uncontaminated (baseline pass) $A$ of the same R1 recording; $n$ counts recordings, split by whether any of their injections overlapped.",
    [("lllcccc", r"Condition & Fault & Subset & Injections & Reached & Recall [95\% CI] & Reached $n$", b1, "(a) Recall by overlap subset"),
     ("lllccc", r"Condition & Marker & Recordings & $n$ & Median $|\Delta A|$ & P95 $|\Delta A|$", b2, "(b) Marker displacement in recordings with or without overlap")]))

# ---------------- S8 amplitude-matched D (b5_*, pw_p8)
M5, T5, X8 = res("b5_matched.csv"), res("b5_transfer.csv"), res("pw_p8_cross_transfer.csv")
VL = {"primary": "primary", "horizontal_pv": "horizontal PV", "centred": "centred events", "band_3_10": r"3--10$^\circ$ band", "band_15_25": r"15--25$^\circ$ band"}
b1 = [f"{r.task} & {VL[r.variant]} & {ni(r.n_pairs)} & {ni(r.n_pid)} & {ci(r.mean_D, r.mean_D_lo, r.mean_D_hi)} & {t3(r.median_D)} & {t3(r.frac_pos)} & "
      f"{ni(r.icc_n) if r.icc_n > 0 else '--'} & {ci(r.icc_a1, r.icc_a1_lo, r.icc_a1_hi)}" for r in M5.itertuples()]
b2 = [f"{r.task} & {VL[r.variant]} & {ni(r.n)} & {ci(r.rho, r.rho_lo, r.rho_hi)} & {'judged' if r.judged else ('inconclusive' if r.inconclusive else 'sensitivity')}" for r in T5.itertuples()]
KL = {"cross": "cross-session", "same": "same-session", "cross_avg_common": "cross-session, mean (common $n$)",
      "same_avg_common": "same-session, mean (common $n$)", "same-cross": "same $-$ cross"}
b3 = [f"{r.task} & {VL[r.variant]} & {KL[r.kind]} & {'S%d' % r.A_session if not pd.isna(r.A_session) else '--'} & {'S%d' % r.D_session if not pd.isna(r.D_session) else '--'} & "
      f"{ni(r.n)} & {ci(r.rho, r.rho_lo, r.rho_hi)}" for r in X8.itertuples()]
made.append(table(8, "Amplitude-matched left/right peak-velocity difference $D$ (participant level) and its association with the reading marker $A$ "
    "(R1, centred mean-speed construction). Within a recording, leftward and rightward saccades are paired greedily on amplitude, mid-point and direction, "
    "and $D = \\tfrac12(l_\\text{leftward} - l_\\text{rightward})$ per pair (positive: relatively faster left eye in leftward saccades). Participant $D$ is the mean "
    "of the two R1 session medians and requires both sessions and $\\ge$ 5 pairs in total; the ICC requires $\\ge$ 5 pairs in each session. (a) $D$ by task and variant, with R1 S1--S2 repeatability where estimable; (b) Spearman $\\rho$ between "
    "participant $D$ and two-session mean $A$; (c) session-specific associations. RAN session-specific estimates rest on 6--11 participants and are "
    "inconclusive; RAN repeatability was not estimable because no participant had $\\ge$ 5 pairs in each session. Pairs: all matched pairs of the task "
    "and variant, including those of participants who did not qualify for participant $D$. Status: an association is judged only for the primary variant "
    "with $n \\ge$ 100 participants (frozen rule); it is inconclusive when $n <$ 100, and other variants are sensitivity analyses.",
    [("llrrcccrc", r"Task & Variant & Pairs & $n$ & Mean $D$ [95\% CI] & Median & $D>0$ & ICC $n$ & ICC(A,1) [95\% CI]", b1, "(a)"),
     ("llrcl", r"Task & Variant & $n$ & $\rho$ [95\% CI] & Status", b2, "(b)"),
     ("lllccrc", r"Task & Variant & Kind & $A$ session & $D$ session & $n$ & $\rho$ [95\% CI]", b3, "(c)")], size=r"\scriptsize"))

# ---------------- S9 B3a full injection grid (b3a_injection, b3a_layers, b3a_agreement)
I9 = res("b3a_injection.csv"); L9 = res("b3a_layers.csv").rename(columns={"Unnamed: 0": "cond"}).set_index("cond"); AG = res("b3a_agreement.csv").iloc[0]
def cond_lab(c):
    k, v = c.split("_", 1) if "_" in c else (c, "")
    return {"F1": f"F1 one-eye step, {v[1:]}$^\\circ$", "F2a": f"F2a one-eye slowing, factor {v[1:]}", "F2b": f"F2b peak offset, {v[1:]}$^\\circ$",
            "F3": "F3 one-eye dropout", "N1": f"N1 mild one-eye scaling, $\\kappa$ = {v[1:]}", "C": f"contamination (F2a), {v[1:]}\\% of events"}[k]
b1 = []
for r in I9.itertuples():
    st = "judged" if r.judged else ("inconclusive" if r.inconclusive else ("specificity" if r.cond.startswith("N1") else ("descriptive" if r.cond.startswith("C_") else "unjudged")))
    val = ci(r.excess, r.excess_lo, r.excess_hi) if r.cond.startswith("N1") else ci(r.recall, r.recall_lo, r.recall_hi)
    b1.append(f"{cond_lab(r.cond)} & {st} & {ni(r.injections)} & {t3(r.admitted_rate)} & {val} & {pct(r.flip_rate, 2)} [{pct(r.flip_hi, 2)}]")
b2 = [f"{cond_lab(c)} & " + " & ".join(t3(L9.loc[c, k]) if not pd.isna(L9.loc[c, k]) else "0" for k in lay_cols) for c in I9.cond]
made.append(table(9, "Synthetic-fault injection grid (B3a). Recall: fraction of injections reaching centring abstained from by the centred rule "
    "(95\\% participant-cluster CI); reached: fraction of all injections that reached the centring step. N1 (specificity control) scales the within-saccade "
    "displacement of one randomly chosen eye by $\\kappa$, a negative-control perturbation with $|\\Delta l| \\approx .15 < \\log 2$; for N1 the column gives excess rejection, "
    "the rejection rate of the same target events reaching centring after scaling minus that in the baseline pass. "
    "Frozen criterion for judged conditions: recall $\\ge$ .80 with lower bound $>$ .70. Flips: percentage of untouched events whose keep/abstain status changed "
    "(upper 95\\% bound in brackets). Panel (b): outcome layer of each injection (fractions of all injections). Most F3 injections were removed at the "
    f"missingness layer ({ni(I9.set_index('cond').loc['F3', 'n_admitted'])} of {ni(I9.set_index('cond').loc['F3', 'injections'])} reached centring), so F3 recall was not interpreted. Agreement between centred rejection and absent per-eye support in uncontaminated data: odds ratio {AG.OR:.1f} "
    f"[{f'{AG.OR_lo:.1f}'}, {f'{AG.OR_hi:.1f}'}], $\\kappa$ = {ci(AG.kappa, AG.kappa_lo, AG.kappa_hi)}.",
    [("llrccc", r"Condition & Status & Injections & Reached & Recall or excess [95\% CI] & Flips, \%", b1, "(a)"),
     ("lcccccc", r"Condition & Retained & Abstained & Not detected & Out of class & Duration & Missingness", b2, "(b)")], size=r"\scriptsize"))

# ---------------- S10 quality covariates (b3a_quality)
Q = res("b3a_quality.csv")
QL = {"d_rms": "left/right difference in fixation RMS-S2S", "d_miss": "left/right difference in missing fraction",
      "d_sd": "left/right difference in fixation SD", "worst_rms": "worse-eye fixation RMS-S2S"}
b1 = [f"{QL[r.covariate]} & {ni(r.n)} & {ni(r.n_pid)} & {ci(r.spearman, r.lo, r.hi)}" for r in Q[Q.covariate.isin(QL)].itertuples()]
w = Q[Q.covariate == "ICC_S1S2_excl_worst_decile"].iloc[0]
made.append(table(10, "Association of per-recording rejection fraction (centred mean-speed construction) with per-eye data-quality covariates "
    f"(Spearman $\\rho$, 95\\% participant-cluster CI). RMS-S2S: root-mean-square sample-to-sample deviation within fixations, median over fixations "
    f"$\\ge$ 60 ms, per eye; worse-eye: the larger of the left- and right-eye values. Excluding participants with any recording above the 90th percentile of "
    f"worse-eye RMS-S2S over all baseline-pass recordings, "
    f"R1 S1--S2 ICC(A,1) of the centred marker (floor 5) was {ci(w.icc_a1, w.icc_a1_lo, w.icc_a1_hi)} ($n$ = {int(w.n)}).",
    [("lrrc", r"Covariate & Recordings & Participants & $\rho$ [95\% CI]", b1)]))

# ---------------- S11 per-eye gain grid (b6_gain)
G = res("b6_gain.csv")
def glab(c):
    if c.startswith("gain_asym"):
        return "1", r"%s.15 / %s.15" % (("+", "$-$") if c.endswith("_p") else ("$-$", "+"))
    q, dd = c[len("gain_q"):].split("_d"); return q, f"{float(dd):.2f}".replace("0.", ".", 1).replace("-", "$-$")
b1 = [" & ".join(glab(r.cond)) + f" & {'judged' if r.judged else 'exploratory'} & {ni(r.cond_records)} & "
      + f"{r.coverage_loss_pp:+.2f}".replace("-", "$-$") +
      f" & {t3(r.med_absdA)} & {t3(r.p95_absdA)} & {ni(r.n_pairs)} & {ci(r.d_icc, r.d_icc_lo, r.d_icc_hi)} & {'yes' if r['pass'] else 'no'}" for _, r in G.iterrows()]
made.append(table(11, "Constant per-eye gain applied to raw gaze before re-extraction: left gain $q e^{\\delta/2}$, right gain $q e^{-\\delta/2}$ "
    "(centred mean-speed construction, floor 5). Rows with $\\delta$ = $\\pm$.15 / $\\mp$.15: asymmetric condition, opposite sign in S1 / S2. Pass: median $|\\Delta A| \\le .01$, P95 $|\\Delta A| \\le .03$, $|\\Delta\\mathrm{ICC}| \\le .05$ "
    "and coverage loss $\\le 5$ percentage points. Conditions with $|\\delta| \\le .15$ were judged; $|\\delta|$ = .30 conditions were exploratory. "
    "Coverage loss: baseline-eligible minus condition-eligible recordings, in percentage points of all recordings (negative values: recordings gained). "
    "$|\\Delta A|$ is computed on recordings eligible in both baseline and condition; $\\Delta$ICC and $n$ use participants with R1 S1--S2 pairs in both.",
    [("rrlrrcccc c".replace(" ", ""), r"$q$ & $\delta$ & Status & Recordings & Coverage loss (pp) & Median $|\Delta A|$ & P95 $|\Delta A|$ & $n$ & $\Delta$ICC [95\% CI] & Pass", b1)],
    size=r"\scriptsize"))

# ---------------- S12 short-term ICCs by construction, joint and own cohorts (pw_p1_joint_icc, floor 5)
J5 = res("pw_p1_joint_icc.csv").query("floor == 5 and comparison == 'R1_S1S2'")
b1 = []
for m in CONS:
    j6 = J5.query("cohort == 'joint6' and method == @m").iloc[0]
    j4 = J5.query("cohort == 'joint4' and method == @m")
    b1.append(f"{CL[m]} & {ci(j6.icc_a1, j6.icc_a1_lo, j6.icc_a1_hi)} & {t3(j6.spearman)} & "
              f"{ci(*j4[['icc_a1', 'icc_a1_lo', 'icc_a1_hi']].iloc[0]) if len(j4) else '--'} & {t3(j6.own_icc)} ({int(j6.own_n)})")
n6 = int(J5.query("cohort == 'joint6'").n.iloc[0]); n4 = int(J5.query("cohort == 'joint4'").n.iloc[0])
made.append(table(12, "Short-term repeatability (R1 S1--S2, floor 5) by construction: tabular form of Fig.~2a. "
    "joint6: participants whose two R1 sessions are eligible under all six constructions (the cohort of Fig.~2a); joint4: a separate cohort eligible under "
    "the four mean-speed and per-eye constructions (admitted constructions excluded); own cohort: each construction's eligible participants "
    "(point estimates; own $n$ in parentheses). Spearman: rank correlation of participants' S1 and S2 values of $A$ in the joint6 cohort.",
    [("lcccc", f"Construction & joint6 ($n$ = {n6}) ICC [95\\% CI] & Spearman & joint4 ($n$ = {n4}) ICC [95\\% CI] & Own cohort", b1)]))

F4 = lambda x: f"{x:.4f}".replace("0.", ".", 1).replace("-", "$-$")
# ---------------- S13 variance components (rv_v3) and Bland-Altman coverage (rv_v2)
V3 = res("rv_v3_variance.csv").query("floor == 5 and comparison == 'R1_S1S2' and cohort == 'joint6'").set_index("method")
V3M = ["mean_u", "adm_u", "mean_c", "pereye_u"]
b1 = [f"{CL[m]} & {int(V3.loc[m].n)} & {t3(V3.loc[m].icc_a1)} & {t3(V3.loc[m].sd_s1)} & {t3(V3.loc[m].sd_s2)} & {t3(V3.loc[m].sd_between)} & "
      f"{t3(V3.loc[m].sd_diff)} & {F4(V3.loc[m].mean_diff)} & {pct(V3.loc[m].top1_share, 0)} & {pct(V3.loc[m].top3_share, 0)} & {pct(V3.loc[m].top10_share, 0)}" for m in V3M]
CV = res("rv_v2_ba_coverage.csv").query("floor == 5")
CVR = list(CV.itertuples())
f4 = lambda x: f"{x:.4f}".replace("0.", ".", 1).replace("-", "$-$")
STAT = [("Recordings", lambda r: ni(r.n_rec)), ("Participants", lambda r: str(int(r.n_pid))),
        ("Bias", lambda r: " ".join([f4(r.bias), f"[{f4(r.bias_lo)}, {f4(r.bias_hi)}]"])), ("Lower limit", lambda r: ci(r.loa_lo, r.loa_lo_lo, r.loa_lo_hi)),
        ("Upper limit", lambda r: ci(r.loa_hi, r.loa_hi_lo, r.loa_hi_hi)),
        ("Inside limits", lambda r: f"{ni(r.n_inside)} ({pct(r.coverage)}\\% [{pct(r.coverage_lo)}, {pct(r.coverage_hi)}])"),
        ("P2.5", lambda r: ci(r.p025, r.p025_lo, r.p025_hi)), ("P97.5", lambda r: ci(r.p975, r.p975_lo, r.p975_hi))]
b2 = [f"{lab} & " + " & ".join(f(r) for r in CVR) for lab, f in STAT]
b2head = "Statistic & " + " & ".join(f"{CL[r.y]} vs {CL[r.x]}" for r in CVR)
made.append(table(13, "Session-difference structure in the joint cohort and empirical coverage of the Bland--Altman limits. (a) R1 S1--S2 (floor 5, "
    "joint cohort): standard deviations of $A$ in each session, between-participant SD (SD of the participants' two-session means), SD and "
    "mean of the session differences S2 $-$ S1, and the percentage of the summed squared deviations of the session differences from their mean "
    "contributed by the largest one, three and ten participants. (b) Agreement of recording-level $A$ (Fig.~4; floor 5, all rounds): bias, limits of agreement (bias $\\pm$ 1.96 SD), recordings inside "
    "the limits, and empirical 2.5th and 97.5th percentiles of the differences. All quantities in (b) are pooled over recordings, so participants with more "
    "eligible recordings carry more weight; brackets are 95\\% participant-cluster bootstrap intervals (5{,}000 draws, participants resampled with all of "
    "their recordings).",
    [("lcccccccccc", r"Construction & $n$ & ICC & SD S1 & SD S2 & SD of participant means & SD diff. & Mean diff. & Top 1 (\%) & Top 3 (\%) & Top 10 (\%)", b1, "(a) Joint-cohort session differences"),
     ("l" + "c" * len(CVR), b2head, b2, "(b) Bland--Altman coverage")], size=r"\scriptsize"))

# ---------------- S14 post hoc peak-velocity ceiling comparator (rv_v7, experiments/r3_revision.py)
CE = res("rv_v7_ceiling_icc.csv"); CD = res("rv_v7_ceiling_delta.csv"); CEV = res("rv_v7_ceiling_events.csv")
TAGS = [(fl, c) for fl in (5, 10) for c in (750, 1000, 1500)]
def cm(m, c): return m.replace("ceiladm", f"ceiladm{c}").replace("ceil_", f"ceil{c}_") if "ceil" in m else m
ICM = [("mean_u", "Mean-speed, unfiltered"), ("adm_u", "Admitted, unfiltered"), ("ceil_u", "Ceiling, unfiltered"),
       ("ceil_c", "Ceiling, centred"), ("ceiladm_u", "Ceiling + admission, unfiltered"), ("ceiladm_c", "Ceiling + admission, centred")]
def ce(fl, c, m): return row(CE, floor=fl, ceiling=c, method=cm(m, c))
b1 = [f"{lab} & " + " & ".join(ci(ce(fl, c, m).icc_a1, ce(fl, c, m).icc_a1_lo, ce(fl, c, m).icc_a1_hi) for fl, c in TAGS) for m, lab in ICM]
b1 = [r"Joint $n$ & " + " & ".join(str(int(ce(fl, c, "mean_u").n)) for fl, c in TAGS)] + b1
DCM = [("ceil_u-mean_u", "Ceiling $-$ mean-speed (u)"), ("adm_u-ceil_u", "Admitted $-$ ceiling (u)"),
       ("ceiladm_u-ceil_u", "Ceiling + admission $-$ ceiling (u)"), ("adm_u-mean_u", "Admitted $-$ mean-speed (u)")]
def cdl(fl, c, k, sub):
    a, b = k.split("-"); r = row(CD, tag=f"floor{fl}_ceil{c}", subset=sub, contrast=f"{cm(a, c)}-{cm(b, c)}"); return ci(r.delta, r.lo, r.hi)
b2 = [f"{lab}, {'all' if sub == 'all' else 'one omitted'} & " + " & ".join(cdl(fl, c, k, sub) for fl, c in TAGS) for sub in ("all", "without_top") for k, lab in DCM]
hd = "& " + " & ".join(f"F{fl}, {c:,}".replace(",", "{,}") for fl, c in TAGS)
ev = lambda c, cls: row(CEV, ceiling=c, cls=cls)
b3 = []
for c in (750, 1000, 1500):
    for cls, lab in (("rs", "Return sweeps"), ("prog", "Progressive"), ("all", "All")):
        r = ev(c, cls)
        b3.append(f"{c:,}".replace(",", "{,}") + f" & {lab} & {ni(r.n_events)} & {ni(r.n_over)} ({pct(r.frac_over, 2)}) & {pct(r.frac_over_without_support)} & "
                  f"{pct(r.frac_no_support_over)} & {ni(r.n_files_with_over)} & {ni(r.n_over_L)} / {ni(r.n_over_R)}")
p1, p2 = ev(1000, "pid198_S1_rs"), ev(1000, "pid198_S2_rs")
made.append(table(14, "Post hoc peak-velocity ceiling comparator (not in the frozen plan). Ceiling constructions drop classified mean-speed events whose larger "
    "per-eye peak velocity exceeds the ceiling (750, 1{,}000 or 1{,}500\\textdegree/s); `ceiling + admission' applies the ceiling and then the admission rule; centred "
    "variants recompute the centring rule within the retained set. (a) R1 S1--S2 ICC(A,1) [95\\% participant-bootstrap CI, 5{,}000 draws] in the joint cohort "
    "eligible under the six principal and four ceiling constructions (F, event floor). (b) Paired ICC differences in the same cohorts, with all participants and with the participant contributing "
    "most to the unfiltered mean-speed session differences omitted. (c) Classified mean-speed events (all rounds) above each ceiling: number (percentage of the class), percentage of "
    "over-ceiling events without a matched left/right per-eye pair, percentage of events without such a pair that exceed the ceiling, recordings with at least one over-ceiling event, and "
    f"over-ceiling events in which the left / right eye exceeded the ceiling (counted per eye; an event with both eyes above the ceiling counts in both). At 1{{,}}000\\textdegree/s, {int(p1.n_over)} of the {int(p1.n_events)} session-1 return sweeps of the most influential participant "
    f"exceeded the ceiling (median $l$ {t3(p1.med_l_before)} before and {t3(p1.med_l_after)} after), against {int(p2.n_over)} of {int(p2.n_events)} in session 2.",
    [("l" + "c" * len(TAGS), "Construction " + hd, b1, "(a) ICC(A,1)"), ("l" + "c" * len(TAGS), "Contrast " + hd, b2, "(b) $\\Delta$ICC"),
     ("llcccccc", r"Ceiling (\textdegree/s) & Class & Events & Over ceiling (\%) & Without pair (\%) & No-pair events over (\%) & Recordings & Left / right", b3, "(c) Events above the ceiling")],
    size=r"\scriptsize"))

# ---------------- S15 typical magnitude of session differences |S2 - S1| of A (rv_v8)
AD = res("rv_v8_absdiff.csv"); ADD = res("rv_v8_absdiff_delta.csv")
AM = [("mean_u", "Mean-speed, unfiltered"), ("adm_u", "Admitted, unfiltered"), ("ceil1000_u", "Ceiling 1{,}000, unfiltered"),
      ("ceiladm1000_u", "Ceiling + admission, unfiltered"), ("mean_c", "Mean-speed, centred"), ("adm_c", "Admitted, centred"),
      ("pereye_u", "Per-eye, unfiltered"), ("pereye_c", "Per-eye, centred")]
def a4(r, k): return f"{F4(r[k])} [{F4(r[k + '_lo'])}, {F4(r[k + '_hi'])}]"
b1 = [f"{lab} & " + " & ".join(a4(row(AD, subset=sub, method=m), k) for sub in ("all", "without_top") for k in ("median", "p90", "mean")) for m, lab in AM]
DM = [("adm_u-mean_u", "Admitted $-$ mean-speed (u)"), ("ceil1000_u-mean_u", "Ceiling $-$ mean-speed (u)"), ("adm_u-ceil1000_u", "Admitted $-$ ceiling (u)"),
      ("ceiladm1000_u-ceil1000_u", "Ceiling + admission $-$ ceiling (u)"), ("mean_c-mean_u", "Centred $-$ unfiltered (mean-speed)"),
      ("adm_c-adm_u", "Centred $-$ unfiltered (admitted)"), ("pereye_u-mean_u", "Per-eye $-$ mean-speed (u)")]
b2 = [f"{lab} & " + " & ".join(a4(row(ADD, subset=sub, contrast=k), q) for sub in ("all", "without_top") for q in ("median", "p90", "mean")) for k, lab in DM]
n0, n1 = int(row(AD, subset="all", method="mean_u").n), int(row(AD, subset="without_top", method="mean_u").n)
hd = f"& Median & P90 & Mean & Median & P90 & Mean"
made.append(table(15, "Typical magnitude of the session difference of $A$ (post hoc). Absolute R1 session difference $|A_{S2} - A_{S1}|$ in the floor-5 joint cohort "
    f"($n$ = {n0}) and with the most influential participant omitted ($n$ = {n1}); brackets, 95\\% participant-bootstrap intervals (5{{,}}000 draws; "
    "differences in (b) use paired draws). $A$ is half a difference of log ratios, so $|\\Delta A| = 0.02$ corresponds to a change of about 4\\% in the ratio of the return-sweep and "
    "progressive median left/right peak-velocity ratios ($e^{2 \\times 0.02} \\approx 1.04$).",
    [("l" + "c" * 6, f"Construction & \\multicolumn{{3}}{{c}}{{All ($n$ = {n0})}} & \\multicolumn{{3}}{{c}}{{One omitted ($n$ = {n1})}} \\\\ " + hd, b1, "(a) $|A_{S2} - A_{S1}|$"),
     ("l" + "c" * 6, f"Difference & \\multicolumn{{3}}{{c}}{{All}} & \\multicolumn{{3}}{{c}}{{One omitted}} \\\\ " + hd, b2, "(b) Paired differences")],
    size=r"\scriptsize"))

# ---------------- S16 detector-setting sensitivity (rv_v9)
SV = res("rv_v9_sens_summary.csv"); RP = res("rv_v9_reproduction.csv").set_index("detector")
VL = {"ref": "Frozen (40, 7, 3)", "vth30": "Threshold 30", "vth50": "Threshold 50", "sg5": "SG window 5", "sg9": "SG window 9", "gap1": "Pairing gap 1", "gap5": "Pairing gap 5"}
b1 = [f"{VL[r.variant]} & {int(r.n_joint)} & {ni(r.n_events_mean)} & {ni(r.n_events_pereye)} & {pct(r.adm_retention)} & {t3(r.icc_mean_u)} & {t3(r.icc_adm_u)} & {t3(r.icc_ceil_u)} & "
      f"{ci(r.delta, r.delta_lo, r.delta_hi)} & {ci(r.delta_wo_top, r.delta_wo_top_lo, r.delta_wo_top_hi)} ({int(r.top_pid)}) & {ci(r.d_adm_ceil, r.d_adm_ceil_lo, r.d_adm_ceil_hi)} & "
      f"{t3(r.d_meanc_meanu)} & {t3(r.d_admc_admu)}" for r in SV.itertuples()]
SI = res("rv_v9_sens_icc.csv")
b2 = [f"{VL[v]} & " + " & ".join(t3(row(SI, variant=v, method=m).icc_a1) for m in CONS) for v in SV.variant]
ok = all(bool(RP.loc[d].same_windows) and RP.loc[d].max_abs_diff_pv_amp == 0 for d in ("mean", "pereye"))
assert ok
made.append(table(16, "Post hoc detector-setting sensitivity. All R1 recordings were re-extracted with one detector setting changed at a time: velocity threshold "
    "(\\textdegree/s), Savitzky--Golay smoothing window (samples) or the maximum left/right onset gap for pairing per-eye events (samples); the pairing gap affects only "
    "per-eye events and therefore admission. Re-extraction with the frozen settings reproduced the analysed event sets exactly "
    f"({ni(RP.loc['mean'].n_ref)} mean-speed and {ni(RP.loc['pereye'].n_ref)} per-eye events; identical windows; maximum absolute difference in peak velocity or amplitude 0). "
    "(a) Per variant, in its own floor-5 joint cohort: classified mean-speed and per-eye events in R1, admission retention (percentage of classified mean-speed events kept), "
    "R1 S1--S2 ICC(A,1) of unfiltered mean-speed, admitted and 1{,}000\\textdegree/s ceiling events (the ceiling ICC and admitted $-$ ceiling in the joint cohort that also requires ceiling eligibility), the admission contrast with all participants and with the most influential participant omitted "
    "(identifier in parentheses), admitted minus ceiling, and the centring contrasts. Brackets, 95\\% paired participant-bootstrap intervals (5{,}000 draws). (b) ICC(A,1) of the six constructions.",
    [("lcccccccccccc", r"Variant & $n$ & Mean-speed events & Per-eye events & Retention (\%) & Mean-speed u & Admitted u & Ceiling u & Adm. $-$ mean (u) & Adm. $-$ mean, one omitted & Adm. $-$ ceiling (u) & Centring, mean-speed & Centring, admitted", b1, "(a) Summary"),
     ("l" + "c" * len(CONS), "Variant & " + " & ".join(CL[c] for c in CONS), b2, "(b) ICC(A,1) by construction")],
    size=r"\scriptsize"))

Path(__file__).with_name("SI_TABLES.tex").write_text("\n".join(rf"\input{{figures/si/TABLE_S{k}.tex}}" for k in made) + "\n")
print("wrote", [f"TABLE_S{k}" for k in made])
