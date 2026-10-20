"""Table 1: longitudinal (R1-R3) repeatability of A for the six constructions, floor 5.
Sources: results/pw_p1_joint_icc.csv (joint6 cohorts, own_n/own_icc), results/b1_reliability(_pereye).csv
(R1_R2 centred; cross-check of own cohorts), results/rv_v1_own_r1r3.csv (own-cohort R1-R3 intervals for all six
constructions, recomputed by experiments/r1_revision.py), results/pw_p1_joint_delta.csv (R1_R3 pereye_c-mean_c)."""
from paper_plot_style import res, f3, CONS, cons_label
from pathlib import Path

FL = 5
J = res("pw_p1_joint_icc.csv").query("floor == @FL and cohort == 'joint6'").set_index(["comparison", "method"])
D = res("pw_p1_joint_delta.csv").query("floor == @FL and cohort == 'joint6'").set_index(["comparison", "contrast"])
V1 = res("rv_v1_own_r1r3.csv").query("floor == @FL and comparison == 'R1_R3' and sample == 'own'").set_index("method")
B = {"base": res("b1_reliability.csv"), "pereye": res("b1_reliability_pereye.csv")}
OWN_SRC = {"mean_u": ("base", "unfiltered"), "mean_c": ("base", "centred"),
           "pereye_u": ("pereye", "unfiltered"), "pereye_c": ("pereye", "centred")}

def t3(x): return f3(x).replace("-", "$-$")
def ci(v, lo, hi): return f"{t3(v)} [{t3(lo)}, {t3(hi)}]"
def b1(cond, comp, filt, sample):
    q = B[cond].query("comparison == @comp and floor == @FL and filt == @filt and sample == @sample")
    assert len(q) == 1; return q.iloc[0]

ns = {c: sorted(J.xs(c).n.unique()) for c in ("R1_R3", "R1S1_R3S1", "R1_S1S2")}
assert ns == {"R1_R3": [53], "R1S1_R3S1": [56], "R1_S1S2": [376]}, ns
rows = []
for m in CONS:
    r3, s1, st = J.loc[("R1_R3", m)], J.loc[("R1S1_R3S1", m)], J.loc[("R1_S1S2", m)]
    own_n = int(r3.own_n); v = V1.loc[m]
    assert int(v.n) == own_n and abs(v.icc_a1 - r3.own_icc) < 1e-9, m
    if m in OWN_SRC:
        o = b1(OWN_SRC[m][0], "R1_R3", OWN_SRC[m][1], "own")
        assert int(o.n) == own_n and max(abs(o.icc_a1 - v.icc_a1), abs(o.icc_a1_lo - v.icc_a1_lo), abs(o.icc_a1_hi - v.icc_a1_hi)) < 1e-9, m
    own = ci(v.icc_a1, v.icc_a1_lo, v.icc_a1_hi)
    rows.append(f"{cons_label(m)} & {ci(r3.icc_a1, r3.icc_a1_lo, r3.icc_a1_hi)} & {own} ({own_n}) & "
                f"{ci(s1.icc_a1, s1.icc_a1_lo, s1.icc_a1_hi)} & {ci(st.icc_a1, st.icc_a1_lo, st.icc_a1_hi)} \\\\")

d = D.loc[("R1_R3", "pereye_c-mean_c")]
r2 = b1("base", "R1_R2", "centred", "intersection")
assert int(d.n) == 53 and int(r2.n) == 29
foot = [
    f"R1--R3 contrast, per-eye centred $-$ mean-speed centred (joint cohort, $n$ = {int(d.n)}): "
    f"$\\Delta$ICC = {ci(d.delta, d.lo, d.hi)}.",
    f"R1--R2 round averages, mean-speed centred, participants eligible both unfiltered and centred ($n$ = {int(r2.n)}): ICC = {ci(r2.icc_a1, r2.icc_a1_lo, r2.icc_a1_hi)}.",
]
tex = r"""\begin{table*}[t]
\centering
\caption{Longitudinal repeatability of $A$ for the six event constructions (eligibility floor 5). Entries are ICC(A,1) with 95\%
participant-bootstrap intervals (5{,}000 resamples). R1--R3: per-participant mean of the two R1 session values versus the mean of the two R3 session values
(unit: round average; all four sessions eligible), in the joint cohort eligible under all six constructions and in each construction's own cohort (own $n$ in parentheses).
R1S1--R3S1: single sessions (unit: one session). R1 S1--S2: the short-term comparison of Fig.~2a (unit: one session), for reference.
The two sessions of a round were about 30 min apart; R3 was recorded about 19--26 months after R1 (dataset description; Methods). Differences between
columns combine changes in interval, aggregation and cohort, so they do not identify an effect of elapsed time.}
\label{tab:longitudinal}
\small
\begin{tabular}{lcccc}
\toprule
 & \multicolumn{2}{c}{R1--R3, round averages} & R1S1--R3S1 & R1 S1--S2 \\
\cmidrule(lr){2-3}
Construction & joint, $n$ = 53 & own cohort ($n$) & joint, $n$ = 56 & joint, $n$ = 376 \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}

\smallskip
\begin{minipage}{0.95\textwidth}\footnotesize
""" + " ".join(foot) + r"""
\end{minipage}
\end{table*}
"""
Path(__file__).with_name("TABLE_1.tex").write_text(tex)
print(tex)
