"""Revision analyses (paper review round 2; CPU only, local). Reads figures/source_data/recording_A.csv; no new extraction.

V4  paired leave-one-participant-out sensitivity of the joint-cohort (R1 S1-S2, floor 5) admission contrast
    mean_u -> adm_u: each omission is applied to both constructions; reports the full-cohort contrast, the contrast
    after omitting the participant with the largest centred squared unfiltered session difference, and the
    leave-one-out range; paired participant bootstrap CI (5000 draws) with and without that participant.

V5  joint-cohort distribution of participant A (mean of the two R1 sessions) per construction.
V6  R1-R3 round-average cohorts: participants in each construction's own cohort but outside the joint cohort.

usage: python experiments/r2_revision.py -> results/rv_v4_loo_admission.csv, results/rv_v4_loo_summary.csv,
       results/rv_v5_A_distribution.csv, results/rv_v6_r1r3_cohorts.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import common as C, stats as S

FL = 5
R = pd.read_csv("figures/source_data/recording_A.csv").query("floor == @FL")
MK = {m: g.drop(columns=["floor", "construction"]) for m, g in R.groupby("construction")}
CONS = ["mean_u", "mean_c", "adm_u", "adm_c", "pereye_u", "pereye_c"]
P = {m: C.s1s2(MK[m]) for m in CONS}
ids = None
for m in CONS:
    ids = P[m].index if ids is None else ids.intersection(P[m].index)
x1, x2 = P["mean_u"].loc[ids, 1].values, P["mean_u"].loc[ids, 2].values
y1, y2 = P["adm_u"].loc[ids, 1].values, P["adm_u"].loc[ids, 2].values
n = len(ids)
full = C.icc(y1, y2)[0] - C.icc(x1, x2)[0]
d = x2 - x1; sq = (d - d.mean()) ** 2; top = int(np.argmax(sq))
rows = []
for i in range(n):
    k = np.arange(n) != i
    ix, iy = C.icc(x1[k], x2[k])[0], C.icc(y1[k], y2[k])[0]
    rows.append(dict(pid=int(ids[i]), icc_mean_u=ix, icc_adm_u=iy, delta=iy - ix, sq_share=sq[i] / sq.sum()))
L = pd.DataFrame(rows); L.to_csv("results/rv_v4_loo_admission.csv", index=False)
t = L.iloc[top]
s = pd.DataFrame([dict(floor=FL, n=n, icc_mean_u=C.icc(x1, x2)[0], icc_adm_u=C.icc(y1, y2)[0], delta_full=full,
                       top_pid=int(t.pid), top_share=t.sq_share, icc_mean_u_wo_top=t.icc_mean_u, icc_adm_u_wo_top=t.icc_adm_u,
                       delta_wo_top=t.delta, delta_min=L.delta.min(), delta_max=L.delta.max(),
                       pid_min=int(L.pid[L.delta.idxmin()]), pid_max=int(L.pid[L.delta.idxmax()]),
                       n_delta_below_0_1=int((L.delta < .1).sum()))])
k = np.arange(n) != top
bt = S.paired_delta(y1[k], y2[k], x1[k], x2[k], C.boot_idx(n - 1, S.B, C.SEED))
s["delta_wo_top_lo"], s["delta_wo_top_hi"] = bt["lo"], bt["hi"]
bf = S.paired_delta(y1, y2, x1, x2, C.boot_idx(n, S.B, C.SEED))
s["delta_full_lo"], s["delta_full_hi"] = bf["lo"], bf["hi"]
s.to_csv("results/rv_v4_loo_summary.csv", index=False)
print(s.round(4).T.to_string()); s.to_csv("results/rv_v4_loo_summary.csv", index=False)
print(L.sort_values("delta").head(5).round(4).to_string())

# ---- V5
rows = []
for m in CONS:
    a = P[m].loc[ids].mean(axis=1)
    rows.append(dict(floor=FL, construction=m, n=len(a), mean=a.mean(), median=a.median(), sd=a.std(ddof=1),
                     p05=a.quantile(.05), p95=a.quantile(.95), share_positive=(a > 0).mean()))
v5 = pd.DataFrame(rows); print(v5.round(4).to_string()); v5.to_csv("results/rv_v5_A_distribution.csv", index=False)

# ---- V6
W = {m: C.rounds(MK[m], 3) for m in CONS}
jj = None
for m in CONS:
    jj = W[m].index if jj is None else jj.intersection(W[m].index)
v6 = pd.DataFrame([dict(floor=FL, construction=m, n_own=len(W[m]), n_joint=len(jj),
                        extra_pids=" ".join(str(int(x)) for x in sorted(set(W[m].index) - set(jj)))) for m in CONS])
print(v6.to_string()); v6.to_csv("results/rv_v6_r1r3_cohorts.csv", index=False)
