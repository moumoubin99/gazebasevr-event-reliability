"""Revision analyses (paper review round 1; CPU only, local). Reads figures/source_data/recording_A.csv (per-recording A
for the six constructions, written by figures/make_source_data.py from the frozen passes); no new extraction.

V1  own-cohort R1-R3 ICC(A,1) with participant bootstrap CI for every construction (Table 1 own-cohort column)
V2  Bland-Altman empirical coverage: share of recordings inside bias +/- 1.96 SD, 2.5th/97.5th percentiles of the
    differences, participant-cluster bootstrap CIs for bias, limits, coverage and percentiles
V3  joint-cohort (R1 S1-S2, floor 5) variance summary, mean_u vs adm_u: ICC components, session means, between-person SD,
    SD of session differences, share of sum of squared session differences from the top contributors

usage: python experiments/r1_revision.py -> results/rv_v1_own_r1r3.csv, results/rv_v2_ba_coverage.csv,
       results/rv_v3_variance.csv, results/rv_v3_top.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import common as C, stats as S

FL = 5
R = pd.read_csv("figures/source_data/recording_A.csv").query("floor == @FL")
MK = {m: g.drop(columns=["floor", "construction"]) for m, g in R.groupby("construction")}
CONS = ["mean_u", "mean_c", "adm_u", "adm_c", "pereye_u", "pereye_c"]

# ---- V1
rows = []
for m in CONS:
    w = C.rounds(MK[m], 3)
    ps = S.pair_stats(w[1].values, w[3].values)
    rows.append(dict(floor=FL, comparison="R1_R3", method=m, sample="own", **{k: ps[k] for k in ("n", "icc_a1", "icc_a1_lo", "icc_a1_hi")}))
v1 = pd.DataFrame(rows); print(v1.round(4).to_string()); v1.to_csv("results/rv_v1_own_r1r3.csv", index=False)

# ---- V2
rows = []
rng_seed = C.SEED
for x, y in (("mean_c", "pereye_c"), ("mean_u", "pereye_u")):
    a = MK[x].set_index(["pid", "round", "session"]).A; b = MK[y].set_index(["pid", "round", "session"]).A
    j = a.index.intersection(b.index); d = (b.loc[j] - a.loc[j]).values
    pid = j.get_level_values("pid").values; up = np.unique(pid); rows_of = [np.where(pid == p)[0] for p in up]

    def summ(dd):
        mu, sd = dd.mean(), dd.std(ddof=1); lo, hi = mu - 1.96 * sd, mu + 1.96 * sd
        return np.array([mu, lo, hi, np.mean((dd >= lo) & (dd <= hi)), np.percentile(dd, 2.5), np.percentile(dd, 97.5)])
    pt = summ(d)
    bi = C.boot_idx(len(up), S.B, rng_seed)
    bs = np.array([summ(d[np.concatenate([rows_of[i] for i in k])]) for k in bi])
    names = ["bias", "loa_lo", "loa_hi", "coverage", "p025", "p975"]
    out = dict(floor=FL, x=x, y=y, n_rec=len(d), n_pid=len(up), n_inside=int(round(pt[3] * len(d))))
    for i, nm in enumerate(names):
        out[nm] = pt[i]; out[nm + "_lo"], out[nm + "_hi"] = C.pct(bs[:, i])
    rows.append(out)
v2 = pd.DataFrame(rows); print(v2.round(4).T.to_string()); v2.to_csv("results/rv_v2_ba_coverage.csv", index=False)

# ---- V3
P = {m: C.s1s2(MK[m]) for m in CONS}
ids = None
for m in CONS:
    ids = P[m].index if ids is None else ids.intersection(P[m].index)
rows, top = [], []
for m in ("mean_u", "adm_u", "mean_c", "pereye_u"):
    o1, o2 = P[m].loc[ids, 1].values, P[m].loc[ids, 2].values
    n = len(ids); rm = (o1 + o2) / 2; d = o2 - o1; gm = rm.mean()
    msr = 2 * ((rm - gm) ** 2).sum() / (n - 1)
    sse = ((o1 - rm - o1.mean() + gm) ** 2 + (o2 - rm - o2.mean() + gm) ** 2).sum(); mse = sse / (n - 1)
    msc = n * ((o1.mean() - gm) ** 2 + (o2.mean() - gm) ** 2)
    sq = (d - d.mean()) ** 2; order = np.argsort(-sq)
    rows.append(dict(floor=FL, comparison="R1_S1S2", cohort="joint6", method=m, n=n, icc_a1=C.icc(o1, o2)[0],
                     mean_s1=o1.mean(), mean_s2=o2.mean(), sd_s1=o1.std(ddof=1), sd_s2=o2.std(ddof=1),
                     sd_between=rm.std(ddof=1), sd_diff=d.std(ddof=1), mean_diff=d.mean(),
                     var_subject=(msr - mse) / 2, var_error=mse, var_session=max((msc - mse) / n, 0),
                     msr=msr, mse=mse, msc=msc,
                     top1_share=sq[order[:1]].sum() / sq.sum(), top3_share=sq[order[:3]].sum() / sq.sum(),
                     top5_share=sq[order[:5]].sum() / sq.sum(), top10_share=sq[order[:10]].sum() / sq.sum()))
    for k, i in enumerate(order[:5]):
        top.append(dict(method=m, rank=k + 1, pid=int(ids[i]), s1=o1[i], s2=o2[i], diff=d[i], share=sq[i] / sq.sum()))
v3 = pd.DataFrame(rows); print(v3.round(4).T.to_string()); v3.to_csv("results/rv_v3_variance.csv", index=False)
t = pd.DataFrame(top); print(t.round(4).to_string()); t.to_csv("results/rv_v3_top.csv", index=False)
