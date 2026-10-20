"""B1 main reliability + B2 selection vs exclusion (and B4 when run with cond=pereye).

usage: python experiments/analyze_b1_b2.py [base|pereye]   -> results/b1_reliability[_pereye].csv, results/b2_*.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import common as C, stats as S

cond = sys.argv[1] if len(sys.argv) > 1 else "base"
tag = "" if cond == "base" else f"_{cond}"
os.makedirs("results", exist_ok=True)
ev, cl = S.load_pass(cond)
FLOORS = (1, 5, 10)
M = {(filt, fl): C.marker(cl, filt=filt, rs_floor=fl) for filt in (False, True) for fl in FLOORS}

# ---- denominators
den = dict(cond=cond, recordings=ev.file.nunique(), events_detected=len(ev), events_pool=len(cl),
           events_kept=int(cl.keep.sum()), retention=cl.keep.mean(),
           rs_pool=int((cl.cls == "rs").sum()), prog_pool=int((cl.cls == "prog").sum()))
for (filt, fl), r in M.items():
    k = f"{'centred' if filt else 'unfilt'}_f{fl}"
    den[f"eligible_rec_{k}"] = len(r); den[f"pairs_R1_{k}"] = len(C.s1s2(r))
pd.DataFrame([den]).to_csv(f"results/b1_denominators{tag}.csv", index=False)

# ---- B1: own-sample and intersection-pair rows
rows = []
comps = {key: S.comparisons(r) for key, r in M.items()}
for fl in FLOORS:
    for name in comps[(True, fl)]:
        cu, cc = comps[(False, fl)].get(name), comps[(True, fl)].get(name)
        if cu is None or cc is None:
            continue
        inter = cu.index.intersection(cc.index)
        idx = C.boot_idx(len(inter)) if len(inter) else None
        for filt, w in ((False, cu), (True, cc)):
            rows.append(dict(cond=cond, comparison=name, floor=fl, filt="centred" if filt else "unfiltered", sample="own",
                             **S.pair_stats(w.o1.values, w.o2.values)))
            wi = w.loc[inter]
            rows.append(dict(cond=cond, comparison=name, floor=fl, filt="centred" if filt else "unfiltered",
                             sample="intersection", **S.pair_stats(wi.o1.values, wi.o2.values, idx)))
        if len(inter) >= 5:
            a, b = cc.loc[inter], cu.loc[inter]
            rows.append(dict(cond=cond, comparison=name, floor=fl, filt="centred-unfiltered", sample="intersection",
                             **{f"d_{k}": v for k, v in S.paired_delta(a.o1.values, a.o2.values, b.o1.values, b.o2.values, idx).items()}))
B1 = pd.DataFrame(rows)
B1.to_csv(f"results/b1_reliability{tag}.csv", index=False)

# returners vs non-returners baseline level (R1, centred, floor 5)
r = M[(True, 5)]; w = C.s1s2(r); ret = set(r[r["round"] > 1].pid)
lvl = pd.DataFrame([dict(group=g, n=len(x), mean_A=x.mean(1).mean(), sd_A=x.mean(1).std())
                    for g, x in (("returners", w[w.index.isin(ret)]), ("non_returners", w[~w.index.isin(ret)]))])
lvl.to_csv(f"results/b1_returner_baseline{tag}.csv", index=False)

if cond != "base":
    print(B1[(B1.comparison == "R1_S1S2") & (B1.floor == 5)][["filt", "sample", "n", "icc_a1", "icc_a1_lo", "icc_a1_hi", "d_delta", "d_lo", "d_hi"]].round(3).to_string())
    sys.exit()

# ---- B2: counts before/after per participant-session-class; flagged participants
cnt = cl.groupby(["pid", "round", "session", "cls"]).keep.agg(before="size", after="sum").reset_index()
cnt.to_csv("results/b2_counts.csv", index=False)
flag = cl[cl.pid.isin([101, 245, 319]) & (cl["round"] == 1)].groupby(["pid", "session", "cls"]).agg(
    n=("l", "size"), kept=("keep", "sum"), med_l=("l", "median"), med_l_kept=("l", lambda x: np.nan))
for (p, s, c), g in cl[cl.pid.isin([101, 245, 319]) & (cl["round"] == 1)].groupby(["pid", "session", "cls"]):
    flag.loc[(p, s, c), "med_l_kept"] = g.l[g.keep].median()
flag.reset_index().to_csv("results/b2_flagged.csv", index=False)

# ---- B2: count-matched random baseline (R1 S1-S2, floors as B1), 200 reps
R1 = cl[cl["round"] == 1]
groups = {k: (g.l.values, int(g.keep.sum())) for k, g in R1.groupby(["pid", "session", "cls"])}
NREP = 200
reps = []
for rep in range(NREP):
    rng = np.random.default_rng(C.SEED + rep)
    med = {k: (np.median(v[rng.permutation(len(v))[:m]]) if m > 0 else np.nan) for k, (v, m) in groups.items()}
    s = pd.Series(med); s.index.names = ["pid", "session", "cls"]
    reps.append(s.unstack("cls"))
rows = []
for fl in FLOORS:
    cc = C.s1s2(M[(True, fl)]); cu = C.s1s2(M[(False, fl)])
    inter = cc.index.intersection(cu.index)          # random keeps centred counts -> centred eligibility
    rnd = np.stack([((0.5 * (x.rs - x.prog)).unstack("session").reindex(inter)[[1, 2]].values) for x in reps])  # (NREP, n, 2)
    idx = C.boot_idx(len(inter))
    a, b = cc.loc[inter], cu.loc[inter]
    rand_pt = C.icc(rnd[:, :, 0], rnd[:, :, 1])[0].mean()
    bc = C.icc(a[1].values[idx], a[2].values[idx])[0]
    bu = C.icc(b[1].values[idx], b[2].values[idx])[0]
    br = np.zeros(len(idx))
    for j in range(0, len(idx), 250):                 # chunked: (NREP, 250, n)
        ii = idx[j:j + 250]
        br[j:j + 250] = C.icc(rnd[:, ii, 0], rnd[:, ii, 1])[0].mean(0)
    c_pt = C.icc(a[1].values, a[2].values)[0]; u_pt = C.icc(b[1].values, b[2].values)[0]
    rows.append(dict(floor=fl, n=len(inter), centred=c_pt, unfiltered=u_pt, random_mean=rand_pt,
                     random_rep_lo=np.percentile(C.icc(rnd[:, :, 0], rnd[:, :, 1])[0], 2.5),
                     random_rep_hi=np.percentile(C.icc(rnd[:, :, 0], rnd[:, :, 1])[0], 97.5),
                     d_c_u=c_pt - u_pt, d_c_u_lo=C.pct(bc - bu)[0], d_c_u_hi=C.pct(bc - bu)[1],
                     d_c_r=c_pt - rand_pt, d_c_r_lo=C.pct(bc - br)[0], d_c_r_hi=C.pct(bc - br)[1]))
B2 = pd.DataFrame(rows)
B2.to_csv("results/b2_selection.csv", index=False)

pd.set_option("display.width", 250)
print(pd.DataFrame([den]).T.to_string())
print(B1[(B1.floor == 5)][["comparison", "filt", "sample", "n", "icc_a1", "icc_a1_lo", "icc_a1_hi", "spearman", "d_delta", "d_lo", "d_hi"]].round(3).to_string())
print(B2.round(3).to_string())
