"""B3a automated validation: synthetic one-eye failure injection + independent per-eye agreement + quality covariates.

usage: python experiments/analyze_b3a.py   (needs passes base, pereye, F1_*, F2a_*, F2b_*, F3, N1_*, C_p*, REC)
"""
import glob, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(__file__))
import common as C, stats as S

os.makedirs("results", exist_ok=True)
GRID = ["F1_a2", "F1_a4", "F1_a8", "F2a_f1.5", "F2a_f2", "F2a_f3", "F2b_b3", "F2b_b6", "F3", "N1_k0.85", "N1_k1.15"]
JUDGED = {"F1_a4", "F1_a8", "F2a_f3", "F2b_b6"}
CONTAM = ["C_p5", "C_p10", "C_p20", "C_p30"]

base_ev, base = S.load_pass("base")
base_key = base.set_index(["file", "s"])


def status(ev, cl, inj):
    """For each injection: best-overlapping re-detected event and the layer that handled it."""
    clk = cl.set_index(["file", "s"])
    out = []
    by = {f: g for f, g in ev.groupby("file")}
    for _, j in inj.iterrows():
        g = by.get(j.file)
        st = dict(file=j.file, pid=j.pid, session=j.session, kind=j.kind, eye=j.eye, param=j.param,
                  target_cls=j.get("target_cls"), detected=False, layer="not_detected", cls=None, keep=None, l=np.nan)
        if g is not None:
            ov = np.minimum(g.e1.values, j.w1) - np.maximum(g.s.values, j.w0)
            if (ov >= 0).any():
                e = g.iloc[int(np.argmax(ov))]
                st["detected"] = True
                if e.nan_raw:
                    st["layer"] = "missingness"
                elif e.dur_ms > C.MAXDUR:
                    st["layer"] = "duration"
                elif (j.file, e.s) not in clk.index:
                    st["layer"] = "out_of_class"
                else:
                    c = clk.loc[(j.file, e.s)]
                    st.update(cls=c.cls, keep=bool(c.keep), l=c.l, layer="retained" if c.keep else "abstained")
        out.append(st)
    return pd.DataFrame(out)


def untouched_flips(ev, cl, inj):
    """Decision flips of events >= 25 samples from every injection window AND every modified sample range
    (SPEC amendment A3), matched to baseline by onset +-1."""
    win = {f: np.c_[np.minimum(g.w0, g.m0), np.maximum(g.w1, g.m1 - 1)] for f, g in inj.groupby("file")}
    cl = cl[cl.file.isin(win)]
    keepf = []
    for f, g in cl.groupby("file"):
        w = win[f]
        far = np.all((g.s.values[:, None] > w[:, 1] + 25) | (g.e1.values[:, None] < w[:, 0] - 25), axis=1)
        keepf.append(g[far])
    g = pd.concat(keepf)
    b = base[base.file.isin(win)][["file", "s", "keep"]].rename(columns={"keep": "keep_base"})
    m = pd.concat([g.merge(b.assign(s=b.s + d), on=["file", "s"]) for d in (0, -1, 1)]).drop_duplicates(["file", "s"])
    m["flip"] = m.keep != m.keep_base
    return m


rows, layers = [], []
for cond in GRID + CONTAM:
    if not os.path.exists(f"{C.PASSES}/{cond}/injections.parquet"):
        print("missing", cond); continue
    ev, cl = S.load_pass(cond)
    inj = pd.read_parquet(f"{C.PASSES}/{cond}/injections.parquet")
    st = status(ev, cl, inj)
    st["admitted"] = st.layer.isin(["retained", "abstained"]); st["rejected"] = st.layer == "abstained"
    lay = st.layer.value_counts(normalize=True).rename(cond); layers.append(lay)
    adm = st[st.admitted]
    rec = S.cluster_rate(adm, "rejected")
    fl = untouched_flips(ev, cl, inj)
    flr = S.cluster_rate(fl, "flip")
    row = dict(cond=cond, injections=len(st), admitted_rate=st.admitted.mean(), recall=rec["rate"], recall_lo=rec["lo"],
               recall_hi=rec["hi"], n_admitted=rec["n"], untouched_n=flr["n"], flip_rate=flr["rate"], flip_hi=flr["hi"],
               judged=cond in JUDGED, inconclusive=st.admitted.mean() < 0.30)
    if cond.startswith("N1"):
        # excess rejection = paired change on the same admitted targets: rejected after N1 minus rejected in the
        # baseline pass (SPEC amendment A4); participant-cluster bootstrap. Admission loss reported separately.
        bk = base_key.loc[list(zip(adm.file, inj.loc[adm.index, "target_s"].astype(int))), "keep"].values
        adm = adm.assign(diff=adm.rejected.astype(float) - (~bk.astype(bool)).astype(float))
        a = adm.groupby("pid")["diff"].agg(["sum", "count"])
        idx = np.random.default_rng(C.SEED).integers(0, len(a), (S.B, len(a)))
        A = a.values
        ex = A[idx, 0].sum(1) / A[idx, 1].sum(1)
        row.update(excess=A[:, 0].sum() / A[:, 1].sum(), excess_lo=C.pct(ex)[0], excess_hi=C.pct(ex)[1],
                   admission_loss=1 - st.admitted.mean())
    if cond.startswith("C_p"):
        # marker displacement vs baseline, filtered vs unfiltered (R1 recordings eligible in both)
        for filt in (True, False):
            mb = C.marker(base[base["round"] == 1], filt=filt).set_index(["pid", "session"]).A
            mi = C.marker(cl, filt=filt).set_index(["pid", "session"]).A
            j = mb.index.intersection(mi.index)
            d = (mi.loc[j] - mb.loc[j]).abs()
            row[f"absdA_{'centred' if filt else 'unfilt'}_med"] = d.median()
            row[f"absdA_{'centred' if filt else 'unfilt'}_p95"] = d.quantile(.95)
    rows.append(row)
    st.to_parquet(f"{C.PASSES}/{cond}/status.parquet")
R = pd.DataFrame(rows)
R.to_csv("results/b3a_injection.csv", index=False)
pd.concat(layers, axis=1).T.to_csv("results/b3a_layers.csv")

# ---- reliability recovery (REC): R1 S1-S2 floor 5, intersection of all four
rec_rows = []
if os.path.exists(f"{C.PASSES}/REC/events.parquet"):
    _, rc = S.load_pass("REC")
    W = {(k, filt): C.s1s2(C.marker(c, filt=filt)) for k, c in (("base", base[base["round"] == 1]), ("rec", rc)) for filt in (True, False)}
    inter = W[("base", True)].index
    for w in W.values():
        inter = inter.intersection(w.index)
    idx = C.boot_idx(len(inter))
    v = {k: (w.loc[inter, 1].values, w.loc[inter, 2].values) for k, w in W.items()}
    for k, (a, b) in v.items():
        rec_rows.append(dict(run=k[0], filt="centred" if k[1] else "unfiltered", **S.pair_stats(a, b, idx)))
    rec_rows.append(dict(run="rec-base", filt="centred", **{f"d_{k}": x for k, x in S.paired_delta(*v[("rec", True)], *v[("base", True)], idx).items()}))
    rec_rows.append(dict(run="rec", filt="centred-unfiltered", **{f"d_{k}": x for k, x in S.paired_delta(*v[("rec", True)], *v[("rec", False)], idx).items()}))
    pd.DataFrame(rec_rows).to_csv("results/b3a_recovery.csv", index=False)

# ---- independent agreement with the per-eye matched detector (real data, no injection)
pe = pd.read_parquet(f"{C.PASSES}/pereye/events.parquet")
pew = {f: g[["s", "e1"]].values for f, g in pe.groupby("file")}
matched = []
for f, g in base.groupby("file"):
    w = pew.get(f, np.zeros((0, 2), int))
    m = np.any((g.s.values[:, None] <= w[:, 1]) & (g.e1.values[:, None] >= w[:, 0]), axis=1) if len(w) else np.zeros(len(g), bool)
    matched.append(pd.Series(m, index=g.index))
base["pe_matched"] = pd.concat(matched)
base["unmatched"] = ~base.pe_matched


def agree(d):
    a = ((~d.keep) & d.unmatched).sum(); b = ((~d.keep) & ~d.unmatched).sum()
    c = (d.keep & d.unmatched).sum(); e = (d.keep & ~d.unmatched).sum()
    orr = (a + .5) * (e + .5) / ((b + .5) * (c + .5))
    n = a + b + c + e; po = (a + e) / n
    pe_ = ((a + b) * (a + c) + (c + e) * (b + e)) / n ** 2
    return orr, (po - pe_) / (1 - pe_), (a, b, c, e)


orr, kap, tab = agree(base)
cells = np.stack([((~base.keep) & base.unmatched), ((~base.keep) & ~base.unmatched),
                  (base.keep & base.unmatched), (base.keep & ~base.unmatched)], 1).astype(float)
P_ = pd.DataFrame(cells, index=base.pid.values).groupby(level=0).sum().values      # per-participant 2x2 counts
idx = np.random.default_rng(C.SEED).integers(0, len(P_), (S.B, len(P_)))
T = P_[idx].sum(1); a_, b_, c_, e_ = T.T
bo = (a_ + .5) * (e_ + .5) / ((b_ + .5) * (c_ + .5))
n_ = T.sum(1); po = (a_ + e_) / n_; pe_ = ((a_ + b_) * (a_ + c_) + (c_ + e_) * (b_ + e_)) / n_ ** 2
bk = (po - pe_) / (1 - pe_)
AG = dict(OR=orr, OR_lo=C.pct(bo)[0], OR_hi=C.pct(bo)[1], kappa=kap, kappa_lo=C.pct(bk)[0], kappa_hi=C.pct(bk)[1],
          rej_unm=tab[0], rej_m=tab[1], kept_unm=tab[2], kept_m=tab[3],
          p_unmatched_rejected=tab[0] / (tab[0] + tab[1]), p_unmatched_kept=tab[2] / (tab[2] + tab[3]))
pd.DataFrame([AG]).to_csv("results/b3a_agreement.csv", index=False)

# ---- per-eye quality covariates
q = pd.read_parquet(f"{C.PASSES}/base/quality.parquet")
rr = base.groupby("file").keep.agg(rej=lambda x: 1 - x.mean(), n="size").reset_index().merge(q, on="file")
rr["d_rms"] = (rr.rms_l - rr.rms_r).abs(); rr["d_miss"] = (rr.miss_l - rr.miss_r).abs()
rr["d_sd"] = (rr.sd_l - rr.sd_r).abs(); rr["worst_rms"] = rr[["rms_l", "rms_r"]].max(1)
rr.to_csv("results/b3a_quality_by_recording.csv", index=False)
qrows = []
QC = ("d_rms", "d_miss", "d_sd", "worst_rms")
pids = rr.pid.unique(); rows_of = {p: np.where(rr.pid.values == p)[0] for p in pids}
pidx = C.boot_idx(len(pids))                       # participant-cluster draws, shared across covariates
bsq = {c: np.empty(S.B) for c in QC}
for k in range(S.B):
    sel = rr.iloc[np.concatenate([rows_of[pids[i]] for i in pidx[k]])]
    for c in QC:
        d = sel[["rej", c]].dropna(); bsq[c][k] = spearmanr(d.rej, d[c])[0]
for c in QC:
    d = rr[["rej", c]].dropna(); rho = spearmanr(d.rej, d[c])[0]
    qrows.append(dict(covariate=c, n=len(d), n_pid=d.index.map(rr.pid).nunique(), spearman=rho,
                      lo=C.pct(bsq[c])[0], hi=C.pct(bsq[c])[1]))
# sensitivity: R1 S1-S2 (centred, floor 5) excluding participants with a recording in the worst-quality decile
cut = rr.worst_rms.quantile(.9); bad = set(rr[rr.worst_rms > cut].pid)
w = C.s1s2(C.marker(base)); w2 = w[~w.index.isin(bad)]
qrows.append(dict(covariate="ICC_S1S2_excl_worst_decile", n=len(w2), **{k: v for k, v in S.pair_stats(w2[1].values, w2[2].values).items() if k.startswith("icc_a1")}))
pd.DataFrame(qrows).to_csv("results/b3a_quality.csv", index=False)

pd.set_option("display.width", 250)
print(R.round(3).to_string()); print(pd.DataFrame(rec_rows).round(3).to_string())
print(pd.DataFrame([AG]).round(3).to_string()); print(pd.DataFrame(qrows).round(3).to_string())
