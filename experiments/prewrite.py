"""Pre-writing analyses P1-P8 (review round 1; review-stage/PRE_WRITING_ANALYSES.md). Exploratory / post-hoc:
existing passes only, no new extraction, no threshold changes.

usage: python experiments/prewrite.py   (needs passes base, pereye, ran, REC, C_p*, and C_p*/status.parquet from B3a)
-> results/pw_p1_joint_icc.csv ... results/pw_p8_cross_transfer.csv, results/logs/prewrite.txt (via stdout)
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
sys.path.insert(0, os.path.dirname(__file__))
import common as C, stats as S
import b5

FLOORS = (1, 5, 10)
os.makedirs("results", exist_ok=True)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)


def show(name, df):
    print(f"\n==== {name}"); print(df.round(4).to_string()); df.to_csv(f"results/pw_{name}.csv", index=False)


# ---------------------------------------------------------------- methods: event tables
_, base = S.load_pass("base")
_, pe = S.load_pass("pereye")

# P3 admission-only intermediate: mean-speed events overlapping any per-eye matched event (same rule as B3a
# agreement), original mean-speed windows / peaks / classes kept; centring recomputed within the admitted set.
pe_all = pd.read_parquet(f"{C.PASSES}/pereye/events.parquet")
pew = {f: g[["s", "e1"]].values for f, g in pe_all.groupby("file")}
sup = []
for f, g in base.groupby("file"):
    w = pew.get(f, np.zeros((0, 2), int))
    m = np.any((g.s.values[:, None] <= w[:, 1]) & (g.e1.values[:, None] >= w[:, 0]), axis=1) if len(w) else np.zeros(len(g), bool)
    sup.append(pd.Series(m, index=g.index))
base["pe_support"] = pd.concat(sup)
adm = base[base.pe_support].copy()
adm["keep"] = (adm.l - adm.groupby("file").l.transform("median")).abs() < C.LOG2

EV = {"mean": base, "adm": adm, "pereye": pe}
METHODS = {f"{k}_{'c' if filt else 'u'}": (k, filt) for k in EV for filt in (False, True)}


def markers(fl):
    return {m: C.marker(EV[k], filt=filt, rs_floor=fl) for m, (k, filt) in METHODS.items()}


def comps(r):
    """R1_S1S2, R1_R3 (two-session means), R1S1_R3S1 -> DataFrame(index=pid, [o1, o2])."""
    c = S.comparisons(r)
    return {k: c[k] for k in ("R1_S1S2", "R1_R3", "R1S1_R3S1") if k in c}


MK = {fl: markers(fl) for fl in FLOORS}
CP = {fl: {m: comps(r) for m, r in MK[fl].items()} for fl in FLOORS}

# ---------------------------------------------------------------- P1 joint-cohort ICC (+ P3 intermediate)
DELTAS = [("mean_c", "mean_u"), ("adm_u", "mean_u"), ("pereye_u", "adm_u"), ("pereye_u", "mean_u"),
          ("pereye_c", "pereye_u"), ("adm_c", "adm_u"), ("pereye_c", "mean_c"), ("adm_c", "mean_c"), ("pereye_u", "mean_c")]
rows, drows = [], []
for fl in FLOORS:
    for comp in ("R1_S1S2", "R1_R3", "R1S1_R3S1"):
        for cohort, ms in (("joint4", ["mean_u", "mean_c", "pereye_u", "pereye_c"]), ("joint6", list(METHODS))):
            ids = None
            for m in ms:
                i = CP[fl][m][comp].index
                ids = i if ids is None else ids.intersection(i)
            if len(ids) < 5:
                rows.append(dict(floor=fl, comparison=comp, cohort=cohort, n=len(ids))); continue
            idx = C.boot_idx(len(ids))
            V = {m: (CP[fl][m][comp].loc[ids, "o1"].values, CP[fl][m][comp].loc[ids, "o2"].values) for m in METHODS
                 if all(ids.isin(CP[fl][m][comp].index))}
            for m in ms:
                ps = S.pair_stats(*V[m], idx)
                rows.append(dict(floor=fl, comparison=comp, cohort=cohort, method=m, n=len(ids),
                                 own_n=len(CP[fl][m][comp]), own_icc=C.icc(CP[fl][m][comp].o1, CP[fl][m][comp].o2)[0],
                                 **{k: ps[k] for k in ("icc_a1", "icc_a1_lo", "icc_a1_hi", "spearman")}))
            for x, y in DELTAS:
                if x in ms and y in ms:
                    drows.append(dict(floor=fl, comparison=comp, cohort=cohort, contrast=f"{x}-{y}", **S.paired_delta(*V[x], *V[y], idx)))
show("p1_joint_icc", pd.DataFrame(rows))
show("p1_joint_delta", pd.DataFrame(drows))

# ---------------------------------------------------------------- P2 recording-level A agreement between methods
PAIRS2 = [("mean_c", "pereye_c"), ("mean_u", "pereye_u"), ("mean_c", "mean_u"), ("pereye_c", "pereye_u"),
          ("mean_u", "adm_u"), ("adm_u", "pereye_u"), ("mean_c", "pereye_u")]
rows = []
for fl in FLOORS:
    R = {m: r.set_index(["pid", "round", "session"]).A for m, r in MK[fl].items()}
    for x, y in PAIRS2:
        j = R[x].index.intersection(R[y].index)
        a, b = R[x].loc[j].values, R[y].loc[j].values
        d = b - a
        # participant-cluster bootstrap for correlations and bias
        pid = j.get_level_values("pid").values; up = np.unique(pid); rows_of = [np.where(pid == p)[0] for p in up]
        bi = C.boot_idx(len(up), 1000)
        bp, bs, bb = [], [], []
        for k in range(len(bi)):
            s = np.concatenate([rows_of[i] for i in bi[k]])
            bp.append(np.corrcoef(a[s], b[s])[0, 1]); bs.append(spearmanr(a[s], b[s])[0]); bb.append(d[s].mean())
        rows.append(dict(floor=fl, x=x, y=y, n_rec=len(j), n_pid=len(up), pearson=pearsonr(a, b)[0], pearson_lo=C.pct(bp)[0],
                         pearson_hi=C.pct(bp)[1], spearman=spearmanr(a, b)[0], spearman_lo=C.pct(bs)[0], spearman_hi=C.pct(bs)[1],
                         bias=d.mean(), bias_lo=C.pct(bb)[0], bias_hi=C.pct(bb)[1], sd_diff=d.std(ddof=1),
                         loa_lo=d.mean() - 1.96 * d.std(ddof=1), loa_hi=d.mean() + 1.96 * d.std(ddof=1),
                         sd_x=a.std(ddof=1), sd_y=b.std(ddof=1),
                         prop_bias_r=pearsonr(0.5 * (a + b), d)[0]))
    # participant-level (R1 two-session mean) agreement, the unit used for transfer
    for x, y in PAIRS2[:2]:
        wx, wy = C.s1s2(MK[fl][x]).mean(1), C.s1s2(MK[fl][y]).mean(1)
        j = wx.index.intersection(wy.index)
        rows.append(dict(floor=fl, x=x + "@pid_R1mean", y=y + "@pid_R1mean", n_pid=len(j), pearson=pearsonr(wx[j], wy[j])[0],
                         spearman=spearmanr(wx[j], wy[j])[0], bias=(wy[j] - wx[j]).mean(), sd_diff=(wy[j] - wx[j]).std(ddof=1)))
show("p2_agreement", pd.DataFrame(rows))

# P3 event-level bookkeeping: what admission removes
ev_rows = []
for nm, e in (("mean_pool", base), ("adm_pool", adm), ("pereye_pool", pe)):
    for cls in ("rs", "prog"):
        x = e[e.cls == cls]
        ev_rows.append(dict(pool=nm, cls=cls, n=len(x), kept_frac=x.keep.mean(), med_absl=x.l.abs().median(),
                            p95_absl=x.l.abs().quantile(.95)))
for cls in ("rs", "prog"):
    x = base[base.cls == cls]
    for sup_, g in x.groupby("pe_support"):
        ev_rows.append(dict(pool=f"mean_pool|pe_support={sup_}", cls=cls, n=len(g), kept_frac=g.keep.mean(),
                            med_absl=g.l.abs().median(), p95_absl=g.l.abs().quantile(.95)))
show("p3_admission_events", pd.DataFrame(ev_rows))

# ---------------------------------------------------------------- P4 leave-one-participant-out influence
def loo(a, b):
    n = len(a); M = np.array([np.r_[0:i, i + 1:n] for i in range(n)])
    return C.icc(a[M], b[M])[0]


rows, top = [], []
SETS = [(fl, "R1_S1S2", "mean_c", "mean_u") for fl in FLOORS] + [(fl, "R1_S1S2", "pereye_c", "pereye_u") for fl in FLOORS] + \
       [(fl, comp, "pereye_u", "mean_u") for fl in (5,) for comp in ("R1_R3", "R1S1_R3S1")] + \
       [(fl, comp, "mean_c", "mean_u") for fl in (5,) for comp in ("R1_R3", "R1S1_R3S1")]
for fl, comp, x, y in SETS:
    for cohort in ("own_x", "own_y", "intersection"):
        cx, cy = CP[fl][x][comp], CP[fl][y][comp]
        if cohort == "intersection":
            ids = cx.index.intersection(cy.index); w, wy = cx.loc[ids], cy.loc[ids]
        else:
            w = cx if cohort == "own_x" else cy; wy = None
        if len(w) < 6:
            continue
        a, b = w.o1.values, w.o2.values
        full = C.icc(a, b)[0]; li = loo(a, b)
        sq = (b - a) ** 2; share = sq / sq.sum()
        dev = ((0.5 * (a + b)) - 0.5 * (a + b).mean()) ** 2
        r = dict(floor=fl, comparison=comp, method=x if cohort != "own_y" else y, cohort=cohort, n=len(a), icc=full,
                 loo_min=li.min(), loo_max=li.max(), loo_range=li.max() - li.min(),
                 n_shift_gt_05=int((np.abs(li - full) > .05).sum()), top1_sqdiff_share=np.sort(share)[-1],
                 top3_sqdiff_share=np.sort(share)[-3:].sum())
        order = np.argsort(-np.abs(li - full))
        for k in (1, 2, 3, 5):                                   # drop the k most influential participants
            keep = np.setdiff1d(np.arange(len(a)), order[:k]); r[f"icc_drop_top{k}"] = C.icc(a[keep], b[keep])[0]
        if wy is not None:                                       # paired delta influence on the intersection cohort
            ay, by = wy.o1.values, wy.o2.values
            dfull = full - C.icc(ay, by)[0]; dl = li - loo(ay, by)
            r.update(dicc=dfull, dloo_min=dl.min(), dloo_max=dl.max(), n_dshift_gt_05=int((np.abs(dl - dfull) > .05).sum()))
            for k in (1, 2, 3, 5):
                o = np.argsort(-np.abs(dl - dfull)); keep = np.setdiff1d(np.arange(len(a)), o[:k])
                r[f"dicc_drop_top{k}"] = C.icc(a[keep], b[keep])[0] - C.icc(ay[keep], by[keep])[0]
        rows.append(r)
        for i in order[:5]:
            top.append(dict(floor=fl, comparison=comp, method=r["method"], cohort=cohort, pid=w.index[i], o1=a[i], o2=b[i],
                            loo_icc=li[i], shift=li[i] - full, sqdiff_share=share[i], mean_dev_sq=dev[i]))
show("p4_loo", pd.DataFrame(rows))
show("p4_loo_top", pd.DataFrame(top))

# ---------------------------------------------------------------- P5 eligibility transitions + retention
rows = []
for nm, e in (("mean", base), ("pereye", pe)):
    g = e.groupby(["pid", "round", "session", "cls"]).keep.agg(n="size", kept="sum").unstack("cls").fillna(0)
    for fl in FLOORS:
        eu = (g[("n", "rs")] >= fl) & (g[("n", "prog")] >= C.FLOOR_PROG)
        ec = (g[("kept", "rs")] >= fl) & (g[("kept", "prog")] >= C.FLOOR_PROG)
        for rnd in (1, 2, 3):
            m = g.index.get_level_values("round") == rnd
            u, c = eu[m], ec[m]
            rows.append(dict(detector=nm, floor=fl, round=rnd, unit="recording", n=int(m.sum()), both=int((u & c).sum()),
                             unfilt_only=int((u & ~c).sum()), centred_only=int((~u & c).sum()), neither=int((~u & ~c).sum()),
                             lost_rs=int((u & ~c & (g[("kept", "rs")][m] < fl)).sum()),
                             lost_prog=int((u & ~c & (g[("kept", "prog")][m] < C.FLOOR_PROG)).sum())))
        # participant-level R1 S1-S2 pair eligibility
        pu = eu[g.index.get_level_values("round") == 1].groupby(level="pid").agg(lambda x: x.sum() == 2 and len(x) == 2)
        pc = ec[g.index.get_level_values("round") == 1].groupby(level="pid").agg(lambda x: x.sum() == 2 and len(x) == 2)
        rows.append(dict(detector=nm, floor=fl, round=1, unit="R1_pair", n=len(pu), both=int((pu & pc).sum()),
                         unfilt_only=int((pu & ~pc).sum()), centred_only=int((~pu & pc).sum()), neither=int((~pu & ~pc).sum())))
show("p5_eligibility", pd.DataFrame(rows))
rows = []
for nm, e in (("mean", base), ("pereye", pe)):
    r = e.groupby(["file", "cls"]).keep.mean().unstack("cls")
    for cls in ("rs", "prog"):
        x = r[cls].dropna()
        rows.append(dict(detector=nm, cls=cls, n_rec=len(x), mean=x.mean(), p05=x.quantile(.05), p25=x.quantile(.25),
                         median=x.median(), p75=x.quantile(.75), p95=x.quantile(.95), frac_rec_below_80=(x < .8).mean(),
                         frac_rec_below_50=(x < .5).mean()))
    x = e.groupby("file").keep.mean()
    rows.append(dict(detector=nm, cls="all", n_rec=len(x), mean=x.mean(), p05=x.quantile(.05), p25=x.quantile(.25), median=x.median(),
                     p75=x.quantile(.75), p95=x.quantile(.95), frac_rec_below_80=(x < .8).mean(), frac_rec_below_50=(x < .5).mean()))
show("p5_retention", pd.DataFrame(rows))


# ---------------------------------------------------------------- P6 + P7 injections
def status(ev, cl, inj):
    """Same rule as analyze_b3a.status: best-overlapping re-detected event and the layer that handled it."""
    clk = cl.set_index(["file", "s"]); by = {f: g for f, g in ev.groupby("file")}; out = []
    for _, j in inj.iterrows():
        g = by.get(j.file); lay = "not_detected"
        if g is not None:
            ov = np.minimum(g.e1.values, j.w1) - np.maximum(g.s.values, j.w0)
            if (ov >= 0).any():
                e = g.iloc[int(np.argmax(ov))]
                lay = ("missingness" if e.nan_raw else "duration" if e.dur_ms > C.MAXDUR else
                       "out_of_class" if (j.file, e.s) not in clk.index else ("retained" if clk.loc[(j.file, e.s)].keep else "abstained"))
        out.append(lay)
    return np.array(out)


def overlaps(inj):
    """Per injection: modified range [m0, m1) overlaps another injection's modified range / evaluation window."""
    om, ow = np.zeros(len(inj), bool), np.zeros(len(inj), bool)
    for f, g in inj.groupby("file"):
        i = g.index.values; m0, m1, w0, w1 = g.m0.values, g.m1.values, g.w0.values, g.w1.values
        mm = (m0[:, None] < m1[None]) & (m0[None] < m1[:, None]); np.fill_diagonal(mm, False)
        mw = (m0[None] <= w1[:, None]) & (w0[:, None] < m1[None]); np.fill_diagonal(mw, False)   # [i] window hit by [j] range
        om[inj.index.get_indexer(i)] = mm.any(1); ow[inj.index.get_indexer(i)] = mw.any(1)
    return om, ow


ev_rec, rc = S.load_pass("REC")
inj = pd.read_parquet(f"{C.PASSES}/REC/injections.parquet").reset_index(drop=True)
inj["layer"] = status(ev_rec, rc, inj)
inj["ov_mod"], inj["ov_win"] = overlaps(inj)
tg = set(zip(inj.pid, inj.session))
rows = []
for fl in FLOORS:
    for filt in (True, False):
        mb = C.marker(base[base["round"] == 1], filt=filt, rs_floor=fl).set_index(["pid", "session"]).A
        mi = C.marker(rc, filt=filt, rs_floor=fl).set_index(["pid", "session"]).A
        T = [k for k in tg]
        eb = [k for k in T if k in mb.index]; ei = [k for k in T if k in mi.index]; j = [k for k in eb if k in mi.index]
        d = (mi.loc[j] - mb.loc[j]) if j else pd.Series(dtype=float)
        rows.append(dict(floor=fl, filt="centred" if filt else "unfiltered", injected_recs=len(T), eligible_base=len(eb),
                         eligible_rec=len(ei), lost=len(set(eb) - set(ei)), gained=len(set(ei) - set(eb)), n_shared=len(j),
                         dA_mean=d.mean(), absdA_med=d.abs().median(), absdA_mean=d.abs().mean(), absdA_max=d.abs().max()))
P6 = pd.DataFrame(rows)
# paired difference of injection effects (|dA| unfiltered - |dA| centred) on recordings eligible in all four, floor 5
M4 = {(k, filt): C.marker(c, filt=filt).set_index(["pid", "session"]).A
      for k, c in (("base", base[base["round"] == 1]), ("rec", rc)) for filt in (True, False)}
j = [k for k in tg if all(k in m.index for m in M4.values())]
du = (M4[("rec", False)].loc[j] - M4[("base", False)].loc[j]).abs().values
dc = (M4[("rec", True)].loc[j] - M4[("base", True)].loc[j]).abs().values
bi = C.boot_idx(len(j)); pdiff = (du - dc)
P6 = pd.concat([P6, pd.DataFrame([dict(floor=5, filt="unfilt-centred |dA|", n_shared=len(j), dA_mean=pdiff.mean(),
                                         absdA_med=np.median(pdiff), dA_lo=C.pct(pdiff[bi].mean(1))[0], dA_hi=C.pct(pdiff[bi].mean(1))[1])])])
# ICC restricted to injected participants (both sessions eligible in base and REC under both filters)
B1_ = base[base["round"] == 1]
W = {(k, filt): C.s1s2(C.marker(cc, filt=filt)) for k, cc in (("base", B1_), ("rec", rc)) for filt in (True, False)}
tp = sorted({p for p, _ in tg}); ids = pd.Index(tp)
for w in W.values():
    ids = ids.intersection(w.index)
if len(ids) >= 5:
    bi = C.boot_idx(len(ids)); V = {k: (w.loc[ids, 1].values, w.loc[ids, 2].values) for k, w in W.items()}
    ex = [dict(floor=5, filt=f"ICC_injected_pids {k[0]}_{'c' if k[1] else 'u'}", n_shared=len(ids),
               **{q: v for q, v in S.pair_stats(*V[k], bi).items() if q.startswith("icc_a1")}) for k in V]
    for filt in (True, False):
        dd = S.paired_delta(*V[("rec", filt)], *V[("base", filt)], bi)
        ex.append(dict(floor=5, filt=f"dICC_injected_pids rec-base {'c' if filt else 'u'}", n_shared=len(ids), icc_a1=dd["delta"],
                       icc_a1_lo=dd["lo"], icc_a1_hi=dd["hi"]))
    P6 = pd.concat([P6, pd.DataFrame(ex)])
show("p6_rec_injected", P6)
lay = inj.groupby("kind").layer.value_counts(normalize=True).unstack().fillna(0)
lay["n"] = inj.groupby("kind").size(); lay["frac_ov_mod"] = inj.groupby("kind").ov_mod.mean(); lay["frac_ov_win"] = inj.groupby("kind").ov_win.mean()
show("p6_rec_layers", lay.reset_index())

rows = []
for cond in ("C_p5", "C_p10", "C_p20", "C_p30", "REC"):
    if cond == "REC":
        ij = inj
    else:
        ij = pd.read_parquet(f"{C.PASSES}/{cond}/injections.parquet").reset_index(drop=True)
        st = pd.read_parquet(f"{C.PASSES}/{cond}/status.parquet").reset_index(drop=True)
        assert len(st) == len(ij) and (st.file.values == ij.file.values).all(), cond
        ij["layer"] = st.layer.values; ij["ov_mod"], ij["ov_win"] = overlaps(ij)
    ij["clean"] = ~(ij.ov_mod | ij.ov_win)
    for kind, g in ij.groupby("kind"):
        for sub, h in (("all", g), ("non_overlapping", g[g.clean]), ("overlapping", g[~g.clean])):
            a = h[h.layer.isin(["retained", "abstained"])].assign(rej=lambda x: x.layer == "abstained")
            cr = S.cluster_rate(a, "rej") if len(a) else dict(rate=np.nan, lo=np.nan, hi=np.nan, n=0)
            rows.append(dict(cond=cond, kind=kind, subset=sub, injections=len(h), frac_of_all=len(h) / len(g),
                             admitted_rate=h.layer.isin(["retained", "abstained"]).mean() if len(h) else np.nan,
                             recall=cr["rate"], recall_lo=cr["lo"], recall_hi=cr["hi"], n_admitted=cr["n"],
                             files_any_overlap=int(h.groupby("file").clean.agg(lambda x: (~x).any()).sum()) if sub == "all" else np.nan))
    if cond.startswith("C_p"):                                   # marker displacement: recordings without any overlap
        _, cl = S.load_pass(cond)
        cln = ij.groupby(["pid", "session"]).clean.all()
        for filt in (True, False):
            mb = C.marker(base[base["round"] == 1], filt=filt).set_index(["pid", "session"])
            mi = C.marker(cl, filt=filt).set_index(["pid", "session"])
            jj = mb.index.intersection(mi.index)
            cf = np.array([bool(cln.get(k, True)) for k in jj])
            for sub, kk in (("rec_no_overlap", jj[cf]), ("rec_with_overlap", jj[~cf])):
                d = (mi.loc[kk, "A"] - mb.loc[kk, "A"]).abs()
                rows.append(dict(cond=cond, kind=f"absdA_{'centred' if filt else 'unfilt'}", subset=sub, injections=len(kk),
                                 absdA_med=d.median(), absdA_p95=d.quantile(.95)))
show("p7_overlap", pd.DataFrame(rows))

# ---------------------------------------------------------------- P8 cross-session transfer (TEX A vs matched D)
rA = C.marker(base[base["round"] == 1])                      # centred, floor 5, R1 session-level A
Aw = rA.pivot(index="pid", columns="session", values="A")
rows = []
for task in ("TEX", "RAN"):
    for variant in ("primary", "centred"):
        p = pd.read_parquet(f"{C.PASSES}/ran/pairs_{task}_{variant}.parquet")
        _, med, cnt = b5.participant_D(p)
        Dw = med.where(cnt >= 5)                                 # session D needs >= 5 pairs in that session (B5 rule)
        for la, ld in ((1, 2), (2, 1), (1, 1), (2, 2)):
            a, d = Aw.get(la), Dw.get(ld)
            if a is None or d is None:
                continue
            j = a.dropna().index.intersection(d.dropna().index)
            r = dict(task=task, variant=variant, A_session=la, D_session=ld, kind="cross" if la != ld else "same", n=len(j))
            if len(j) >= 5:
                bi = C.boot_idx(len(j)); x, y = a.loc[j].values, d.loc[j].values
                bs = S.spearman_v(x[bi], y[bi])
                r.update(rho=spearmanr(x, y)[0], rho_lo=C.pct(bs)[0], rho_hi=C.pct(bs)[1])
            rows.append(r)
        # both cross directions on one common cohort, averaged
        j = Aw.dropna().index.intersection(Dw.dropna().index)
        if len(j) >= 5:
            bi = C.boot_idx(len(j))
            f = lambda x1, x2, y1, y2: 0.5 * (S.spearman_v(x1, y2) + S.spearman_v(x2, y1))
            g = lambda x1, x2, y1, y2: 0.5 * (S.spearman_v(x1, y1) + S.spearman_v(x2, y2))
            v = [Aw.loc[j, 1].values, Aw.loc[j, 2].values, Dw.loc[j, 1].values, Dw.loc[j, 2].values]
            for kind, fn in (("cross_avg_common", f), ("same_avg_common", g)):
                bs = fn(*(x[bi] for x in v))
                rows.append(dict(task=task, variant=variant, kind=kind, n=len(j), rho=fn(*(x[None] for x in v))[0],
                                 rho_lo=C.pct(bs)[0], rho_hi=C.pct(bs)[1]))
            bd = g(*(x[bi] for x in v)) - f(*(x[bi] for x in v))
            rows.append(dict(task=task, variant=variant, kind="same-cross", n=len(j),
                             rho=(g(*(x[None] for x in v)) - f(*(x[None] for x in v)))[0], rho_lo=C.pct(bd)[0], rho_hi=C.pct(bd)[1]))
show("p8_cross_transfer", pd.DataFrame(rows))
print("\nDONE")
