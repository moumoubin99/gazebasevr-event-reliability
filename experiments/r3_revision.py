"""Revision analyses (independent review 2026-10-20; post hoc, CPU only, local). Frozen settings are unchanged; every
analysis here is a declared sensitivity analysis added after inspection of the results.

V7  peak-velocity ceiling comparator. Mean-speed events whose per-eye peak speed exceeds a ceiling in either eye are
    dropped after classification (ceiling 1,000 deg/s; 750 and 1,500 as sensitivity values); centring is recomputed
    within the remaining events, as for admission. Constructions ceil_u/c (ceiling only) and ceiladm_u/c (ceiling, then
    admission). Joint cohort = participants eligible (R1 S1-S2) under the six frozen constructions and the four ceiling
    constructions; paired participant-bootstrap contrasts with and without the participant with the largest centred
    squared mean_u session difference; event bookkeeping and recording eligibility reported separately.
V8  size of session differences: median, 90th percentile and mean of |A(S2) - A(S1)| per construction on the joint
    cohort, with paired participant-bootstrap differences between constructions.
V9  detector sensitivity: re-extraction of all R1 recordings under the frozen settings (reference; must reproduce the
    frozen passes) and with one setting changed at a time: speed threshold 30 or 50 deg/s (frozen 40), Savitzky-Golay
    window 5 or 9 samples (frozen 7), per-eye pairing tolerance 1 or 5 samples (frozen 3). Per variant: joint-cohort
    ICC(A,1) for the six constructions, the admission contrast with and without the most influential participant,
    eligibility and admission retention.

usage: python experiments/r3_revision.py [v7] [v8] [v9]   (default all; WORKERS env for V9 extraction, default 16)
-> results/rv_v7_ceiling_icc.csv, rv_v7_ceiling_delta.csv, rv_v7_ceiling_events.csv, rv_v7_ceiling_elig.csv,
   results/rv_v8_absdiff.csv, rv_v8_absdiff_delta.csv,
   results/rv_v9_sens_icc.csv, rv_v9_sens_summary.csv, results/passes_sens/<variant>/{mean,pereye}.parquet
"""
import os, sys
from functools import partial
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.signal import savgol_filter
sys.path.insert(0, os.path.dirname(__file__))
import common as C, stats as S

FL = 5
CEILS = (1000.0, 750.0, 1500.0)
SENS = os.environ.get("SENS_DIR", "results/passes_sens")
FROZEN = ["mean_u", "mean_c", "adm_u", "adm_c", "pereye_u", "pereye_c"]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)


def show(name, df):
    print(f"\n==== {name}"); print(df.round(4).to_string()); df.to_csv(f"results/{name}.csv", index=False)


def admit(base, pe_raw):
    """P3 admission rule (verbatim from prewrite.py): mean-speed event overlapping any per-eye matched event."""
    pew = {f: g[["s", "e1"]].values for f, g in pe_raw.groupby("file")}
    sup = []
    for f, g in base.groupby("file"):
        w = pew.get(f, np.zeros((0, 2), int))
        m = np.any((g.s.values[:, None] <= w[:, 1]) & (g.e1.values[:, None] >= w[:, 0]), axis=1) if len(w) else np.zeros(len(g), bool)
        sup.append(pd.Series(m, index=g.index))
    base = base.copy(); base["pe_support"] = pd.concat(sup)
    return base


def recentre(ev):
    ev = ev.copy(); ev["keep"] = (ev.l - ev.groupby("file").l.transform("median")).abs() < C.LOG2
    return ev


def pools(base, pe):
    """Frozen event pools: mean (classified base), adm (admitted, centring recomputed), pereye (classified)."""
    return {"mean": base, "adm": recentre(base[base.pe_support]), "pereye": pe}


def s1s2_all(EV, fl=FL):
    out = {}
    for k, ev in EV.items():
        for filt in (False, True):
            out[f"{k}_{'c' if filt else 'u'}"] = C.s1s2(C.marker(ev, filt=filt, rs_floor=fl))
    return out


def joint(P, names):
    ids = None
    for m in names:
        ids = P[m].index if ids is None else ids.intersection(P[m].index)
    return ids


def top_pid(P, ids, m="mean_u"):
    d = (P[m].loc[ids, 2] - P[m].loc[ids, 1]).values; sq = (d - d.mean()) ** 2
    return ids[int(np.argmax(sq))], sq.max() / sq.sum()


def contrasts(P, ids, pairs, tag):
    """Paired ICC(A,1) contrasts x - y on ids, with and without the top mean_u contributor."""
    rows = []
    tp, share = top_pid(P, ids)
    for sub, jj in (("all", ids), ("without_top", ids[ids != tp])):
        bi = C.boot_idx(len(jj))
        for x, y in pairs:
            d = S.paired_delta(P[x].loc[jj, 1].values, P[x].loc[jj, 2].values, P[y].loc[jj, 1].values, P[y].loc[jj, 2].values, bi)
            rows.append(dict(tag=tag, subset=sub, contrast=f"{x}-{y}", top_pid=int(tp), top_share=share, **d))
    return rows


def iccs(P, ids, names, tag):
    rows = []
    bi = C.boot_idx(len(ids))
    for m in names:
        ps = S.pair_stats(P[m].loc[ids, 1].values, P[m].loc[ids, 2].values, bi)
        rows.append(dict(tag=tag, method=m, n=len(ids), own_n=len(P[m]), **{k: ps[k] for k in ("icc_a1", "icc_a1_lo", "icc_a1_hi", "spearman")}))
    return rows


def frozen_pools():
    _, base = S.load_pass("base")
    _, pe = S.load_pass("pereye")
    base = admit(base, pd.read_parquet(f"{C.PASSES}/pereye/events.parquet"))
    return base, pe


# ------------------------------------------------------------------------------------------------ V7 + V8
def v7_v8(do7=True, do8=True):
    base, pe = frozen_pools()
    EV = pools(base, pe)
    pk = np.maximum(base.pvL, base.pvR)
    cons = {}
    for ce in CEILS:
        cb = base[pk <= ce]
        cons[f"{int(ce)}"] = {f"ceil{int(ce)}": recentre(cb), f"ceiladm{int(ce)}": recentre(cb[cb.pe_support])}
    P, P10 = s1s2_all(EV), s1s2_all(EV, 10)
    for EVc in cons.values():
        P.update(s1s2_all(EVc)); P10.update(s1s2_all(EVc, 10))
    ids6 = joint(P, FROZEN)
    assert len(ids6) == 376, len(ids6)
    assert abs(C.icc(P["mean_u"].loc[ids6, 1], P["mean_u"].loc[ids6, 2])[0] - .3738) < 5e-4

    if do7:
        rows_i, rows_d, rows_e, rows_el = [], [], [], []
        for fl, Pf in ((FL, P), (10, P10)):
            for t in cons:
                new = [f"ceil{t}_u", f"ceil{t}_c", f"ceiladm{t}_u", f"ceiladm{t}_c"]
                ids = joint(Pf, FROZEN + new)
                tag = f"floor{fl}_ceil{t}"
                rows_i += [dict(floor=fl, ceiling=int(t), **r) for r in iccs(Pf, ids, FROZEN + new, tag)]
                pairs = [(f"ceil{t}_u", "mean_u"), ("adm_u", f"ceil{t}_u"), (f"ceiladm{t}_u", f"ceil{t}_u"),
                         (f"ceiladm{t}_u", "adm_u"), ("adm_u", "mean_u"), (f"ceil{t}_c", "mean_c"), ("mean_c", f"ceil{t}_u")]
                rows_d += [dict(floor=fl, ceiling=int(t), **r) for r in contrasts(Pf, ids, pairs, tag)]
                if fl == FL and t == "1000":                     # leave-one-out range of adm_u - ceil_u
                    x, y = Pf["adm_u"].loc[ids], Pf[f"ceil{t}_u"].loc[ids]
                    n = len(ids); dl = []
                    for i in range(n):
                        k = np.arange(n) != i
                        dl.append(C.icc(x[1].values[k], x[2].values[k])[0] - C.icc(y[1].values[k], y[2].values[k])[0])
                    dl = np.array(dl)
                    rows_d.append(dict(floor=fl, ceiling=1000, tag=tag, subset="loo_range", contrast=f"adm_u-ceil{t}_u",
                                       delta=np.nan, lo=dl.min(), hi=dl.max(), n=n,
                                       top_pid=int(ids[int(np.argmin(dl))]), top_share=np.nan))
            # eligibility: own-cohort R1 S1-S2 n per construction and how many frozen-joint participants each keeps
            for t in cons:
                for m in (f"ceil{t}_u", f"ceil{t}_c", f"ceiladm{t}_u", f"ceiladm{t}_c"):
                    rows_el.append(dict(floor=fl, ceiling=int(t), construction=m, own_n=len(Pf[m]),
                                        frozen_joint_n=len(joint(Pf, FROZEN)), in_frozen_joint=len(joint(Pf, FROZEN).intersection(Pf[m].index))))
        # event bookkeeping on the classified mean-speed pool (all rounds)
        for ce in CEILS:
            over = pk > ce
            for cls in ("rs", "prog", "all"):
                m = (base.cls == cls) if cls != "all" else np.ones(len(base), bool)
                b = base[m]; o = over[m]
                rows_e.append(dict(ceiling=int(ce), cls=cls, n_events=len(b), n_over=int(o.sum()), frac_over=o.mean(),
                                   frac_over_without_support=(~b.pe_support[o]).mean() if o.any() else np.nan,
                                   n_no_support=int((~b.pe_support).sum()),
                                   frac_no_support_over=(o[~b.pe_support.values]).mean(),
                                   n_files_with_over=int(b[o].file.nunique()),
                                   n_over_L=int((b.pvL[o] > ce).sum()), n_over_R=int((b.pvR[o] > ce).sum())))
        # the influential session under the ceiling (frozen-joint top participant)
        tp, _ = top_pid(P, ids6)
        for ce in CEILS:
            g = base[(base.pid == tp) & (base["round"] == 1)]
            for s, gs in g.groupby("session"):
                rs = gs[gs.cls == "rs"]; rk = rs[np.maximum(rs.pvL, rs.pvR) <= ce]
                rows_e.append(dict(ceiling=int(ce), cls=f"pid{tp}_S{s}_rs", n_events=len(rs), n_over=len(rs) - len(rk),
                                   med_l_before=rs.l.median(), med_l_after=rk.l.median() if len(rk) else np.nan))
        # Fig 2b: every path-specific contrast on the floor-5 ceiling-1000 cohort, all / without the top participant,
        # with the single-participant leave-one-out range
        ids = joint(P, FROZEN + ["ceil1000_u", "ceil1000_c", "ceiladm1000_u", "ceiladm1000_c"])
        PATHS = [("adm_u", "mean_u"), ("ceil1000_u", "mean_u"), ("adm_u", "ceil1000_u"), ("pereye_u", "adm_u"), ("adm_c", "adm_u"),
                 ("mean_c", "mean_u"), ("pereye_c", "pereye_u"), ("pereye_c", "mean_c")]
        pr = pd.DataFrame(contrasts(P, ids, PATHS, "floor5_ceil1000_paths"))
        n = len(ids); lo_, hi_ = [], []
        for x, y in PATHS:
            X, Y = P[x].loc[ids], P[y].loc[ids]; dl = []
            for i in range(n):
                k = np.arange(n) != i
                dl.append(C.icc(X[1].values[k], X[2].values[k])[0] - C.icc(Y[1].values[k], Y[2].values[k])[0])
            lo_.append(min(dl)); hi_.append(max(dl))
        loo = pd.DataFrame(dict(contrast=[f"{x}-{y}" for x, y in PATHS], loo_min=lo_, loo_max=hi_))
        show("rv_v7_paths", pr.merge(loo, on="contrast"))
        show("rv_v7_ceiling_icc", pd.DataFrame(rows_i))
        show("rv_v7_ceiling_delta", pd.DataFrame(rows_d))
        show("rv_v7_ceiling_events", pd.DataFrame(rows_e))
        show("rv_v7_ceiling_elig", pd.DataFrame(rows_el))

    if do8:
        names = FROZEN + ["ceil1000_u", "ceiladm1000_u"]
        ids = joint(P, FROZEN + ["ceil1000_u", "ceil1000_c", "ceiladm1000_u", "ceiladm1000_c"])
        tp, _ = top_pid(P, ids)
        rows, drows = [], []

        def metr(d):                                             # d: (..., n) absolute differences
            return np.stack([np.median(d, -1), np.percentile(d, 90, axis=-1), d.mean(-1)], -1)
        MN = ["median", "p90", "mean"]
        for sub, jj in (("all", ids), ("without_top", ids[ids != tp])):
            bi = C.boot_idx(len(jj))
            D = {m: np.abs(P[m].loc[jj, 2].values - P[m].loc[jj, 1].values) for m in names}
            B_ = {m: metr(D[m][bi]) for m in names}
            for m in names:
                pt = metr(D[m])
                r = dict(subset=sub, method=m, n=len(jj), top_pid=int(tp))
                for i, nm in enumerate(MN):
                    r[nm] = pt[i]; r[nm + "_lo"], r[nm + "_hi"] = C.pct(B_[m][:, i])
                rows.append(r)
            for x, y in (("adm_u", "mean_u"), ("ceil1000_u", "mean_u"), ("adm_u", "ceil1000_u"), ("ceiladm1000_u", "ceil1000_u"),
                         ("mean_c", "mean_u"), ("adm_c", "adm_u"), ("pereye_u", "mean_u")):
                pt = metr(D[x]) - metr(D[y]); bd = B_[x] - B_[y]
                r = dict(subset=sub, contrast=f"{x}-{y}", n=len(jj))
                for i, nm in enumerate(MN):
                    r[nm] = pt[i]; r[nm + "_lo"], r[nm + "_hi"] = C.pct(bd[:, i])
                drows.append(r)
        show("rv_v8_absdiff", pd.DataFrame(rows))
        show("rv_v8_absdiff_delta", pd.DataFrame(drows))


# ------------------------------------------------------------------------------------------------ V9
VARIANTS = {"ref": (40.0, 7, 3), "vth30": (30.0, 7, 3), "vth50": (50.0, 7, 3), "sg5": (40.0, 5, 3), "sg9": (40.0, 9, 3),
            "gap1": (40.0, 7, 1), "gap5": (40.0, 7, 5)}
_MATCH = C.match_pereye


def _svel(a, w):
    return savgol_filter(np.where(np.isnan(a), 0.0, a), w, 2, deriv=1, delta=C.DT)


def _one(args):
    f, (vth, win, gap), need_mean = args
    C.VTH = vth; C.svel = partial(_svel, w=win); C.match_pereye = partial(_MATCH, max_gap=gap)
    m = C.meta(f); raw = C.load(f); out = []
    for det in (("mean", "pereye") if need_mean else ("pereye",)):
        ev = C.extract(raw, det)
        for k, v in m.items():
            ev[k] = v
        out.append(ev)
    return out if need_mean else [None] + out


def v9_extract(workers):
    files = C.tex_files((1,))
    for name, par in VARIANTS.items():
        out = f"{SENS}/{name}"
        if os.path.exists(f"{out}/pereye.parquet"):
            continue
        os.makedirs(out, exist_ok=True)
        need_mean = not name.startswith("gap")                     # pairing tolerance only changes per-eye events
        with Pool(workers) as p:
            res = p.map(_one, [(f, par, need_mean) for f in files], chunksize=4)
        if need_mean:
            pd.concat([r[0] for r in res], ignore_index=True).to_parquet(f"{out}/mean.parquet")
        pd.concat([r[1] for r in res], ignore_index=True).to_parquet(f"{out}/pereye.parquet")
        print(name, par, len(files), "files", flush=True)


def v9_analyse():
    ref_mean = pd.read_parquet(f"{SENS}/ref/mean.parquet"); ref_pe = pd.read_parquet(f"{SENS}/ref/pereye.parquet")
    # reproduction check against the frozen passes (R1 only)
    fb = pd.read_parquet(f"{C.PASSES}/base/events.parquet").query("round == 1")
    fp = pd.read_parquet(f"{C.PASSES}/pereye/events.parquet").query("round == 1")
    chk = []
    for nm, a, b in (("mean", ref_mean, fb), ("pereye", ref_pe, fp)):
        a = a.sort_values(["file", "s"]).reset_index(drop=True); b = b.sort_values(["file", "s"]).reset_index(drop=True)
        same_n = len(a) == len(b)
        maxdiff = float(np.nanmax(np.abs(a[["pvL", "pvR", "ampx"]].values - b[["pvL", "pvR", "ampx"]].values))) if same_n else np.nan
        chk.append(dict(detector=nm, n_ref=len(a), n_frozen=len(b), same_windows=same_n and (a[["file", "s", "e1"]].values == b[["file", "s", "e1"]].values).all(),
                        max_abs_diff_pv_amp=maxdiff))
    show("rv_v9_reproduction", pd.DataFrame(chk))
    rows, srows = [], []
    for name, par in VARIANTS.items():
        mean = ref_mean if name.startswith("gap") else pd.read_parquet(f"{SENS}/{name}/mean.parquet")
        pe_raw = pd.read_parquet(f"{SENS}/{name}/pereye.parquet")
        base = admit(C.classify(mean), pe_raw)
        EV = pools(base, C.classify(pe_raw))
        P = s1s2_all(EV)
        ids = joint(P, FROZEN)
        rows += [dict(variant=name, vth=par[0], sg_window=par[1], pair_gap=par[2], **r) for r in iccs(P, ids, FROZEN, name)]
        d = {(r["subset"], r["contrast"]): r for r in contrasts(P, ids, [("adm_u", "mean_u"), ("mean_c", "mean_u"), ("adm_c", "adm_u")], name)}
        a, w = d[("all", "adm_u-mean_u")], d[("without_top", "adm_u-mean_u")]
        tp, share = top_pid(P, ids)
        x, y = P["adm_u"].loc[ids], P["mean_u"].loc[ids]; n = len(ids); dl = []
        for i in range(n):
            k = np.arange(n) != i
            dl.append(C.icc(x[1].values[k], x[2].values[k])[0] - C.icc(y[1].values[k], y[2].values[k])[0])
        cb = base[np.maximum(base.pvL, base.pvR) <= CEILS[0]]
        P.update(s1s2_all({"ceil1000": recentre(cb)}))
        idc = joint(P, FROZEN + ["ceil1000_u"])
        bc = C.boot_idx(len(idc))
        dc = S.paired_delta(*(P[m].loc[idc, k].values for m in ("adm_u", "ceil1000_u") for k in (1, 2)), bc)
        Pm = P["mean_u"]; a198 = [Pm.loc[198, k] if 198 in Pm.index else np.nan for k in (1, 2)]
        srows.append(dict(variant=name, n_joint_ceil=len(idc), icc_ceil_u=C.icc(P["ceil1000_u"].loc[idc, 1], P["ceil1000_u"].loc[idc, 2])[0],
                          d_adm_ceil=dc["delta"], d_adm_ceil_lo=dc["lo"], d_adm_ceil_hi=dc["hi"],
                          A198_meanu_S1=a198[0], A198_meanu_S2=a198[1], vth=par[0], sg_window=par[1], pair_gap=par[2], n_joint=n,
                          own_n_mean_u=len(P["mean_u"]), own_n_adm_u=len(P["adm_u"]),
                          n_events_mean=len(base), n_events_pereye=len(EV["pereye"]), adm_retention=base.pe_support.mean(),
                          icc_mean_u=C.icc(y[1], y[2])[0], icc_adm_u=C.icc(x[1], x[2])[0],
                          delta=a["delta"], delta_lo=a["lo"], delta_hi=a["hi"], top_pid=int(tp), top_share=share,
                          delta_wo_top=w["delta"], delta_wo_top_lo=w["lo"], delta_wo_top_hi=w["hi"],
                          loo_min=min(dl), loo_max=max(dl),
                          d_meanc_meanu=d[("all", "mean_c-mean_u")]["delta"], d_admc_admu=d[("all", "adm_c-adm_u")]["delta"]))
    show("rv_v9_sens_icc", pd.DataFrame(rows))
    show("rv_v9_sens_summary", pd.DataFrame(srows))


if __name__ == "__main__":
    todo = sys.argv[1:] or ["v7", "v8", "v9"]
    if "v7" in todo or "v8" in todo:
        v7_v8("v7" in todo, "v8" in todo)
    if "v9" in todo:
        v9_extract(int(os.environ.get("WORKERS", 16)))
        v9_analyse()
    print("\nDONE")
