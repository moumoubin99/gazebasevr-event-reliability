"""B5 matched-direction controls (TEX and RAN) + transfer to TEX A_r (SPEC: B5, frozen before RAN unzip).

usage: python experiments/b5.py    (needs results/passes/base; RAN files *_RAN.csv under GBVR_SRC; WORKERS env)
-> results/passes/ran/events.parquet, results/b5_*.csv
D_j = 1/2 (l_left - l_right) on matched leftward/rightward pairs; l = log(PV_L / PV_R). Abduction predicts D > 0.
"""
import glob, os, re, sys
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(__file__))
import common as C, stats as S

RAN = re.compile(r"S_(\d)(\d{3})_S(\d)_(\d)_RAN\.csv")
W0, W1 = 20, 150                                   # 80-600 ms after the target jump (samples)


def ran_one(f):
    m = RAN.search(os.path.basename(f))
    meta = dict(round=int(m[1]), pid=int(m[2]), session=int(m[3]), file=os.path.basename(f))
    df = pd.read_csv(f, usecols=list(C.CH) + ["xT", "yT"])
    raw = {c: df[c].values.astype(float) for c in C.CH}
    ev = C.extract(raw)
    xT, yT = df.xT.values, df.yT.values
    jumps = np.where((np.diff(xT) != 0) | (np.diff(yT) != 0))[0] + 1
    rows = []
    for t in jumps:
        k = np.searchsorted(ev.s.values, t + W0)              # first detected event with onset in window
        if k >= len(ev) or ev.s.values[k] > t + W1:
            continue
        e = ev.iloc[k].to_dict()
        e.update(t=int(t), tdx=xT[t] - xT[t - 1], tdy=yT[t] - yT[t - 1],
                 tmidx=0.5 * (xT[t] + xT[t - 1]), tmidy=0.5 * (yT[t] + yT[t - 1]), lat_ms=(e["s"] - t) * 4.0)
        rows.append(e)
    out = pd.DataFrame(rows)
    if len(out):                                          # one physical event linked to >1 jump is ambiguous: drop
        out = out[~out.duplicated("s", keep=False)]
    for k_, v in meta.items():
        out[k_] = v
    return out


def angle(ax, ay):
    return np.degrees(np.arctan2(np.abs(ay), np.abs(ax)))


def match(g, cal, midx, midy):
    """Greedy 1:1 left<->right matching within a recording; cal = (d|amp|, dmidx, dmidy, dangle)."""
    L, R = g[g.ampx < 0], g[g.ampx > 0]
    if len(L) == 0 or len(R) == 0:
        return []
    aL, aR = np.abs(L.ampx.values)[:, None], np.abs(R.ampx.values)[None]
    dA = np.abs(aL - aR)
    dX = np.abs(L[midx].values[:, None] - R[midx].values[None])
    dY = np.abs(L[midy].values[:, None] - R[midy].values[None])
    angL = np.degrees(np.arctan2(-L.ampy.values, -L.ampx.values))[:, None]   # leftward reversed
    angR = np.degrees(np.arctan2(R.ampy.values, R.ampx.values))[None]
    dG = np.abs((angL - angR + 180) % 360 - 180)
    ok = (dA <= cal[0]) & (dX <= cal[1]) & (dY <= cal[2]) & (dG <= cal[3])
    i, j = np.where(ok)
    cost = dA[i, j] / cal[0] + dX[i, j] / cal[1] + dY[i, j] / cal[2] + dG[i, j] / cal[3]
    out, ul, ur = [], set(), set()
    for q in np.argsort(cost, kind="stable"):
        if i[q] in ul or j[q] in ur:
            continue
        ul.add(i[q]); ur.add(j[q]); out.append((L.index[i[q]], R.index[j[q]]))
    return out


def pairs(pool, cal, midx, midy, col_l):
    """Matched pairs; only events with a finite ratio enter (undefined ratios counted separately by the caller)."""
    pool = pool[np.isfinite(pool[col_l])]
    rows = []
    for f, g in pool.groupby("file"):
        for a, b in match(g, cal, midx, midy):
            rows.append(dict(file=f, pid=g.pid.iat[0], session=g.session.iat[0], ampL=abs(g.ampx[a]), ampR=abs(g.ampx[b]),
                             D=0.5 * (g[col_l][a] - g[col_l][b])))
    return pd.DataFrame(rows, columns=["file", "pid", "session", "ampL", "ampR", "D"])


def centred(pool, col_l):
    return pool[(pool[col_l] - pool.groupby("file")[col_l].transform("median")).abs() < C.LOG2]


def participant_D(p):
    """Session median D, participant D = mean of the two R1 session medians; >= 5 pairs total.
    Returns (D series, session medians (pid x [1, 2]), session counts (pid x [1, 2]))."""
    g = p.groupby(["pid", "session"]).D
    med = g.median().unstack("session").reindex(columns=[1, 2]) if len(p) else pd.DataFrame(columns=[1, 2])
    cnt = g.size().unstack("session").reindex(columns=[1, 2]).fillna(0) if len(p) else pd.DataFrame(columns=[1, 2])
    both = med.notna().all(1)
    D = med[both].mean(1)[cnt[both].sum(1) >= 5]
    return D, med, cnt


# shared participant-ID bootstrap draws (SPEC: shared indices within a table); cohorts are subsets of PIDS
PIDS, WB = None, None


def init_boot(pids):
    global PIDS, WB
    PIDS = np.array(sorted(pids))
    idx = C.boot_idx(len(PIDS))
    WB = np.stack([np.bincount(r, minlength=len(PIDS)) for r in idx])          # (B, |PIDS|) draw multiplicities


def boot(stat, ids, *arrays, min_n=1):
    """Statistic over the shared draws restricted to cohort `ids` (participants resampled by identity);
    draws with fewer than min_n cohort members are NaN (count reported by ci())."""
    w = WB[:, np.searchsorted(PIDS, np.asarray(ids))]
    out = np.full(len(w), np.nan)
    for k in range(len(w)):
        if w[k].sum() >= min_n:
            out[k] = stat(*(np.repeat(a, w[k]) for a in arrays))
    return out


def ci(b):
    """Percentile CI over finite draws, plus number of omitted draws."""
    f = b[np.isfinite(b)]
    return (*(C.pct(f) if len(f) else (np.nan, np.nan)), int(len(b) - len(f)))


def summarize(task, variant, p, n_undef):
    D, med, cnt = participant_D(p)
    row = dict(task=task, variant=variant, n_pairs=len(p), n_undefined_ratio_events=n_undef, n_pid=len(D))
    if len(D) >= 5:
        lo, hi, om = ci(boot(np.mean, D.index, D.values, min_n=1))
        row.update(mean_D=D.mean(), mean_D_lo=lo, mean_D_hi=hi, mean_D_omitted_draws=om, median_D=D.median(), frac_pos=(D > 0).mean())
    rep = med[(cnt >= 5).all(1) & med.notna().all(1)]   # repeatability sample: >= 5 pairs in each session
    if len(rep) >= 5:
        a, b = rep[1].values, rep[2].values
        lo, hi, om = ci(boot(lambda x, y: C.icc(x, y)[0], rep.index, a, b, min_n=2))
        row.update(icc_n=len(rep), icc_a1=C.icc(a, b)[0], icc_a1_lo=lo, icc_a1_hi=hi, icc_omitted_draws=om)
    else:
        row.update(icc_n=len(rep))
    return row, D


if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    out = f"{C.PASSES}/ran"; os.makedirs(out, exist_ok=True)
    files = [f for f in sorted(glob.glob(f"{C.SRC}/*_RAN.csv")) if RAN.search(f)[1] == "1"]
    with Pool(int(os.environ.get("WORKERS", 30))) as pl:
        ran = pd.concat(pl.map(ran_one, files, chunksize=2), ignore_index=True)
    ran.to_parquet(f"{out}/events.parquet")

    # ---- pools (same duration / missingness exclusions)
    ran = ran[(ran.dur_ms <= C.MAXDUR) & ~ran.nan_raw]
    ran = ran[(ran.ampx.abs() >= 3) & (ran.ampx.abs() <= 25) & (angle(ran.ampx, ran.ampy) <= 30) & (np.sign(ran.ampx) == np.sign(ran.tdx))].copy()
    tex, _ = S.load_pass("base")
    tex = tex[(tex["round"] == 1) & (tex.dur_ms <= C.MAXDUR) & ~tex.nan_raw]
    tex = tex[(tex.ampx.abs() >= 1) & (tex.ampx.abs() <= 10) & (angle(tex.ampx, tex.ampy) <= 15)].copy()
    for d in (ran, tex):
        d["l"] = np.log(d.pvL / d.pvR); d["lh"] = np.log(d.pvLh / d.pvRh)
    CAL_T, CAL_R = (1, 2, 1, 10), (2, 3, 2, 10)

    init_boot(set(tex.pid) | set(ran.pid))
    rows, Ds = [], {}
    for task, pool, cal, mx, my in (("TEX", tex, CAL_T, "midx", "midy"), ("RAN", ran, CAL_R, "tmidx", "tmidy")):
        fin = pool[np.isfinite(pool.l)]                  # undefined ratios removed before centring
        for variant, col, src in (("primary", "l", pool), ("horizontal_pv", "lh", pool), ("centred", "l", centred(fin, "l"))):
            nu = int((~np.isfinite(pool[col])).sum())
            p = pairs(src, cal, mx, my, col)
            p.to_parquet(f"{out}/pairs_{task}_{variant}.parquet")
            r, D = summarize(task, variant, p, nu); rows.append(r); Ds[(task, variant)] = D
            if task == "RAN" and variant == "primary":
                for lo, hi in ((3, 10), (15, 25)):
                    pb = p[p.ampL.between(lo, hi) & p.ampR.between(lo, hi)]
                    r, D = summarize(task, f"band_{lo}_{hi}", pb, nu); rows.append(r); Ds[(task, f"band_{lo}_{hi}")] = D
    R = pd.DataFrame(rows)
    R.to_csv("results/b5_matched.csv", index=False)       # D > 0 and RAN repeatability, written before transfer

    # ---- transfer: Spearman rho(participant D, two-session-mean TEX A_r; R1, centred, floor 5)
    A = C.s1s2(C.marker(S.load_pass("base")[1])).mean(1)
    tr = []
    for (task, variant), D in Ds.items():
        j = D.index.intersection(A.index)
        row = dict(task=task, variant=variant, n=len(j), judged=variant == "primary" and len(j) >= 100, inconclusive=len(j) < 100)
        if len(j) >= 5:
            lo, hi, om = ci(boot(lambda x, y: spearmanr(x, y)[0], j, D.loc[j].values, A.loc[j].values, min_n=2))
            row.update(rho=spearmanr(D.loc[j], A.loc[j])[0], rho_lo=lo, rho_hi=hi, rho_omitted_draws=om)
        tr.append(row)
    T = pd.DataFrame(tr)
    T.to_csv("results/b5_transfer.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.round(3).to_string()); print(T.round(3).to_string())
