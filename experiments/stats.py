"""Reliability statistics shared by the analysis scripts (SPEC: statistics)."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
import common as C

B = 5000


def load_pass(cond):
    ev = pd.read_parquet(f"{C.PASSES}/{cond}/events.parquet")
    return ev, C.classify(ev)


def spearman_v(a, b):
    """Spearman over last axis (vectorised via ranks)."""
    ra = rankdata(a, axis=-1); rb = rankdata(b, axis=-1)
    ra -= ra.mean(-1, keepdims=True); rb -= rb.mean(-1, keepdims=True)
    return (ra * rb).sum(-1) / np.sqrt((ra ** 2).sum(-1) * (rb ** 2).sum(-1))


def pair_stats(a, b, idx=None):
    """Point estimates + bootstrap CIs for paired vectors a (occasion 1), b (occasion 2)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = len(a)
    if n < 5:
        return dict(n=n)
    idx = C.boot_idx(n, B) if idx is None else idx
    i2, i3 = C.icc(a, b); bi2, bi3 = C.icc(a[idx], b[idx])
    d = b - a; bd = d[idx]
    sp = spearmanr(a, b)[0]; bsp = spearman_v(a[idx], b[idx])
    out = dict(n=n, icc_a1=i2, icc_a1_lo=C.pct(bi2)[0], icc_a1_hi=C.pct(bi2)[1],
               icc_c1=i3, icc_c1_lo=C.pct(bi3)[0], icc_c1_hi=C.pct(bi3)[1],
               spearman=sp, spearman_lo=C.pct(bsp)[0], spearman_hi=C.pct(bsp)[1],
               mean_diff=d.mean(), mean_diff_lo=C.pct(bd.mean(1))[0], mean_diff_hi=C.pct(bd.mean(1))[1],
               sd_diff=d.std(ddof=1), sd_diff_lo=C.pct(bd.std(1, ddof=1))[0], sd_diff_hi=C.pct(bd.std(1, ddof=1))[1])
    out["loa_lo"], out["loa_hi"] = d.mean() - 1.96 * d.std(ddof=1), d.mean() + 1.96 * d.std(ddof=1)
    return out


def comparisons(r):
    """Named paired comparisons from a per-recording marker table -> {name: DataFrame(index=pid, cols=[o1, o2])}."""
    out = {"R1_S1S2": C.s1s2(r).rename(columns={1: "o1", 2: "o2"})}
    for rr in (2, 3):
        w = C.rounds(r, rr)
        out[f"R1_R{rr}"] = w.rename(columns={1: "o1", rr: "o2"})
        ret = set(w.index)
        s = C.s1s2(r); out[f"R1_S1S2_R{rr}cohort"] = s[s.index.isin(ret)].rename(columns={1: "o1", 2: "o2"})
        p = r.pivot_table(index="pid", columns=["round", "session"], values="A")
        if (1, 1) in p.columns and (rr, 1) in p.columns:
            out[f"R1S1_R{rr}S1"] = p[[(1, 1), (rr, 1)]].dropna().set_axis(["o1", "o2"], axis=1)
    return out


def paired_delta(x1, x2, y1, y2, idx):
    """ICC(A,1)(x) - ICC(A,1)(y) on the same pairs and bootstrap draws."""
    dx = C.icc(x1, x2)[0] - C.icc(y1, y2)[0]
    bd = C.icc(x1[idx], x2[idx])[0] - C.icc(y1[idx], y2[idx])[0]
    lo, hi = C.pct(bd)
    return dict(delta=dx, lo=lo, hi=hi, n=len(x1))


def cluster_rate(df, col, by="pid", B_=B, seed=C.SEED):
    """Proportion of df[col] (bool) with participant-cluster bootstrap CI."""
    g = df.groupby(by)[col].agg(["sum", "count"])
    if len(g) == 0:
        return dict(k=0, n=0, rate=np.nan, lo=np.nan, hi=np.nan)
    s, c = g["sum"].values, g["count"].values
    idx = np.random.default_rng(seed).integers(0, len(g), (B_, len(g)))
    bs = s[idx].sum(1) / np.maximum(c[idx].sum(1), 1)
    lo, hi = C.pct(bs)
    return dict(k=int(s.sum()), n=int(c.sum()), rate=s.sum() / c.sum(), lo=lo, hi=hi)
