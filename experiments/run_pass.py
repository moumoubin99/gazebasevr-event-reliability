"""Run one full extraction pass under a named condition (SPEC: B3a / B4 / B6 / baseline).

usage: python experiments/run_pass.py <cond> [<cond> ...]   (from project root; WORKERS env, default 30)
conditions: base, pereye, gain_q<q>_d<delta>, gain_asym_p, gain_asym_n,
            F1_a<a>, F2a_f<f>, F2b_b<b>, F3, N1_k<k>, C_p<pct>, REC
Writes results/passes/<cond>/events.parquet (+ injections.parquet for B3a, quality.parquet for base).
"""
import os, sys
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
import common as C

TAU = 0.15
_BASE = None


def base_events():
    global _BASE
    if _BASE is None:
        ev = pd.read_parquet(f"{C.PASSES}/base/events.parquet")
        cl = C.classify(ev)
        ev = ev.merge(cl[["file", "s", "cls", "keep"]], on=["file", "s"], how="left")
        _BASE = {f: g.sort_values("s").reset_index(drop=True) for f, g in ev.groupby("file")}
    return _BASE


def rec_targets():
    ev = base_events()
    r1 = sorted({int(f[3:6]) for f in ev if f.startswith("S_1")})
    rng = np.random.default_rng(C.SEED)
    pids = rng.choice(r1, int(round(0.05 * len(r1))), replace=False)
    return {int(p): int(rng.integers(1, 3)) for p in sorted(pids)}


NDEC = int(5 * TAU / C.DT)


def decay(n, start, amp):
    """Exponential return of an offset amp starting at sample start, for 5 tau."""
    k = np.arange(0, max(min(NDEC, n - start), 0))
    return start, amp * np.exp(-k * C.DT / TAU)


def spaced(cands, k, rng, gap=125, taken=None):
    """Up to k candidates >= gap apart from each other and from `taken` (shared list, updated in place)."""
    out, used = [], ([] if taken is None else taken)
    for i in rng.permutation(len(cands)):
        c = cands[i]
        if all(abs(c - u) >= gap for u in used):
            out.append(i); used.append(c)
        if len(out) >= k:
            break
    return out


def f1(raw, eye, c, a, log):
    x = raw[eye + "x"]; n = len(x)
    if c + 5 >= n:
        return
    k = np.arange(1, 6)
    x[c:c + 5] += a * (1 - np.cos(np.pi * k / 5)) / 2
    st, d = decay(n, c + 5, a)
    x[st:st + len(d)] += d
    log.append(dict(kind="F1", eye=eye, w0=c, w1=c + 8, m0=c, m1=min(st + len(d), n), param=a))


def f2a(raw, eye, s, e1, fct, log):
    n = len(raw["lx"]); s0 = max(s - 1, 0); s1 = min(e1 + 1, n - 1)
    L = s1 - s0 + 1; Lp = min(int(np.ceil(L * fct)), n - s0)
    t = s0 + np.arange(Lp) / fct
    for ch in (eye + "x", eye + "y"):
        x = raw[ch]; seg = x[s0:s1 + 1].copy()
        x[s0:s0 + Lp] = np.interp(t, np.arange(s0, s1 + 1), seg)   # NaN in seg propagates locally
    log.append(dict(kind="F2a", eye=eye, w0=s - 2, w1=e1 + 2, m0=s0, m1=s0 + Lp, param=fct))


def f2b(raw, eye, s, e1, b, sp, log):
    pk = s + int(np.nanargmax(sp[s:e1 + 1]))
    raw[eye + "x"][max(pk - 1, 0):pk + 2] += b
    log.append(dict(kind="F2b", eye=eye, w0=s - 2, w1=e1 + 2, m0=max(pk - 1, 0), m1=pk + 2, param=b))


def f3(raw, eye, s, e1, rng, log):
    k = int(rng.integers(3, 13)); st = max(s + int(rng.integers(-10, 11)), 0)
    raw[eye + "x"][st:st + k] = np.nan; raw[eye + "y"][st:st + k] = np.nan
    log.append(dict(kind="F3", eye=eye, w0=s - 2, w1=e1 + 2, m0=st, m1=st + k, param=k))


def n1(raw, eye, s, e1, kap, log):
    n = len(raw["lx"]); pre = max(s - 1, 0)
    for ch in (eye + "x", eye + "y"):
        x = raw[ch]; x0 = x[pre]
        off = (kap - 1) * (x[e1] - x0)
        x[s:e1 + 1] = x[s:e1 + 1] + (kap - 1) * (x[s:e1 + 1] - x0)
        if np.isfinite(off):
            st, d = decay(n, e1 + 1, off); x[st:st + len(d)] += d
    log.append(dict(kind="N1", eye=eye, w0=s - 2, w1=e1 + 2, m0=s, m1=min(e1 + 1 + NDEC, n), param=kap))


def inject(cond, f, raw, rec):
    """Apply B3a condition to raw (in place). Returns injection log."""
    b = base_events().get(os.path.basename(f))
    log = []
    if b is None or len(b) == 0:
        return log
    rng = np.random.default_rng(C.seed_for(cond, f))
    m = C.meta(f)
    elig = b[b.keep.fillna(False).astype(bool)]
    # fixation intervals >= 300 ms between consecutive detected events with no raw NaN (tracking loss excluded)
    nan = np.zeros(len(raw["lx"]), bool)
    for ch in C.CH:
        nan |= np.isnan(raw[ch])
    cum = np.r_[0, np.cumsum(nan)]
    f0, f1_ = b.e1.values[:-1], b.s.values[1:]          # fixation starts at e1 = first sub-threshold sample (A1)
    valid = ((f1_ - f0) >= 75) & (cum[np.minimum(f1_, len(nan))] - cum[np.minimum(f0, len(nan))] == 0)
    fix_c = (f0 + 20)[valid]
    eye = lambda: "l" if rng.random() < 0.5 else "r"
    _, _, _, sl, sr = C.prep(raw) if cond.startswith("F2b") else (None,) * 5

    def pick_events(k_rs, k_prog, frac=None):
        idx, taken = [], []
        for cl, k in (("rs", k_rs), ("prog", k_prog)):
            e = elig[elig.cls == cl]
            if frac is not None:
                kk = int(round(frac * len(e))); idx += list(e.index[rng.permutation(len(e))[:kk]])
            else:
                idx += [e.index[i] for i in spaced(e.s.values, k, rng, taken=taken)]
        return b.loc[sorted(idx)]

    if cond.startswith("F1"):
        a = float(cond.split("_a")[1])
        for i in spaced(fix_c, 6, rng):
            f1(raw, eye(), int(fix_c[i]), a, log)
    elif cond == "REC":
        if m["round"] != 1 or rec.get(m["pid"]) != m["session"]:
            return log
        tg = pick_events(0, 0, frac=0.30)
        for _, r in tg.iterrows():
            f2a(raw, eye(), int(r.s), int(r.e1), 3.0, log)
        far = [c for c in fix_c if all(abs(c - s) >= 25 for s in tg.s.values)]
        for c in rng.permutation(far)[:int(round(0.10 * len(fix_c)))]:
            f1(raw, eye(), int(c), 4.0, log)
    else:
        if cond.startswith("C_p"):
            tg = pick_events(0, 0, frac=int(cond[3:]) / 100)
        else:
            tg = pick_events(3, 3)
        for _, r in tg.iterrows():
            s, e1, ey = int(r.s), int(r.e1), eye()
            if cond.startswith("F2a") or cond.startswith("C_p"):
                f2a(raw, ey, s, e1, 3.0 if cond.startswith("C_p") else float(cond.split("_f")[1]), log)
            elif cond.startswith("F2b"):
                f2b(raw, ey, s, e1, float(cond.split("_b")[1]), sl if ey == "l" else sr, log)
            elif cond == "F3":
                f3(raw, ey, s, e1, rng, log)
            elif cond.startswith("N1"):
                n1(raw, ey, s, e1, float(cond.split("_k")[1]), log)
            log[-1].update(target_s=s, target_cls=r.cls)
    return log


def quality(raw, ev):
    """Per-eye data quality: missing fraction, fixation RMS-S2S and SD (median over fixations >= 60 ms)."""
    q = {}
    s, e1 = ev.s.values, ev.e1.values
    for eye in ("l", "r"):
        x, y = raw[eye + "x"], raw[eye + "y"]
        q[f"miss_{eye}"] = float(np.mean(np.isnan(x) | np.isnan(y)))
        rms, sd = [], []
        for a, bb in zip(e1[:-1] + 5, s[1:] - 1):
            if bb - a < 15:
                continue
            xx, yy = x[a:bb], y[a:bb]
            dd = np.diff(xx) ** 2 + np.diff(yy) ** 2
            if np.isfinite(dd).sum() >= 10:
                rms.append(np.sqrt(np.nanmean(dd))); sd.append(np.sqrt(np.nanvar(xx) + np.nanvar(yy)))
        q[f"rms_{eye}"] = float(np.median(rms)) if rms else np.nan
        q[f"sd_{eye}"] = float(np.median(sd)) if sd else np.nan
    return q


def gains(cond, m):
    if cond.startswith("gain_asym"):
        d = 0.15 if (cond.endswith("_p") == (m["session"] == 1)) else -0.15
        return np.exp(d / 2), np.exp(-d / 2)
    q, d = cond[len("gain_q"):].split("_d")
    q, d = float(q), float(d)
    return q * np.exp(d / 2), q * np.exp(-d / 2)


def one(args):
    f, cond, rec = args
    m = C.meta(f)
    raw = C.load(f)
    log, q = [], None
    if cond.startswith("gain"):
        raw = C.gain(raw, *gains(cond, m))
    elif cond not in ("base", "pereye"):
        log = inject(cond, f, raw, rec)
    ev = C.extract(raw, "pereye" if cond == "pereye" else "mean")
    for k, v in m.items():
        ev[k] = v
    if cond == "base":
        q = dict(**m, **quality(C.load(f), ev))
    for d in log:
        d.update(m)
    return ev, log, q


def run(cond, workers):
    files = C.tex_files((1,) if cond.startswith(("F", "N1", "C_p", "REC")) else (1, 2, 3))
    step = int(os.environ.get("FILES_STEP", 1))                # smoke tests only: every step-th participant
    files = [f for f in files if C.meta(f)["pid"] % step == 0]
    rec = rec_targets() if cond == "REC" else {}
    out = f"{C.PASSES}/{cond}"; os.makedirs(out, exist_ok=True)
    with Pool(workers) as p:
        res = p.map(one, [(f, cond, rec) for f in files], chunksize=2)
    pd.concat([r[0] for r in res], ignore_index=True).to_parquet(f"{out}/events.parquet")
    inj = [d for r in res for d in r[1]]
    pd.DataFrame([C.meta(f) for f in files]).to_parquet(f"{out}/manifest.parquet")   # input universe of the pass
    if inj:
        pd.DataFrame(inj).to_parquet(f"{out}/injections.parquet")
    if cond == "base":
        pd.DataFrame([r[2] for r in res]).to_parquet(f"{out}/quality.parquet")
    print(cond, len(files), "files", len(inj), "injections", flush=True)


if __name__ == "__main__":
    W = int(os.environ.get("WORKERS", 30))
    for cond in sys.argv[1:]:
        run(cond, W)
