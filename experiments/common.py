"""Shared frozen functions (B0). See experiments/SPEC.md. Run scripts from the project root."""
import glob, os, re, zlib
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

SRC = os.environ.get("GBVR_SRC", "data/collected/gazebasevr/full")
PASSES = os.environ.get("PASS_DIR", "results/passes")
FS, DT = 250.0, 0.004
VTH, MIN_SAC, INTERP = 40.0, 3, 12
MAXDUR, RS_AMP, PROG = 150.0, -15.0, (0.5, 10.0)
LOG2 = np.log(2.0)
FLOOR_PROG = 20
SEED = 20261018
NAME = re.compile(r"S_(\d)(\d{3})_S(\d)_(\d)_TEX\.csv")
CH = ("lx", "ly", "rx", "ry")


def tex_files(rounds=(1, 2, 3)):
    return [f for f in sorted(glob.glob(f"{SRC}/*_TEX.csv")) if int(NAME.search(f)[1]) in rounds]


def meta(f):
    m = NAME.search(os.path.basename(f))
    return dict(round=int(m[1]), pid=int(m[2]), session=int(m[3]), file=os.path.basename(f))


def load(f):
    return {c: v.astype(float) for c, v in pd.read_csv(f, usecols=list(CH)).items()}


def runs(mask):
    e = np.diff(np.r_[0, mask.astype(np.int8), 0])
    return np.where(e == 1)[0], np.where(e == -1)[0]


def svel(a):
    return savgol_filter(np.where(np.isnan(a), 0.0, a), 7, 2, deriv=1, delta=DT)


def interp(a):
    return pd.Series(a).interpolate(limit=INTERP, limit_area="inside").values


def prep(raw):
    """raw dict of channels -> (interpolated channels, raw NaN mask, finite mask, per-eye speeds)."""
    nan = np.zeros(len(raw["lx"]), bool)
    for c in CH:
        nan |= np.isnan(raw[c])
    d = {c: interp(raw[c]) for c in CH}
    ok = np.isfinite(d["lx"] + d["ly"] + d["rx"] + d["ry"])
    sl = np.hypot(svel(d["lx"]), svel(d["ly"]))
    sr = np.hypot(svel(d["rx"]), svel(d["ry"]))
    return d, nan, ok, sl, sr


def eye_ok(d, eye):
    """Own-eye finite mask after interpolation (B4)."""
    return np.isfinite(d[eye + "x"] + d[eye + "y"])


def detect(speed, ok):
    st, en = runs((speed > VTH) & ok)
    k = (en - st) >= MIN_SAC
    return st[k], en[k]


def event_rows(d, nan, sl, sr, st, en, iv_l=None, iv_r=None):
    """One row per event. st/en = union interval (en exclusive; e1 = en, the first sub-threshold sample, as in the
    pilot -- SPEC amendment A1). iv_l/iv_r = per-eye own intervals [start, end) for PV (B4)."""
    n = len(nan)
    lx, ly, rx, ry = d["lx"], d["ly"], d["rx"], d["ry"]
    hl, hr = np.abs(svel(lx)), np.abs(svel(rx))          # horizontal-component speeds (B5 sensitivity)
    out = []
    for i, (s, e) in enumerate(zip(st, en)):
        e1 = min(e, n - 1)
        pre, post = max(s - 1, 0), min(e1 + 20, n - 1)
        a, b = iv_l[i] if iv_l is not None else (s, e1 + 1)
        c, dd = iv_r[i] if iv_r is not None else (s, e1 + 1)
        aLx, aRx = lx[e1] - lx[pre], rx[e1] - rx[pre]
        out.append((s, e1, (e1 - s) * 4.0, aLx, aRx, ly[e1] - ly[pre], ry[e1] - ry[pre],
                    np.nanmax(sl[a:b]), np.nanmax(sr[c:dd]), bool(nan[max(pre - 5, 0):post + 1].any()),
                    0.5 * (lx[s] + rx[s] + lx[e1] + rx[e1]) / 2, 0.5 * (ly[s] + ry[s] + ly[e1] + ry[e1]) / 2,
                    np.nanmax(hl[a:b]), np.nanmax(hr[c:dd])))
    cols = ["s", "e1", "dur_ms", "ampLx", "ampRx", "ampLy", "ampRy", "pvL", "pvR", "nan_raw", "midx", "midy", "pvLh", "pvRh"]
    ev = pd.DataFrame(out, columns=cols)
    ev["ampx"] = 0.5 * (ev.ampLx + ev.ampRx)
    ev["ampy"] = 0.5 * (ev.ampLy + ev.ampRy)
    return ev


def extract(raw, detector="mean"):
    """Primary (mean-speed) or B4 per-eye matched detector -> event table."""
    d, nan, ok, sl, sr = prep(raw)
    if detector == "mean":
        st, en = detect(0.5 * (sl + sr), ok)
        return event_rows(d, nan, sl, sr, st, en)
    pairs = match_pereye(detect(sl, eye_ok(d, "l")), detect(sr, eye_ok(d, "r")))
    if not pairs:
        return event_rows(d, nan, sl, sr, np.array([], int), np.array([], int))
    st = np.array([min(p[0][0], p[1][0]) for p in pairs]); en = np.array([max(p[0][1], p[1][1]) for p in pairs])
    return event_rows(d, nan, sl, sr, st, en, [p[0] for p in pairs], [p[1] for p in pairs])


def match_pereye(L, R, max_gap=3):
    """1:1 matching of per-eye intervals: overlap AND |onset gap| <= max_gap; greatest overlap, then smallest gap;
    exact ties (same best score for a candidate on either side) discarded."""
    (ls, le), (rs, re_) = L, R
    cand = []
    j0 = 0
    for i in range(len(ls)):
        while j0 < len(rs) and re_[j0] <= ls[i]:
            j0 += 1
        j = j0
        while j < len(rs) and rs[j] < le[i]:
            ov = min(le[i], re_[j]) - max(ls[i], rs[j]); gap = abs(int(ls[i]) - int(rs[j]))
            if ov > 0 and gap <= max_gap:
                cand.append((-ov, gap, i, j))
            j += 1
    cand.sort()
    score = {}
    for ov, gap, i, j in cand:
        score.setdefault(("L", i), []).append((ov, gap)); score.setdefault(("R", j), []).append((ov, gap))
    tied = {k for k, v in score.items() if len(v) > 1 and sorted(v)[0] == sorted(v)[1]}
    used_l, used_r, out = set(), set(), []
    for ov, gap, i, j in cand:
        if i in used_l or j in used_r or ("L", i) in tied or ("R", j) in tied:
            continue
        used_l.add(i); used_r.add(j)
        out.append(((ls[i], le[i]), (rs[j], re_[j])))
    out.sort(key=lambda p: min(p[0][0], p[1][0]))
    return out


def classify(ev):
    """Apply exclusions, classes, l, recording-centred abstention. ev must have 'file'."""
    ev = ev[(ev.dur_ms <= MAXDUR) & ~ev.nan_raw].copy()
    ev["cls"] = np.where(ev.ampx < RS_AMP, "rs", np.where((ev.ampx > PROG[0]) & (ev.ampx < PROG[1]), "prog", None))
    ev = ev[ev.cls.notna()].copy()
    ev["l"] = np.log(ev.pvL / ev.pvR)
    ev["keep"] = (ev.l - ev.groupby("file").l.transform("median")).abs() < LOG2
    return ev


def marker(ev, filt=True, rs_floor=5, prog_floor=FLOOR_PROG, col="l"):
    """Per-recording A_r with counts; ineligible recordings dropped (abstain)."""
    e = ev[ev.keep] if filt else ev
    g = e.groupby(["pid", "round", "session", "cls"])[col]
    r = g.median().unstack("cls"); n = g.size().unstack("cls").fillna(0)
    r = r.reindex(columns=["rs", "prog"]); n = n.reindex(columns=["rs", "prog"]).fillna(0)
    r["n_rs"], r["n_prog"] = n.rs, n.prog
    r = r[(r.n_rs >= rs_floor) & (r.n_prog >= prog_floor)].dropna(subset=["rs", "prog"])
    r["A"] = 0.5 * (r.rs - r.prog)
    return r.reset_index()


def s1s2(r, rnd=1):
    return r[r["round"] == rnd].pivot(index="pid", columns="session", values="A").reindex(columns=[1, 2]).dropna()


def rounds(r, rr):
    w = r.pivot_table(index="pid", columns=["round", "session"], values="A")
    need = [(1, 1), (1, 2), (rr, 1), (rr, 2)]
    if not all(c in w.columns for c in need):
        return pd.DataFrame(columns=[1, rr])
    w = w[need].dropna()
    return pd.DataFrame({1: w[[(1, 1), (1, 2)]].mean(1), rr: w[[(rr, 1), (rr, 2)]].mean(1)})


def icc(a, b):
    """ICC(A,1) and ICC(C,1) for paired vectors; vectorised over leading axes (last axis = subjects)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = a.shape[-1]
    gm = 0.5 * (a.mean(-1) + b.mean(-1))
    rm = 0.5 * (a + b)
    msr = 2 * ((rm - gm[..., None]) ** 2).sum(-1) / (n - 1)
    msc = n * ((a.mean(-1) - gm) ** 2 + (b.mean(-1) - gm) ** 2)
    sse = ((a - rm - a.mean(-1)[..., None] + gm[..., None]) ** 2 + (b - rm - b.mean(-1)[..., None] + gm[..., None]) ** 2).sum(-1)
    mse = sse / (n - 1)
    return (msr - mse) / (msr + mse + 2 * (msc - mse) / n), (msr - mse) / (msr + mse)


def boot_idx(n, B=5000, seed=SEED):
    return np.random.default_rng(seed).integers(0, n, (B, n))


def pct(x, q=(2.5, 97.5)):
    return np.percentile(x, q)


def gain(raw, gl, gr):
    """B6: scale each eye's H/V about that channel's recording median."""
    out = dict(raw)
    for c, g in (("lx", gl), ("ly", gl), ("rx", gr), ("ry", gr)):
        med = np.nanmedian(raw[c]); out[c] = med + g * (raw[c] - med)
    return out


def seed_for(cond, f):
    return zlib.crc32(f"{cond}|{os.path.basename(f)}".encode())
