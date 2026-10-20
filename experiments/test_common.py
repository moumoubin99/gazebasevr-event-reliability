"""B0 unit tests. python experiments/test_common.py (from project root)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import common as C


def fake_events(rng, n_rec=6):
    rows = []
    for r in range(n_rec):
        for k in range(80):
            cls = "rs" if k < 15 else "prog"
            ampx = -rng.uniform(16, 30) if cls == "rs" else rng.uniform(1, 9)
            pv = rng.uniform(150, 500)
            rows.append(dict(file=f"f{r}", pid=r // 2, round=1, session=r % 2 + 1, dur_ms=40.0, nan_raw=False,
                             ampx=ampx, pvL=pv * np.exp(rng.normal(0, .3)), pvR=pv))
    return pd.DataFrame(rows)


def test_gain_cancellation():
    rng = np.random.default_rng(0)
    ev = fake_events(rng)
    base = C.marker(C.classify(ev))
    for c in (-0.4, 0.13, 0.9):
        e2 = ev.copy(); e2["pvL"] = e2.pvL * np.exp(c)       # constant per-eye gain ratio on a fixed catalogue
        k1, k2 = C.classify(ev), C.classify(e2)
        assert (k1.keep.values == k2.keep.values).all()
        assert np.max(np.abs(C.marker(k2).A.values - base.A.values)) < 1e-12


def test_icc_matches_reference():
    rng = np.random.default_rng(1)
    a = rng.normal(size=50); b = a + rng.normal(0.3, 0.5, 50)
    X = np.c_[a, b]; n, k = X.shape; gm = X.mean()
    msr = k * ((X.mean(1) - gm) ** 2).sum() / (n - 1); msc = n * ((X.mean(0) - gm) ** 2).sum() / (k - 1)
    mse = ((X - X.mean(1, keepdims=True) - X.mean(0, keepdims=True) + gm) ** 2).sum() / ((n - 1) * (k - 1))
    ref2 = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n); ref3 = (msr - mse) / (msr + (k - 1) * mse)
    i2, i3 = C.icc(a, b)
    assert abs(i2 - ref2) < 1e-12 and abs(i3 - ref3) < 1e-12
    I2, I3 = C.icc(np.stack([a, a]), np.stack([b, b]))
    assert np.allclose(I2, ref2) and np.allclose(I3, ref3)


def test_match_pereye():
    L = (np.array([10, 50, 90]), np.array([20, 58, 95]))
    R = (np.array([11, 60, 89]), np.array([21, 66, 96]))
    m = C.match_pereye(L, R)
    assert [(a[0], b[0]) for a, b in m] == [(10, 11), (90, 89)]      # 50 vs 60: gap 10 > 3 -> unmatched


def test_floor():
    rng = np.random.default_rng(2)
    ev = fake_events(rng)
    ev.loc[(ev.file == "f0") & (ev.ampx < -15), "ampx"] = 5.0       # f0 loses all RS
    r = C.marker(C.classify(ev))
    assert "f0" not in set(r.pid.astype(str) + "") and len(r) == 5


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v(); print("PASS", k)
