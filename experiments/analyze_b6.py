"""B6 gain re-extraction grid: marker and reliability change vs baseline on shared records/pairs.

usage: python experiments/analyze_b6.py   -> results/b6_gain.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import common as C, stats as S

CONDS = [f"gain_q{q}_d{d}" for q in ("0.9", "1", "1.1") for d in ("-0.3", "-0.15", "0", "0.15", "0.3")] + ["gain_asym_p", "gain_asym_n"]
_, base = S.load_pass("base")
MAN = set(pd.read_parquet(f"{C.PASSES}/base/manifest.parquet").file)
U = len(MAN)                                      # fixed universe: all input recordings (SPEC amendment A5)
mb = C.marker(base).set_index(["pid", "round", "session"])
wb = C.s1s2(mb.reset_index())
rows = []
for cond in CONDS:
    if not os.path.exists(f"{C.PASSES}/{cond}/events.parquet"):
        print("missing", cond); continue
    _, cl = S.load_pass(cond)
    assert set(pd.read_parquet(f"{C.PASSES}/{cond}/manifest.parquet").file) == MAN, cond
    mc = C.marker(cl).set_index(["pid", "round", "session"])
    sh = mb.index.intersection(mc.index)
    dA = (mc.loc[sh, "A"] - mb.loc[sh, "A"]).abs()
    wc = C.s1s2(mc.reset_index()); pr = wb.index.intersection(wc.index)
    idx = C.boot_idx(len(pr))
    d = S.paired_delta(wc.loc[pr, 1].values, wc.loc[pr, 2].values, wb.loc[pr, 1].values, wb.loc[pr, 2].values, idx)
    dlt = 0.15 if cond.startswith("gain_asym") else abs(float(cond.split("_d")[1]))
    row = dict(cond=cond, abs_delta=dlt, universe=U, base_records=len(mb), cond_records=len(mc), shared_records=len(sh),
               coverage_loss_pp=100 * (len(mb) - len(mc)) / U, lost_base_eligible=len(mb) - len(sh),
               med_absdA=dA.median(), p95_absdA=dA.quantile(.95), retention=cl.keep.mean(), base_retention=base.keep.mean(),
               n_pairs=len(pr), icc_cond=C.icc(wc.loc[pr, 1].values, wc.loc[pr, 2].values)[0],
               icc_base=C.icc(wb.loc[pr, 1].values, wb.loc[pr, 2].values)[0], d_icc=d["delta"], d_icc_lo=d["lo"], d_icc_hi=d["hi"])
    row["pass"] = (row["med_absdA"] <= .01) and (row["p95_absdA"] <= .03) and (abs(row["d_icc"]) <= .05) and (row["coverage_loss_pp"] <= 5)
    row["judged"] = dlt <= 0.15
    rows.append(row)
R = pd.DataFrame(rows)
os.makedirs("results", exist_ok=True)
R.to_csv("results/b6_gain.csv", index=False)
pd.set_option("display.width", 250)
print(R.round(4).to_string())
JUDGED = [c for c in CONDS if c.startswith("gain_asym") or abs(float(c.split("_d")[1])) <= 0.15]
done = set(R.cond) if len(R) else set()
c4 = "INCOMPLETE" if not set(JUDGED) <= done else ("PASS" if R[R.judged]["pass"].all() else "FAIL")
print(f"C4 ({len(JUDGED)} judged conditions, |delta|<=.15):", c4)
