"""Export recording-/participant-level source data for the paper figures from the existing event passes.
No new extraction, no threshold changes: reuses experiments/common.py, stats.py, b5.py and the P3 admission rule
of experiments/prewrite.py verbatim.

usage (project root): python figures/make_source_data.py
-> figures/source_data/recording_A.csv     (6 constructions x floors 5,10; one row per eligible recording)
-> figures/source_data/participant_D.csv   (B5 participant D, TEX/RAN primary, with transfer A = R1 two-session-mean mean_c)
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "experiments")
import common as C, stats as S, b5

OUT = "figures/source_data"; os.makedirs(OUT, exist_ok=True)

_, base = S.load_pass("base")
_, pe = S.load_pass("pereye")
# admission (P3 rule, copied from prewrite.py): mean-speed event overlapping any per-eye matched event
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

rows = []
for fl in (5, 10):
    for k, ev in EV.items():
        for filt in (False, True):
            r = C.marker(ev, filt=filt, rs_floor=fl)
            r["construction"] = f"{k}_{'c' if filt else 'u'}"; r["floor"] = fl
            rows.append(r)
R = pd.concat(rows, ignore_index=True)[["floor", "construction", "pid", "round", "session", "A", "rs", "prog", "n_rs", "n_prog"]]
R.to_csv(f"{OUT}/recording_A.csv", index=False)
print("recording_A", R.groupby(["floor", "construction"]).size().to_string())

# B5 participant D (primary) + transfer A (as in b5.py: R1 two-session mean of centred mean-speed A, floor 5)
A = C.s1s2(C.marker(base)).mean(1).rename("A_tex_mean_c")
drows = []
for task in ("TEX", "RAN"):
    p = pd.read_parquet(f"{C.PASSES}/ran/pairs_{task}_primary.parquet")
    D, med, cnt = b5.participant_D(p)
    d = pd.DataFrame({"D": D}); d["task"] = task
    d = d.join(A, how="left")
    drows.append(d.reset_index().rename(columns={"index": "pid"}))
P = pd.concat(drows, ignore_index=True)
P.to_csv(f"{OUT}/participant_D.csv", index=False)
print(P.groupby("task").agg(n=("D", "size"), meanD=("D", "mean"), nA=("A_tex_mean_c", "count")))
