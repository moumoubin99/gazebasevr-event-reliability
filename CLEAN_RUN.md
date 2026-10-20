# Clean-run check — 20 October 2026

**What was rerun.** In an empty directory (scratch copy of `experiments/` and `figures/make_source_data.py`; frozen
`results/b*.csv` copied; `PASS_DIR` pointing read-only at `results/passes`; fresh `SENS_DIR`), the `analyses` stage was
run end to end with the versions in `requirements.txt`, `WORKERS=8`, on a CPU workstation:
`prewrite.py` → `make_source_data.py` → `r1_revision.py` → `r2_revision.py` → `r3_revision.py v7 v8 v9`
(about 2 minutes, including the V9 re-extraction of all 814 R1 recordings under seven detector settings).

**Result.** All 32 analysis CSV files (`results/pw_*.csv`, `results/rv_*.csv`, `figures/source_data/*.csv`) matched the
shipped files exactly (same shape and columns; maximum absolute numeric difference 0; non-numeric columns identical).
The V9 reference re-extraction reproduced the frozen R1 passes exactly (203,307 mean-speed and 157,381 per-eye events;
identical windows; maximum absolute difference in peak velocity or amplitude 0).

**Not rerun here.** The full `passes` stage (all rounds, injection and gain conditions; `experiments/run_all.sh`) was run on
the rented server and was not repeated; the V9 reference covers the R1 extraction code path for the base and per-eye
detectors. The `paper` stage is rerun at every build.
