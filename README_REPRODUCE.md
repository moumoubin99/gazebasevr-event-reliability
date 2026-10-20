# Reproducing "Event construction shapes the test–retest reliability of a binocular saccade-velocity contrast in head-mounted display reading"

## Data
GazeBaseVR (Lohr et al., 2023), reading task (TEX) and random-saccade task (RAN) recordings, unzipped into
`data/collected/gazebasevr/full` (or set `GBVR_SRC`). The original TEX/RAN subset is supplied as `texran.zip` in
[release v1.0.0](https://github.com/moumoubin99/gazebasevr-event-reliability/releases/tag/v1.0.0),
under the upstream CC BY 4.0 license. Full derived outputs are in `derived_results.tar.gz`.
Verify the release assets with `sha256sum -c SHA256SUMS`. The dataset contains public,
anonymized human eye-movement recordings. This secondary analysis collected no new
recordings and recruited no new participants; original ethics procedures are documented
in the data descriptor.
The original data descriptor is https://doi.org/10.1038/s41597-023-02075-5;
the upstream Version 2 data record is https://doi.org/10.6084/m9.figshare.21308391.v2.
See `DATA_LICENSE.md` for complete attribution and license information.
`GBVR_SRC` must directly contain the task CSV files (the loader is not recursive).

## Specification and amendments
- `experiments/SPEC.md` — frozen specification (18 October 2026) with amendments A1–A6 (made after code review, before any full-data pass).
- `experiments/SPEC_20261018_frozen_v0.md` — pre-amendment text.
- Post hoc analyses added after the results were seen are declared in each script's docstring and in the Methods
  ("Transparency, specification and amendments"; "Post hoc comparators"):
  `prewrite.py` (admitted constructions, joint cohorts), `r1_revision.py`, `r2_revision.py` (leave-one-participant-out),
  `r3_revision.py` (V7 peak-velocity ceiling, V8 absolute session differences, V9 detector settings).

## Environment
`requirements.txt` (Python 3.13.11 for the local analyses). TeX Live 2023 with `latexmk` for the PDF.
`paper/references.bib` is the checked bibliography; use the shipped file.

## Order
```
bash reproduce.sh passes     # experiments/run_all.sh: unit tests, extraction passes (base, pereye, injections, gain grid, REC, RAN) and frozen analyses B1–B6
bash reproduce.sh analyses   # prewrite.py -> figures/make_source_data.py -> r1/r2/r3_revision.py (CPU)
bash reproduce.sh paper      # figures/gen_*.py -> gen_si_tables.py -> gen_latex_includes.py -> figures/split_display.py -> latexmk
```
| Stage | Inputs | Outputs |
|---|---|---|
| passes | GazeBaseVR TEX/RAN | `results/passes/<condition>/events.parquet`, `results/b*.csv` |
| analyses | passes | `results/pw_*.csv`, `figures/source_data/*.csv`, `results/rv_*.csv`, `results/passes_sens/<variant>/` |
| paper | `results/*.csv`, source data | `figures/*.pdf`, `figures/latex_includes.tex`, `figures/SI_TABLES.tex`, `paper/display/*.tex`, `paper/main.pdf`, `paper/supplementary_main.pdf` |

Every number in figures, captions and SI tables is read from these CSV files by the scripts; captions are generated
by `figures/gen_latex_includes.py` (edit it, not `paper/display/*.tex`).

## Compute
The full passes were run on a rented multi-core Linux server (`WORKERS` processes). Everything after the passes,
including the V9 re-extraction of all R1 recordings, runs on a CPU workstation (32 threads, `WORKERS=8`); no GPU is used.

## Checks
- `experiments/test_common.py` — unit tests of the detector, classification and marker.
- V9 reference: re-extraction of all R1 recordings with the frozen settings reproduces the frozen event passes exactly
  (203,307 mean-speed and 157,381 per-eye events; identical windows; maximum absolute difference 0;
  `results/rv_v9_reproduction.csv`).
- Clean run (20 October 2026): the `analyses` stage, rerun from the existing passes in an empty directory, reproduced all
  32 analysis CSV files exactly (`CLEAN_RUN.md`).

## Archive layout and journal formatting

The derived archive preserves the source provenance directories `server/proj`, `server/pw`
and `analysis/results`; inspect its listing to locate the desired passes. Set `PASS_DIR`
to the extracted directory containing the pass condition folders. The frozen pass directory is `server/proj/results/passes`. Original raw task CSV
files are separate in `texran.zip`, under its flat `data/` directory. The release also
includes `DERIVED_FILE_MANIFEST.csv` with per-file SHA-256 checksums and provenance. Duplicate-content files in the derived tar archive
may be stored as standard hard links, preserving all retained paths on extraction.

The paper stage rebuilds the figure/table includes from result CSV files and can replace
editorial changes to caption or supplementary-table formatting. Review the resulting
PDFs before submission. It builds both the main manuscript and supplementary document.
This release preparation does not constitute a new full-data extraction run.
