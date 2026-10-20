# GazeBaseVR event-construction reliability

Data and code for **Event construction shapes the test–retest reliability of a binocular saccade-velocity contrast in head-mounted display reading**.

Jibin Zhang and Hui Gui contributed equally. Correspondence: Hui Gui, guihui0221@whut.edu.cn.

The study examines how constructing binocular saccade events changes the test–retest reliability of a reading-related velocity contrast. The manuscript has six figures and Supplementary Tables S1–S16. This is a secondary analysis of public, anonymized human eye-movement recordings from GazeBaseVR. No new participants were recruited and no new recordings were collected. The original data descriptor documents the source study and its ethics procedures.

## Interpretation and transparency

The frozen specification, subsequent amendments and post hoc analyses are distinguished in `experiments/SPEC.md`, its original frozen version and the manuscript. The admission gate, joint-cohort comparisons, influence analysis, peak-velocity ceiling, absolute session differences and detector-setting sensitivity include exploratory analyses added after inspecting results. These comparisons diagnose event construction within this dataset; they do not establish external validity, a validated clinical biomarker or a prospectively confirmed benefit. See `README_REPRODUCE.md` for the script-to-analysis map.

## Repository and release assets

Small result CSV files, figure source tables, source code, the analysis specification and LaTeX manuscript sources are versioned here. The large data assets are attached to [release v1.0.0](https://github.com/moumoubin99/gazebasevr-event-reliability/releases/tag/v1.0.0):

| Asset | Contents |
| --- | --- |
| `texran.zip` | Original GazeBaseVR TEX (reading) and RAN (random-saccade) CSV files used for the study. |
| `derived_results.tar.gz` | Full derived output archive, preserving `server/proj`, `server/pw` and `analysis/results` paths; original raw CSV files are supplied separately in `texran.zip`. Repeated identical content can use standard tar hard links; extraction restores the retained paths. |
| `DERIVED_FILE_MANIFEST.csv` | Per-file archive-relative paths, byte counts, SHA-256 digests and provenance for the derived archive. |
| `SHA256SUMS` | SHA-256 digests for the release assets. |

The original data were created by Dillon Lohr, Samantha Aziz, Lee Friedman and Oleg Komogortsev and are available from [Figshare, Version 2](https://doi.org/10.6084/m9.figshare.21308391.v2). The data descriptor is Lohr et al. (2023), [GazeBaseVR, a large-scale, longitudinal, binocular eye-tracking dataset collected in virtual reality](https://doi.org/10.1038/s41597-023-02075-5). Upstream data are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The precise Figshare record and attribution are recorded in `DATA_LICENSE.md`.

Download the four release assets into the repository root, then verify and unpack them into a new directory:

```bash
sha256sum -c SHA256SUMS
mkdir -p release_data
unzip texran.zip -d release_data
tar -xzf derived_results.tar.gz -C release_data
export GBVR_SRC="$PWD/release_data/data"
export PASS_DIR="$PWD/release_data/server/proj/results/passes"
```

Use an archive utility with standard tar hard-link support. Extract to a new directory to preserve the archive paths. Do not extract over an existing analysis run without a backup. The raw ZIP contains a flat `data/` directory with 2,008 original CSV files (1,004 TEX and 1,004 RAN); their contents were retained without modification. The Python loader expects `GBVR_SRC` to point directly to this directory. The archived frozen passes are in `server/proj/results/passes`. `DERIVED_FILE_MANIFEST.csv` provides per-file checksums and provenance inside the derived archive.

`FILE_MANIFEST.sha256` verifies the tracked source and small-result files for this release. Run `sha256sum -c FILE_MANIFEST.sha256` from the repository root. It excludes itself and the large release assets.

## Reproduction

Use Linux, Bash and a Python environment installed from `requirements.txt`. The recorded CPU analysis environment used Python 3.13.11. Full data extraction can be expensive; the archived derived results allow analysis of existing passes. Set `GBVR_SRC` to the flat directory of extracted raw CSV files, `PASS_DIR` to the extracted `results/passes` directory when it is outside the repository, and `SENS_DIR` to a new writable sensitivity-output directory if desired. `WORKERS` controls parallelism and `PY` selects the Python executable.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
export PY=python
export WORKERS=8
# Retain the GBVR_SRC and PASS_DIR values set after extraction above.
PASS_DIR="$PWD/results/passes" bash reproduce.sh passes  # optional full rerun using new local outputs
bash reproduce.sh analyses  # post hoc analyses from existing passes
bash reproduce.sh paper     # regenerate six figures, captions, tables and both PDFs
```

`passes` writes event passes under `PASS_DIR` and summaries under `results/`; `analyses` writes its summaries under `results/`. The optional command above redirects regenerated passes away from the downloaded archive. Run in a scratch clone to preserve the shipped small CSV files for comparison. The paper stage additionally requires TeX Live with `latexmk`, `apacite`, `newtx`, `lineno` and the packages listed in `paper/jov_preamble.tex`. Figure 1 uses CairoSVG and DejaVu Sans fonts. Rebuilding captions and SI tables can replace editorial formatting in the journal-target manuscript; review the rebuilt documents before submission. The released paper source retains the currently prepared journal format.

The historical clean-run check in `CLEAN_RUN.md` reproduced all 32 analysis/source-data CSV files exactly from the existing passes. That check is distinct from re-running the complete original extraction stage. See the file for its inputs and limits.

## File map

| Path | Purpose |
| --- | --- |
| `experiments/` | Frozen extraction, analysis code, detector tests, specification and declared amendments. |
| `results/*.csv` | Small frozen and revision-analysis summaries. |
| `figures/*.py` | Six-figure, table, caption and source-table generation. |
| `figures/source_data/*.csv` | Recording-level and participant-level tables used for the figures. |
| `paper/` | LaTeX manuscript and supplementary source, author information and bibliography. |
| `reproduce.sh`, `README_REPRODUCE.md` | Stage entry point and detailed reproduction instructions. |

## Licensing

GazeBaseVR is third-party material under CC BY 4.0; cite its original creators and data descriptor when using it. This permission does not automatically license the authors' original analysis code. No additional code license is granted here; all rights reserved, and contact the authors for licensing. See `DATA_LICENSE.md` for the data-source notice. The repository contains no server credentials or machine access configuration.
