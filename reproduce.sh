#!/usr/bin/env bash
# Reproduce the paper from the GazeBaseVR TEX/RAN recordings. Run from the project root.
#   bash reproduce.sh passes     # full extraction passes + frozen analyses (server: ~1,000 recordings x 37 conditions)
#   bash reproduce.sh analyses   # post hoc analyses from existing passes (CPU; local is fine)
#   bash reproduce.sh paper      # figures, captions, SI tables and PDF from results/*.csv
#   bash reproduce.sh all        # all three in order
# Env: GBVR_SRC (unzipped TEX/RAN directory; default data/collected/gazebasevr/full), WORKERS (default 8), PY (default python3)
set -euo pipefail
PY=${PY:-python3}
export WORKERS=${WORKERS:-8}
log() { echo "[$(date +%H:%M:%S)] $*"; }
mkdir -p results/logs

passes() {
  PY=$PY bash experiments/run_all.sh                       # tests, base/pereye/injection/gain passes, B1-B6
}
analyses() {
  log prewrite;  $PY experiments/prewrite.py        > results/logs/prewrite.txt
  log source;    $PY figures/make_source_data.py
  log r1;        $PY experiments/r1_revision.py     > results/logs/rv_r1.txt
  log r2;        $PY experiments/r2_revision.py     > results/logs/rv_r2.txt
  log r3;        $PY experiments/r3_revision.py v7 v8 v9 > results/logs/r3_all.txt   # V9 re-extracts R1 (~minutes x WORKERS)
}
paper() {
  for s in figures/gen_fig1_schematic.py figures/gen_fig2_icc.py figures/gen_fig3_influence.py \
           figures/gen_fig4_bland_altman.py figures/gen_fig5_injection.py figures/gen_fig6_matched.py \
           figures/gen_table1_longitudinal.py figures/gen_si_tables.py figures/gen_latex_includes.py; do
    log "$s"; $PY "$s"
  done
  $PY figures/split_display.py
  (cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex supplementary_main.tex)
}
case "${1:-all}" in
  passes) passes ;; analyses) analyses ;; paper) paper ;;
  all) passes; analyses; paper ;;
  *) echo "usage: bash reproduce.sh [passes|analyses|paper|all]"; exit 2 ;;
esac
