#!/usr/bin/env bash
# Server driver: full passes + analyses (SPEC v0 + amendments A1-A6). Run from project root on the server.
# usage: bash experiments/run_all.sh  (GBVR_SRC must point at the unzipped TEX/RAN directory)
set -euo pipefail
PY=${PY:-python3}
export WORKERS=${WORKERS:-120}
mkdir -p results/logs
log() { echo "[$(date +%H:%M:%S)] $*"; }

log tests;  $PY experiments/test_common.py
log base;   $PY experiments/run_pass.py base
GRID="pereye F1_a2 F1_a4 F1_a8 F2a_f1.5 F2a_f2 F2a_f3 F2b_b3 F2b_b6 F3 N1_k0.85 N1_k1.15 C_p5 C_p10 C_p20 C_p30 REC"
GAIN=""
for q in 0.9 1 1.1; do for d in -0.3 -0.15 0 0.15 0.3; do GAIN="$GAIN gain_q${q}_d${d}"; done; done
GAIN="$GAIN gain_asym_p gain_asym_n"
log passes; $PY experiments/run_pass.py $GRID $GAIN
log b1b2;   $PY experiments/analyze_b1_b2.py base   > results/logs/b1_b2.txt
log b4;     $PY experiments/analyze_b1_b2.py pereye > results/logs/b4.txt
log b3a;    $PY experiments/analyze_b3a.py          > results/logs/b3a.txt
log b6;     $PY experiments/analyze_b6.py           > results/logs/b6.txt
log b5;     $PY experiments/b5.py                   > results/logs/b5.txt
log done
