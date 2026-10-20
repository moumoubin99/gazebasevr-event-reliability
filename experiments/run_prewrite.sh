#!/usr/bin/env bash
# Server driver for pre-writing analyses P1-P8 (existing passes only). Run from project root on the server.
set -euo pipefail
PY=${PY:-python3}
mkdir -p results/logs
echo "[$(date +%H:%M:%S)] prewrite"
$PY experiments/prewrite.py > results/logs/prewrite.txt 2>&1
echo "[$(date +%H:%M:%S)] done"
