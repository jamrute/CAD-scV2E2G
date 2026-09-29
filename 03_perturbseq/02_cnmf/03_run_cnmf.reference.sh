#!/usr/bin/env bash
# =============================================================================
# 03_perturbseq/06_cnmf/03_run_cnmf.sh
# Purpose : cNMF v1.7.1 — k = 5..10, 100 iterations, seed 42; select k = 8;
#           local density threshold 0.01 (99/800 outlier spectra removed)
# Paper   : ED Fig. 9b–c
# Env     : envs/py-singlecell.yml
# Status  : REFERENCE IMPLEMENTATION — confirm --numgenes (not stated in Methods)
# =============================================================================
set -euo pipefail
WORK=${WORK:?results/perturbseq/06_cnmf}
NAME=${NAME:-hcasmc_perturbseq}
H5AD=${WORK}/input/counts.h5ad
NUMGENES=${NUMGENES:-2000}      # TODO confirm
WORKERS=${WORKERS:-8}

cnmf prepare --output-dir "${WORK}" --name "${NAME}" -c "${H5AD}" \
  -k 5 6 7 8 9 10 --n-iter 100 --seed 42 --numgenes "${NUMGENES}"

for w in $(seq 0 $((WORKERS - 1))); do
  cnmf factorize --output-dir "${WORK}" --name "${NAME}" \
    --worker-index "${w}" --total-workers "${WORKERS}" &
done
wait

cnmf combine --output-dir "${WORK}" --name "${NAME}"
cnmf k_selection_plot --output-dir "${WORK}" --name "${NAME}"      # ED 9c
cnmf consensus --output-dir "${WORK}" --name "${NAME}" \
  --components 8 --local-density-threshold 0.01 --show-clustering  # ED 9b
