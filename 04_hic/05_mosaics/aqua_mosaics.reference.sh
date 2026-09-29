#!/usr/bin/env bash
# =============================================================================
# 04_hic/03_aqua_mosaics.sh
# Purpose : Cohort 3D structures with AQuA Tools (Chakraborty et al. 2025)
#   1. extract_bedpe : per-sample loops/structures (BEDPE) at 5 kb
#   2. union_bedpe   : stack + merge coordinates across 15 samples (common scaffold)
#   3. cluster_bedpe : assign mosaic IDs (column 7) by shared anchor overlap
#   4. query_bedpe   : per-sample contact values (inherent and CPM) on the mosaic BEDPE
# Paper   : Fig. 4a; ED 10a–d; 5,130 mosaics / 11,494 protein-coding genes
# Status  : SCAFFOLD — ARGUMENTS BELOW ARE PLACEHOLDERS; replace with the exact
#           AQuA Tools calls (and version) used for the paper.
# =============================================================================
set -euo pipefail
AQUA=${AQUA:?path to AQuA Tools bin}
IN=${IN:?dir with per-sample .hic files}
OUT=${OUT:?results/hic/aqua}
RES=5000
mkdir -p "${OUT}"/{per_sample,query}

for hic in "${IN}"/*.hic; do
  s=$(basename "${hic}" .hic)
  "${AQUA}/extract_bedpe" TODO_ARGS "${hic}" "${RES}" > "${OUT}/per_sample/${s}.bedpe"
done

"${AQUA}/union_bedpe"   TODO_ARGS "${OUT}"/per_sample/*.bedpe > "${OUT}/cohort_union.bedpe"
"${AQUA}/cluster_bedpe" TODO_ARGS "${OUT}/cohort_union.bedpe"  > "${OUT}/cohort_mosaics.bedpe"

for hic in "${IN}"/*.hic; do
  s=$(basename "${hic}" .hic)
  for norm in inherent cpm; do
    "${AQUA}/query_bedpe" TODO_ARGS "${hic}" "${OUT}/cohort_mosaics.bedpe" "${norm}" \
      > "${OUT}/query/${s}.${norm}.bedpe"
  done
done
