#!/usr/bin/env bash
# =============================================================================
# 02_qtl_discovery/01_genotypes/05_rfmix_local_ancestry.sh
# Purpose : Local ancestry with RFMix v2.03-r0 (defaults);
#           references 1000G YRI (n=186) = AFR, CEU (n=183) = EUR
# Paper   : ED Fig. 3b
# Status  : SCAFFOLD — supply reference VCF, sample map, genetic map
# =============================================================================
set -euo pipefail
QUERY=${1:?phased query vcf.gz}
REF=${2:?phased 1000G ref vcf.gz}
SAMPLE_MAP=${3:?ref sample -> AFR/EUR map}
GMAP=${4:?genetic map (chr pos cM)}
OUT=${5:?output dir}
mkdir -p "${OUT}"
for chr in $(seq 1 22); do
  rfmix -f "${QUERY}" -r "${REF}" -m "${SAMPLE_MAP}" -g "${GMAP}" \
        -o "${OUT}/chr${chr}" --chromosome="chr${chr}"
done
