#!/usr/bin/env bash
# =============================================================================
# 02_qtl_discovery/01_genotypes/04_genotype_pcs.sh
# Purpose : Genotype principal components (plink --pca); first 4 used as covariates
# Paper   : ED Fig. 3a, 3c
# Status  : REFERENCE IMPLEMENTATION — confirm LD-pruning parameters
# =============================================================================
set -euo pipefail
VCF=${1:?merged vcf.gz}
OUT=${2:?output prefix}
plink --vcf "${VCF}" --double-id --maf 0.05 --indep-pairwise 50 5 0.2 --out "${OUT}.prune"   # TODO confirm
plink --vcf "${VCF}" --double-id --extract "${OUT}.prune.prune.in" --pca 10 --out "${OUT}"
