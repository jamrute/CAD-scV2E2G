#!/usr/bin/env bash
# =============================================================================
# 02_qtl_discovery/01_genotypes/01_array_qc.sh
# Purpose : Illumina GSA-24v3 variant QC (after GenomeStudio v2.0.5 export)
#           Exclude variants with (MAF > 5% AND call rate < 95%)
#                              or (MAF < 5% AND call rate < 99%)
#           No individuals excluded at the <95% sample call-rate filter.
# Env     : envs/genomics-cli.yml (PLINK 1.9)
# Status  : REFERENCE IMPLEMENTATION from Methods
# =============================================================================
set -euo pipefail
IN=${1:?PLINK bfile prefix exported from GenomeStudio}
OUT=${2:?output prefix}

plink --bfile "${IN}" --freq    --out "${OUT}.stats"
plink --bfile "${IN}" --missing --out "${OUT}.stats"

# .frq: CHR SNP A1 A2 MAF NCHROBS ; .lmiss: CHR SNP N_MISS N_GENO F_MISS
awk 'NR==FNR { if (FNR>1) maf[$2]=$5; next }
     FNR>1 {
       cr = 1 - $5; m = maf[$2]
       if ((m >  0.05 && cr < 0.95) || (m <= 0.05 && cr < 0.99)) print $2
     }' "${OUT}.stats.frq" "${OUT}.stats.lmiss" > "${OUT}.exclude_snps.txt"

echo "Excluding $(wc -l < "${OUT}.exclude_snps.txt") variants"
plink --bfile "${IN}" --exclude "${OUT}.exclude_snps.txt" --mind 0.05 \
      --make-bed --out "${OUT}"
# Next: convert to per-chromosome VCF (hg38) and upload to the Michigan
# Imputation Server v1.7.1 (TOPMed panel, minimac4 v1.0.2, Eagle v2.4.1 phasing).
