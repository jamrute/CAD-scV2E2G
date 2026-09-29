#!/usr/bin/env bash
# =============================================================================
# 02_qtl_discovery/01_genotypes/03_merge_turner.sh
# Purpose : Merge WashU imputed genotypes with Turner et al. (low-coverage WGS)
#           keeping only shared SNPs -> 5,229,397 variants (Methods)
# Status  : REFERENCE IMPLEMENTATION — confirm allele harmonisation used
# =============================================================================
set -euo pipefail
WASHU=${1:?filtered WashU vcf.gz}
TURNER=${2:?Turner vcf.gz}
OUT=${3:?merged vcf.gz}

tmp=$(mktemp -d)
# Shared sites (same CHROM/POS/REF/ALT)
bcftools isec -n=2 -c none -w1 "${WASHU}" "${TURNER}" -Oz -o "${tmp}/washu.shared.vcf.gz"
bcftools isec -n=2 -c none -w2 "${WASHU}" "${TURNER}" -Oz -o "${tmp}/turner.shared.vcf.gz"
tabix -p vcf "${tmp}/washu.shared.vcf.gz"; tabix -p vcf "${tmp}/turner.shared.vcf.gz"
bcftools merge "${tmp}/washu.shared.vcf.gz" "${tmp}/turner.shared.vcf.gz" -Oz -o "${OUT}"
tabix -p vcf "${OUT}"
echo "Merged variants: $(bcftools index -n "${OUT}")"
rm -rf "${tmp}"
