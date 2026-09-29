#!/usr/bin/env bash
# =============================================================================
# 02_qtl_discovery/01_genotypes/02_post_imputation_filter.sh
# Purpose : Filter TOPMed-imputed VCFs for caQTL mapping
#           Methods: final set = 7,250,405 variants.
#           NOTE: Methods text reads "variants with imputation R2 < 0.3 and
#           MAF > 0.05 were further filtered"; implemented here as KEEP
#           R2 >= 0.3 AND MAF > 0.05 (ED Fig. 3a legend). Confirm.
# Status  : REFERENCE IMPLEMENTATION
# =============================================================================
set -euo pipefail
IMPUTED_DIR=${1:?dir with chr*.dose.vcf.gz}
OUT=${2:?output vcf.gz}

tmp=$(mktemp -d)
for chr in $(seq 1 22); do
  bcftools view -i 'INFO/R2>=0.3 && INFO/MAF>0.05' -m2 -M2 -v snps \
    "${IMPUTED_DIR}/chr${chr}.dose.vcf.gz" -Oz -o "${tmp}/chr${chr}.vcf.gz"
done
bcftools concat "${tmp}"/chr{1..22}.vcf.gz -Oz -o "${OUT}"
tabix -p vcf "${OUT}"
echo "Variants retained: $(bcftools index -n "${OUT}")"
rm -rf "${tmp}"
