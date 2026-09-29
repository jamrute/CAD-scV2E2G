#!/usr/bin/env bash
# =============================================================================
# tests/run_tests.sh — simulation tests for the reference implementations.
# Requires R with data.table, optparse, yaml, lme4, qvalue, GenomicRanges,
# rtracklayer, R.utils.   Usage: bash tests/run_tests.sh
# =============================================================================
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
export CAD_REPO=${REPO}
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
cd "${T}"
pass() { echo "PASS  $1"; }
fail() { echo "FAIL  $1"; exit 1; }

# 1. PME static + dynamic caQTL: recovers simulated effects
Rscript "${REPO}/tests/sim_pme.R"
Rscript "${REPO}/02_qtl_discovery/07_dynamic_caqtl/reference/pme_caqtl.reference.R" --counts counts.rds --cell-meta cell_meta.tsv \
  --donor-meta donor_meta.tsv --dosage dosage.tsv --pairs pairs.tsv --out out/static.tsv 2>/dev/null
Rscript "${REPO}/02_qtl_discovery/07_dynamic_caqtl/reference/pme_caqtl.reference.R" --counts counts.rds --cell-meta cell_meta.tsv \
  --donor-meta donor_meta.tsv --dosage dosage.tsv --pairs pairs.tsv --out out/dyn.chunk1.tsv --dynamic 2>/dev/null
Rscript "${REPO}/02_qtl_discovery/07_dynamic_caqtl/reference/merge_pme.reference.R" --chunks-glob 'out/dyn.chunk*.tsv' --out out/dynamic.tsv 2>/dev/null
Rscript -e '
library(data.table)
s <- fread("out/static.tsv"); d <- fread("out/dynamic.tsv")
stopifnot(abs(s[variant_id == "rs_static", beta_G] - 0.4) < 0.1,
          s[variant_id == "rs_static", lrt_p] < 1e-4,
          abs(d[variant_id == "rs_dynamic", beta_GxFMC] + 0.35) < 0.1,
          d[dynamic_caqtl == TRUE, variant_id] == "rs_dynamic")' && pass "PME caQTL" || fail "PME caQTL"

# 2. RASQUAL genome-wide empirical FDR is calibrated
Rscript "${REPO}/tests/sim_rasqual.R"
CAD_CONFIG=${T}/paths_test.yaml Rscript "${REPO}/02_qtl_discovery/06_multiple_testing/reference/rasqual_fdr.reference.R" --cell-type SMC >/dev/null
Rscript -e '
library(data.table); x <- fread("rq/qtl/rasqual/SMC/caqtl_leads_fdr.tsv")
x[, true := as.integer(sub("rs", "", rsid)) <= 400]
s <- x[fdr_genomewide <= 0.05]; stopifnot(nrow(s) >= 380, mean(!s$true) < 0.10)' \
  && pass "RASQUAL FDR" || fail "RASQUAL FDR"

# 3. CTCF circular-permutation enrichment detects planted enrichment
Rscript "${REPO}/tests/sim_ctcf.R"
Rscript "${REPO}/04_hic/06_ctcf/ctcf_cad_enrichment.reference.R" --variants cad_vars.tsv --peaks ctcf.bed \
  --chrom-sizes chrom.sizes --n-perm 200 --out ctcf/enrich.tsv >/dev/null 2>&1
Rscript -e 'x <- data.table::fread("ctcf/enrich.tsv"); stopifnot(x$fold_enrichment > 1.3, x$p_empirical < 0.01)' \
  && pass "CTCF enrichment" || fail "CTCF enrichment"

# 4. scV2E2G map counts on a toy locus set
mkdir -p v2g/multiome/scE2G_summary v2g/ctpeaks
Rscript "${REPO}/tests/sim_scv2e2g.R"
CAD_CONFIG=${T}/v2g/paths.yaml Rscript "${REPO}/02_qtl_discovery/08_scV2E2G_map/scV2E2G_map.reference.R" \
  --celltype-peaks-dir v2g/ctpeaks --h3k27ac v2g/h3k27ac.bed >/dev/null 2>&1
Rscript -e '
x <- data.table::fread("v2g/qtl/scV2E2G/Fig3a_scV2E2G_counts.tsv")
stopifnot(x[cell_type == "SMC", c(n_cad_peaks_in_celltype, n_cad_enhancers, n_scE2G_genes, n_caqtl_supported_genes)] == c(2, 1, 2, 2),
          x[cell_type == "Endothelium", n_caqtl_supported_genes] == 0)' \
  && pass "scV2E2G map" || fail "scV2E2G map"

echo "All tests passed."
