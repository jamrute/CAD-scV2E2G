#!/usr/bin/env Rscript
# =============================================================================
# 02_qtl_discovery/06_rasqual_fdr.R
# Purpose : Two-step multiple-testing correction (RASQUAL recommendation):
#   1. Locus level: BH q-value across SNPs within each peak; lead SNP = min q
#      (RASQUAL -t output already reports lead SNP with log10 BH q).
#   2. Genome-wide: empirical null from 4 genotype-permuted RASQUAL runs,
#      averaged; q-value cut-offs for 1%, 5%, 10% genome-wide FDR.
# Paper   : ED Fig. 3g–h; 11,269 caQTLs at 10% FDR; Supp. Tables 15–26
# Inputs  : results/qtl/rasqual/<CellType>/{real,perm1..perm4}/chunk*.txt
# Outputs : results/qtl/rasqual/<CellType>/caqtl_leads_fdr.tsv, fdr_cutoffs.tsv
# Status  : REFERENCE IMPLEMENTATION from Methods — verify column order below
#           against the RASQUAL v1.1 README.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config(); pq <- cfg$params$caqtl

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cell-type", type = "character")
)))
ct  <- opt$`cell-type`
dir <- file.path(cfg_path(cfg, "qtl", "out"), "rasqual", ct)

RASQUAL_COLS <- c(
  "feature", "rsid", "chr", "pos", "ref", "alt", "af", "hwe_chisq", "imp_quality",
  "log10_q", "chisq", "pi", "delta", "phi", "overdisp", "snp_idx", "n_fsnp",
  "n_rsnp", "iter_null", "iter_alt", "ties", "loglik_null", "converged",
  "r2_fsnp", "r2_rsnp"
)

read_rasqual <- function(d) {
  files <- list.files(d, pattern = "^chunk.*\\.txt$", full.names = TRUE)
  if (!length(files)) stop("No RASQUAL output in ", d)
  x <- rbindlist(lapply(files, fread, header = FALSE))
  setnames(x, seq_along(RASQUAL_COLS), RASQUAL_COLS)
  x[, p := pchisq(chisq, df = 1, lower.tail = FALSE)]
  x[, q_locus := 10^log10_q]
  # one lead per feature (guard against duplicates from re-runs)
  x[order(q_locus, p)][!duplicated(feature)]
}

real  <- read_rasqual(file.path(dir, "real"))
perms <- lapply(seq_len(pq$n_permutations), function(i) read_rasqual(file.path(dir, paste0("perm", i))))

#' Genome-wide FDR for a locus-level q threshold t:
#'   FDR(t) = mean_over_perms( #{q_perm <= t} ) / #{q_real <= t}
#' Returns, for each real lead, the minimum FDR achievable at thresholds >= its q
#' (monotone, like a q-value).
empirical_fdr <- function(q_real, q_perm_list) {
  thr <- sort(unique(q_real))
  n_real <- findInterval(thr, sort(q_real))
  n_null <- rowMeans(sapply(q_perm_list, function(qp) findInterval(thr, sort(qp))))
  fdr <- pmin(1, n_null / pmax(n_real, 1))
  fdr <- rev(cummin(rev(fdr)))                         # enforce monotonicity
  fdr[match(q_real, thr)]
}

real[, fdr_genomewide := empirical_fdr(q_locus, lapply(perms, `[[`, "q_locus"))]

cutoffs <- rbindlist(lapply(pq$fdr_levels, function(a) {
  sig <- real[fdr_genomewide <= a]
  data.table(fdr = a,
             q_locus_cutoff = if (nrow(sig)) max(sig$q_locus) else NA_real_,
             n_caqtl = nrow(sig))
}))
cutoffs[, cell_type := ct]
print(cutoffs)

fwrite(real[order(fdr_genomewide, p)], file.path(dir, "caqtl_leads_fdr.tsv"), sep = "\t")
fwrite(cutoffs, file.path(dir, "fdr_cutoffs.tsv"), sep = "\t")
save_session_info(dir, "06_rasqual_fdr.R")
