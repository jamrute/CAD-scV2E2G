#!/usr/bin/env Rscript
# =============================================================================
# 02_qtl_discovery/08_pme_caqtl.R
# Purpose : Single-cell Poisson mixed-effects (PME) caQTL model in SMCs,
#           adapted from Nathan et al. (single-cell eQTL), fit with lme4::glmer.
#
#   Static model (per lead variant–peak pair from pseudobulk RASQUAL):
#     log E[A] = θ + βG·G + βage·age + βsex·sex + βnFrags·log(nFrags)
#                + βTSS·TSSEnrich + Σk βaccPCk·accPCk + Σk βgenoPCk·genoPCk
#                + φ_donor + ε_site
#     Test βG by LRT (full vs. model without G).
#
#   Dynamic model (--dynamic): add FMC score and G×FMC interaction:
#     ... + βFMC·FMC + βG×FMC·(G×FMC)
#     Test βG×FMC by LRT (full vs. model without the interaction),
#     then Storey q-value across all tested variants (q <= 0.05).
#
#   All quantitative covariates are centred and scaled. Genotype is kept as
#   0/1/2 dosage (TODO: confirm whether G was also scaled in the paper).
#   Static-model results are not FDR-corrected (see Methods); Z-values are
#   compared with pseudobulk Z instead.
#
# Paper   : ED Fig. 6a; Fig. 2i–k; 712 dynamic caQTLs (q <= 0.05); ST 29–30
#
# Inputs  : --counts    RDS dgCMatrix, peaks x SMC nuclei (meta-map), ATAC fragment counts
#           --cell-meta TSV: cell, donor, site, nFrags, TSSEnrichment, accPC1..4, FMC_score
#           --donor-meta TSV: donor, age, sex, genoPC1..4
#           --dosage    TSV: variant_id, <donor1>, <donor2>, ...  (0/1/2)
#           --pairs     TSV: peak_id, variant_id  (lead variants, pseudobulk FDR < 5%)
# Outputs : --out TSV, one row per pair
#
# Usage   : Rscript 08_pme_caqtl.R --counts ... --pairs ... --out static.tsv
#           Rscript 08_pme_caqtl.R ... --dynamic --out dynamic.tsv
#           Array jobs: add --chunk i --n-chunks N, then 08b merge + q-values.
#
# Status  : REFERENCE IMPLEMENTATION from Methods
# =============================================================================
suppressPackageStartupMessages({
  library(lme4); library(Matrix); library(data.table)
  library(optparse); library(parallel)
})
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config(); pp <- cfg$params$pme

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts"), make_option("--cell-meta"), make_option("--donor-meta"),
  make_option("--dosage"), make_option("--pairs"), make_option("--out"),
  make_option("--dynamic", action = "store_true", default = FALSE),
  make_option("--chunk", type = "integer", default = 1L),
  make_option("--n-chunks", type = "integer", default = 1L),
  make_option("--cores", type = "integer", default = 1L),
  make_option("--nAGQ", type = "integer", default = 1L,
              help = "1 = Laplace (default); 0 = faster, less accurate")
)))

# ---- Load & assemble per-cell design ---------------------------------------
counts <- readRDS(opt$counts)
cm <- fread(opt$`cell-meta`)
dm <- fread(opt$`donor-meta`)
dos <- fread(opt$dosage)
pairs <- fread(opt$pairs)

pairs <- pairs[(seq_len(.N) - 1L) %% opt$`n-chunks` == (opt$chunk - 1L)]
log_msg("Chunk ", opt$chunk, "/", opt$`n-chunks`, ": ", nrow(pairs), " pairs; dynamic = ", opt$dynamic)

cm <- cm[cell %in% colnames(counts)]
cm <- merge(cm, dm, by = "donor", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(cm$age))
cm <- cm[match(colnames(counts), cell)]
cm <- cm[!is.na(cell)]
counts <- counts[, cm$cell, drop = FALSE]

acc_pcs  <- paste0("accPC",  seq_len(pp$n_acc_pcs))
geno_pcs <- paste0("genoPC", seq_len(pp$n_geno_pcs))

z <- function(x) as.numeric(scale(x))
design <- data.table(
  donor = factor(cm$donor),
  site  = factor(cm$site),
  age   = z(cm$age),
  sex   = factor(cm$sex),
  log_nFrags = z(log(cm$nFrags)),
  TSSEnrich  = z(cm$TSSEnrichment)
)
for (v in c(acc_pcs, geno_pcs)) set(design, j = v, value = z(cm[[v]]))
if (opt$dynamic) design[, FMC := z(cm$FMC_score)]

dos_mat <- as.matrix(dos[, -1, with = FALSE]); rownames(dos_mat) <- dos$variant_id
donor_idx <- match(as.character(design$donor), colnames(dos_mat))
if (anyNA(donor_idx)) stop("Donors in cell-meta missing from dosage file")

# ---- Formulae --------------------------------------------------------------
covars <- c("age", "sex", "log_nFrags", "TSSEnrich", acc_pcs, geno_pcs)
re     <- "(1 | donor) + (1 | site)"
if (!opt$dynamic) {
  f_full <- reformulate(c("G", covars, re), response = "y")
  f_null <- reformulate(c(covars, re),      response = "y")
  test_term <- "G"
} else {
  f_full <- reformulate(c("G", covars, "FMC", "G:FMC", re), response = "y")
  f_null <- reformulate(c("G", covars, "FMC", re),          response = "y")
  test_term <- "G:FMC"
}

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
                     calc.derivs = FALSE)

fit_pair <- function(i) {
  pk <- pairs$peak_id[i]; vr <- pairs$variant_id[i]
  res <- list(peak_id = pk, variant_id = vr, n_cells = NA_integer_, beta_G = NA_real_,
              se_G = NA_real_, z_G = NA_real_, beta_FMC = NA_real_, beta_GxFMC = NA_real_,
              se_GxFMC = NA_real_, z_GxFMC = NA_real_, lrt_chisq = NA_real_,
              lrt_p = NA_real_, converged = NA, error = NA_character_)
  tryCatch({
    if (!pk %in% rownames(counts)) stop("peak not in counts")
    if (!vr %in% rownames(dos_mat)) stop("variant not in dosage")
    d <- copy(design)
    d[, y := as.numeric(counts[pk, ])]
    d[, G := as.numeric(dos_mat[vr, donor_idx])]
    d <- d[!is.na(G)]
    res$n_cells <- nrow(d)
    m1 <- glmer(f_full, data = d, family = poisson, control = ctrl, nAGQ = opt$nAGQ)
    m0 <- glmer(f_null, data = d, family = poisson, control = ctrl, nAGQ = opt$nAGQ)
    cf <- summary(m1)$coefficients
    res$beta_G <- cf["G", "Estimate"]; res$se_G <- cf["G", "Std. Error"]
    res$z_G    <- cf["G", "z value"]
    if (opt$dynamic) {
      res$beta_FMC   <- cf["FMC", "Estimate"]
      res$beta_GxFMC <- cf["G:FMC", "Estimate"]
      res$se_GxFMC   <- cf["G:FMC", "Std. Error"]
      res$z_GxFMC    <- cf["G:FMC", "z value"]
    }
    a <- anova(m0, m1, test = "Chisq")
    res$lrt_chisq <- a$Chisq[2]
    res$lrt_p     <- a$`Pr(>Chisq)`[2]
    res$converged <- is.null(m1@optinfo$conv$lme4$messages)
    res
  }, error = function(e) { res$error <- conditionMessage(e); res })
}

out <- rbindlist(mclapply(seq_len(nrow(pairs)), fit_pair, mc.cores = opt$cores))
out[, `:=`(model = if (opt$dynamic) "dynamic" else "static", tested_term = test_term)]
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(out, opt$out, sep = "\t")
log_msg("Done: ", sum(is.na(out$error)), " fitted, ", sum(!is.na(out$error)), " failed")
