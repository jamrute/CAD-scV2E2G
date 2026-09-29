#!/usr/bin/env Rscript
# =============================================================================
# 03_perturbseq/crispri_validation/01_qpcr_ddct.R
# Purpose : Relative expression (2^-ΔΔCt, reference gene UBC) for CRISPRi
#           enhancer/TSS knockdowns vs non-targeting control
# Paper   : Fig. 3c (MYO9B), Fig. 5c (AMOTL2), ED 7c–e (LOXL1, NAV1, PALLD)
# Inputs  : TSV: experiment, cell_type, condition, replicate, gene, ct
# Status  : REFERENCE IMPLEMENTATION — confirm test choice per panel
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--ct-table"), make_option("--out"),
  make_option("--ref-gene", default = "UBC"), make_option("--control", default = "NTC")
)))
ct <- fread(opt$`ct-table`)
ct <- ct[, .(ct = mean(ct)), by = .(experiment, cell_type, condition, replicate, gene)]
ref <- ct[gene == opt$`ref-gene`, .(experiment, cell_type, condition, replicate, ct_ref = ct)]
d <- merge(ct[gene != opt$`ref-gene`], ref, by = c("experiment", "cell_type", "condition", "replicate"))
d[, dct := ct - ct_ref]
d[, ddct := dct - mean(dct[condition == opt$control]), by = .(experiment, cell_type, gene)]
d[, rel_expr := 2^(-ddct)]

tests <- d[, {
  conds <- setdiff(unique(condition), opt$control)
  if (length(conds) == 1) {
    tt <- t.test(rel_expr[condition == conds], rel_expr[condition == opt$control])
    .(comparison = paste(conds, "vs", opt$control), test = "two-sided t-test", p = tt$p.value)
  } else {
    if (!requireNamespace("multcomp", quietly = TRUE)) stop("install multcomp for Dunnett")
    f <- factor(condition, levels = c(opt$control, conds))
    fit <- aov(rel_expr ~ f)
    dn <- summary(multcomp::glht(fit, linfct = multcomp::mcp(f = "Dunnett")))
    .(comparison = names(dn$test$coefficients), test = "ANOVA + Dunnett",
      p = as.numeric(dn$test$pvalues))
  }
}, by = .(experiment, cell_type, gene)]

fwrite(d, sub("\\.tsv$", "_values.tsv", opt$out), sep = "\t")
fwrite(tests, opt$out, sep = "\t")
