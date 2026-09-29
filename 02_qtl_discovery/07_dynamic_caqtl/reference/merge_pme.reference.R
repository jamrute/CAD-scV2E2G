#!/usr/bin/env Rscript
# =============================================================================
# 02_qtl_discovery/08b_merge_pme.R
# Purpose : Merge PME chunks. Dynamic mode: Storey q-values on interaction
#           LRT p-values (q <= 0.05 -> dynamic caQTL). Static mode: join with
#           pseudobulk RASQUAL results and compare Z-values.
# Paper   : ST 29–30; 712 dynamic caQTLs
# Status  : REFERENCE IMPLEMENTATION from Methods
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(qvalue); library(optparse) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config()

opt <- parse_args(OptionParser(option_list = list(
  make_option("--chunks-glob"), make_option("--out"),
  make_option("--rasqual", default = NA, help = "caqtl_leads_fdr.tsv (static mode)")
)))
x <- rbindlist(lapply(Sys.glob(opt$`chunks-glob`), fread))
mode <- unique(x$model); stopifnot(length(mode) == 1)

if (mode == "dynamic") {
  ok <- !is.na(x$lrt_p)
  storey_q <- function(p) tryCatch(qvalue(p)$qvalues, error = function(e) {
    warning("qvalue pi0 estimation failed (", conditionMessage(e), "); using pi0 = 1")
    qvalue(p, pi0 = 1)$qvalues
  })
  x[ok, q_interaction := storey_q(lrt_p)]
  x[, dynamic_caqtl := !is.na(q_interaction) & q_interaction <= cfg$params$pme$dynamic_qvalue]
  log_msg("Dynamic caQTLs (q <= ", cfg$params$pme$dynamic_qvalue, "): ", sum(x$dynamic_caqtl))
} else if (!is.na(opt$rasqual)) {
  rq <- fread(opt$rasqual)[, .(peak_id = feature, variant_id = rsid, pi, chisq)]
  # RASQUAL pi is centred at 0.5; sign(pi - 0.5) * sqrt(chisq) gives a signed Z
  rq[, z_pseudobulk := sign(pi - 0.5) * sqrt(chisq)]
  x <- merge(x, rq, by = c("peak_id", "variant_id"), all.x = TRUE)
  log_msg("Spearman(z_PME, z_pseudobulk) = ",
          round(cor(x$z_G, x$z_pseudobulk, method = "spearman", use = "complete.obs"), 3))
}
fwrite(x, opt$out, sep = "\t")
