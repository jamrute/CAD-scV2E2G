#!/usr/bin/env Rscript
# =============================================================================
# 01_multiome/06_scE2G/02_summarise_scE2G.R
# Purpose : Threshold scE2G predictions; shared-enhancer heatmap; link counts
# Paper   : Fig. 1d (shared enhancers), Fig. 1e (links per cell type); ST 3–14
# Status  : REFERENCE IMPLEMENTATION from Methods / Fig. 1 legend — verify
#           column names against your sc-E2G output before use
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(pheatmap) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config()
thr <- cfg$params$scE2G$score_threshold
pred_dir <- file.path(cfg_path(cfg, "multiome", "out"), "scE2G_predictions")
out <- file.path(cfg_path(cfg, "multiome", "out"), "scE2G_summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Expected: one file per cell type with columns chr,start,end,TargetGene,Score,isPromoterElement
files <- list.files(pred_dir, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
names(files) <- sub("\\.tsv(\\.gz)?$", "", basename(files))

links <- rbindlist(lapply(names(files), function(ct) {
  x <- fread(files[[ct]])
  x <- x[Score > thr & !as.logical(isPromoterElement)]   # ED Fig. 2 omits promoters
  x[, CellType := ct]
}))
fwrite(links, file.path(out, "scE2G_links_thresholded.tsv.gz"), sep = "\t")

# Fig. 1e: number of links per cell type
fwrite(links[, .N, by = CellType], file.path(out, "Fig1e_links_per_celltype.tsv"), sep = "\t")

# Fig. 1d: fraction of enhancers in A with >50% of their width overlapping enhancers in B
enh <- lapply(split(links, by = "CellType"), function(x)
  reduce(GRanges(unique(x[, .(chr, start, end)]))))
cts <- names(enh)
shared <- matrix(NA_real_, length(cts), length(cts), dimnames = list(cts, cts))
for (a in cts) for (b in cts) {
  ov <- findOverlaps(enh[[a]], enh[[b]])
  if (length(ov) == 0) { shared[a, b] <- 0; next }
  w  <- width(pintersect(enh[[a]][queryHits(ov)], enh[[b]][subjectHits(ov)]))
  ov_bp <- tapply(w, queryHits(ov), sum)
  frac  <- ov_bp / width(enh[[a]])[as.integer(names(ov_bp))]
  shared[a, b] <- sum(frac > 0.5) / length(enh[[a]])
}
fwrite(as.data.table(shared, keep.rownames = "CellType"),
       file.path(out, "Fig1d_shared_enhancers.tsv"), sep = "\t")
pheatmap(shared, filename = file.path(out, "Fig1d_shared_enhancers.pdf"),
         width = 4, height = 4)
