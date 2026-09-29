#!/usr/bin/env Rscript
# =============================================================================
# 02_qtl_discovery/10_scV2E2G_map.R
# Purpose : Build the disease-relevant scV2E2G map (Methods, "scE2G model predictions"):
#   1. ATAC peaks (4 major cell types) containing a CAD GWAS variant   (n = 1,360)
#   2. Cell-type-specific enhancers = those peaks overlapping coronary
#      artery H3K27ac peaks          (SMC 883, EC 729, Fib 784, Mye 760)
#   3. scE2G-linked genes            (SMC 325, EC 295, Fib 295, Mye 325)
#   4. Links also supported by a caQTL in the same cell type
#                                    (SMC 56, EC 8, Fib 13, Mye 27 genes)
# Paper   : Fig. 3a; ST 33–36
# Inputs  : results/multiome/cad_variant_peaks.tsv       (01_multiome/09)
#           per-cell-type peak sets (BED)                 (--celltype-peaks-dir)
#           coronary artery H3K27ac peaks (ENCODE, BED)   (--h3k27ac)
#           results/multiome/scE2G_summary/scE2G_links_thresholded.tsv.gz
#           results/qtl/rasqual/<CT>/caqtl_leads_fdr.tsv
# Status  : REFERENCE IMPLEMENTATION — verify the exact overlap rules used
#           (e.g. whether caQTL support requires the CAD variant itself to be
#           the caQTL, or any caQTL on the same peak) against the counts above.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(optparse) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config()

opt <- parse_args(OptionParser(option_list = list(
  make_option("--celltype-peaks-dir", help = "dir with <CellType>.bed"),
  make_option("--h3k27ac", help = "coronary artery H3K27ac peaks BED"),
  make_option("--caqtl-fdr", type = "double", default = 0.10),
  make_option("--caqtl-variant-must-be-cad", action = "store_true", default = FALSE)
)))

cts <- cfg$params$caqtl$cell_types
mo  <- cfg_path(cfg, "multiome", "out")
out <- file.path(cfg_path(cfg, "qtl", "out"), "scV2E2G"); dir.create(out, recursive = TRUE, showWarnings = FALSE)

pk_gr <- function(dt, chr = "chr", s = "start", e = "end")
  GRanges(dt[[chr]], IRanges(dt[[s]], dt[[e]]))

cadp  <- fread(file.path(mo, "cad_variant_peaks.tsv"))
cadv  <- unique(cadp[, .(rsid, peak_chr, peak_start, peak_end)])
h3    <- rtracklayer::import(opt$h3k27ac)
links <- fread(file.path(mo, "scE2G_summary", "scE2G_links_thresholded.tsv.gz"))

summary_rows <- list(); map_rows <- list()
for (ct in cts) {
  ctp <- rtracklayer::import(file.path(opt$`celltype-peaks-dir`, paste0(ct, ".bed")))
  v_gr <- GRanges(cadv$peak_chr, IRanges(cadv$peak_start, cadv$peak_end), rsid = cadv$rsid)

  # (1) CAD-variant peaks present in this cell type
  in_ct <- v_gr[overlapsAny(v_gr, ctp)]
  # (2) ... that are H3K27ac-marked enhancers
  enh <- in_ct[overlapsAny(in_ct, h3)]
  # (3) scE2G links from these enhancers
  lk  <- links[CellType == ct]
  lk_gr <- pk_gr(lk)
  ov <- findOverlaps(enh, lk_gr)
  m  <- data.table(cell_type = rep(ct, length(ov)), rsid = enh$rsid[queryHits(ov)],
                   enh_chr = as.character(seqnames(enh))[queryHits(ov)],
                   enh_start = start(enh)[queryHits(ov)], enh_end = end(enh)[queryHits(ov)],
                   gene = lk$TargetGene[subjectHits(ov)], scE2G_score = lk$Score[subjectHits(ov)])
  # (4) caQTL support
  cq <- fread(file.path(cfg_path(cfg, "qtl", "out"), "rasqual", ct, "caqtl_leads_fdr.tsv"))
  cq <- cq[fdr_genomewide <= opt$`caqtl-fdr`]
  # feature IDs assumed "chr:start-end"
  cq[, `:=`(fchr   = sub(":.*$", "", feature),
            fstart = as.integer(sub("^.*:([0-9]+)-.*$", "\\1", feature)),
            fend   = as.integer(sub("^.*-", "", feature)))]
  cq_gr <- GRanges(cq$fchr, IRanges(cq$fstart, cq$fend), rsid = cq$rsid)
  if (nrow(m)) {
    m_gr <- GRanges(m$enh_chr, IRanges(m$enh_start, m$enh_end))
    o2 <- findOverlaps(m_gr, cq_gr)
    sup <- data.table(i = queryHits(o2), caqtl_rsid = cq_gr$rsid[subjectHits(o2)])
    if (opt$`caqtl-variant-must-be-cad`) sup <- sup[caqtl_rsid == m$rsid[i]]
    sup_agg <- sup[, .(caqtl_rsids = paste(unique(caqtl_rsid), collapse = ",")), by = i]
    m[, `:=`(row_i = .I, caqtl_rsids = NA_character_)]
    m[sup_agg, caqtl_rsids := i.caqtl_rsids, on = .(row_i = i)]
    m[, `:=`(caqtl_supported = !is.na(caqtl_rsids), row_i = NULL)]
  } else {
    m[, `:=`(caqtl_rsids = character(), caqtl_supported = logical())]
  }

  map_rows[[ct]] <- m
  summary_rows[[ct]] <- data.table(
    cell_type = ct,
    n_cad_peaks_in_celltype = length(unique(in_ct)),
    n_cad_enhancers = length(unique(enh)),
    n_scE2G_genes = uniqueN(m$gene),
    n_caqtl_supported_genes = uniqueN(m[caqtl_supported == TRUE, gene])
  )
}
summ <- rbindlist(summary_rows)
summ[, n_cad_peaks_total := nrow(unique(cadv[, .(peak_chr, peak_start, peak_end)]))]
print(summ)
fwrite(summ, file.path(out, "Fig3a_scV2E2G_counts.tsv"), sep = "\t")
fwrite(rbindlist(map_rows, fill = TRUE), file.path(out, "scV2E2G_map.tsv"), sep = "\t")
