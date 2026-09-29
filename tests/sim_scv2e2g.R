library(data.table)
# 4 CAD-variant peaks; peak1+2 in SMC, peak2 has no H3K27ac; peak3 in EC
fwrite(data.table(rsid = c("rs1","rs2","rs3","rs4","rs1b"),
                  peak_chr = "chr1", peak_start = c(1001, 5001, 9001, 20001, 1001),
                  peak_end = c(1500, 5500, 9500, 20500, 1500), peak_group = "x"),
       "v2g/multiome/cad_variant_peaks.tsv", sep = "\t")
bed <- function(s, e, f) fwrite(data.table("chr1", s - 1, e), f, sep = "\t", col.names = FALSE)
bed(c(1001, 5001), c(1500, 5500), "v2g/ctpeaks/SMC.bed")
bed(c(9001), c(9500), "v2g/ctpeaks/Endothelium.bed")
bed(c(50000), c(50100), "v2g/ctpeaks/Myeloid.bed")
bed(c(50000), c(50100), "v2g/ctpeaks/Fibroblast.bed")
bed(c(1100, 9100), c(1300, 9300), "v2g/h3k27ac.bed")
fwrite(data.table(chr = "chr1", start = c(1001, 1001, 9001), end = c(1500, 1500, 9500),
                  TargetGene = c("GENEA", "GENEB", "GENEC"), Score = c(.5, .3, .9),
                  isPromoterElement = FALSE, CellType = c("SMC", "SMC", "Endothelium")),
       "v2g/multiome/scE2G_summary/scE2G_links_thresholded.tsv.gz", sep = "\t")
for (ct in c("SMC","Endothelium","Myeloid","Fibroblast")) {
  d <- file.path("v2g/qtl/rasqual", ct); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  x <- if (ct == "SMC") data.table(feature = "chr1:1001-1500", rsid = "rs_other", fdr_genomewide = 0.01)
       else data.table(feature = "chr1:70000-70100", rsid = "rsZ", fdr_genomewide = 0.5)
  fwrite(x, file.path(d, "caqtl_leads_fdr.tsv"), sep = "\t")
}
writeLines(c(paste0("root: ", normalizePath("v2g")), "multiome:", "  out: multiome", "qtl:", "  out: qtl"), "v2g/paths.yaml")
