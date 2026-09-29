#!/usr/bin/env Rscript
# =============================================================================
# 04_hic/04_mosaic_annotation.R
# Purpose : Intersect 2D mosaic anchors with 1D features (protein-coding gene
#           TSS, snATAC SMC peaks harbouring CAD variants); summarise mosaics;
#           compare loop counts in CAD vs non-CAD mosaics (Mann–Whitney U).
# Paper   : Fig. 4b (25 vs 14 loops, P = 1.7e-17), Fig. 4c; ST 46–49
# Inputs  : results/hic/aqua/cohort_mosaics.bedpe  (chr1 s1 e1 chr2 s2 e2 mosaic_id ...)
#           GENCODE v38 GTF; SMC CAD-variant peaks (from 02_qtl_discovery/10)
# Status  : REFERENCE IMPLEMENTATION — confirm feature definitions (TSS vs gene body)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(rtracklayer) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config()
aqua <- cfg_path(cfg, "hic", "aqua_dir"); out <- cfg_path(cfg, "hic", "out", create_dir = TRUE)

bp <- fread(file.path(aqua, "cohort_mosaics.bedpe"), header = FALSE)
setnames(bp, 1:7, c("chr1", "s1", "e1", "chr2", "s2", "e2", "mosaic"))
bp[, loop_id := .I]
anchors <- rbind(bp[, .(loop_id, mosaic, chr = chr1, start = s1 + 1L, end = e1)],
                 bp[, .(loop_id, mosaic, chr = chr2, start = s2 + 1L, end = e2)])
a_gr <- GRanges(anchors$chr, IRanges(anchors$start, anchors$end))

gtf <- import(cfg_path(cfg, "refs", "gencode_gtf"))
genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
tss <- resize(genes, 1, fix = "start")

smc <- fread(file.path(cfg_path(cfg, "qtl", "out"), "scV2E2G", "scV2E2G_map.tsv"))[cell_type == "SMC"]
smc_gr <- unique(GRanges(smc$enh_chr, IRanges(smc$enh_start, smc$enh_end)))
cadv <- fread(file.path(cfg_path(cfg, "multiome", "out"), "cad_variant_peaks.tsv"))
cad_gr <- GRanges(cadv$peak_chr, IRanges(cadv$peak_start, cadv$peak_end), rsid = cadv$rsid)

hit <- function(q, s) as.data.table(findOverlaps(q, s))
g <- hit(a_gr, tss);    g[, `:=`(mosaic = anchors$mosaic[queryHits], gene = tss$gene_name[subjectHits])]
e <- hit(a_gr, smc_gr); e[, `:=`(mosaic = anchors$mosaic[queryHits], enh = subjectHits)]
v <- hit(a_gr, cad_gr); v[, `:=`(mosaic = anchors$mosaic[queryHits], rsid = cad_gr$rsid[subjectHits])]

summ <- bp[, .(n_loops = .N), by = mosaic]
summ <- merge(summ, g[, .(n_genes = uniqueN(gene), genes = paste(sort(unique(gene)), collapse = ",")), by = mosaic], all.x = TRUE)
summ <- merge(summ, e[, .(n_smc_cad_enhancers = uniqueN(enh)), by = mosaic], all.x = TRUE)
summ <- merge(summ, v[, .(n_cad_variants = uniqueN(rsid), cad_variants = paste(sort(unique(rsid)), collapse = ",")), by = mosaic], all.x = TRUE)
for (col in c("n_genes", "n_smc_cad_enhancers", "n_cad_variants")) set(summ, which(is.na(summ[[col]])), col, 0L)
summ[, has_cad := n_cad_variants > 0]

log_msg("Mosaics: ", nrow(summ), "; containing protein-coding genes: ", uniqueN(g$gene))
wt <- wilcox.test(n_loops ~ has_cad, data = summ)
fwrite(data.table(mean_loops_cad = summ[has_cad == TRUE, mean(n_loops)],
                  mean_loops_noncad = summ[has_cad == FALSE, mean(n_loops)],
                  p_mann_whitney = wt$p.value), file.path(out, "Fig4b_loops_cad_vs_noncad.tsv"), sep = "\t")
fwrite(summ[order(-n_smc_cad_enhancers, -n_genes)], file.path(out, "Fig4c_mosaic_summary.tsv"), sep = "\t")
