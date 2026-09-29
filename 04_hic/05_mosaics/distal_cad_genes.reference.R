#!/usr/bin/env Rscript
# =============================================================================
# 04_hic/06_distal_cad_genes.R
# Purpose : Link CAD variants to genes by direct Hi-C contact for loops present
#           in all 15 samples (one anchor overlaps a CAD variant, the other a
#           protein-coding TSS); flag non-nearest genes.
# Paper   : Fig. 5a (11 genes: ANTXR1, ATMC12, CNNM2, FKBP5, CAPG, SH2D6,
#           AMOTL2, NUAK2, RGS1, PLCL1, FOXC2)
# Status  : REFERENCE IMPLEMENTATION — confirm "present" definition (score > 0?)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(rtracklayer) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config(); n_samples <- cfg$params$hic$n_samples
aqua <- cfg_path(cfg, "hic", "aqua_dir"); out <- cfg_path(cfg, "hic", "out")

qf <- list.files(file.path(aqua, "query"), "\\.inherent\\.bedpe$", full.names = TRUE)
stopifnot(length(qf) == n_samples)
q <- rbindlist(lapply(qf, function(f) fread(f, header = FALSE)[, sample := basename(f)]))
setnames(q, 1:6, c("chr1", "s1", "e1", "chr2", "s2", "e2"))
presence_col <- "V8"   # TODO: confirm score column and presence threshold
pres <- q[get(presence_col) > 0, .(n_present = uniqueN(sample)), by = .(chr1, s1, e1, chr2, s2, e2)]
fwrite(pres, file.path(out, "Fig5a_loop_presence.tsv"), sep = "\t")
core <- pres[n_present == n_samples]

v <- fread(cfg_path(cfg, "gwas", "cad_variant_superset"))
v[, chr := ifelse(startsWith(as.character(chr), "chr"), chr, paste0("chr", chr))]
v_gr <- GRanges(v$chr, IRanges(v$pos, width = 1), rsid = v$rsid)
gtf <- import(cfg_path(cfg, "refs", "gencode_gtf"))
genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
tss <- resize(genes, 1, fix = "start")

link <- function(vc, vs, ve, gc, gs, ge) {
  A <- GRanges(vc, IRanges(vs + 1L, ve)); B <- GRanges(gc, IRanges(gs + 1L, ge))
  ov_v <- as.data.table(findOverlaps(A, v_gr)); ov_g <- as.data.table(findOverlaps(B, tss))
  m <- merge(ov_v[, .(loop = queryHits, rsid = v_gr$rsid[subjectHits])],
             ov_g[, .(loop = queryHits, gene = tss$gene_name[subjectHits])], by = "loop", allow.cartesian = TRUE)
  m
}
res <- unique(rbind(link(core$chr1, core$s1, core$e1, core$chr2, core$s2, core$e2),
                    link(core$chr2, core$s2, core$e2, core$chr1, core$s1, core$e1)))

nearest_gene <- tss$gene_name[nearest(v_gr, tss, ignore.strand = TRUE)]
res[, nearest := nearest_gene[match(rsid, v_gr$rsid)]]
res[, is_nearest := gene == nearest]
fwrite(res, file.path(out, "Fig5a_cad_variant_distal_genes.tsv"), sep = "\t")
log_msg("Genes linked in all ", n_samples, " samples: ", paste(sort(unique(res$gene)), collapse = ", "))
