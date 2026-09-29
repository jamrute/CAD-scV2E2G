#!/usr/bin/env Rscript
# =============================================================================
# 03_perturbseq/06_cnmf/04_annotate_programs.R
# Purpose : Import cNMF usages + gene-specificity scores; annotate programs by
#           ORA (clusterProfiler::enricher v4.20.0) of top 100 genes vs MSigDB
#           Hallmark + GO:BP; retain 6/8 programs (2 unannotated excluded);
#           add usages to Seurat metadata; program z-score by perturbation.
# Paper   : Fig. 3f; ED 9d (left); ST 43–44
# Status  : REFERENCE IMPLEMENTATION — program labels must be set manually
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(clusterProfiler); library(msigdbr); library(data.table)
})
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config(); pc <- cfg$params$perturbseq$cnmf
base <- cfg_path(cfg, "perturbseq", "out"); wd <- file.path(base, "06_cnmf")
name <- Sys.getenv("CNMF_NAME", "hcasmc_perturbseq")
tag  <- sprintf("k_%d.dt_%s", pc$k_selected, gsub("\\.", "_", as.character(pc$local_density_threshold)))

usage  <- fread(file.path(wd, name, sprintf("%s.usages.%s.consensus.txt", name, tag)))
setnames(usage, 1, "cell")
spectra <- fread(file.path(wd, name, sprintf("%s.gene_spectra_score.%s.txt", name, tag)))
spec <- t(as.matrix(spectra[, -1, with = FALSE])); colnames(spec) <- paste0("P", spectra[[1]])
fwrite(as.data.table(spec, keep.rownames = "gene"), file.path(wd, "program_gene_scores.tsv"), sep = "\t")

t2g <- rbind(
  as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))[, .(gs_name, gene_symbol)],
  as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))[, .(gs_name, gene_symbol)]
)
ora <- rbindlist(lapply(colnames(spec), function(p) {
  top <- names(sort(spec[, p], decreasing = TRUE))[seq_len(pc$top_genes_for_annotation)]
  e <- enricher(top, TERM2GENE = t2g, universe = rownames(spec))
  if (is.null(e)) return(NULL)
  cbind(program = p, as.data.table(e@result))
}))
fwrite(ora, file.path(wd, "program_ora.tsv"), sep = "\t")

# TODO: fill from ORA review (6 retained programs)
program_labels <- c(
  # P1 = "Inflammation", P2 = "ECM remodeling", P3 = "Vascular remodeling",
  # P4 = "Interferon response", P5 = "Proliferation (S)", P6 = "Proliferation (M)"
)

obj <- readRDS(file.path(base, "02_seurat_mixscape.rds"))
u <- as.data.frame(usage); rownames(u) <- u$cell; u$cell <- NULL
colnames(u) <- paste0("cNMF_P", colnames(u))
u <- u / rowSums(u)                                  # normalise usages per cell
obj <- AddMetaData(obj, u)
saveRDS(obj, file.path(base, "06_seurat_cnmf.rds"))

# Fig. 3f: mean usage per perturbation, z-scored across perturbations
md <- as.data.table(obj@meta.data)[mixscape_class.global %in% c("KD", "NT")]
mu <- md[, lapply(.SD, mean), by = gene, .SDcols = colnames(u)]
z  <- mu[, lapply(.SD, function(v) as.numeric(scale(v))), .SDcols = colnames(u)][, gene := mu$gene]
fwrite(z, file.path(wd, "Fig3f_program_zscores.tsv"), sep = "\t")
