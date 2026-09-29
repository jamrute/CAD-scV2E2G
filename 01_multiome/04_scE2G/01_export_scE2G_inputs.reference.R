#!/usr/bin/env Rscript
# =============================================================================
# 01_multiome/06_scE2G/01_export_scE2G_inputs.R
# Purpose : Split fragments and RNA counts by cell type for sc-E2G
# Outputs : results/multiome/scE2G_inputs/<CellType>/...
# Status  : SCAFFOLD — confirm exact input format for the sc-E2G version used
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config()
rna <- readRDS(cfg_path(cfg, "multiome", "seurat_obj"))
out <- file.path(cfg_path(cfg, "multiome", "out"), "scE2G_inputs")

for (ct in unique(rna$CellType)) {
  d <- file.path(out, ct); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  cells <- colnames(rna)[rna$CellType == ct]
  writeLines(cells, file.path(d, "barcodes.txt"))
  # TODO: filter each sample's atac_fragments.tsv.gz to `cells`
  #       (e.g. sinto filterbarcodes), bgzip + tabix;
  #       write raw RNA counts for `cells` in the format sc-E2G expects.
}
