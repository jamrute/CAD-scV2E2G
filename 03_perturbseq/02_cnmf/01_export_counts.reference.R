#!/usr/bin/env Rscript
# =============================================================================
# 03_perturbseq/06_cnmf/01_export_counts.R
# Purpose : Export raw counts of Mixscape KD + NT cells for cNMF
# Outputs : results/perturbseq/06_cnmf/input/{matrix.mtx, genes.tsv, cells.tsv, obs.tsv}
# Status  : REFERENCE IMPLEMENTATION
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
cfg <- load_config()
base <- cfg_path(cfg, "perturbseq", "out")
out <- file.path(base, "06_cnmf", "input"); dir.create(out, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(base, "02_seurat_mixscape.rds"))
obj <- subset(obj, mixscape_class.global %in% c("KD", "NT"))
m <- GetAssayData(obj, assay = "RNA", slot = "counts")
writeMM(m, file.path(out, "matrix.mtx"))                 # genes x cells
writeLines(rownames(m), file.path(out, "genes.tsv"))
writeLines(colnames(m), file.path(out, "cells.tsv"))
fwrite(as.data.table(obj@meta.data[, c("gene", "prtb_type", "mixscape_class.global")],
                     keep.rownames = "cell"), file.path(out, "obs.tsv"), sep = "\t")
