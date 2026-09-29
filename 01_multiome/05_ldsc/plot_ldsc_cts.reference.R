#!/usr/bin/env Rscript
# =============================================================================
# 01_multiome/07_ldsc/04_plot_ldsc.R
# Purpose : Plot -log10 P of cell-type enrichment per trait (Fig. 1f)
# Status  : REFERENCE IMPLEMENTATION
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
source(file.path(Sys.getenv("CAD_REPO", "."), "R", "plot_theme.R"))
cfg <- load_config()
wd  <- file.path(cfg_path(cfg, "multiome", "out"), "ldsc", "results")
res <- rbindlist(lapply(list.files(wd, "\\.cell_type_results\\.txt$", full.names = TRUE), function(f) {
  x <- fread(f); x[, Trait := sub("\\.cell_type_results\\.txt$", "", basename(f))]
}))
p <- ggplot(res, aes(Name, -log10(Coefficient_P_value), fill = Name)) +
  geom_col() + facet_wrap(~Trait, nrow = 1) +
  geom_hline(yintercept = -log10(cfg$params$ldsc$bonferroni_p), linetype = 2) +
  scale_fill_manual(values = pal_celltype) + labs(x = NULL) + theme_cad() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_fig(p, file.path(wd, "Fig1f_ldsc.pdf"), width_mm = 180, height_mm = 60)
fwrite(res, file.path(wd, "Fig1f_ldsc_all_traits.tsv"), sep = "\t")
