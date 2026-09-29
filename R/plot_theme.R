# =============================================================================
# R/plot_theme.R — shared ggplot theme and colour palettes for all figures.
# TODO: replace palette values with those used in the published figures.
# =============================================================================

theme_cad <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_blank(),
      legend.key.size  = grid::unit(3, "mm")
    )
}

# Cell-type palette (11 Multiome cell types; Fig. 1b)
pal_celltype <- c(
  SMC = "#E64B35", Endothelium = "#4DBBD5", Fibroblast = "#00A087",
  Myeloid = "#3C5488", TNKCell = "#F39B7F", BCells = "#8491B4",
  PlasmaCells = "#91D1C2", Adipocyte = "#DC0000", Neuron = "#7E6148",
  Lymphatic = "#B09C85", Epicardial = "#FFDC91"
)

pal_perturb_type <- c(TSS = "#4575B4", Enhancer = "#F4A582", NT = "grey60")

save_fig <- function(p, file, width_mm = 89, height_mm = 60, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, p, width = width_mm, height = height_mm, units = "mm",
                  device = grDevices::cairo_pdf, ...)
}
