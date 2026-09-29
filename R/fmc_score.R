# =============================================================================
# R/fmc_score.R — fibromyocyte (FMC) module score, shared across modules.
#
# Methods: FMC identity score = module score over the top 100 modulated-SMC
# (fibromyocyte) marker genes from Wirka et al. 2019 (Nat Med), computed with
# Seurat::AddModuleScore on RNA (Perturb-seq, Multiome) or on ArchR gene
# activity scores (snATAC meta-map, for dynamic caQTL modelling).
#
# Gene list: resources/gene_sets/fmc_wirka2019_top100.txt
# =============================================================================

read_fmc_genes <- function(cfg, n = cfg$params$pme$fmc_top_n_genes) {
  f <- repo_file(cfg, "resources", "gene_sets", "fmc_wirka2019_top100.txt")
  g <- readLines(f)
  g <- g[nzchar(g) & !startsWith(g, "#")]
  if (length(g) < n) warning("FMC gene list has ", length(g), " genes; expected ", n)
  head(g, n)
}

#' Add an FMC score column to a Seurat object.
#' @param assay "RNA" for expression, or the assay holding gene-activity scores.
add_fmc_score <- function(obj, genes, assay = Seurat::DefaultAssay(obj),
                          name = "FMC_score", seed = 1) {
  genes_present <- intersect(genes, rownames(obj[[assay]]))
  if (length(genes_present) < 0.8 * length(genes))
    warning(sprintf("Only %d/%d FMC genes found in assay '%s'",
                    length(genes_present), length(genes), assay))
  obj <- Seurat::AddModuleScore(obj, features = list(genes_present), assay = assay,
                                name = name, seed = seed)
  # AddModuleScore appends "1" to the name
  obj[[name]] <- obj[[paste0(name, "1")]]
  obj[[paste0(name, "1")]] <- NULL
  obj
}

#' Bin a continuous FMC score into tertiles (ED Fig. 6c; Fig. 2k).
fmc_tertile <- function(x) {
  cut(x, breaks = stats::quantile(x, c(0, 1/3, 2/3, 1), na.rm = TRUE),
      labels = c("Low", "Medium", "High"), include.lowest = TRUE)
}
