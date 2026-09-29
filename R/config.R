# =============================================================================
# R/config.R — load paths + parameters for every R script in the repo.
#
# Usage (from any script):
#   source(file.path(Sys.getenv("CAD_REPO", "."), "R", "config.R"))
#   cfg <- load_config()
#   p   <- cfg_path(cfg, "multiome", "seurat_obj")
#   thr <- cfg$params$multiome$qc$min_tss_enrichment
# =============================================================================

suppressPackageStartupMessages(library(yaml))

repo_root <- function() {
  env <- Sys.getenv("CAD_REPO", unset = NA)
  if (!is.na(env)) return(normalizePath(env))
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "config", "params.yaml"))) {
    parent <- dirname(d)
    if (identical(parent, d)) stop("Cannot locate repo root; set CAD_REPO.")
    d <- parent
  }
  d
}

load_config <- function(paths_file = Sys.getenv("CAD_CONFIG", unset = NA)) {
  root <- repo_root()
  if (is.na(paths_file)) {
    local <- file.path(root, "config", "paths.local.yaml")
    paths_file <- if (file.exists(local)) local else file.path(root, "config", "paths.yaml")
  }
  list(
    repo   = root,
    paths  = yaml::read_yaml(paths_file),
    params = yaml::read_yaml(file.path(root, "config", "params.yaml"))
  )
}

#' Resolve a path from paths.yaml, relative to `root` unless already absolute.
cfg_path <- function(cfg, ..., create_dir = FALSE) {
  rel <- cfg$paths
  for (k in c(...)) {
    rel <- rel[[k]]
    if (is.null(rel)) stop("Missing path key: ", paste(c(...), collapse = "$"))
  }
  out <- if (startsWith(rel, "/")) rel else file.path(cfg$paths$root, rel)
  if (create_dir) dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

#' Resolve a repo-internal resource (e.g. resources/gene_sets/...).
repo_file <- function(cfg, ...) file.path(cfg$repo, ...)

log_msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)

#' Write sessionInfo next to outputs so every result is traceable.
save_session_info <- function(out_dir, script_name) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()),
             file.path(out_dir, paste0(tools::file_path_sans_ext(script_name), ".sessionInfo.txt")))
}
