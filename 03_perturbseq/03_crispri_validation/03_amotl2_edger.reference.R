#!/usr/bin/env Rscript
# =============================================================================
# 03_perturbseq/crispri_validation/03_amotl2_edger.R
# Purpose : Differential expression, AMOTL2 TSS KD vs NTC (stressed HCASMCs)
#           edgeR 3.42.4: negative-binomial GLM (glmFit + glmLRT);
#           significant: |logFC| > 0.5 and FDR < 0.05. CPM for plotting.
#           GO BP enrichment (one-sided Fisher / BH) on up and down genes.
# Paper   : Fig. 5e (heatmap), 5f (GO BP), 5g (contractile genes)
# Inputs  : dir of HTSeq count files named <sample>.txt; samples TSV (sample, condition)
# Status  : REFERENCE IMPLEMENTATION from Methods
# =============================================================================
suppressPackageStartupMessages({ library(edgeR); library(data.table); library(optparse) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts-dir"), make_option("--samples"), make_option("--out"),
  make_option("--gene-map", default = NA, help = "TSV gene_id -> gene_name")
)))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

ss <- fread(opt$samples)
cnt <- Reduce(function(a, b) merge(a, b, by = "gene_id"), lapply(ss$sample, function(s) {
  x <- fread(file.path(opt$`counts-dir`, paste0(s, ".txt")), header = FALSE,
             col.names = c("gene_id", s))
  x[!startsWith(gene_id, "__")]
}))
M <- as.matrix(cnt, rownames = "gene_id")

grp <- factor(ss$condition, levels = c("NTC", "AMOTL2_KD"))
y <- DGEList(M, group = grp)
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ grp)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
tt <- as.data.table(topTags(lrt, n = Inf)$table, keep.rownames = "gene_id")
tt[, sig := abs(logFC) > 0.5 & FDR < 0.05]
if (!is.na(opt$`gene-map`)) tt <- merge(fread(opt$`gene-map`), tt, by = "gene_id", all.y = TRUE)
fwrite(tt, file.path(opt$out, "AMOTL2_KD_vs_NTC_edgeR.tsv"), sep = "\t")
fwrite(as.data.table(cpm(y), keep.rownames = "gene_id"), file.path(opt$out, "cpm.tsv"), sep = "\t")
message("Up: ", tt[sig & logFC > 0, .N], "  Down: ", tt[sig & logFC < 0, .N])

if (requireNamespace("clusterProfiler", quietly = TRUE) && "gene_name" %in% names(tt)) {
  for (d in c("up", "down")) {
    g <- tt[sig & (if (d == "up") logFC > 0 else logFC < 0), gene_name]
    ego <- clusterProfiler::enrichGO(g, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
                                     ont = "BP", universe = tt$gene_name, pAdjustMethod = "BH")
    if (!is.null(ego)) fwrite(as.data.table(ego@result), file.path(opt$out, paste0("Fig5f_GO_BP_", d, ".tsv")), sep = "\t")
  }
}
