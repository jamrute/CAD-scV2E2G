library(Seurat)
library(dplyr)
library(harmony)

merged <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/RNA/RNA_merged_postQC_postScrublet.rds")

DefaultAssay(merged) <- 'RNA'
merged <- SCTransform(merged, vars.to.regress = c("percent.mt", "nCount_RNA"))
merged <- RunPCA(merged, npcs=100, verbose=TRUE)
merged <- RunHarmony(merged, c("sample"), reduction = "pca", reduction.save = "harmony", assay.use = "SCT")
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:50)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:50)
merged <- FindClusters(merged, graph.name = "SCT_snn", algorithm = 3, resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = TRUE)

saveRDS(merged, "/data/Junedh/Coronary_Multiomics/analysis/RNA/integrated_v1.rds")
