library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/", sep=""))

scRNA <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/RNA/globalObjectConstruction/annotated/coronary_RNA_annotated.rds")
proj1 <- loadArchRProject(path = "/data/Junedh/Coronary_Multiomics/analysis/Save-proj1")

markersPeaks <- readRDS("./Save-proj1/ArchR_output/markerPeak_CellType_Seurat_RNA_15000.rds")