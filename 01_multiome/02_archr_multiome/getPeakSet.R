library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/", sep=""))
proj1 <- loadArchRProject(path = "/data/Junedh/Coronary_Multiomics/analysis/Save-proj1")

####### Co-accessibility analysis
peakSet <- proj1@peakSet
saveRDS(peakSet, "./Save-proj1/ArchR_output/peakSet.rds")