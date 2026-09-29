library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/", sep=""))

proj1 <- loadArchRProject(path = "/data/Junedh/Coronary_Multiomics/analysis/Save-proj1")

getGroupBW(
  ArchRProj = proj1,
  groupBy = "CellType_Seurat_RNA",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 30000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
