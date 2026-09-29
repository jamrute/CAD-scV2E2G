library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

proj1 <- loadArchRProject(path = "/data/Junedh/Coronary_Multiomics/analysis/Save-proj1")

# LMOD1
p <- plotBrowserTrack(
    ArchRProj = proj1, 
    groupBy = "CellType_Seurat_RNA",
    useGroups = c("SMC","Fibroblast","Myeloid","Endothelium"),
    geneSymbol = c("PALLD"),
    upstream = 500000,
    downstream = 500000,
    loops = getPeak2GeneLinks(proj1)
)

plotPDF(plotList = p, 
    name = "PALLD.pdf", 
    ArchRProj = proj1, 
    addDOC = FALSE, width = 5, height = 2.5)












