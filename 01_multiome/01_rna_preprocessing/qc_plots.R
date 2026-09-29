library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/", sep=""))

proj1 <- loadArchRProject(path = "/data/Junedh/Coronary_Multiomics/analysis/Save-proj1")

# p1 <- plotGroups(
#     ArchRProj = proj1, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "TSSEnrichment",
#     plotAs = "violin",
#     alpha = 0.4,
#     addBoxPlot = TRUE
#    )

# p2 <- plotGroups(
#     ArchRProj = proj1, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "log10(nFrags)",
#     plotAs = "violin",
#     alpha = 0.4,
#     addBoxPlot = TRUE
#    )

p3 <- plotFragmentSizes(ArchRProj = proj1, groupBy = "Sample")
p4 <- plotTSSEnrichment(ArchRProj = proj1, groupBy = "Sample")

# p5 <- plotGroups(
#     ArchRProj = proj1, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "Gex_nUMI",
#     plotAs = "violin",
#     alpha = 0.4,
#     addBoxPlot = TRUE
#    )

# p6 <- plotGroups(
#     ArchRProj = proj1, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "Gex_nGenes",
#     plotAs = "violin",
#     alpha = 0.4,
#     addBoxPlot = TRUE
#    )

# p7 <- plotGroups(
#     ArchRProj = proj1, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "Gex_MitoRatio",
#     plotAs = "violin",
#     alpha = 0.4,
#     addBoxPlot = TRUE
#    )

# p8 <- plotGroups(
#     ArchRProj = proj1, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "Gex_RiboRatio",
#     plotAs = "violin",
#     alpha = 0.4,
#     addBoxPlot = TRUE
#    )

plotPDF(p3,p4, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj1, addDOC = FALSE, width = 60, height = 10)
