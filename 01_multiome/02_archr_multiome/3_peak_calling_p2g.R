library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/", sep=""))

scRNA <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/RNA/globalObjectConstruction/annotated/coronary_RNA_annotated.rds")
proj1 <- loadArchRProject(path = "/data/Junedh/Coronary_Multiomics/analysis/Save-proj1")

############################## Peak Calling with MACS2 ############################################################################

# Pseudobulk ATAC and Call Peaks
proj2 <- addGroupCoverages(ArchRProj = proj1, groupBy = "CellType_Seurat_RNA")
pathToMacs2 <- findMacs2()
proj1 <- addReproduciblePeakSet(
    ArchRProj = proj1, 
    groupBy = "CellType_Seurat_RNA", 
    pathToMacs2 = pathToMacs2,
    maxPeaks = 50000000
)

proj1 <- addPeakMatrix(proj1)

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)

getAvailableMatrices(proj1)

#Peak Matrix Marker Peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj1, 
    useMatrix = "PeakMatrix", 
    useGroups = c("Myeloid","Endothelium","SMC","Fibroblast","Adipocyte"),
    groupBy = "CellType_Seurat_RNA",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    testMethod = "wilcoxon",
    maxCells = 5000
)

saveRDS(markersPeaks, file="./Save-proj1/ArchR_output/markerPeak_CellType_Seurat_RNA_5000.rds")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.25")

# Plot heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.25",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker_noMax-Heatmap_celltypeRNA_5000.pdf", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)

#Create accessibility list organized by cell type
celltypes <- unique(as.character(scRNA$cell.type.l1))

cellPeaksBed <- NULL
for (type in celltypes) {
    list <- readRDS(paste0("./Save-proj1/PeakCalls/", type, "-reproduciblePeaks.gr.rds"))
    chr <- as.vector(list %>% seqnames(.))
    start <- list %>% start(.)
    end <- list %>% end(.)
    type <- rep(type, each=length(list))    
    typeBed <- data.frame(chr, start, end, type)
    cellPeaksBed <- rbind(cellPeaksBed, typeBed)
}
write.table(cellPeaksBed, "./Save-proj1/ArchR_output/cellPeaks.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

############################################### Peak 2 Gene ###########################################################

proj1 <- addPeak2GeneLinks(ArchRProj = proj1, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")

p <- plotPeak2GeneHeatmap(ArchRProj = proj1, groupBy = "CellType_Seurat_RNA", corCutOff = 0.3)
plotPDF(p, name = "plotPeak2GeneHeatmap.pdf", width = 4, height = 8, ArchRProj = proj1, addDOC = FALSE)

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)

##########################################################  GPC Analysis  ############################################################
#GSM/GEX correlation
corr_vals <- correlateMatrices(proj1, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "GeneExpressionMatrix",reducedDims = "HAR_Combined") #replace with GeneExpressionMatrix for multiome
corr_vals <- as.data.frame(corr_vals)

#remove insignificant correlations
corr_vals <- corr_vals[!is.na(corr_vals$cor),]
corr_vals <- corr_vals[corr_vals$padj < 0.05,]

#get peak 2 gene links
p2g <- getPeak2GeneLinks(proj1, returnLoops = FALSE)
x <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/Save-proj1/Peak2GeneLinks/seRNA-Group-KNN.rds")
p2g$idxRNA <- x@rowRanges$name[p2g$idxRNA]

gpc <- as.data.frame(table(p2g$idxRNA))
rownames(gpc) <- gpc$Var1

gpc <- merge(gpc,corr_vals[c('GeneScoreMatrix_name','cor')], by.x='Var1', by.y='GeneScoreMatrix_name')

#ggscatter(x='Freq',y='cor',data = gpc)

write.csv(gpc,"./Save-proj1/ArchR_output/gpc_pilot_table.csv", quote = FALSE)
write.csv(corr_vals,"./Save-proj1/ArchR_output/gex_gsm_corr_pilot_table.csv", quote = FALSE)

######### ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 

# Approach 1: Create peak2gene bed file
p2g <- getPeak2GeneLinks(
    ArchRProj = proj1,
    corCutOff = 0.3,
    resolution = 1,
    returnLoops = FALSE
)

saveRDS(p2g, "./Save-proj1/ArchR_output/p2g_withCorr.rds")

# p2gChr <- (metadata(p2g)$peakSet %>% seqnames(.))[p2g$idxATAC]
# p2gStart <- (metadata(p2g)$peakSet %>% start(.))[p2g$idxATAC]
# p2gEnd <- (metadata(p2g)$peakSet %>% end(.))[p2g$idxATAC]
# p2gGene <- (metadata(p2g)$geneSet)$name[p2g$idxRNA]

# p2gCorrelation <- 

# p2gBed <- data.frame(p2gChr, p2gStart, p2gEnd, p2gGene, p2gCorrelation)
# write.table(p2gBed, "./Save-proj1/ArchR_output/p2g_withCorr.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# # Approach 2
# p2g <- getPeak2GeneLinks(ArchRProj = proj1, returnLoops = TRUE, resolution = 1, corCutOff = 0.3)
# write.table(p2g$Peak2GeneLinks,file="./Save-proj1/ArchR_output/p2g.txt",quote=F,sep="\t",row.names=F)

# ####### Co-accessibility analysis
# proj1 <- addCoAccessibility(
#     ArchRProj = proj1,
#     reducedDims = "LSI_Combined"
# )

# # Approach 1 to obtain promoter peaks
# cA <- getCoAccessibility(ArchRProj = proj1, returnLoops = FALSE, corCutOff = 0.3, resolution = 1)

# cA_Chr <- (metadata(cA)$peakSet %>% seqnames(.))[cA$queryHits]
# cA_Start <- (metadata(cA)$peakSet %>% start(.))[cA$queryHits]
# cA_End <- (metadata(cA)$peakSet %>% end(.))[cA$queryHits]
# cA_SubjectStart <- (metadata(cA)$peakSet %>% start(.))[cA$subjectHits]
# id <- paste(cA_Chr, cA_SubjectStart)
# cA_Bed <- data.frame(cA_Chr, cA_Start, cA_End, id)

# pro_Chr <- proj1@peakSet@seqnames    #Obtain promoter peaks
# pro_Start <- proj1@peakSet@ranges@start
# id <- paste(pro_Chr, pro_Start)
# peakType <- proj1@peakSet$peakType
# PromoterPeaks <- data.frame(id, peakType)
# PromoterPeaks <- PromoterPeaks[peakType == "Promoter",]

# Promoter_cA <- inner_join(cA_Bed, PromoterPeaks, by="id")   #Intersect promoter peaks with subject peaks
# Promoter_cA <- Promoter_cA[, -c(5)]
# write.table(Promoter_cA, "./Save-proj1/ArchR_output/Promoter_cA.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# # Approach 2
# cA <- getCoAccessibility(ArchRProj = proj1, returnLoops = TRUE, corCutOff = 0.3, resolution = 1)
# write.table(cA$CoAccessibility,file="./Save-proj1/ArchR_output/CoAccessibility_final.txt",quote=F,sep="\t",row.names=F)

# ###################################################### Motif annotation
# markersPeaks <- readRDS("./Save-proj1/ArchR_output/markerPeak_CellType_Seurat_RNA_15000.rds")

# proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "homer", name = "Motif",force = TRUE)

# enrichMotifs <- peakAnnoEnrichment(
#     seMarker = markersPeaks,
#     ArchRProj = proj1,
#     peakAnnotation = "Motif",
#     cutOff = "FDR <= 0.1 & Log2FC >= 0.25"
#   )

# heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapEM, name = "enrichMotifs_homer_celltype_Heatmap_final.pdf", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)

# ##########################################################################################################################################
# chromVAR
proj1 <- addBgdPeaks(proj1)
proj1 <- addDeviationsMatrix(
  ArchRProj = proj1, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj1, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "chromVAR_VariableMotifDeviationScores_final.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)








