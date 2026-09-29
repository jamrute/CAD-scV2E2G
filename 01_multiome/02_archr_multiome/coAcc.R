library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/", sep=""))

scRNA <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/RNA/globalObjectConstruction/annotated/coronary_RNA_annotated.rds")
proj1 <- loadArchRProject(path = "/data/Junedh/Coronary_Multiomics/analysis/Save-proj1")


####### Co-accessibility analysis
proj1 <- addCoAccessibility(
    ArchRProj = proj1,
    reducedDims = "LSI_Combined",
    maxDist = 5e+06
)

# Approach 1 to obtain promoter peaks
cA <- getCoAccessibility(ArchRProj = proj1, returnLoops = FALSE, corCutOff = 0.3, resolution = 1)

cA_Chr <- (metadata(cA)$peakSet %>% seqnames(.))[cA$queryHits]
cA_Start <- (metadata(cA)$peakSet %>% start(.))[cA$queryHits]
cA_End <- (metadata(cA)$peakSet %>% end(.))[cA$queryHits]
cA_SubjectStart <- (metadata(cA)$peakSet %>% start(.))[cA$subjectHits]
id <- paste(cA_Chr, cA_SubjectStart)
cA_Bed <- data.frame(cA_Chr, cA_Start, cA_End, id)

pro_Chr <- proj1@peakSet@seqnames    #Obtain promoter peaks
pro_Start <- proj1@peakSet@ranges@start
id <- paste(pro_Chr, pro_Start)
peakType <- proj1@peakSet$peakType
PromoterPeaks <- data.frame(id, peakType)
PromoterPeaks <- PromoterPeaks[peakType == "Promoter",]

Promoter_cA <- inner_join(cA_Bed, PromoterPeaks, by="id")   #Intersect promoter peaks with subject peaks
Promoter_cA <- Promoter_cA[, -c(5)]
write.table(Promoter_cA, "./Save-proj1/ArchR_output/Promoter_cA_5Mb.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# Approach 2
cA <- getCoAccessibility(ArchRProj = proj1, returnLoops = TRUE, corCutOff = 0.3, resolution = 1)
write.table(cA$CoAccessibility,file="./Save-proj1/ArchR_output/CoAccessibility_final_5Mb.txt",quote=F,sep="\t",row.names=F)




