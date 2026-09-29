library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/", sep=""))
proj1 <- loadArchRProject(path = "/data/Junedh/Coronary_Multiomics/analysis/Save-proj1")

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

write.csv(gpc,"./gpc_pilot_table.csv", quote = FALSE)
write.csv(corr_vals,"./gex_gsm_corr_pilot_table.csv", quote = FALSE)