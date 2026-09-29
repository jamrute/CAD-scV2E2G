library(ArchR)
addArchRGenome("hg38")
addArchRThreads(16)

# samples removed: "MGI3172_MGI3298_TWKH-T1054", "MGI3172_MGI3298_TWKH-T1004", "MGI3172_MGI3298_TWKH-T1037"

samples <- c("MGI2013_MGI2036_TWAP-T1182","MGI2078_MGI2080_TWAP-T1106","MGI2129_MGI2308_MGI2311_TWAP-T1069L",
			 "MGI2129_MGI2308_MGI2311_TWAP-T1106R","MGI2131_MGI2132_TWAP-T1114R","MGI2578_MGI2581_TWKH-D1127",
			 "MGI2578_MGI2581_TWKH-D1144","MGI2578_MGI2581_TWKH-D1148","MGI2578_MGI2581_TWKH-T1103",
			 "MGI2578_MGI2581_TWKH-T1133","MGI2578_MGI2581_TWKH-T1134","MGI2578_MGI2581_TWKH-T1135",
			 "MGI2578_MGI2581_TWKH-T1137","MGI2578_MGI2581_TWKH-T1145","MGI2578_MGI2581_TWKH-T1149",
			 "MGI2578_MGI2581_TWKH-T1163","MGI2701_MGI2706_TWKH-D1173","MGI2701_MGI2706_TWKH-T1096",
			 "MGI2701_MGI2706_TWKH-T1139","MGI2701_MGI2706_TWKH-T1171","MGI2904_MGI2917_TWKH-JA_T1179",
			 "MGI2904_MGI2917_TWKH-JA_T1182","MGI2904_MGI2917_TWKH-JA_T1183","MGI2904_MGI2917_TWKH-JA_T1184",
			 "MGI2904_MGI2917_TWKH-JA_T1185","MGI2904_MGI2917_TWKH-JA_T1186","MGI2904_MGI2917_TWKH-JA_T1187",
			 "MGI2904_MGI2917_TWKH-JA_T1189","MGI3172_MGI3298_TWKH-T1036",
			 "MGI3172_MGI3298_TWKH-T1038","MGI3172_MGI3298_TWKH-T1040",
			 "MGI3172_MGI3298_TWKH-T1090","MGI3172_MGI3298_TWKH-T1091",
			 "MGI3339_MGI3341_TWKH-T1022","MGI3339_MGI3341_TWKH-T1055","MGI3339_MGI3341_TWKH-T1123",
			 "MGI3339_MGI3341_TWKH-T1127","MGI3339_MGI3341_TWKH-T1167","MGI3339_MGI3341_TWKH-T1169",
			 "MGI3339_MGI3341_TWKH-T1190","MGI3339_MGI3341_TWKH-T1201","MGI3339_MGI3341_TWKH-T1210",
			 "MGI3339_MGI3341_TWKH-T1211","MGI3339_MGI3341_TWKH-T1214","MGI3339_MGI3341_TWKH-T1215")

for (sample in samples) {
dir.create(paste("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/", sample, sep=""), showWarnings = FALSE)

#Set working directory
setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/", sample, sep=""))

inputFiles <- paste("/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/", sample, "/outs/atac_fragments.tsv.gz", sep="")
names(inputFiles) <- sample

#Create Arrow Files and ArchR project
createArrowFiles(inputFiles, force = TRUE) 
ArrowFiles <- paste(sample, ".arrow", sep="") 
proj <- ArchRProject(ArrowFiles)

#Import scRNA data
seRNA <- import10xFeatureMatrix(
    input = c(paste("/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/", sample, "/outs/filtered_feature_bc_matrix.h5", sep="")),
    names = c(sample)
)

proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

#Quality control filtering
origProj <- proj
proj <- proj[proj$TSSEnrichment > 2 & proj$nFrags > 1000 & !is.na(proj$Gex_nUMI)]
proj <- proj[proj$Gex_nUMI > 200 & proj$Gex_nUMI < 50000 & proj$Gex_MitoRatio < 0.05]

#Doublet filtration
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj)

#Obtain filtered cell list for Seurat
cells <- proj$cellNames
cells <- gsub(paste(sample, "#", sep=""), "", cells)  

#Save cells
saveRDS(cells, file=paste("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/", sample, "_cells", sep=""))

}