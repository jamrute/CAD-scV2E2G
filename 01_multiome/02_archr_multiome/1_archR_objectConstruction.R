library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

################ QC and Generate the ArchR project after filtering and add gene expression matrix from paired scRNAseq data

#Get input fragment files for each sample
setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/", sep=""))

samples <- c("MGI2013_MGI2036_TWAP-T1182","MGI2078_MGI2080_TWAP-T1106","MGI2129_MGI2308_MGI2311_TWAP-T1069L",
       "MGI2129_MGI2308_MGI2311_TWAP-T1106R","MGI2131_MGI2132_TWAP-T1114R","MGI2578_MGI2581_TWKH-D1127",
       "MGI2578_MGI2581_TWKH-D1144","MGI2578_MGI2581_TWKH-D1148","MGI2578_MGI2581_TWKH-T1103",
       "MGI2578_MGI2581_TWKH-T1133","MGI2578_MGI2581_TWKH-T1134","MGI2578_MGI2581_TWKH-T1135",
       "MGI2578_MGI2581_TWKH-T1137","MGI2578_MGI2581_TWKH-T1145","MGI2578_MGI2581_TWKH-T1149",
       "MGI2578_MGI2581_TWKH-T1163","MGI2701_MGI2706_TWKH-D1173","MGI2701_MGI2706_TWKH-T1096",
       "MGI2701_MGI2706_TWKH-T1139","MGI2701_MGI2706_TWKH-T1171","MGI2904_MGI2917_TWKH-JA_T1179",
       "MGI2904_MGI2917_TWKH-JA_T1182","MGI2904_MGI2917_TWKH-JA_T1183",
       "MGI2904_MGI2917_TWKH-JA_T1185","MGI2904_MGI2917_TWKH-JA_T1186","MGI2904_MGI2917_TWKH-JA_T1187",
       "MGI2904_MGI2917_TWKH-JA_T1189","MGI3172_MGI3298_TWKH-T1036",
       "MGI3172_MGI3298_TWKH-T1038","MGI3172_MGI3298_TWKH-T1040",
       "MGI3172_MGI3298_TWKH-T1090","MGI3172_MGI3298_TWKH-T1091",
       "MGI3339_MGI3341_TWKH-T1022","MGI3339_MGI3341_TWKH-T1055","MGI3339_MGI3341_TWKH-T1123",
       "MGI3339_MGI3341_TWKH-T1127","MGI3339_MGI3341_TWKH-T1167","MGI3339_MGI3341_TWKH-T1169",
       "MGI3339_MGI3341_TWKH-T1190","MGI3339_MGI3341_TWKH-T1201","MGI3339_MGI3341_TWKH-T1210",
       "MGI3339_MGI3341_TWKH-T1211","MGI3339_MGI3341_TWKH-T1214","MGI3339_MGI3341_TWKH-T1215")

ArrowFiles <- c()

for (s in samples) {
  ArrowFiles <- c(ArrowFiles, paste(s, "/", s, ".arrow", sep=""))
}

proj1 <- ArchRProject(ArrowFiles, copyArrows = TRUE)

scRNA <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/RNA/globalObjectConstruction/annotated/coronary_RNA_annotated.rds")
col <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/RNA/globalObjectConstruction/annotated/coronary_multiome_cellList")

#Isolate cells used in Seurat
proj1 <- subsetCells(ArchRProj = proj1, cellNames = col)

# GEX Matrix
gene_matrix_files <- c()
counts_location <- "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/"

for (s in samples) {
  gene_matrix_files <- c(gene_matrix_files, paste(counts_location, s, "/outs/filtered_feature_bc_matrix.h5", sep=""))
}

seRNA <- import10xFeatureMatrix(input = gene_matrix_files, names = samples)

seRNAcombined<-cbind(assay(seRNA[[1]]), assay(seRNA[[2]]), assay(seRNA[[3]]), assay(seRNA[[4]]), assay(seRNA[[5]]), assay(seRNA[[6]]), assay(seRNA[[7]]),
                     assay(seRNA[[8]]), assay(seRNA[[9]]), assay(seRNA[[10]]), assay(seRNA[[11]]), assay(seRNA[[12]]), assay(seRNA[[13]]), assay(seRNA[[14]]),
                     assay(seRNA[[15]]), assay(seRNA[[16]]), assay(seRNA[[17]]), assay(seRNA[[18]]), assay(seRNA[[19]]), assay(seRNA[[20]]), assay(seRNA[[21]]),
                     assay(seRNA[[22]]), assay(seRNA[[23]]), assay(seRNA[[24]]), assay(seRNA[[25]]), assay(seRNA[[26]]), assay(seRNA[[27]]), assay(seRNA[[28]]),
                     assay(seRNA[[29]]), assay(seRNA[[30]]), assay(seRNA[[31]]), assay(seRNA[[32]]), assay(seRNA[[33]]), assay(seRNA[[34]]), assay(seRNA[[35]]),
                     assay(seRNA[[36]]), assay(seRNA[[37]]), assay(seRNA[[38]]), assay(seRNA[[39]]), assay(seRNA[[40]]), assay(seRNA[[41]]), assay(seRNA[[42]]),
                     assay(seRNA[[43]]), assay(seRNA[[44]]))

seRNA2<-SummarizedExperiment(assays=list(counts=seRNAcombined), rowRanges= rowRanges(seRNA[[1]]))

proj1 <- addGeneExpressionMatrix(
  input = proj1,
  seRNA = seRNA2,
  chromSizes = getChromSizes(proj1),
  excludeChr = c("chrM", "chrY"),
  scaleTo = 10000,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addGeneExpressionMatrix")
)

# Transfer cluster and identities from scRNA
clusters <- as.character(scRNA$cell.type)
stent <- as.character(scRNA$stent)
age <- as.character(scRNA$age)
sex <- as.character(scRNA$sex)
race <- as.character(scRNA$race)
HF <- as.character(scRNA$HF)

proj1$celltypeRNA <- clusters
proj1$stent <- stent
proj1$age <- age
proj1$sex <- sex
proj1$race <- race
proj1$HF <- HF

setwd(paste("/data/Junedh/Coronary_Multiomics/analysis/", sep=""))

proj1 <- saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)