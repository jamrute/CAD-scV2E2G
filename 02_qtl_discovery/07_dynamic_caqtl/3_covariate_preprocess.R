library(data.table)
library(DESeq2)
library(dplyr)
library(Matrix)
library(singlecellmethods)

access_raw <- readRDS("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/3_Dynamic_QTL/counts/SMC_counts_processed/SMC_total_sparse_matrix.rds")
peaks <- fread(paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/3_Dynamic_QTL/counts/peaks.txt"),header=F)
rownames(access_raw) <- peaks$V1

#Donor level covariates
donor <- data.frame(donor=sapply(strsplit(colnames(access_raw),"_"), function(x) x[1]))
#Genotyping, Age, Sex, Accessibility Covariates
SMC_cvrt <- fread("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/QTL_covariates/cvrts_no_counts.txt")
colnames(SMC_cvrt) <- c("gPC1","gPC2","gPC3","gPC4","Age","Sex")
#Get SMC cell ID
ID_map <- fread("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/QTL_covariates/ID_map.csv", sep=',', header=TRUE)
ID_map <- ID_map[!(ID_map$s %in% c("D1144","T1171","T1069")),]
SMC_cvrt$ID <- ID_map$s
#Expand out
cvrt <- merge(donor, SMC_cvrt, by.x="donor", by.y="ID", all.x=TRUE)

#Cell level covariates
atacNorm <- normalizeData(access_raw, method = "log")
var_peaks <- vargenes_vst(access_raw, 3000)
atacNorm_top <- atacNorm[var_peaks[1:3000,]$symbol,]
pca_atac <- irlba::prcomp_irlba(t(atacNorm_top), 4)


nUMI <- scale(log(meta$nUMI))
MT <- meta$percent_mito
PC <- pcs[as.character(meta$donor),1:5]
expPC <- pca_res[,1:5]


