library(lme4)
library(Matrix)
library(data.table)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

peak <- args[1]
snp <- args[2]

# Load raw UMI counts matrix and metadata (from GEO: GSE158769) and peak expression PCs

chrom=1
geno <- fread(paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/RASQUAL_ATAC/genotype_matrix/36_samples/36_genotype_matrix_",chrom,".txt"))
CB <- fread(paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/CB_list.txt"),header=F)

#Raw SMC count processing, read from RDS object in the future.
#peaks <- fread(paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/counts/peaks.txt"),header=F)
#access_raw <- fread("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/counts/SMC_counts_processed/SMC_total_count.txt")
#colnames(access_raw) <- CB$V1
#rownames(access_raw) <- peaks$V1
#mat = Matrix(as.matrix(access_raw),sparse=TRUE)
#saveRDS(mat, "/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/counts/SMC_counts_processed/SMC_total_sparse_matrix.rds")


meta <- fread("/path/GSE158769_meta_data.txt.gz")
pca_res <- readRDS("/path/pca_irlba.rds")
pca_res <- pca_res$x

# Load genotype dosages (from dbGaP: phs002025) and PCs
geno <- fread("/path/geno_imputed_dosage_maf05_hg38_QC.vcf.gz")
pcs <- read.delim("/path/geno_PCs.txt", sep = " ", header = F)
CB <- 
row.names(pcs) <- pcs$V2
pcs <- pcs[,3:ncol(pcs)]

# Make data frame of variables for model
E <- as.numeric(access_raw[peak,]) 
G <- t(subset(geno, ID == snp)[,as.character(meta$donor)])[,1]
IND <- factor(meta$donor)
B <- meta$batch
AGE <- scale(meta$tbru_age)
SEX <- meta$Sex
nUMI <- scale(log(meta$nUMI))
MT <- meta$percent_mito
PC <- pcs[as.character(meta$donor),1:5]
expPC <- pca_res[,1:5]

# Load interaction covariate (e.g., batch-corrected CV from CCA)
cov <- readRDS("/path/cca_donor_batch.rds")[,1]

data <- data.frame(E,G,IND,B,AGE,SEX,nUMI,MT,
                   PC1 = PC[,1],PC2 = PC[,2],PC3 = PC[,3],PC4 = PC[,4],PC5 = PC[,5],
                   expPC1 = expPC[,1], expPC2 = expPC[,2], expPC3 = expPC[,3], expPC4 = expPC[,4], expPC5 = expPC[,5], cov)
data$G <- as.numeric(as.character(data$G))
    
full_model <- lme4::glmer(formula = E ~ G + (1 | B) + (1 | IND) + AGE + SEX + nUMI + MT + PC1 + PC2 + PC3 + PC4 + PC5 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + cov + G*cov, 
                          family = "poisson", nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap"))
null_model <- lme4::glmer(formula = E ~ G + (1 | B) + (1 | IND) + AGE + SEX + nUMI + MT + PC1 + PC2 + PC3 + PC4 + PC5 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + cov, 
                          family = "poisson", nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap"))
model_lrt <- anova(null_model, full_model)

out <- summary(test)$coefficients
colnames(out) <- c("Estimate","Std.Error","zvalue","pval")
out <- data.frame(peak=peak,snp=snp,term=row.names(x),x, lrt_pval=model_lrt$`Pr(>Chisq)`[2])