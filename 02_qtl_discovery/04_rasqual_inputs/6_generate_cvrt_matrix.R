library(dplyr)
library(plyr)
library(data.table)
library(rasqualTools)
library(ggplot2)


randomize <- function(x,g=NULL){
  if(is.null(g)){
    n=ncol(x);
    t(apply(x,1,function(xx){xx[order(runif(n))]}))
  }else{
    for(i in unique(g)){
      x[,g==i]=randomize(x[,g==i,drop=F])
    }
    x
  }
}

rasqualMakeCovariates <- function(counts, size_factors) {
  
  #Map parameters to Natsuhiko's variables
  Y = counts
  K = size_factors
  n=ncol(Y)
  
  # fpm calculation
  fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9
  
  # Singular value decomposition
  fpkm.svd   = svd((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd))
  fpkm.svd.r = svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))
  
  # Covariate selection
  sf=log(apply(Y,2,sum))
  covs=fpkm.svd$v[,1:sum(fpkm.svd$d[-n]>fpkm.svd.r$d[-n])]
  if(cor(sf,covs[,1])^2<0.9){covs=covs}
  
  # Write covariates
  return(covs)
}

wustl_ID_meta <- read.table('/ccdg/Active/analysis/l.paullee/multiome_qtl/common_files/metadata/patient_metadata/wustl_id_map_new.tsv', header=TRUE)

wustl_ID_meta <- wustl_ID_meta[, c("g", "p","Age", "Sex_Recoded")]
names(wustl_ID_meta)[names(wustl_ID_meta) == 'p'] <- 's'
wustl_ID_meta$source <- 0

#Add metadata from Clint's study, in order
snatac_meta <- read.table('/ccdg/Active/analysis/l.paullee/multiome_qtl/common_files/metadata/patient_metadata/UVA_patient_meta.txt', sep='\t', header=TRUE)
snatac_meta$Sex_Recoded <- ifelse(snatac_meta$Sex  == "M", 1, 0)
snatac_meta$g <- snatac_meta$scID
snatac_meta$s <- snatac_meta$g
snatac_df <- snatac_meta[, c("g", "s","Age", "Sex_Recoded")]
snatac_df$source <- 1

#concatenate
merge_df <- rbind(wustl_ID_meta, snatac_df)

#Note 79 samples, one will be lost after merge with ancestry

# Read ancestry file, select 4 ancestry PCs
#Note that conversion from VCF to plink inexplicably dropped the number prefixes in front of the WUSTL samples
#New file plink_reid.eigenvec generated.
ancestry_pcs <- read.table("/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/plink/plink_reid.eigenvec", sep=' ', header=TRUE)
ancestry_pcs$g <- c(ancestry_pcs$FID_IID[1:38],ancestry_pcs$FID[39:79])
ancestry_pcs <- ancestry_pcs[, c('g','PC1', 'PC2', 'PC3', 'PC4')]
colnames(ancestry_pcs) <- c("g","ancestry_PC1", "ancestry_PC2", "ancestry_PC3", "ancestry_PC4")

# Merge with ancestry_pcs and drop unnecessary columns
result <- join(ancestry_pcs, merge_df, by="g")


#Generate covariates with and without counts
ct.list = c("SMC","Myeloid", "Endothelium", "Fibroblast")
for(ct in ct.list){
print(ct)
ct_subset <- read.table(paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/",ct,"_sample_list.tsv"))
colnames(ct_subset) <- c("g")
out <- join(ct_subset, result, by="g")
#result <- result[, -5]
out <- out[,-which(names(out) %in% c("g","s"))]
#Without counts
write.table(out, paste0("/ccdg/Active/analysis//l.paullee/multiome_qtl/qtl/data/covariates/",ct,"_cvrts_no_counts.txt"), sep='\t', quote=FALSE, row.names=FALSE, col.names = FALSE)
#Create count covariates using rasualMakeCovariates
counts = read.table(paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/",ct,"_count_mtx.tsv"))
offset = read.table(paste0("/ccdg/Active/analysis//l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/",ct,"_offset.txt"))
Y=as.matrix(counts)
K=as.matrix(offset)
count_cov <- rasqualMakeCovariates(Y,K)
out_counts <- cbind(out,count_cov)
write.table(out_counts, paste0("/ccdg/Active/analysis//l.paullee/multiome_qtl/qtl/data/covariates/",ct,"_cvrts_counts.txt"), sep='\t', quote=FALSE, row.names=FALSE, col.names = FALSE)
}

#Add count covariates (Old)
#ct.list = c("SMC","Myeloid", "Endothelium", "Fibroblast")
#for(ct in ct.list){
#bamlist = read.table("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/data/VCF/atac_bam_list/bamlist.txt")
#counts = read.table(paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/0b_Signac_RASQUAL/counts/",ct,"_total_count.txt"))
#mtx <- data.matrix(counts[,c(-1)])
#colnames(mtx) <- unlist(bamlist["V1"])
#rownames(mtx) = counts$V1
#type <- rep(1,36)
#SampleName <- seq(1:36)
#colData <- data.frame(SampleName, type)
#dds <- DESeqDataSetFromMatrix(countData=mtx, colData=colData, design = ~1)
#vsd <- vst(dds)
#vst_counts <- assay(vsd)
##Count PCs
#pca <- prcomp(t(assay(vsd)))
#result[,c("PC1","PC2","PC3","PC4")] <- pca$x[,1:4]
#plt <- plotPCA(vsd, intgroup = "type", ntop = 500, returnData = FALSE) + geom_label(aes(label = name))
#ggsave(paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/data/QTL_covariates/atac_covariates_signac/",ct,"_pca_plt.png"))
#write.table(result, paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/data/QTL_covariates/atac_covariates_signac/",ct,"_cvrt_with_counts.txt"), sep='\t', quote=FALSE, row.names=FALSE, col.names = FALSE)
#}


#Quick PCA plot of the meta-analyzed genotypes
# Read the PCA result file generated by PLINK
#pca_result <- read.table("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/data/Meta_QTL_covariates/plink/multiome_snatac_combined_filtered_chr.ancestry.pca.eigenvec", sep=' ')
#pca_result <- pca_result[, c(1, 3, 4, 5, 6, 7)]

# Rename the columns for better readability (assuming the default column names)
#colnames(pca_result) <- c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5")

# Create a scatter plot of PC1 and PC2
plt <- ggplot(pca_result, aes(x = PC1, y = PC2)) +
  geom_point() +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PLINK PCA Results")