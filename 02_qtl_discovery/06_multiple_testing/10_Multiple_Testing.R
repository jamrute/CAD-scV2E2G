#12/30/23: Exploratory script to test some testing conditions for caQTL to maximize power

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(dplyr)
library(tidyverse)


getFDR <- function(q1, q0, alpha=0.1, z=NULL, subset=NULL){
	if(is.null(z)){
		a=0
		for(itr in 1:10){
			a=getFDR(q1,q0,alpha,rev(a+0:100/100^itr),subset)
		}
		a
	}else{
		if(!is.null(subset)){
			q1=q1[subset]
			q0=q0[subset]
		}
		q1=q1[!is.na(q1)]
		q0=q0[!is.na(q0)]
		x=NULL;
		for(i in z){
			x=c(x,sum(q0<i)/length(q0)/(sum(q1<i)/length(q1)))
		};
		max(c(0,z[x<alpha]),na.rm=T)
	}
}


find_caqtl <- function(celltype,subset,covs){
   original <- read.table(paste0('/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/',celltype,'/QTL_Top_cov7/caQTL_',celltype,'_original_cov_full_test.txt'), sep='\t', header=FALSE, fill=TRUE)
   colnames(original) <- c("fID", "rsID", "chrom", "pos", "ref", "alt", "freq", "hwe_x", "IA", "log10Q", "X2", "B", "Delta", "Phi", "overdisperse", "SNP_ID", "no_f_SNP", "no_t_SNP", "null_iter", "alt_iter", "rand_location", "log_likelihood", "convergence_stat", "R2_prio_post_fSNP", "R2_prio_post_rSNP")
   #Drop duplicate,skipped
   original <- original[!duplicated(original), ]
   original <- original[original$rsID != 'SKIPPED',]
   original$pval <- pchisq(original$X2,1,lower=F)
   original$BF_P <- with(original, ifelse(no_t_SNP * pval > 1, 1, no_t_SNP * pval))
   original$Q <- 10^original$log10Q
   original <- original[order(original$Q),]

   #Process Permutation Runs
   perm <- read.table(paste0('/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/',celltype,'/QTL_Perm_cov7/caQTL_',celltype,'_perm_cov_full.txt'), sep='\t', header=FALSE, fill=TRUE)
   colnames(perm) <- c("fID", "rsID", "chrom", "pos", "ref", "alt", "freq", "hwe_x", "IA", "log10Q", "X2", "B", "Delta", "Phi", "overdisperse", "SNP_ID", "no_f_SNP", "no_t_SNP", "null_iter", "alt_iter", "rand_location", "log_likelihood", "convergence_stat", "R2_prio_post_fSNP", "R2_prio_post_rSNP", "index")
   perm <- perm[!duplicated(perm), ]
   perm <- perm[perm$rsID != 'SKIPPED',]

   if(!missing(subset)) {
      original <- original[original$fID %in% subset,]
      perm <- perm[perm$fID %in% subset,]
   }

   #Add P val based on X2 statistic
   perm$pval <- pchisq(perm$X2,1,lower=F)
   perm$Q <- 10^(perm$log10Q)
   perm$BF_P <- with(perm, ifelse(no_t_SNP * pval > 1, 1, no_t_SNP * pval))
   #perm <- perm[,c("index","Q")]
   perm <- perm[order(perm$pval),]

   #perm_df <- data.frame(perm[perm$index==1,]$Q, perm[perm$index==2,]$Q, perm[perm$index==3,]$Q, perm[perm$index==4,]$Q)
   perm_df <- data.frame(perm[perm$index==1,]$BF_P, perm[perm$index==2,]$BF_P, perm[perm$index==3,]$BF_P, perm[perm$index==4,]$BF_P)
   colnames(perm_df) <- c("run1", "run2", "run3", "run4")
   perm_df <- perm_df %>% mutate(average_perm_q = rowMeans(select(.,c(1:4))))
   
   #BH_Perm_FDR
   
   FDR_010_cut <- getFDR(original$BF_P, perm_df$average_perm_q, 0.1)

   FDR_010_cut <- getFDR(original$Q, perm_df$average_perm_q, 0.1)
   FDR_005_cut <- getFDR(original$Q, perm_df$average_perm_q, 0.05)
   FDR_001_cut <- getFDR(original$Q, perm_df$average_perm_q, 0.01)
   print(FDR_010_cut)
   print(FDR_005_cut)
   print(FDR_001_cut)


   print(summary(perm$Q))

   print(paste0("Covs: ",covs))

   BH_Perm_FDR_010 <- original[original$Q < FDR_010_cut,]
   BH_Perm_FDR_005 <- original[original$Q < FDR_005_cut,]
   BH_Perm_FDR_001 <- original[original$Q < FDR_001_cut,]
   print(paste0("1% ", dim(BH_Perm_FDR_001)[[1]]))
   print(paste0("5% ", dim(BH_Perm_FDR_005)[[1]]))
   print(paste0("10% ", dim(BH_Perm_FDR_010)[[1]]))

   return(list(BH_Perm_FDR_001, BH_Perm_FDR_005, BH_Perm_FDR_010))
}

SMC <- find_caqtl("SMC",covs=7)
Endo <- find_caqtl("Endothelium",covs=7)
Myelo <- find_caqtl("Myeloid",covs=7)
Fibro <- find_caqtl("Fibroblast",covs=7)


for (ind in c(1,2,3,4))
{
	celltype <- c("SMC","Endothelium", "Myeloid", "Fibroblast")[ind]
	df <- c(SMC, Endo, Myelo, Fibro)
	print(celltype)
	write.table(df[ind*3], paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_FDR_010_Perm5.txt"), sep='\t', quote=FALSE, row.names=FALSE)
	write.table(df[(ind*3)-1], paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_FDR_005_Perm5.txt"), sep='\t', quote=FALSE, row.names=FALSE)
	write.table(df[(ind*3)-2], paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_FDR_001_Perm5.txt"), sep='\t', quote=FALSE, row.names=FALSE)
}

#Full dataframe at 5% fdr
SMC_FDR <- SMC[[2]]
Endothelium_FDR <- Endo[[2]]
Myeloid_FDR <- Myelo[[2]]
Fibroblast_FDR <- Fibro[[2]]

SMC_FDR$caQTL_Cell_Type <- "SMC"
Endothelium_FDR$caQTL_Cell_Type <- "Endothelium"
Myeloid_FDR$caQTL_Cell_Type <- "Myeloid"
Fibroblast_FDR$caQTL_Cell_Type <- "Fibroblast"
FDR_full <- rbind(SMC_FDR, Endothelium_FDR, Myeloid_FDR, Fibroblast_FDR)
write.table(FDR_full, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/All_caQTL_FDR_005_Perm5.txt"), sep='\t', quote=FALSE, row.names=FALSE)



# Multiple Testing + Extract Full Locus Results

find_caqtl <- function(celltype,subset,covs){
   original <- read.table(paste0('/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/',celltype,'/QTL_Top_cov7/caQTL_',celltype,'_original_cov_full_test.txt'), sep='\t', header=FALSE, fill=TRUE)
   colnames(original) <- c("fID", "rsID", "chrom", "pos", "ref", "alt", "freq", "hwe_x", "IA", "log10Q", "X2", "B", "Delta", "Phi", "overdisperse", "SNP_ID", "no_f_SNP", "no_t_SNP", "null_iter", "alt_iter", "rand_location", "log_likelihood", "convergence_stat", "R2_prio_post_fSNP", "R2_prio_post_rSNP")
   #Drop duplicate,skipped
   original <- original[!duplicated(original), ]
   original <- original[original$rsID != 'SKIPPED',]
   original$pval <- pchisq(original$X2,1,lower=F)
   original$Q <- 10^original$log10Q
   original <- original[order(original$Q),]

   #Process Permutation Runs
   perm <- read.table(paste0('/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/',celltype,'/QTL_Perm_cov7/caQTL_',celltype,'_perm_cov_full.txt'), sep='\t', header=FALSE, fill=TRUE)
   colnames(perm) <- c("fID", "rsID", "chrom", "pos", "ref", "alt", "freq", "hwe_x", "IA", "log10Q", "X2", "B", "Delta", "Phi", "overdisperse", "SNP_ID", "no_f_SNP", "no_t_SNP", "null_iter", "alt_iter", "rand_location", "log_likelihood", "convergence_stat", "R2_prio_post_fSNP", "R2_prio_post_rSNP", "index")
   perm <- perm[!duplicated(perm), ]
   perm <- perm[perm$rsID != 'SKIPPED',]

   if(!missing(subset)) {
      original <- original[original$fID %in% subset,]
      perm <- perm[perm$fID %in% subset,]
   }

   #Add P val based on X2 statistic
   perm$pval <- pchisq(perm$X2,1,lower=F)
   perm$Q <- 10^(perm$log10Q)
   #perm <- perm[,c("index","Q")]
   perm <- perm[order(perm$pval),]

   perm_df <- data.frame(perm[perm$index==1,]$Q, perm[perm$index==2,]$Q, perm[perm$index==3,]$Q, perm[perm$index==4,]$Q)
   colnames(perm_df) <- c("run1", "run2", "run3", "run4")
   perm_df <- perm_df %>% mutate(average_perm_q = rowMeans(select(.,c(1:4))))
   
   #BH_Perm_FDR
   FDR_010_cut <- getFDR(original$Q, perm_df$average_perm_q, 0.1)
   FDR_005_cut <- getFDR(original$Q, perm_df$average_perm_q, 0.05)
   FDR_001_cut <- getFDR(original$Q, perm_df$average_perm_q, 0.01)
   print(FDR_010_cut)
   print(FDR_005_cut)
   print(FDR_001_cut)


   print(summary(perm$Q))

   print(paste0("Covs: ",covs))


   return(list(FDR_010_cut, FDR_005_cut, FDR_001_cut))
}

SMC <- find_caqtl("SMC",covs=7)
Endo <- find_caqtl("Endothelium",covs=7)
Myelo <- find_caqtl("Myeloid",covs=7)
Fibro <- find_caqtl("Fibroblast",covs=7)

get_sig <- function(celltype,cutoffs,index){
	full_locus <- read.table(paste0('/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/',celltype,'/',celltype,'_010_full.txt'), sep='\t', header=FALSE, fill=TRUE)
	colnames(full_locus) <- c("fID", "rsID", "chrom", "pos", "ref", "alt", "freq", "hwe_x", "IA", "log10Q", "X2", "B", "Delta", "Phi", "overdisperse", "SNP_ID", "no_f_SNP", "no_t_SNP", "null_iter", "alt_iter", "rand_location", "log_likelihood", "convergence_stat", "R2_prio_post_fSNP", "R2_prio_post_rSNP")
	full_locus <- full_locus[full_locus$rsID != 'SKIPPED',]
	full_locus$pval <- pchisq(full_locus$X2,1,lower=F)
	full_locus$Q <- 10^full_locus$log10Q
	full_locus <- full_locus[order(full_locus$Q),]
	sig_locus <- full_locus[full_locus$Q< cutoffs[[index]],]
	sig_locus
}

celltype="SMC"
sig_SMC <- get_sig("SMC",SMC,1)
write.table(sig_SMC, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_full_locus_FDR_010.txt"), sep='\t', quote=FALSE, row.names=FALSE)
sig_SMC <- get_sig("SMC",SMC,2)
write.table(sig_SMC, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_full_locus_FDR_005.txt"), sep='\t', quote=FALSE, row.names=FALSE)

celltype="Endothelium"
sig_Endo <- get_sig("Endothelium",Endo,1)
write.table(sig_Endo, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_full_locus_FDR_010.txt"), sep='\t', quote=FALSE, row.names=FALSE)
sig_Endo <- get_sig("Endothelium",Endo,2)
write.table(sig_Endo, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_full_locus_FDR_005.txt"), sep='\t', quote=FALSE, row.names=FALSE)

celltype="Myeloid"
sig_Myelo <- get_sig("Myeloid",Myelo,1)
write.table(sig_Myelo, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_full_locus_FDR_010.txt"), sep='\t', quote=FALSE, row.names=FALSE)
sig_Myelo <- get_sig("Myeloid",Myelo,2)
write.table(sig_Myelo, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_full_locus_FDR_005.txt"), sep='\t', quote=FALSE, row.names=FALSE)

celltype="Fibroblast"
sig_Fibro <- get_sig("Fibroblast",Fibro,1)
write.table(sig_Fibro, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_full_locus_FDR_010.txt"), sep='\t', quote=FALSE, row.names=FALSE)
sig_Fibro <- get_sig("Fibroblast",Fibro,2)
write.table(sig_Fibro, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_full_locus_FDR_005.txt"), sep='\t', quote=FALSE, row.names=FALSE)
