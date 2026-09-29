#Script to Generate Genotype Matrix from Imputed VCFs
#https://www.biostars.org/p/241751/
#Also setting a MAF threshold of 0.05, imputation quality threshold of 0.9 (Rsq>0.9)

library(dplyr)
library(VariantAnnotation)
library(snpStats)

args = commandArgs(trailingOnly=TRUE)

row_na_cut = 0.2
fname = paste("/common/multiome_coronary_genotyping/TOPMED_imputed_VCFs/chr",chrom,".dose.vcf.gz",sep="")
tab <- TabixFile(fname, yieldSize=500000)
param <- ScanVcfParam(fixed="ALT", geno=c("GT"), info=c("MAF","R2"))
open(tab)

geno_list = c()
info = c()

while (nrow(vcf_yield <- readVcf(tab, "hg38", param=param))){
	mat <- genotypeToSnpMatrix(filtered)
	nums0 = t(as(mat$genotype, "numeric")) %>% as.data.frame
	row_rem_NA = nums0[rowSums(is.na(nums0)) < row_na_cut*ncol(nums0),]
	if(nrow(row_rem_NA)*ncol(row_rem_NA) > 0){
	geno_list = rbind( geno_list, row_rem_NA )
	}
}
close(tab)

write.table(geno_list,paste("/common/multiome_coronary_genotyping/genotype_matrix/genotype_matrix_",chrom,".txt",sep=""),quote=FALSE)