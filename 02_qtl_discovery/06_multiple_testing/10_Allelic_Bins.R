#Script to assess allelic frequency bins for discovery

library(dplyr)
library(stringr)

# Minor allele frequency files from VCFtools
allele_freqs <- read.table("/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/genotyping_combined/genotype_allele_freq.frq", fill=TRUE, row.names = NULL)
colnames(allele_freqs) <- c('chrom', 'pos', 'n_alleles', 'n_chr', 'allele_freq1', 'allele_freq2')
allele_freqs$MAF <- pmin(as.numeric(sub(".*:", "", allele_freqs$allele_freq1)), as.numeric(sub(".*:", "", allele_freqs$allele_freq2)))
allele_freqs$coord <- paste0(allele_freqs$chrom,":", allele_freqs$pos)

# 
cts <- c("SMC", "Fibroblast", "Endothelium","Myeloid")
for (i in seq(1,4))
{
celltype=cts[i]
caqtls <- read.table(paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_FDR_005_Perm5.txt"), sep='\t', header=TRUE)
caqtls$coord <- paste0("chr",caqtls$chrom,":",caqtls$pos)

out <- merge(caqtls, allele_freqs, by='coord', all.x=TRUE)

#Note that this query sadly messes up the multiallelic SNPs
query_df <- get_qtls((caqtls$query))
query_df_sub <- query_df[((query_df$source=="GTEx v8") & ((query_df$tissue=="Artery - Aorta") | (query_df$tissue=="Artery - Coronary") | (query_df$tissue=="Artery - Tibial"))),]
query_df_sub$query <- paste0("GRCh38:",query_df_sub$chr,":",query_df_sub$var_pos_hg38)



write.table(out, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_qtlizer.txt"), sep='\t', quote=FALSE, row.names=FALSE)

proportion=length(unique(query_df_sub$query_term))/length(unique(caqtls$query))
props_vector[i] <- proportion
cts[i] <- celltype
}


# I added the Qtlizer annotations to the peak file. Let's do a quick overlap to find some interesting hits.
peaks_gwas <- read.table("/ccdg/Active/analysis/l.paullee/multiome_qtl/overlap_annotate/results/gwas_peaks.txt",header=TRUE)
SMC_qtlizer <- read.table(paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/SMC/SMC_qtlizer.txt"), sep='\t',header=TRUE)
Endothelium_qtlizer <- read.table(paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/Endothelium/Endothelium_qtlizer.txt"), sep='\t',header=TRUE)
Fibroblast_qtlizer <- read.table(paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/Fibroblast/Fibroblast_qtlizer.txt"), sep='\t',header=TRUE)
Myeloid_qtlizer <- read.table(paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/Myeloid/Myeloid_qtlizer.txt"), sep='\t',header=TRUE)

qtlizer <- list()
qtlizer[["SMC"]] <- SMC_qtlizer
qtlizer[["Endothelium"]] <- Endothelium_qtlizer
qtlizer[["Fibroblast"]] <- Fibroblast_qtlizer
qtlizer[["Myeloid"]] <- Myeloid_qtlizer

for (ct in cts){
	peaks_qtl <- peaks_gwas[peaks_gwas[paste0(ct,"_qtlizer")]==TRUE,]
	peaks_qtl_overlap <- merge(qtlizer[[ct]],peaks_qtl,by='fID')
	write.table(peaks_qtl_overlap, paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/qtlizer/",ct,"_qtlizer_gwas_peaks_overlap.tsv"), sep='\t', quote=FALSE, row.names=FALSE)
}
