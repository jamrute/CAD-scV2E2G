#Script to check for GTEx eQTL caQTL overlap
library(Qtlizer)
library(dplyr)
library(stringr)

props_vector <- c()
cts <- c()

# Add QTLizer Result

cts <- c("SMC", "Fibroblast", "Endothelium","Myeloid")

for (i in seq(1,4))
{
celltype=cts[i]
caqtls <- read.table(paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/",celltype,"/",celltype,"_caQTL_FDR_005_Perm5.txt"), sep='\t', header=TRUE)
caqtls$query <- paste0("GRCh38:",caqtls$chrom,":",caqtls$pos)
#Note that this query sadly messes up the multiallelic SNPs
query_df <- get_qtls((caqtls$query))
query_df_sub <- query_df[((query_df$source=="GTEx v8") & ((query_df$tissue=="Artery - Aorta") | (query_df$tissue=="Artery - Coronary") | (query_df$tissue=="Artery - Tibial"))),]
query_df_sub$query <- paste0("GRCh38:",query_df_sub$chr,":",query_df_sub$var_pos_hg38)

out <- merge(query_df_sub, caqtls, by='query', all.x=TRUE)

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
