library(rasqualTools)
options(scipen=999)

peak_info_path="/ccdg/Active/analysis/l.paullee/multiome_qtl/common_files/peaks_metadata/peaks.bed"
snp_info_path="/ccdg/Active/analysis/l.paullee/multiome_qtl/common_files/snp_metadata/snp_coords_meta.txt"

snps = read.table(snp_info_path)
peaks = read.table(peak_info_path)

colnames(snps) <- c("chr","pos","snp_id")
colnames(peaks) <- c("chr","exon_starts","exon_ends")

#Based on the strsplit function, the exon_starts and exon_ends need to be strings. Lmao
peaks$exon_starts <- as.character(peaks$exon_starts) 
peaks$exon_ends <- as.character(peaks$exon_ends) 

peaks$gene_id <- paste(peaks$chr, "_", peaks$exon_starts, "_", peaks$exon_ends, sep="")
peaks$strand <- "*"

#Test region within 10KB, as done by Miller group. Could consider going up to 50KB, if desired.
#snp_counts = countSnpsOverlapingExons(peaks, snps, cis_window = 10000)
snp_counts = countSnpsOverlapingExons(peaks, snps, cis_window = 10000)

#snp_counts = countSnpsOverlapingExons(test_peaks, test_snps, cis_window = 5e
out <- dplyr::select(snp_counts, chromosome_name, gene_id, exon_starts, exon_ends, feature_snp_count, cis_snp_count)
out$chromosome_name <- sub("^chr", "", out$chromosome_name)
out$region_starts <- (as.numeric(out$exon_starts))-10000
out$region_ends <- (as.numeric(out$exon_ends))+10000
out$region <- paste(out$chromosome_name, ":", out$region_starts, "-", out$region_ends, sep="")
write.table(data.frame(out), paste0("/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/params/atac_param_10kb.txt"), sep='\t', quote=FALSE, row.names=FALSE)