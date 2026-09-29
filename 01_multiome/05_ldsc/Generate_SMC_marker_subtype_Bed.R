library(ArchR)


#Read Marker Peaks File
markers <- readRDS("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/coronary_multiome_atac_RDS/Save-proj1/ArchR_output/markerPeak_CellType_Seurat_RNA_15000.rds")
SMC <- data.frame(readRDS("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/coronary_multiome_atac_RDS/Save-proj1/PeakCalls/SMC-reproduciblePeaks.gr.rds"))
#markerList <- getMarkers(markers, cutOff = “FDR <= 0.1 & Log2FC >= 0.25”, returnGR = TRUE)

ct.list = c("SMC","Myeloid", "Endothelium", "Fibroblast", "Adipocyte")

for(ct in ct.list){

FC <- assays(markers)$Log2FC[ct]
FDR <- assays(markers)$FDR[ct]
#coords <- rowData(markers)[as.vector((abs(FC) >= 0.25) & (FDR <= 0.1)),]

#FDR <= 0.01 & Log2FC >= 1 like Miller group and Kellis group
coords <- rowData(markers)[as.vector((abs(FC) >= 1) & (FDR <= 0.01)),]

#bed <- data.frame(`#Chr`=gsub("chr","",coords$seqnames),
bed <- data.frame(`#Chr`=coords$seqnames,
start=as.numeric(coords$start),
end=as.numeric(coords$end),
peakid=rownames(coords), 
check.names=FALSE)

#write.table(bed,paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/LDSC/peaks/",ct,".peaks.bed"),row.names=F,sep="\t",quote=F)
write.table(bed,paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/LDSC/peaks/",ct,".peaks.FDR001.bed"),row.names=F,sep="\t",quote=F)
}