#Script to obtain cell labels for snATAC based on the multiome reference mapping
library(Signac)
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

#Read ArchR Object
path <- file.path("/ccdg/Active/analysis/l.paullee/Miller_SC/ArchRObj/Save-Proj-postQC")
snatac_obj <- loadArchRProject(path, force = TRUE, showLogo = TRUE)
snatac_df <- getCellColData(ArchRProj = snatac_obj)   

#Read Seurat object
snatac_seurat_ref_mapped <- readRDS("/ccdg/Active/analysis/l.paullee/Miller_SC/SignacObj/final/snATAC_Turner_multiome_reference_mapped.rds")
meta <- snatac_seurat_ref_mapped@meta.data
setwd("/scratch1/fs1/ccdg/l.paullee/1_ATAC_Preprocess/")

#Subset cell type assignment, save
rownames(meta)
for (sample in unique(meta$sample)){
  sub_df <- meta[grep(paste0(sample,"_"), rownames(meta)), ]
  ids <- sub(".*_", "", rownames(sub_df))
  celltype <- sub_df$predicted.id
  df <- data.frame(ids = ids, type=celltype)
  print(sample)
  print(dim(df[df$type=="SMC",]))
  print(dim(df[df$type=="Myeloid",]))
  print(dim(df[df$type=="Endothelium",]))
  print(dim(df[df$type=="Fibroblast",]))
  write.table(df,paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/1c_snATAC_Preprocess/pseudobulk_bams/CB_annotations/1k_peaks/",sample,".txt"),col.names=F,row.names=F,sep="\t",quote=F)
}

#Ground truth annotation from Miller et al.
path <- file.path("/ccdg/Active/analysis/l.paullee/coronary_multiome/coronary_ATAC_Turner_et_al/Turner_2022_ATAC_NEW")
true_obj <- loadArchRProject(path, force = TRUE, showLogo = TRUE)
true_df <- getCellColData(ArchRProj = true_obj)
true_sub_df <- true_df[grep(paste0(sample,"#"), rownames(true_df)), ]

#Looks reasonable based on spot check.