#Script to obtain cell labels for snATAC based on the multiome reference mapping
library(Signac)
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

#Read meta dataframe with celltype predictions.
meta <- read.table(paste0("/ccdg/Active/analysis/l.paullee/Miller_SC/ArchRObj/integrated_labels.tsv"),sep="\t", head=T)

#Subset cell type assignment, save
rownames(meta)
for (sample in unique(meta$sample)){
  print(sample)
  sub_df <- meta[grep(paste0("\\b",sample,"\\b"), meta$sample), ]
  ids <- sub_df$sample_cell
  celltype <- sub_df$predicted.id
  df <- data.frame(ids = ids, type=celltype)
  smc_count = dim(df[df$type=="SMC",])
  myelo_count=dim(df[df$type=="Myeloid",])
  endo_count=dim(df[df$type=="Endothelium",])
  fibro_count=dim(df[df$type=="Fibroblast",])
  write.table(df,paste0("/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/1c_snATAC_Preprocess/pseudobulk_bams/CB_annotations/1k_peaks/",sample,".txt"),col.names=F,row.names=F,sep="\t",quote=F)
}

#Ground truth annotation from Miller et al.
path <- file.path("/ccdg/Active/analysis/l.paullee/coronary_multiome/coronary_ATAC_Turner_et_al/Turner_2022_ATAC_NEW")
true_obj <- loadArchRProject(path, force = TRUE, showLogo = TRUE)
true_df <- getCellColData(ArchRProj = true_obj)
true_sub_df <- true_df[grep(paste0(sample,"#"), rownames(true_df)), ]

#Looks reasonable based on spot check.