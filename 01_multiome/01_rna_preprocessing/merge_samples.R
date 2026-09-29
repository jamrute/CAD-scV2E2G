library(Seurat)
library(dplyr)
library(harmony)
#############################################
#Load scRNA files
#############################################

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2013_MGI2036_TWAP-T1182/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2013_MGI2036_TWAP-T1182_cells")
s1 <- subset(scRNA, cells = cells)
s1$sample <- "MGI2013_MGI2036_TWAP-T1182"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2078_MGI2080_TWAP-T1106/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2078_MGI2080_TWAP-T1106_cells")
s2 <- subset(scRNA, cells = cells)
s2$sample <- "MGI2078_MGI2080_TWAP-T1106"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2129_MGI2308_MGI2311_TWAP-T1069L/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2129_MGI2308_MGI2311_TWAP-T1069L_cells")
s3 <- subset(scRNA, cells = cells)
s3$sample <- "MGI2129_MGI2308_MGI2311_TWAP-T1069L"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2129_MGI2308_MGI2311_TWAP-T1106R/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2129_MGI2308_MGI2311_TWAP-T1106R_cells")
s4 <- subset(scRNA, cells = cells)
s4$sample <- "MGI2129_MGI2308_MGI2311_TWAP-T1106R"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2131_MGI2132_TWAP-T1114R/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2131_MGI2132_TWAP-T1114R_cells")
s5 <- subset(scRNA, cells = cells)
s5$sample <- "MGI2131_MGI2132_TWAP-T1114R"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-D1127/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-D1127_cells")
s6 <- subset(scRNA, cells = cells)
s6$sample <- "MGI2578_MGI2581_TWKH-D1127"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-D1144/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-D1144_cells")
s7 <- subset(scRNA, cells = cells)
s7$sample <- "MGI2578_MGI2581_TWKH-D1144"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-D1148/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-D1148_cells")
s8 <- subset(scRNA, cells = cells)
s8$sample <- "MGI2578_MGI2581_TWKH-D1148"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-T1103/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-T1103_cells")
s9 <- subset(scRNA, cells = cells)
s9$sample <- "MGI2578_MGI2581_TWKH-T1103"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-T1133/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-T1133_cells")
s10 <- subset(scRNA, cells = cells)
s10$sample <- "MGI2578_MGI2581_TWKH-T1133"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-T1134/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-T1134_cells")
s11 <- subset(scRNA, cells = cells)
s11$sample <- "MGI2578_MGI2581_TWKH-T1134"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-T1135/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-T1135_cells")
s12 <- subset(scRNA, cells = cells)
s12$sample <- "MGI2578_MGI2581_TWKH-T1135"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-T1137/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-T1137_cells")
s13 <- subset(scRNA, cells = cells)
s13$sample <- "MGI2578_MGI2581_TWKH-T1137"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-T1145/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-T1145_cells")
s14 <- subset(scRNA, cells = cells)
s14$sample <- "MGI2578_MGI2581_TWKH-T1145"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-T1149/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-T1149_cells")
s15 <- subset(scRNA, cells = cells)
s15$sample <- "MGI2578_MGI2581_TWKH-T1149"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2578_MGI2581_TWKH-T1163/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2578_MGI2581_TWKH-T1163_cells")
s16 <- subset(scRNA, cells = cells)
s16$sample <- "MGI2578_MGI2581_TWKH-T1163"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2701_MGI2706_TWKH-D1173/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2701_MGI2706_TWKH-D1173_cells")
s17 <- subset(scRNA, cells = cells)
s17$sample <- "MGI2701_MGI2706_TWKH-D1173"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2701_MGI2706_TWKH-T1096/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2701_MGI2706_TWKH-T1096_cells")
s18 <- subset(scRNA, cells = cells)
s18$sample <- "MGI2701_MGI2706_TWKH-T1096"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2701_MGI2706_TWKH-T1139/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2701_MGI2706_TWKH-T1139_cells")
s19 <- subset(scRNA, cells = cells)
s19$sample <- "MGI2701_MGI2706_TWKH-T1139"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2701_MGI2706_TWKH-T1171/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2701_MGI2706_TWKH-T1171_cells")
s20 <- subset(scRNA, cells = cells)
s20$sample <- "MGI2701_MGI2706_TWKH-T1171"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2904_MGI2917_TWKH-JA_T1179/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2904_MGI2917_TWKH-JA_T1179_cells")
s21 <- subset(scRNA, cells = cells)
s21$sample <- "MGI2904_MGI2917_TWKH-JA_T1179"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2904_MGI2917_TWKH-JA_T1182/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2904_MGI2917_TWKH-JA_T1182_cells")
s22 <- subset(scRNA, cells = cells)
s22$sample <- "MGI2904_MGI2917_TWKH-JA_T1182"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2904_MGI2917_TWKH-JA_T1183/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2904_MGI2917_TWKH-JA_T1183_cells")
s23 <- subset(scRNA, cells = cells)
s23$sample <- "MGI2904_MGI2917_TWKH-JA_T1183"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2904_MGI2917_TWKH-JA_T1184/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2904_MGI2917_TWKH-JA_T1184_cells")
s24 <- subset(scRNA, cells = cells)
s24$sample <- "MGI2904_MGI2917_TWKH-JA_T1184"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2904_MGI2917_TWKH-JA_T1185/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2904_MGI2917_TWKH-JA_T1185_cells")
s25 <- subset(scRNA, cells = cells)
s25$sample <- "MGI2904_MGI2917_TWKH-JA_T1185"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2904_MGI2917_TWKH-JA_T1186/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2904_MGI2917_TWKH-JA_T1186_cells")
s26 <- subset(scRNA, cells = cells)
s26$sample <- "MGI2904_MGI2917_TWKH-JA_T1186"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2904_MGI2917_TWKH-JA_T1187/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2904_MGI2917_TWKH-JA_T1187_cells")
s27 <- subset(scRNA, cells = cells)
s27$sample <- "MGI2904_MGI2917_TWKH-JA_T1187"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI2904_MGI2917_TWKH-JA_T1189/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI2904_MGI2917_TWKH-JA_T1189_cells")
s28 <- subset(scRNA, cells = cells)
s28$sample <- "MGI2904_MGI2917_TWKH-JA_T1189"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3172_MGI3298_TWKH-T1036/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3172_MGI3298_TWKH-T1036_cells")
s29 <- subset(scRNA, cells = cells)
s29$sample <- "MGI3172_MGI3298_TWKH-T1036"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3172_MGI3298_TWKH-T1038/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3172_MGI3298_TWKH-T1038_cells")
s30 <- subset(scRNA, cells = cells)
s30$sample <- "MGI3172_MGI3298_TWKH-T1038"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3172_MGI3298_TWKH-T1040/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3172_MGI3298_TWKH-T1040_cells")
s31 <- subset(scRNA, cells = cells)
s31$sample <- "MGI3172_MGI3298_TWKH-T1040"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3172_MGI3298_TWKH-T1090/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3172_MGI3298_TWKH-T1090_cells")
s32 <- subset(scRNA, cells = cells)
s32$sample <- "MGI3172_MGI3298_TWKH-T1090"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3172_MGI3298_TWKH-T1091/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3172_MGI3298_TWKH-T1091_cells")
s33 <- subset(scRNA, cells = cells)
s33$sample <- "MGI3172_MGI3298_TWKH-T1091"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1022/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1022_cells")
s34 <- subset(scRNA, cells = cells)
s34$sample <- "MGI3339_MGI3341_TWKH-T1022"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1055/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1055_cells")
s35 <- subset(scRNA, cells = cells)
s35$sample <- "MGI3339_MGI3341_TWKH-T1055"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1123/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1123_cells")
s36 <- subset(scRNA, cells = cells)
s36$sample <- "MGI3339_MGI3341_TWKH-T1123"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1127/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1127_cells")
s37 <- subset(scRNA, cells = cells)
s37$sample <- "MGI3339_MGI3341_TWKH-T1127"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1167/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1167_cells")
s38 <- subset(scRNA, cells = cells)
s38$sample <- "MGI3339_MGI3341_TWKH-T1167"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1169/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1169_cells")
s39 <- subset(scRNA, cells = cells)
s39$sample <- "MGI3339_MGI3341_TWKH-T1169"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1190/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1190_cells")
s40 <- subset(scRNA, cells = cells)
s40$sample <- "MGI3339_MGI3341_TWKH-T1190"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1201/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1201_cells")
s41 <- subset(scRNA, cells = cells)
s41$sample <- "MGI3339_MGI3341_TWKH-T1201"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1210/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1210_cells")
s42 <- subset(scRNA, cells = cells)
s42$sample <- "MGI3339_MGI3341_TWKH-T1210"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1211/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1211_cells")
s43 <- subset(scRNA, cells = cells)
s43$sample <- "MGI3339_MGI3341_TWKH-T1211"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1214/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1214_cells")
s44 <- subset(scRNA, cells = cells)
s44$sample <- "MGI3339_MGI3341_TWKH-T1214"

matrix <- Read10X(data.dir = "/data/Junedh/Coronary_Multiomics/Multiome_Data/counts/MGI3339_MGI3341_TWKH-T1215/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "coronaryMultiome")
cells <- readRDS("/data/Junedh/Coronary_Multiomics/analysis/samplePostQC/MGI3339_MGI3341_TWKH-T1215_cells")
s45 <- subset(scRNA, cells = cells)
s45$sample <- "MGI3339_MGI3341_TWKH-T1215"


rm(matrix)
rm(scRNA)
rm(cells)

#Create merged Seurat object
merged <- merge(s1, y = c(s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,
						  s28,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45), 
                add.cell.ids = c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12","s13","s14","s15",
                				 "s16","s17","s18","s19","s20","s21","s22","s23","s24","s25","s26","s27","s28","s29","s30",
                				 "s31","s32","s33","s34","s35","s36","s37","s38","s39","s40","s41","s42","s43","s44","s45"), project = "coronaryMultiome")

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
saveRDS(merged, "/data/Junedh/Coronary_Multiomics/analysis/RNA/RNA_merged_postQC.rds")



