library(data.table); set.seed(11)
cs <- data.table(chr = paste0("chr", 1:3), len = c(5e6, 4e6, 3e6)); fwrite(cs, "chrom.sizes", sep="\t", col.names=FALSE)
pk <- cs[, .(start0 = sort(sample(len - 1000, 400))), by = chr][, end := start0 + 400]
fwrite(pk, "ctcf.bed", sep = "\t", col.names = FALSE)
# 2000 random variants + 60 placed inside peaks (enrichment); include one variant 1bp from chr end
v <- cs[, .(pos = sample(len, 2000 * len / sum(cs$len))), by = chr]
inpk <- pk[sample(.N, 60), .(chr, pos = start0 + 200)]
fwrite(rbind(v, inpk, data.table(chr = "chr3", pos = 3e6)), "cad_vars.tsv", sep = "\t")
