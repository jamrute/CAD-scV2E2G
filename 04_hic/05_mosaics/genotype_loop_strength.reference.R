#!/usr/bin/env Rscript
# =============================================================================
# 04_hic/05_genotype_loop_strength.R
# Purpose : Genotype-dependent loop strength for a variant-anchored loop
#   (default: rs4887091, mosaic-1664 800 kb boundary loop)
#   - inherent score per sample grouped by genotype (T|T, T|C / C|T, C|C)   Fig. 4e
#   - phased allele-specific inherent score, T vs C (Mann–Whitney)          Fig. 4f
# Inputs  : --loop       chr1:s1-e1,chr2:s2-e2 (mosaic BEDPE coordinates)
#           --query-dir  per-sample query_bedpe outputs (<sample>.inherent.bedpe)
#           --phased-dir per-haplotype query outputs (<sample>.hap1/.hap2.inherent.bedpe)
#           --genotypes  TSV: sample, gt (phased, e.g. "T|C")
# Status  : REFERENCE IMPLEMENTATION — confirm score column in query_bedpe output
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--loop"), make_option("--query-dir"), make_option("--phased-dir", default = NA),
  make_option("--genotypes"), make_option("--ref", default = "T"), make_option("--alt", default = "C"),
  make_option("--score-col", type = "integer", default = 8L), make_option("--out")
)))
parse_loop <- function(s) {
  p <- strsplit(s, "[,:-]")[[1]]
  list(chr1 = p[1], s1 = as.integer(p[2]), e1 = as.integer(p[3]),
       chr2 = p[4], s2 = as.integer(p[5]), e2 = as.integer(p[6]))
}
L <- parse_loop(opt$loop)
get_score <- function(f) {
  x <- fread(f, header = FALSE)
  hit <- x[V1 == L$chr1 & V2 == L$s1 & V3 == L$e1 & V4 == L$chr2 & V5 == L$s2 & V6 == L$e2]
  if (nrow(hit) == 0) NA_real_ else hit[[opt$`score-col`]][1]
}
geno <- fread(opt$genotypes)
geno[, score := vapply(sample, function(s) get_score(file.path(opt$`query-dir`, paste0(s, ".inherent.bedpe"))), 0)]
geno[, dosage_alt := lengths(regmatches(gt, gregexpr(opt$alt, gt, fixed = TRUE)))]
geno[, genotype := factor(c(paste0(opt$ref, "|", opt$ref), "het", paste0(opt$alt, "|", opt$alt))[dosage_alt + 1],
                        levels = c(paste0(opt$ref, "|", opt$ref), "het", paste0(opt$alt, "|", opt$alt)))]
trend <- cor.test(geno$dosage_alt, geno$score, method = "spearman", exact = FALSE)
fwrite(geno, paste0(opt$out, "_by_genotype.tsv"), sep = "\t")
message("Spearman(dosage, score): rho = ", round(trend$estimate, 3), ", P = ", signif(trend$p.value, 3))

if (!is.na(opt$`phased-dir`)) {
  hap <- rbindlist(lapply(seq_len(nrow(geno)), function(i) {
    a <- strsplit(geno$gt[i], "|", fixed = TRUE)[[1]]
    data.table(sample = geno$sample[i], hap = 1:2, allele = a,
               score = c(get_score(file.path(opt$`phased-dir`, paste0(geno$sample[i], ".hap1.inherent.bedpe"))),
                         get_score(file.path(opt$`phased-dir`, paste0(geno$sample[i], ".hap2.inherent.bedpe")))))
  }))
  mw <- wilcox.test(score ~ allele, data = hap)
  fwrite(hap, paste0(opt$out, "_phased.tsv"), sep = "\t")
  message("Phased ", opt$ref, " vs ", opt$alt, " Mann–Whitney P = ", signif(mw$p.value, 3))
}
