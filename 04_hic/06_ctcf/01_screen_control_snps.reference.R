#!/usr/bin/env Rscript
# =============================================================================
# 04_hic/ctcf_allelic/01_screen_control_snps.R
# Purpose : In silico screen for endogenous allelic control SNPs in CTCF peaks
#   - SNPs in CTCF peaks with complete genotype gradient across the 4 lines
#     (hom-ref, het, hom-alt all observed)
#   - CTCF motif (JASPAR MA0139.1) log-odds score for REF and ALT: best
#     window covering the SNP, both strands
#   - positive: |ΔScore| >= 6.0 and max(score) > 8.0
#   - negative: |ΔScore| <= 0.2 and both scores > 8.0
# Status  : REFERENCE IMPLEMENTATION — confirm PWM pseudocounts / log base
# =============================================================================
suppressPackageStartupMessages({
  library(VariantAnnotation); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38); library(TFBSTools); library(JASPAR2020)
  library(data.table); library(optparse)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--vcf", help = "merged WGS VCF of the 4 HCASMC lines"),
  make_option("--peaks"), make_option("--out")
)))
genome <- BSgenome.Hsapiens.UCSC.hg38

pfm <- getMatrixByID(JASPAR2020, ID = "MA0139.1")
pwm <- Matrix(toPWM(pfm, type = "log2probratio", pseudocounts = 0.8))
w <- ncol(pwm)

peaks <- rtracklayer::import(opt$peaks)
vcf <- readVcf(opt$vcf, "hg38")
vcf <- vcf[isSNV(vcf) & overlapsAny(rowRanges(vcf), peaks)]
gt  <- geno(vcf)$GT
dos <- apply(gt, 2, function(g) ifelse(g %in% c("0/0", "0|0"), 0L,
                                ifelse(g %in% c("0/1", "1/0", "0|1", "1|0"), 1L,
                                ifelse(g %in% c("1/1", "1|1"), 2L, NA_integer_))))
dos <- matrix(dos, nrow = nrow(gt), dimnames = dimnames(gt))
gradient <- apply(dos, 1, function(d) all(c(0L, 1L, 2L) %in% d))
vcf <- vcf[gradient]; dos <- dos[gradient, , drop = FALSE]

rr  <- rowRanges(vcf)
ctx <- getSeq(genome, resize(rr, 2 * w - 1, fix = "center"))
ref <- as.character(rr$REF); alt <- vapply(rr$ALT, function(a) as.character(a[1]), "")

best_score <- function(s) {
  f <- max(PWMscoreStartingAt(pwm, s, starting.at = 1:w))
  r <- max(PWMscoreStartingAt(pwm, reverseComplement(s), starting.at = 1:w))
  max(f, r)
}
ref_seq <- ctx; alt_seq <- ctx
subseq(alt_seq, start = w, width = 1) <- DNAStringSet(alt)
sc_ref <- vapply(seq_along(ref_seq), function(i) best_score(ref_seq[[i]]), 0)
sc_alt <- vapply(seq_along(alt_seq), function(i) best_score(alt_seq[[i]]), 0)

snp_ids <- if (is.null(names(rr))) paste0(seqnames(rr), ":", start(rr)) else names(rr)
res <- data.table(snp = snp_ids, chr = as.character(seqnames(rr)), pos = start(rr),
                  ref = ref, alt = alt, score_ref = sc_ref, score_alt = sc_alt,
                  delta = sc_alt - sc_ref, dos)
res[, class := fifelse(abs(delta) >= 6.0 & pmax(score_ref, score_alt) > 8.0, "positive",
               fifelse(abs(delta) <= 0.2 & pmin(score_ref, score_alt) > 8.0, "negative", NA_character_))]
fwrite(res[order(class, -abs(delta))], opt$out, sep = "\t")
print(res[, .N, by = class])
