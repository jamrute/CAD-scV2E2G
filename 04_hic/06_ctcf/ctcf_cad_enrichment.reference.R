#!/usr/bin/env Rscript
# =============================================================================
# 04_hic/chipseq/02_ctcf_cad_enrichment.R
# Purpose : Enrichment of CAD variants in HCASMC CTCF ChIP-seq peaks
#   - variants as 1-bp intervals (hg38); observed = # variants overlapping >= 1 peak
#   - null: chromosome-wise circular permutation of peak coordinates
#       offset ~ U{1, ..., L-1}; wraparound at chromosome end; widths and
#       chromosome assignment preserved; B = 1,000
#   - fold = (observed / N) / (mean(null) / N)
#   - empirical one-sided P = (#{null >= observed} + 1) / (B + 1)
#   - parametric: p0 = mean(null) / N; exact binomial P(X >= k | N, p0)
# Paper   : ED Fig. 10g (fold = 1.25; empirical P = 0.019; binomial P = 8e-7);
#           498 CAD variants in CTCF peaks
# Inputs  : --variants  TSV with chr, pos (hg38, 1-based)
#           --peaks     CTCF peak BED
#           --chrom-sizes  UCSC chrom.sizes (hg38)
# Status  : REFERENCE IMPLEMENTATION from Methods (complete)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--variants"), make_option("--peaks"), make_option("--chrom-sizes"),
  make_option("--n-perm", type = "integer", default = 1000L),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--out")
)))
set.seed(opt$seed)

cs <- fread(opt$`chrom-sizes`, header = FALSE, col.names = c("chr", "len"))
chr_len <- setNames(as.numeric(cs$len), cs$chr)

v <- fread(opt$variants)
v[, chr := ifelse(startsWith(as.character(chr), "chr"), as.character(chr), paste0("chr", chr))]
v <- unique(v[chr %in% names(chr_len), .(chr, pos)])
v_gr <- GRanges(v$chr, IRanges(v$pos, width = 1))
N <- length(v_gr)

pk <- fread(opt$peaks, header = FALSE, select = 1:3, col.names = c("chr", "start0", "end"))
pk <- pk[chr %in% names(chr_len)]
pk[, `:=`(start = start0 + 1L, width = end - start0)]

count_overlapping_variants <- function(peaks_dt) {
  gr <- GRanges(peaks_dt$chr, IRanges(peaks_dt$start, width = peaks_dt$width))
  sum(overlapsAny(v_gr, gr))
}

#' Circularly shift peaks within each chromosome; split any peak that wraps.
circular_shift <- function(p) {
  offs <- vapply(names(chr_len), function(ch) sample.int(chr_len[[ch]] - 1L, 1L), numeric(1))
  p <- copy(p)
  L <- chr_len[p$chr]
  p[, start := ((start - 1 + offs[chr]) %% L) + 1]
  p[, end_new := start + width - 1]
  wrap <- p$end_new > L
  if (any(wrap)) {
    tail_part <- p[wrap, .(chr, start = 1, width = end_new - L[wrap])]
    p[wrap, width := L[wrap] - start + 1]
    p <- rbind(p[, .(chr, start, width)], tail_part)
  } else p <- p[, .(chr, start, width)]
  p
}

observed <- count_overlapping_variants(pk)
null <- vapply(seq_len(opt$`n-perm`), function(i) count_overlapping_variants(circular_shift(pk)), numeric(1))

fold  <- (observed / N) / (mean(null) / N)
p_emp <- (sum(null >= observed) + 1) / (opt$`n-perm` + 1)
p0    <- mean(null) / N
p_bin <- pbinom(observed - 1, N, p0, lower.tail = FALSE)

res <- data.table(n_variants = N, observed = observed, null_mean = mean(null),
                  null_sd = sd(null), fold_enrichment = fold,
                  p_empirical = p_emp, p_binomial = p_bin, n_perm = opt$`n-perm`, seed = opt$seed)
print(res)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, opt$out, sep = "\t")
fwrite(data.table(null_overlaps = null), sub("\\.tsv$", "_null.tsv", opt$out), sep = "\t")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  p <- ggplot(data.table(null = null), aes(null)) + geom_histogram(bins = 40) +
    geom_vline(xintercept = observed, linetype = 2) +
    labs(x = "CAD variants overlapping permuted CTCF peaks", y = "Permutations",
         subtitle = sprintf("fold = %.2f; empirical P = %.3g; binomial P = %.2g", fold, p_emp, p_bin)) +
    theme_classic(base_size = 8)
  ggsave(sub("\\.tsv$", ".pdf", opt$out), p, width = 3.5, height = 2.5)
}
