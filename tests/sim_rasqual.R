library(data.table); set.seed(3)
root <- "rq/qtl/rasqual/SMC"
mk <- function(n, n_sig, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  chisq <- c(rchisq(n_sig, 1, ncp = 40), rchisq(n - n_sig, 1))
  p <- pchisq(chisq, 1, lower.tail = FALSE)
  q <- pmin(1, p * 30)  # crude within-locus BH for ~30 SNPs/peak
  x <- data.table(sprintf("chr1:%d-%d", 1:n*1000, 1:n*1000+500), paste0("rs", 1:n), "chr1", 1:n*1000, "A","G",
                  0.3, 0.1, 1, log10(q), chisq, runif(n, .2, .8), 0.01, 0.5, 0.1, 1, 2, 30, 5, 5, 0, -100, 0, 1, 1)
  half <- split(x, rep(1:2, length.out = n))
  for (i in 1:2) fwrite(half[[i]], file.path(dir, sprintf("chunk%d.txt", i)), sep = "\t", col.names = FALSE)
}
mk(3000, 400, file.path(root, "real"))
for (i in 1:4) mk(3000, 0, file.path(root, paste0("perm", i)))
writeLines(c(paste0("root: ", getwd()), "qtl:", "  out: rq/qtl"), "paths_test.yaml")
