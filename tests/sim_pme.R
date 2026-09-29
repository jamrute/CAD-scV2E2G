suppressPackageStartupMessages({library(Matrix); library(data.table)})
set.seed(7)
nd <- 24; cpd <- 120
donors <- sprintf("D%02d", 1:nd)
dm <- data.table(donor = donors, age = rnorm(nd, 55, 10), sex = sample(c("M","F"), nd, TRUE),
                 genoPC1 = rnorm(nd), genoPC2 = rnorm(nd), genoPC3 = rnorm(nd), genoPC4 = rnorm(nd))
cm <- data.table(cell = paste0("c", 1:(nd*cpd)), donor = rep(donors, each = cpd))
cm[, site := ifelse(donor %in% donors[1:12], "WashU", "UVA")]
cm[, `:=`(nFrags = rpois(.N, 8000), TSSEnrichment = runif(.N, 5, 15), FMC_score = rnorm(.N))]
for (k in 1:4) cm[[paste0("accPC", k)]] <- rnorm(nrow(cm))
G <- matrix(sample(0:2, 3*nd, TRUE, prob = c(.4,.4,.2)), 3, dimnames = list(c("rs_static","rs_dynamic","rs_null"), donors))
re <- setNames(rnorm(nd, 0, 0.2), donors)
g <- function(v) G[v, cm$donor]
fmc <- scale(cm$FMC_score)[,1]
base <- -1 + 0.5*scale(log(cm$nFrags))[,1] + re[cm$donor]
mu <- rbind(peakA = exp(base + 0.4*g("rs_static")),
            peakB = exp(base + 0.1*g("rs_dynamic") - 0.35*g("rs_dynamic")*fmc + 0.2*fmc),
            peakC = exp(base))
cnt <- Matrix(matrix(rpois(length(mu), mu), nrow(mu), dimnames = list(rownames(mu), cm$cell)), sparse = TRUE)
saveRDS(cnt, "counts.rds"); fwrite(cm, "cell_meta.tsv", sep="\t"); fwrite(dm, "donor_meta.tsv", sep="\t")
fwrite(data.table(variant_id = rownames(G), G), "dosage.tsv", sep="\t")
fwrite(data.table(peak_id = c("peakA","peakB","peakC"), variant_id = c("rs_static","rs_dynamic","rs_null")), "pairs.tsv", sep="\t")
