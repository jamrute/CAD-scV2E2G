# Module 2 — Single-cell caQTL discovery and the scV2E2G map

Pseudobulk caQTLs in SMC, endothelium, fibroblast and myeloid cells across 88 donors (245,562 nuclei),
single-cell dynamic caQTLs that depend on SMC state, and the integration of CAD variants, enhancers,
scE2G links and caQTLs into the scV2E2G map.

Produces: Fig. 2, Fig. 3a; ED Fig. 3, 4, 6, 7.

## Order of execution

| Step | Directory | Status | What it does |
|---|---|---|---|
| 1 | `01_genotypes/` | mixed | `GenotypeMatrix.R` (original) converts TOPMed-imputed VCFs to a dosage matrix (MAF > 0.05, R² > 0.9, ≤20% missing). The array QC, post-imputation filtering, cohort merge, genotype PC and RFMix local-ancestry steps are `*.reference.sh` implementations from the Methods — replace with the originals |
| 2 | `02_pseudobulk_bams/` | original | Per-sample, per-cell-type BAM subsetting and indexing for RASQUAL |
| 3 | `03_asvcf/` | original | Allele-specific VCFs (`createASVCF`), sample subsetting, chromosome renaming |
| 4 | `04_rasqual_inputs/` | original | Offsets (`4_makeOffsets.R`), feature/SNP parameter files (`5_createParam*.R`, ±10 kb and ±100 kb windows), covariate matrix (`6_generate_cvrt_matrix.R`: age, sex, site, 4 genotype PCs), text→binary conversion (`7_txt2bin.R`) |
| 5 | `05_rasqual_run/` | original | RASQUAL v1.1 jobs: full model, population-only, allele-specific-only, each with matched permutation runs |
| 6 | `06_multiple_testing/` | original | `10_Multiple_Testing.R` calibrates the genome-wide empirical FDR against 4 permutation runs (`getFDR`) and reports caQTL counts at 1/5/10% FDR; `10_Allelic_Bins.R` and `Qtlizer_GWAS.R` for allelic-imbalance bins and GWAS annotation |
| 7 | `07_dynamic_caqtl/` | mixed | Per-cell peak counts (featureCounts by cell barcode), count merging, covariate preprocessing (accessibility PCs, genotype PCs, age, sex) and the Poisson mixed-effects model with a genotype × FMC interaction. See the warning below |
| 8 | `08_scV2E2G_map/` | reference | Intersects the CAD variant superset with cell-type peaks, H3K27ac and scE2G links, and flags caQTL support (Fig. 3a counts) |
| 9 | `09_annotation/` | external | seq2PRINT / PRINT transcription-factor binding prediction at caQTL variants (Fig. 2f, ED Fig. 3l) |

## Important: the dynamic caQTL script

`07_dynamic_caqtl/4_Univariate_Poisson.R` as supplied is a working draft adapted from another study. It
still refers to GSE158769, dbGaP phs002025, `tbru_age`, `percent_mito` and a CCA covariate, loads `geno`
twice, has an empty `CB <-` assignment, and its final block references undefined objects (`test`, `x`), so
it will not run as written. `3_covariate_preprocess.R` likewise ends with a trailing block referencing an
undefined `meta`.

The Methods describe a different model: Poisson GLMM with random effects for donor and site, fixed effects
for age, sex, nFrags, TSS enrichment, 4 accessibility PCs and 4 genotype PCs, plus the FMC score and its
interaction with genotype, tested by LRT with Storey q ≤ 0.05 (712 dynamic caQTLs).
`reference/pme_caqtl.reference.R` implements that model and `reference/merge_pme.reference.R` does the
q-value step; both are covered by `tests/run_tests.sh`. **Please send the script that actually produced the
712 dynamic caQTLs** — it should replace both the draft and the reference implementation before release.

## caQTL counts and FDR

The paper reports 11,269 caQTLs. `10_Multiple_Testing.R` prints counts at 1%, 5% and 10% FDR; record in the
Supplementary Tables which threshold each table uses.
