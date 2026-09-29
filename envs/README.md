# Environments

The paper used several distinct software stacks (e.g. Seurat v5 for the tissue atlas but Seurat v4 for Perturb-seq; R 4.3 for the HCASMC/HiChIP validation work). Each stack gets its own environment so versions can match the Methods exactly.

| File | Used by | Key pinned versions (from Methods) |
|---|---|---|
| `r-multiome.yml` | 01_multiome, 02_qtl_discovery (R steps) | ArchR 1.0.3, Seurat 5.5.1, sctransform 0.4.3, harmony 2.0.5, Signac 1.17.1, lme4 2.0-6, MungeSumstats 1.20.0, IRanges 2.46.0, slam 0.1-56, Qtlizer 1.22.0 |
| `r-perturbseq.yml` | 03_perturbseq (R steps) | Seurat v4, DESeq2 1.52.0, fgsea 1.38.0, msigdbr 26.1.1, clusterProfiler 4.20.0, decoupleR 2.17.0, RobustRankAggreg 1.2.1, pheatmap 1.0.13 |
| `r-validation.yml` | 03_perturbseq/crispri_validation, 04_hic/hichip | R 4.3.1, edgeR 3.42.4, locuszoomr 0.3.5 |
| `py-singlecell.yml` | scrublet, cNMF, Palantir | scrublet 0.2.3, cnmf 1.7.1 |
| `genomics-cli.yml` | genotype QC, RASQUAL inputs, ChIP-seq, HiChIP | PLINK 1.9, samtools 1.13, bcftools, MACS2 2.2.9.1, BWA, Picard 3.4.0, HiCUP 0.7.4, Primer3 2.6.1 |
| `hic.yml` | 04_hic | Nextflow + nf-core/hic 2.0.0, cooltools 0.7.1, FitHiC 2.0.8, Juicer tools |
| `ldsc.yml` | 01_multiome/07_ldsc | LDSC 1.0.1 (Python 2.7 fork or Python 3 port) |

Tools installed outside conda (see each module README): Cell Ranger ARC 6.1, Cell Ranger 8.0.0 / 9.0.1, Parse split-pipe 1.6.3, RASQUAL 1.1, RFMix 2.03-r0, OnTAD 1.4, HiCRes 2.0, FitHiChIP 11.0, DFilter, AQuA Tools, scE2G, seq2PRINT.

> **TODO before release:** replace these hand-written specs with exports from the environments that actually produced the figures (`conda env export --no-builds`, and `renv::snapshot()` or `sessionInfo()` output for R). Store the R session info under `docs/session_info/`.
