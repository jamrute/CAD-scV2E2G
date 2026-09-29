# Module 4 — Coronary artery Hi-C, regulatory mosaics and CTCF

Genome-wide Hi-C in 15 coronary arteries (Arima, ^GATC and G^ANTC), processed with nf-core/hic, analysed
for loops, TADs and regulatory "mosaics", and used to connect CAD variants to distal genes.

Produces: Fig. 4, Fig. 5a–b; ED Fig. 10.

| Step | Directory | Status | What it does |
|---|---|---|---|
| 1 | `01_nfcore_hic/` | external | [nf-core/hic](https://github.com/nf-core/hic) 2.0.0 run configuration (`nfcore_hic_nextflow.config`, `conf/`); 5 kb matrices, cooltools 0.7.1 |
| 2 | `02_resolution_inference/` | original | Nextflow workflow that merges per-sample BAMs and runs HiCRes 2.0 resolution inference |
| 3 | `03_loop_calling/cooltools/` | original | `Loop_Calls.ipynb`, `Loop_Calls_Tissues_10kb.ipynb` — cooltools dot calling |
| 3 | `03_loop_calling/fithic2/` | original | FitHiC2 container definition, array jobs, HiC-Pro → FitHiC conversion, nearby-interaction merging and filtering |
| 4 | `04_qc_and_comparisons/` | original | `Coronary_QC_figures.ipynb`, `HiC_meta_comparisons.Rmd`, `TADs_Loops_Sizes_Counts.Rmd`, `intersect_plotting.R` — loop/TAD counts and sizes, power-law QC, FitHiC2-vs-HiCCUPS overlap (ED Fig. 10b–f) |
| 5 | `05_mosaics/` | **missing** | Mosaic construction with AQuA Tools, mosaic annotation, genotype-stratified loop strength, distal gene mapping. Only reference skeletons are present |
| 6 | `06_ctcf/` | mixed | CTCF ChIP-seq processing, CAD-variant/CTCF-peak circular-permutation enrichment (reference, tested), allelic control-SNP screen and masked primer design for haplotype ChIP-qPCR |
| 7 | `07_hichip/` | reference | HiCUP 0.7.4 → Juicer pre → FitHiChIP v11.0 peak-to-all at 5 kb; locuszoomr panel for rs9876658 (Fig. 5b) |

## What is missing here

The mosaic analysis is the core of Fig. 4a–f and ED Fig. 10a and no code for it was supplied. To complete
this module I need:

1. The exact **AQuA Tools** calls (`extract_bedpe`, `union_bedpe`, `cluster_bedpe`, `query_bedpe`) with
   their arguments — `05_mosaics/aqua_mosaics.reference.sh` has `TODO_ARGS` placeholders.
2. The script that produced the 5,130 mosaics / 11,494 genes and the CAD-variant mosaic enrichment
   (25 vs 14 loops, P = 1.7 × 10⁻¹⁷).
3. The inherent-normalisation implementation used for Fig. 4d–e and the phased/allele-specific loop
   quantification (Axiotl Turnkey, HiC-Pro 3.1.0) behind Fig. 4f.
4. OnTAD 1.4 calls (25 kb, nested score 1) if TADs in ED Fig. 10e–f came from OnTAD rather than the
   nf-core/hic TAD step.

The upstream nf-core/hic source tree that shipped with the scripts is not committed here; only your run
configuration is. If you modified the pipeline locally, send the diff and I will vendor the modified
modules.
