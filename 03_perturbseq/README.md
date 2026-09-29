# Module 3 — Arrayed CRISPRi Perturb-seq in HCASMCs

CRISPRi against 8 CAD loci, each targeted at the gene TSS and at the scV2E2G-linked enhancer, profiled by
Parse Evercode scRNA-seq (76,179 cells) alongside live-cell proliferation imaging.

Produces: Fig. 3d–h; ED Fig. 8, ED Fig. 9.

| Step | Script | Status | What it does |
|---|---|---|---|
| 1 | `01_preprocessing_and_mixscape/splitpipe.reference.sh` | reference | Parse `split-pipe` 1.6.3 demultiplexing to DGE_filtered |
| 2 | `01_preprocessing_and_mixscape/arrayed_perturbseq.Rmd` | original | The full downstream analysis: Seurat QC (200–7,500 genes, < 25% MT), perturbation signature (`CalcPerturbSig`, k = 20, 30 PCs), Mixscape knockdown classification, LDA, distance-to-NT in PRTB PCA space, pseudobulk DESeq2 with 3 pseudo-replicates, fgsea on Hallmark and GO BP, cross-perturbation similarity, TF activity (decoupleR + DoRothEA), FMC and stimulation scores, and the convergent gene program (RobustRankAggreg) |
| 3 | `02_cnmf/` | external | Consensus NMF with [cNMF](https://github.com/dylkot/cNMF) 1.7.1 (python workflow) |
| 4 | `03_crispri_validation/` | reference | qPCR ΔΔCt, AMOTL2 knockdown bulk RNA-seq (STAR + HTSeq + edgeR), immunofluorescence quantification, proliferation curves |

## cNMF

Programs were derived with the upstream cNMF python workflow, not with code here. Settings used:
k = 5–10 scanned, **k = 8 selected**, 100 iterations, seed 42, local density threshold 0.01,
99/800 outlier components filtered, and **6 of the 8 programs retained** and annotated (interferon
responsive, proliferating S phase, proliferating M phase, inflammation, ECM remodeling, vascular
remodeling).

`02_cnmf/01_export_counts.R` and `02_make_h5ad.py` prepare the post-Mixscape counts for cNMF,
`03_run_cnmf.sh` records the `cnmf prepare/factorize/combine/consensus` calls, and
`04_annotate_programs.reference.R` maps program gene scores to the six labels and makes the heatmaps.

## Still needed

The guide-to-target barcode map (which guide sequence corresponds to which TSS/enhancer target) and the
Incucyte export used for Fig. 3e / ED Fig. 8a are not in this repository.
