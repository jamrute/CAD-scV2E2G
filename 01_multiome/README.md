# Module 1 — Coronary artery Multiome atlas

Paired snRNA-seq + snATAC-seq (10x Multiome) of 126,804 nuclei from 44 coronary arteries, plus the
integrated meta-map with the UVA scATAC-seq cohort (Turner *et al.*) that underlies caQTL discovery.

Produces: Fig. 1b–f; ED Fig. 1, ED Fig. 2, ED Fig. 5.

## Order of execution

| Step | Directory / script | Status | What it does |
|---|---|---|---|
| 1.1 | `01_rna_preprocessing/sampleQC.R` | original | Per-sample ArchR arrow creation and QC (TSS > 4, nFrags > 1,000); records excluded samples |
| 1.2 | `01_rna_preprocessing/merge_samples.R` | original | Loads all 44 filtered feature matrices, merges into one Seurat object, SCTransform + Harmony |
| 1.3 | `01_rna_preprocessing/run_scrublet.ipynb` | original | Doublet detection (scrublet 0.2.3) |
| 1.4 | `01_rna_preprocessing/RNA_integration.R` | original | Integration and clustering of the merged RNA object |
| 1.5 | `01_rna_preprocessing/qc_plots.R` | original | QC summary plots (ED Fig. 1) |
| 2.1 | `02_archr_multiome/1_archR_objectConstruction.R` | original | Builds the ArchR project from 44 arrows; imports the paired GEX matrix; subsets to Seurat-retained cells |
| 2.2 | `02_archr_multiome/2_multiomic_analysis.R` | original | Iterative LSI on RNA and ATAC, combined embedding, Harmony, clustering, cell-type transfer |
| 2.3 | `02_archr_multiome/3_peak_calling_p2g.R` | original | Pseudobulk coverages, MACS2 peak calling, marker peaks, motif enrichment, peak-to-gene links |
| 2.4 | `02_archr_multiome/getPeakSet.R`, `peaks_heatmap.R`, `coAcc.R`, `gpc_analysis.R`, `bigwig.R`, `visualization.R` | original | Union peak set export, marker-peak heatmap (ED Fig. 1e), co-accessibility, gene score/expression correlation, bigwig tracks, atlas UMAPs |
| 3 | `03_snatac_pseudobulk/` | original | Reference-mapped labels for the UVA snATAC cohort and per-sample, per-cell-type pseudobulk BAMs (sinto) used by RASQUAL |
| 4 | `04_scE2G/` | external | Enhancer–gene links with [scE2G](https://github.com/EngreitzLab/scE2G) (Multiome model, 5 Mb window, threshold 0.171) |
| 5 | `05_ldsc/` | original | Cell-type marker peaks → LDSC cell-type-specific partitioned heritability (Fig. 1f) |
| 6 | `06_cell_states/FMC_characterization.Rmd` | original | SMC/EC substates, the FMC de-differentiation program and its orthogonal validation (ED Fig. 5) |

## Notes

- `03_snatac_pseudobulk/ID_map.csv` maps internal specimen IDs to sequencing batches; it is needed to
  reproduce the covariate tables. Confirm before release that publishing these internal IDs is acceptable
  under the dbGaP consent (they are not the dbGaP identifiers).
- The 10x ARC ATAC↔GEX barcode translation list (`atac_gex_barcodes.txt`, 24 MB) is **not** committed.
  It ships with CellRanger ARC at `lib/python/atac/barcodes/737K-arc-v1.txt`.
- `full_bamlist.txt` (absolute paths to per-sample pseudobulk BAMs on the WashU cluster) is not committed;
  regenerate it locally with `find <pseudobulk_bams> -name '*.bam'`.
