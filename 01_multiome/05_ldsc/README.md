# LDSC cell-type-specific partitioned heritability

Fig. 1f. Marker peaks per cell type are converted to BED annotations and tested for enrichment of
cardiometabolic trait heritability with [LDSC](https://github.com/bulik/ldsc) v1.0.1 `--h2-cts`.

| Script | Status | What it does |
|---|---|---|
| `Generate_Cell_Specific_Bed.R` | original | Marker peaks (FDR ≤ 0.01, \|log2FC\| ≥ 1) per cell type → BED |
| `Generate_SMC_marker_subtype_Bed.R` | original | Same, for SMC substates |
| `make_annot.py` | original | BED → per-chromosome LDSC annotation files |
| `LDSC_analysis.sh` | original | munge_sumstats, annotation, LD scores, `--h2-cts` across the six traits |
| `plot_ldsc_cts.reference.R` | reference | Plots the enrichment results; significance line at Bonferroni P < 8.33 × 10⁻³ (6 traits) |

Traits: CAD (Aragam *et al.*, GCST90132314), myocardial infarction, systolic and diastolic blood pressure,
stroke, aneurysm, type 2 diabetes. Summary statistics are lifted to the LDSC hg19 reference before munging.
