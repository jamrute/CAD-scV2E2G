# Figure map

Panel → script. `reference` marks a panel whose original script has not yet been supplied (see
`docs/STATUS.md`).

## Figure 1 — Coronary artery multiome atlas

| Panel | Script |
|---|---|
| 1a | Schematic (BioRender) |
| 1b | `01_multiome/02_archr_multiome/visualization.R` |
| 1c | `01_multiome/01_rna_preprocessing/qc_plots.R` |
| 1d–e | `01_multiome/04_scE2G/02_summarise_scE2G.reference.R` |
| 1f | `01_multiome/05_ldsc/LDSC_analysis.sh` → `plot_ldsc_cts.reference.R` |

## Figure 2 — Single-cell caQTL discovery

| Panel | Script |
|---|---|
| 2b | `02_qtl_discovery/06_multiple_testing/10_Multiple_Testing.R` |
| 2c–d | `02_qtl_discovery/06_multiple_testing/Qtlizer_GWAS.R`; tracks from `01_multiome/02_archr_multiome/bigwig.R` |
| 2e | GTEx v8 portal data |
| 2f | `02_qtl_discovery/09_annotation/seq2print_tfbs.md` (external) |
| 2g | Motif analysis (JASPAR AP-1) |
| 2h | EMSA, wet-lab |
| 2i–k | `02_qtl_discovery/07_dynamic_caqtl/` (original script pending; `reference/pme_caqtl.reference.R`) |

## Figure 3 — scV2E2G map and Perturb-seq

| Panel | Script |
|---|---|
| 3a | `02_qtl_discovery/08_scV2E2G_map/scV2E2G_map.reference.R` (reference) |
| 3b | Track plot; scE2G links from `01_multiome/04_scE2G` |
| 3c | `03_perturbseq/03_crispri_validation/01_qpcr_ddct.reference.R` (reference) |
| 3d | Schematic |
| 3e | `03_perturbseq/03_crispri_validation/05_amotl2_proliferation.reference.R` (reference; Incucyte export needed) |
| 3f | `03_perturbseq/02_cnmf/04_annotate_programs.reference.R` |
| 3g–h | `03_perturbseq/01_preprocessing_and_mixscape/arrayed_perturbseq.Rmd` |

## Figure 4 — Disease chromatin network

| Panel | Script |
|---|---|
| 4a | Schematic + `04_hic/01_nfcore_hic/` |
| 4b–c | `04_hic/05_mosaics/` (**missing** — AQuA mosaic construction) |
| 4d–e | `04_hic/05_mosaics/genotype_loop_strength.reference.R` (inherent score implementation needed) |
| 4f | `04_hic/05_mosaics/phased_hic.md` (phased loop quantification needed) |
| 4g | `04_hic/06_ctcf/h3k27ac_chipseq.reference.sh` |
| 4h | CTCF motif (MA0139.2) |
| 4i | `04_hic/06_ctcf/03_chip_qpcr_percent_input.reference.R` (reference) |

## Figure 5 — Distal gene mapping and AMOTL2

| Panel | Script |
|---|---|
| 5a | `04_hic/05_mosaics/distal_cad_genes.reference.R` |
| 5b | `04_hic/07_hichip/04_locuszoom_rs9876658.reference.R` |
| 5c | `03_perturbseq/03_crispri_validation/01_qpcr_ddct.reference.R` |
| 5d | Immunofluorescence, `03_crispri_validation/04_if_quantification.reference.R` |
| 5e–g | `03_perturbseq/03_crispri_validation/03_amotl2_edger.reference.R` |
| 5h | `03_perturbseq/03_crispri_validation/05_amotl2_proliferation.reference.R` |

## Extended Data

| Panel | Script |
|---|---|
| ED 1a–e | `01_multiome/01_rna_preprocessing/qc_plots.R`, `02_archr_multiome/peaks_heatmap.R`, `3_peak_calling_p2g.R` |
| ED 2a–g | `01_multiome/04_scE2G/02_summarise_scE2G.reference.R` + scE2G benchmarking output |
| ED 3a–c | `02_qtl_discovery/01_genotypes/` (reference) |
| ED 3d–f | `01_multiome/03_snatac_pseudobulk/1_Re_QC.R`, `1_Re_QC_ArchR.r` |
| ED 3g–k, 3m | `02_qtl_discovery/06_multiple_testing/10_Multiple_Testing.R`, `10_Allelic_Bins.R` |
| ED 3l | seq2PRINT (external) |
| ED 4a–d | `02_qtl_discovery/06_multiple_testing/10_Multiple_Testing.R` (genotype boxplots) |
| ED 4e–j | `01_multiome/06_cell_states/FMC_characterization.Rmd` (EC panels: confirm source) |
| ED 4k | EMSA supershift, wet-lab |
| ED 5a–k | `01_multiome/06_cell_states/FMC_characterization.Rmd` |
| ED 6a–k | `02_qtl_discovery/07_dynamic_caqtl/` |
| ED 7a–e | `02_qtl_discovery/08_scV2E2G_map/` + CRISPRi qPCR |
| ED 8a–j | `03_perturbseq/01_preprocessing_and_mixscape/arrayed_perturbseq.Rmd` |
| ED 9a–i | same Rmd + `03_perturbseq/02_cnmf/` |
| ED 10a | `04_hic/05_mosaics/` (missing) |
| ED 10b–f | `04_hic/04_qc_and_comparisons/` |
| ED 10g | `04_hic/06_ctcf/ctcf_cad_enrichment.reference.R` |
| ED 10h–j | `04_hic/06_ctcf/01_screen_control_snps.reference.R`, `03_chip_qpcr_percent_input.reference.R` |
