# Script status

`original` = the script as run for the paper (committed unmodified).
`reference` = written from the Methods because the original has not been supplied; replace before release.
`external config` = run settings for an upstream pipeline.

| Directory | File | Status |
|---|---|---|
| `01_multiome/01_rna_preprocessing` | `RNA_integration.R` | original |
| `01_multiome/01_rna_preprocessing` | `merge_samples.R` | original |
| `01_multiome/01_rna_preprocessing` | `qc_plots.R` | original |
| `01_multiome/01_rna_preprocessing` | `run_scrublet.ipynb` | original |
| `01_multiome/01_rna_preprocessing` | `sampleQC.R` | original |
| `01_multiome/02_archr_multiome` | `1_archR_objectConstruction.R` | original |
| `01_multiome/02_archr_multiome` | `2_multiomic_analysis.R` | original |
| `01_multiome/02_archr_multiome` | `3_peak_calling_p2g.R` | original |
| `01_multiome/02_archr_multiome` | `bigwig.R` | original |
| `01_multiome/02_archr_multiome` | `coAcc.R` | original |
| `01_multiome/02_archr_multiome` | `getPeakSet.R` | original |
| `01_multiome/02_archr_multiome` | `gpc_analysis.R` | original |
| `01_multiome/02_archr_multiome` | `peaks_heatmap.R` | original |
| `01_multiome/02_archr_multiome` | `visualization.R` | original |
| `01_multiome/03_snatac_pseudobulk` | `1_Re_QC.R` | original |
| `01_multiome/03_snatac_pseudobulk` | `1_Re_QC_ArchR.r` | original |
| `01_multiome/03_snatac_pseudobulk` | `2_AddBaiSub.sh` | original |
| `01_multiome/03_snatac_pseudobulk` | `2_subset_bam.sh` | original |
| `01_multiome/03_snatac_pseudobulk` | `2b_bam_subset_test.sh` | original |
| `01_multiome/03_snatac_pseudobulk` | `3_rename_bam.sh` | original |
| `01_multiome/03_snatac_pseudobulk` | `ID_map.csv` | original |
| `01_multiome/03_snatac_pseudobulk` | `sender.sh` | original |
| `01_multiome/03_snatac_pseudobulk` | `sender_bai.sh` | original |
| `01_multiome/04_scE2G` | `01_export_scE2G_inputs.reference.R` | reference |
| `01_multiome/04_scE2G` | `02_summarise_scE2G.reference.R` | reference |
| `01_multiome/04_scE2G` | `config_template.yaml` | original |
| `01_multiome/05_ldsc` | `Generate_Cell_Specific_Bed.R` | original |
| `01_multiome/05_ldsc` | `Generate_SMC_marker_subtype_Bed.R` | original |
| `01_multiome/05_ldsc` | `LDSC_analysis.sh` | original |
| `01_multiome/05_ldsc` | `make_annot.py` | original |
| `01_multiome/05_ldsc` | `plot_ldsc_cts.reference.R` | reference |
| `01_multiome/06_cell_states` | `FMC_characterization.Rmd` | original |
| `02_qtl_discovery/01_genotypes` | `01_array_qc.reference.sh` | reference |
| `02_qtl_discovery/01_genotypes` | `02_post_imputation_filter.reference.sh` | reference |
| `02_qtl_discovery/01_genotypes` | `03_merge_turner.reference.sh` | reference |
| `02_qtl_discovery/01_genotypes` | `04_genotype_pcs.reference.sh` | reference |
| `02_qtl_discovery/01_genotypes` | `05_rfmix_local_ancestry.reference.sh` | reference |
| `02_qtl_discovery/01_genotypes` | `GenotypeMatrix.R` | original |
| `02_qtl_discovery/02_pseudobulk_bams` | `1a_subset_bam.sh` | original |
| `02_qtl_discovery/02_pseudobulk_bams` | `1b_AddBaiSub.sh` | original |
| `02_qtl_discovery/03_asvcf` | `2a_samplesubsetASVCF.sh` | original |
| `02_qtl_discovery/03_asvcf` | `2b_coronary_multiome_ASVCF.sh` | original |
| `02_qtl_discovery/03_asvcf` | `2c_rename_chr.sh` | original |
| `02_qtl_discovery/03_asvcf` | `createASVCF.sh` | original |
| `02_qtl_discovery/03_asvcf` | `createASVCF_new.sh` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `4_create_offsets.sh` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `4_makeOffsets.R` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `5_createParam.R` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `5_createParamSub.sh` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `5_createParam_100kb.R` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `6_generate_cvrt_matrix.R` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `6_generatecvrtSub.sh` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `7_convertBinSub.sh` | original |
| `02_qtl_discovery/04_rasqual_inputs` | `7_txt2bin.R` | original |
| `02_qtl_discovery/05_rasqual_run` | `rasqual_sub_top.sh` | original |
| `02_qtl_discovery/05_rasqual_run` | `rasqual_sub_top_as_only.sh` | original |
| `02_qtl_discovery/05_rasqual_run` | `rasqual_sub_top_perm.sh` | original |
| `02_qtl_discovery/05_rasqual_run` | `rasqual_sub_top_perm_as_only.sh` | original |
| `02_qtl_discovery/05_rasqual_run` | `rasqual_sub_top_perm_pop_only.sh` | original |
| `02_qtl_discovery/05_rasqual_run` | `rasqual_sub_top_pop_only.sh` | original |
| `02_qtl_discovery/06_multiple_testing` | `10_Allelic_Bins.R` | original |
| `02_qtl_discovery/06_multiple_testing` | `10_Multiple_Testing.R` | original |
| `02_qtl_discovery/06_multiple_testing` | `Qtlizer_GWAS.R` | original |
| `02_qtl_discovery/06_multiple_testing/reference` | `rasqual_fdr.reference.R` | reference |
| `02_qtl_discovery/07_dynamic_caqtl` | `3_covariate_preprocess.R` | original |
| `02_qtl_discovery/07_dynamic_caqtl` | `4_Univariate_Poisson.R` | original |
| `02_qtl_discovery/07_dynamic_caqtl` | `FeatureCountsCBSub.sh` | original |
| `02_qtl_discovery/07_dynamic_caqtl` | `mergecounts.sh` | original |
| `02_qtl_discovery/07_dynamic_caqtl` | `sender.sh` | original |
| `02_qtl_discovery/07_dynamic_caqtl/reference` | `merge_pme.reference.R` | reference |
| `02_qtl_discovery/07_dynamic_caqtl/reference` | `pme_caqtl.reference.R` | reference |
| `02_qtl_discovery/08_scV2E2G_map` | `scV2E2G_map.reference.R` | reference |
| `02_qtl_discovery/09_annotation` | `seq2print_tfbs.md` | original |
| `03_perturbseq/01_preprocessing_and_mixscape` | `arrayed_perturbseq.Rmd` | original |
| `03_perturbseq/01_preprocessing_and_mixscape` | `splitpipe.reference.sh` | reference |
| `03_perturbseq/02_cnmf` | `01_export_counts.reference.R` | reference |
| `03_perturbseq/02_cnmf` | `02_make_h5ad.reference.py` | reference |
| `03_perturbseq/02_cnmf` | `03_run_cnmf.reference.sh` | reference |
| `03_perturbseq/02_cnmf` | `04_annotate_programs.reference.R` | reference |
| `03_perturbseq/03_crispri_validation` | `01_qpcr_ddct.reference.R` | reference |
| `03_perturbseq/03_crispri_validation` | `02_amotl2_rnaseq_align.reference.sh` | reference |
| `03_perturbseq/03_crispri_validation` | `03_amotl2_edger.reference.R` | reference |
| `03_perturbseq/03_crispri_validation` | `04_if_quantification.reference.R` | reference |
| `03_perturbseq/03_crispri_validation` | `05_amotl2_proliferation.reference.R` | reference |
| `04_hic/01_nfcore_hic` | `nfcore_hic_nextflow.config` | external config |
| `04_hic/01_nfcore_hic/conf` | `base.config` | external config |
| `04_hic/01_nfcore_hic/conf` | `igenomes.config` | external config |
| `04_hic/01_nfcore_hic/conf` | `modules.config` | external config |
| `04_hic/01_nfcore_hic/conf` | `test.config` | external config |
| `04_hic/01_nfcore_hic/conf` | `test_full.config` | external config |
| `04_hic/02_resolution_inference` | `main.nf` | original |
| `04_hic/02_resolution_inference` | `nextflow.config` | external config |
| `04_hic/03_loop_calling/cooltools` | `Loop_Calls.ipynb` | original |
| `04_hic/03_loop_calling/cooltools` | `Loop_Calls_Tissues_10kb.ipynb` | original |
| `04_hic/03_loop_calling/fithic2` | `CombineNearbyInteraction.py` | original |
| `04_hic/03_loop_calling/fithic2` | `HiCPro2FitHiC.py` | original |
| `04_hic/03_loop_calling/fithic2` | `fithic.def` | original |
| `04_hic/03_loop_calling/fithic2` | `fithic_array.slurm` | original |
| `04_hic/03_loop_calling/fithic2` | `merge-filter-parallelized.sh` | original |
| `04_hic/03_loop_calling/fithic2` | `merge_all.sbatch` | original |
| `04_hic/04_qc_and_comparisons` | `Coronary_QC_figures.ipynb` | original |
| `04_hic/04_qc_and_comparisons` | `HiC_meta_comparisons.Rmd` | original |
| `04_hic/04_qc_and_comparisons` | `TADs_Loops_Sizes_Counts.Rmd` | original |
| `04_hic/04_qc_and_comparisons` | `intersect_plotting.R` | original |
| `04_hic/05_mosaics` | `aqua_mosaics.reference.sh` | reference |
| `04_hic/05_mosaics` | `distal_cad_genes.reference.R` | reference |
| `04_hic/05_mosaics` | `genotype_loop_strength.reference.R` | reference |
| `04_hic/05_mosaics` | `mosaic_annotation.reference.R` | reference |
| `04_hic/05_mosaics` | `phased_hic.md` | original |
| `04_hic/06_ctcf` | `01_screen_control_snps.reference.R` | reference |
| `04_hic/06_ctcf` | `02_design_masked_primers.reference.py` | reference |
| `04_hic/06_ctcf` | `03_chip_qpcr_percent_input.reference.R` | reference |
| `04_hic/06_ctcf` | `ctcf_cad_enrichment.reference.R` | reference |
| `04_hic/06_ctcf` | `h3k27ac_chipseq.reference.sh` | reference |
| `04_hic/07_hichip` | `01_hicup.reference.sh` | reference |
| `04_hic/07_hichip` | `02_juicer_pre.reference.sh` | reference |
| `04_hic/07_hichip` | `03_fithichip.reference.sh` | reference |
| `04_hic/07_hichip` | `04_locuszoom_rs9876658.reference.R` | reference |
| `04_hic/07_hichip` | `fithichip_config_template.txt` | original |

## Still missing

| Analysis | Figures | What is needed |
|---|---|---|
| Regulatory mosaics (AQuA Tools) | Fig. 4a–e, ED 10a | Exact AQuA calls and the mosaic construction / enrichment script |
| Inherent normalisation and phased loops | Fig. 4d–f | Implementation used for inherent score and allele-specific loop quantification |
| Dynamic caQTL model as run | Fig. 2i–k, Fig. 4, ED 6 | Script that produced the 712 dynamic caQTLs |
| Genotype array QC → imputation | ED 3a–c | Original array QC, imputation filtering, cohort merge, PCs, RFMix |
| scV2E2G map assembly | Fig. 3a | Original script behind the 1,360 peaks / per-cell-type enhancer and gene counts |
| CTCF ChIP-seq and haplotype ChIP-qPCR | Fig. 4g–i, ED 10g–j | Alignment/peak calling script and the %input quantification |
| EMSA and CRISPRi qPCR | Fig. 2h, 3c, 5c; ED 4k | Quantification scripts, if any beyond the plots |
| AMOTL2 knockdown RNA-seq | Fig. 5e–g | STAR/HTSeq/edgeR scripts as run |
| H3K27ac ChIP-seq and HiChIP | Fig. 5b | HiCUP/FitHiChIP run scripts as used |
| Perturb-seq guide map and Incucyte export | Fig. 3e, ED 8a | Guide-to-target table and proliferation source data |
| Seq2PRINT / PRINT | Fig. 2f, ED 3l | Commands and model version |

## Release checklist

- [ ] Replace every `*.reference.*` script with the original, or state in the README that it is a
      reimplementation.
- [ ] Fix the Code Availability URL in the manuscript (currently `https://github.com/jamrute/ CAD-scV2E2G`
      with a stray space).
- [ ] Add a Zenodo DOI for the tagged release and cite it in Code Availability.
- [ ] Confirm the sample identifiers in `ID_map.csv` may be published.
- [ ] Fill `resources/gene_sets/fmc_wirka2019_top100.txt` and the Hi-C sample sheet.
- [ ] Add `docs/session_info/` with `sessionInfo()` and `conda env export` for each environment.
- [ ] Deposit the Perturb-seq data and replace "GEO (pending)" with the accession.
- [ ] Confirm which FDR threshold each caQTL Supplementary Table uses.
