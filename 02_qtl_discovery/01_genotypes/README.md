# Genotype processing

`GenotypeMatrix.R` (original) reads the TOPMed-imputed VCFs per chromosome and writes a dosage matrix,
keeping variants with MAF > 0.05, imputation R² > 0.9 and < 20% missing genotypes.

The remaining steps are reference implementations written from the Methods and are not the original code:

| Script | Step |
|---|---|
| `01_array_qc.reference.sh` | Illumina GSA-24v3 QC: call rate (95% common / 99% rare), HWE, sex check |
| `02_post_imputation_filter.reference.sh` | TOPMed imputation output filtering → 7,250,405 variants |
| `03_merge_turner.reference.sh` | Merge with the Turner *et al.* cohort → 5,229,397 shared variants |
| `04_genotype_pcs.reference.sh` | Genotype PCs used as caQTL covariates (ED Fig. 3a, 3c) |
| `05_rfmix_local_ancestry.reference.sh` | RFMix 2.03-r0 local ancestry, YRI (n=186) / CEU (n=183) references (ED Fig. 3b) |

Imputation was done with the TOPMed server (minimac4 1.0.2, Eagle 2.4.1). Note the Methods sentence about
"imputation R2 < 0.3 and MAF > 0.05" reads as excluding common variants; these scripts assume the intended
filter, keep R² ≥ 0.3 and MAF > 0.05.
