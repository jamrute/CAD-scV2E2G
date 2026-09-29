# Data

## Generated in this study

| Dataset | Platform | n | Accession |
|---|---|---|---|
| Coronary artery Multiome (snRNA + snATAC) | 10x Multiome, CellRanger ARC 6.1, GRCh38 | 44 arteries / 126,804 nuclei | GEO GSE328800 |
| Coronary artery Hi-C | Arima, nf-core/hic 2.0.0, 5 kb | 15 arteries | GEO GSE328970 |
| Genotypes | Illumina GSA-24v3, TOPMed imputation | 88 donors (meta-map) | dbGaP phs004118.v1.p1 |
| Arrayed CRISPRi Perturb-seq | Parse Evercode WT v3, split-pipe 1.6.3 | 76,179 cells | GEO (pending) |
| HCASMC CTCF ChIP-seq | Illumina | — | GEO GSE319778 |
| H3K27ac ChIP-seq / HiChIP (HCASMC, HCAEC) | Illumina | — | GEO GSE282556, GSE282557 |
| AMOTL2 knockdown bulk RNA-seq | Illumina, STAR 2.7.0, GENCODE v38 | — | GEO GSE329169 |

## External data used

| Dataset | Use |
|---|---|
| Turner *et al.* scATAC-seq (UVA cohort) | Integrated meta-map for caQTL discovery |
| Wirka *et al.* 2019 | FMC gene program (top 100 modulated-SMC genes) |
| Schnitzler *et al.* | V2G2P gene sets, program 8 (ED Fig. 4i–j) |
| Aragam *et al.* CAD GWAS (GCST90132314) | Lead and conditional variants, LDSC |
| GTEx v8 | Coronary and tibial artery eQTLs |
| 1000 Genomes Phase 3 (YRI n=186, CEU n=183) | RFMix reference panels, LD |
| GSE131780, GSE175621, GSE188422 | Spatial, CITE-seq and carotid validation of the FMC program |

## Reference files

GRCh38 / hg38 throughout, except LDSC, which runs on the hg19 reference; marker peaks are lifted over
before annotation. Gene annotation: GENCODE v38 for bulk RNA-seq, the CellRanger ARC 2020-A reference for
Multiome.
