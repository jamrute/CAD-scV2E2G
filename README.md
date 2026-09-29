# Single-cell variant-to-enhancer-to-gene map for coronary artery disease

Code accompanying Amrute, Lee, Eres, Yamawaki *et al.*, *Nature* (2026).

This repository holds the analysis code for a single-cell variant-to-enhancer-to-gene (scV2E2G) map of
the human coronary artery: paired snRNA/snATAC Multiome profiling of 126,804 nuclei from 44 coronary
arteries, single-cell caQTL discovery in 88 donors, arrayed CRISPRi Perturb-seq in HCASMCs, and
genome-wide Hi-C in 15 coronary arteries.

## Repository organization

| Module | Contents | Main figures |
|---|---|---|
| [`01_multiome/`](01_multiome) | Multiome QC, ArchR/Seurat integration, peak calling, scE2G enhancer–gene links, LDSC heritability, SMC/EC cell states and the FMC program | Fig. 1; ED 1, 2, 5 |
| [`02_qtl_discovery/`](02_qtl_discovery) | Genotype processing, pseudobulk RASQUAL caQTLs, multiple-testing calibration, single-cell dynamic (FMC-interaction) caQTLs, scV2E2G map | Fig. 2, 3a; ED 3, 4, 6, 7 |
| [`03_perturbseq/`](03_perturbseq) | Arrayed CRISPRi Perturb-seq: Mixscape, pseudobulk DE, fgsea, cNMF programs, TF activity, convergent program, CRISPRi validation | Fig. 3; ED 8, 9 |
| [`04_hic/`](04_hic) | nf-core/hic processing, resolution inference, loop calling, QC, regulatory mosaics, CTCF allelic analyses, HiChIP | Fig. 4, 5; ED 10 |

Shared infrastructure:

```
config/      paths.yaml + params.yaml  (all thresholds used in the paper, in one place)
R/           small helpers (config loader, FMC score, plot theme)
envs/        conda environment files, one per software stack
resources/   primer/probe tables, gene sets, sample sheets
docs/        STATUS.md, DATA.md, PUBLISHING.md
figures/     FIGURE_MAP.md — every panel mapped to the script that makes it
tests/       simulation tests for the statistical reference implementations
```

## How to read this repository

Every script carries one of three labels, listed per file in [`docs/STATUS.md`](docs/STATUS.md):

- **original** — the script as run for the paper, unmodified. Paths are the authors' cluster paths and are
  kept verbatim so the record matches what was executed.
- **`*.reference.*`** — a clean implementation written from the Methods, for steps whose original script is
  not yet in this repository. These reproduce the described procedure but are **not** the code that produced
  the published numbers; they are replaced as original scripts arrive.
- **external** — analysis run with a published workflow; the directory README records the version, config and
  parameters rather than duplicating upstream code
  ([scE2G](https://github.com/EngreitzLab/scE2G), [cNMF](https://github.com/dylkot/cNMF),
  [nf-core/hic](https://github.com/nf-core/hic), [RASQUAL](https://github.com/natsuhiko/rasqual),
  [LDSC](https://github.com/bulik/ldsc), [FitHiChIP](https://github.com/ay-lab/FitHiChIP)).

Original scripts are left exactly as run rather than edited, so the committed code matches what was
executed for the paper.

## Getting started

```bash
git clone https://github.com/jamrute/CAD-scV2E2G.git
cd CAD-scV2E2G
cp config/paths.yaml config/paths.local.yaml   # edit to point at your data
export CAD_CONFIG=$PWD/config/paths.local.yaml
conda env create -f envs/r-multiome.yml          # or whichever stack you need
```

`config/params.yaml` holds every analytic threshold in the paper (QC cut-offs, FDR levels, window sizes,
model formulas), so a parameter can be checked or changed in one place.

## Tests

`bash tests/run_tests.sh` checks, on simulated data, that the dynamic caQTL mixed model recovers known
genotype and genotype×FMC effects, that the RASQUAL genome-wide empirical FDR is calibrated, that the CTCF
circular-permutation test detects planted enrichment, and that the scV2E2G map counts are correct.

## Data availability

| Data | Accession |
|---|---|
| Coronary artery Multiome (snRNA + snATAC) | GEO GSE328800 |
| Coronary artery Hi-C | GEO GSE328970 |
| Genotypes | dbGaP phs004118.v1.p1 |
| Arrayed Perturb-seq | GEO GSE347104 |
| HCASMC CTCF ChIP-seq | GEO GSE319778 |
| H3K27ac ChIP-seq and HiChIP | GEO GSE282556, GSE282557 |
| AMOTL2 knockdown bulk RNA-seq | GEO GSE329169 |

External datasets used: GSE131780, GSE175621, GSE188422, GTEx v8, 1000 Genomes Phase 3.

## Citation

See [`CITATION.cff`](CITATION.cff).

Correspondence: junedh.amrute@mdc-berlin.de and jamrute@wustl.edu

## License

MIT (see [`LICENSE`](LICENSE)). Original scripts remain the work of their authors.
