# scE2G enhancer–gene links

Enhancer–gene links were generated with the Engreitz lab [scE2G](https://github.com/EngreitzLab/scE2G)
pipeline, not with code in this repository. Clone that repository at the released version and run it with
the settings recorded here.

- Model: `scE2G_multiome` (`multiome_powerlaw_v3`), the Multiome version of the model.
- Candidate region window: 5 Mb around each TSS.
- Score threshold: **0.171**.
- Inputs per cell type: the fragment file and the ATAC union peak set from `01_multiome/02_archr_multiome`,
  plus the paired RNA count matrix.
- Run per cell type (SMC, Endothelium, Fibroblast, Myeloid and the remaining atlas cell types).

`config_template.yaml` is the scE2G config skeleton with these values filled in; point its input paths at
your own outputs. `01_export_scE2G_inputs.reference.R` writes the per-cell-type fragment/peak/RNA inputs from the
ArchR project, and `02_summarise_scE2G.reference.R` collects the thresholded links into the per-cell-type link counts
and shared-enhancer correlation matrix in Fig. 1d–e and ED Fig. 2.

The scE2G-versus-scATAC power-law comparison in ED Fig. 2g comes from the pipeline's own benchmarking
output (`multiome_powerlaw_v3` vs `scATAC_powerlaw_v3`).
