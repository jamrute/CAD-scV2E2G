# SMC / EC cell states and the FMC program

`FMC_characterization.Rmd` (original) defines the fibromyocyte (FMC) de-differentiation program and
validates it across external datasets: CITE-seq reference mapping, spatial transcriptomics of early and
late coronary lesions, carotid plaque stenosis, Palantir pseudotime, the in vitro serum-starve/replenish
HCASMC model, and the in vivo SMC lineage-tracing time course (ED Fig. 5).

The FMC score is the top 100 modulated-SMC genes from Wirka *et al.* 2019; place that list at
`resources/gene_sets/fmc_wirka2019_top100.txt` so `R/fmc_score.R` and the dynamic caQTL model use the same
definition. The score drives the genotype × FMC interaction term in `02_qtl_discovery/07_dynamic_caqtl`.

Still to add: the endothelial reference-mapping and EndoMT panels of ED Fig. 4e–j, if they were run from a
separate script.
