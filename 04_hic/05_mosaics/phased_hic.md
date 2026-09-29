# Phased (haplotype-resolved) Hi-C

Phased contact matrices were generated with sample-specific genotypes on Axiotl's Turnkey platform (tinker.axiotl.com), which runs HiC-Pro v3.1.0 in allele-specific mode on cloud compute.

TODO (Axiotl): document
* phasing input (VCF source, phasing tool),
* HiC-Pro allele-specific configuration (`ALLELE_SPECIFIC_SNP`, N-masked genome),
* how per-haplotype `.hic` files were produced and then queried with `query_bedpe` (input to `05_genotype_loop_strength.R --phased-dir`).
