#BSUB -G compute-ccdg
#BSUB -q ccdg
#BSUB -a "docker(pl92297/bcftools)"
#BSUB -J "l.paullee generate_ASVCF"
#BSUB -e error/subset_ASVCF_cell.err
#BSUB -o out/subset_ASVCF_cell.out
#BSUB -M 150GB
#BSUB -R 'select[mem>150GB] rusage[mem=150GB] span[hosts=1]'

# Script to subset base VCF to samples meeting QC criteria for each cell type.
ct=CELL
echo $ct
vcf_file="/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/genotyping_combined/multiome_snatac_combined_filtered.vcf.gz"
sample_list="/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/${ct}_sample_list.tsv"

# If running further sample subset, then change sample_list
subset=SUBSET
if $subset != FULL
then
    echo "Using sample list: $sample_list"
else
    echo "Using default sample list: $sample_list"
fi

# Subset based on samples passing QC
bcftools view -S $sample_list $vcf_file -O v -o "/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/genotyping_combined/${ct}_combined_filtered.vcf"

# bgzip and tabix
bgzip "/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/genotyping_combined/${ct}_combined_filtered.vcf"
tabix -p vcf "/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/genotyping_combined/${ct}_combined_filtered.vcf.gz"