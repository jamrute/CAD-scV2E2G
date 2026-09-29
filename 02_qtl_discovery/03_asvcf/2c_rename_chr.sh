#BSUB -G compute-ccdg
#BSUB -q ccdg
#BSUB -a "docker(pl92297/bcftools)"
#BSUB -J "l.paullee rename ASVCF"
#BSUB -e error/subset_ASVCF_cell.err
#BSUB -o out/subset_ASVCF_cell.out
#BSUB -M 150GB
#BSUB -R 'select[mem>150GB] rusage[mem=150GB] span[hosts=1]'

#Script to rename ASVCF chromsome ids (required for RASQUAL)
dir="/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/ASVCFs/"
for ct in SMC Myeloid Fibroblast Endothelium
do
echo $ct
zcat "${dir}/${ct}.corrected.all.dose.filter.vcf.gz" | sed 's/^chr//' | bgzip > "${dir}/${ct}_combined_filtered_no_chr.vcf.gz"
tabix -p vcf ${dir}/${ct}_combined_filtered_no_chr.vcf.gz
done

