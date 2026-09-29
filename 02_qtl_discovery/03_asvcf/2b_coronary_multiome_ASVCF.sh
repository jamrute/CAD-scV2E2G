#BSUB -G compute-ccdg
#BSUB -q general
#BSUB -n 2
#BSUB -a "docker(ndatth/rasqual:v0.0.0)"
#BSUB -J "l.paullee generate_ASVCF"
#BSUB -e error/generate_peakinfo_ASVCF_CELL.err
#BSUB -o out/generate_peakinfo_ASVCF_CELL.out
#BSUB -M 100GB
#BSUB -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]'

peaks=peakinfo
export RASQUALDIR="/rasqual"
echo $RASQUALDIR

bash /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/scripts/2_create_ASVCFs/createASVCF.sh \
paired_end \
"/ccdg/Active/analysis/l.paullee/multiome_qtl/common_files/bam_file_list/CELL_bam_filelist.txt" \
/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/genotyping_combined/CELL_combined_filtered.vcf.gz \
"/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/ASVCFs/CELL.corrected.all.dose.filter.vcf.gz" atac