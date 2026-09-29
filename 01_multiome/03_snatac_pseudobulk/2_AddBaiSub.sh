#BSUB -G compute-ccdg
#BSUB -q ccdg 
#BSUB -n 1
#BSUB -a "docker(ndatth/rasqual:v0.0.0)"
#BSUB -J "l.paullee add bai"  
#BSUB -e addbai.err 
#BSUB -o addbai.out 
#BSUB -M 100GB
#BSUB -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]'

#Script to add bai to all bam files
#for bam in $(ls /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/1_ATAC_Preprocess/pseudobulk_bams/*/*/*/*.bam); do samtools index $bam; done
samtools index bamfile