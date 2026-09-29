#BSUB -G compute-ccdg
#BSUB -q general
#BSUB -n 1
#BSUB -a "docker(ndatth/rasqual:v0.0.0)"
#BSUB -J "l.paullee add bai"  
#BSUB -e addbai.err 
#BSUB -o addbai.out 
#BSUB -M 5GB
#BSUB -R 'select[mem>5GB] rusage[mem=5GB] span[hosts=1]'

#Script to add bai to all bam files
cd /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/pseudobulk_bams/sampleid/
for celltype in SMC Myeloid Endothelium Fibroblast
do
bamfile="${celltype}.bam"
mv $bamfile sampleid_${bamfile}
samtools index sampleid_$bamfile
done