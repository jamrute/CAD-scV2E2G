#BSUB -G compute-ccdg
#BSUB -q general
#BSUB -n 3
#BSUB -a "docker(pegi3s/feature-counts:latest)"
#BSUB -J "l.paullee FeatureCounts"  
#BSUB -e error/FeatureCounts_samp_celltype.err 
#BSUB -o out/FeatureCounts_samp_celltype.out
#BSUB -M 150GB
#BSUB -R 'select[mem>150GB] rusage[mem=150GB] span[hosts=1]'

for i in $(seq NUMSTART 1 NUMEND)
#for i in $(seq 1 1 3)
do
echo $i
#SAMPCB=$(cat /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/CB_list.txt | sed -n ${i}p)
SAMPCB=$(cat /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/Missing.txt | sed -n ${i}p)

featureCounts -p -O -F SAF --ignoreDup \
-a /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/RASQUAL_ATAC/counts/multiome_peaks.SAF \
-o /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/counts/SMC/${SAMPCB}.txt -T 6 \
/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/CB_bams/CB_bams/SMC/${SAMPCB}.bam
done