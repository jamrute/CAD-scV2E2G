#BSUB -G compute-ccdg
#BSUB -q general
#BSUB -a "docker(ndatth/rasqual:v0.0.0)"
#BSUB -J "l.paullee generate_Offset"
#BSUB -e error/generate_peaksinfo_CELLTYPE_Offset.err
#BSUB -o out/generate_peaksinfo_CELLTYPE_Offset.out
#BSUB -M 100GB
#BSUB -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]'

export RASQUALDIR="/rasqual"
echo $RASQUALDIR
count="/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx"
peakinfo="/ccdg/Active/analysis/l.paullee/multiome_qtl/common_files/peaks_metadata"
echo $count

R --vanilla --quiet --args ${count}/CELLTYPE_count_mtx.tsv ${peakinfo}/GC.txt ${count}/CELLTYPE_offset.txt < 4_create_offsets/4_makeOffsets.R > log