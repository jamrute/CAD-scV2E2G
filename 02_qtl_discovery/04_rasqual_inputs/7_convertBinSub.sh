#BSUB -G compute-ccdg
#BSUB -q ccdg
#BSUB -a "docker(ndatth/rasqual:v0.0.0)"
#BSUB -J "l.paullee generate_ASVCF"
#BSUB -e error/convert_bin.err
#BSUB -o out/convert_bin.out
#BSUB -M 100GB
#BSUB -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]'

#Script to convert to Binary
#Note that RASQUAL txt2bin assumes that both y and k have the first column as names
#So I rewrote txt2bin.R so that it takes read count files and covariate files with no first columns.
export RASQUALDIR="/rasqual"
echo $RASQUALDIR

for celltype in SMC Fibroblast Endothelium Myeloid
do
R --vanilla --quiet --args /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/${celltype}_count_mtx.txt /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/${celltype}_offset.txt  /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/covariates/${celltype}_cvrts_no_counts.txt < 7_convert_binary/7_txt2bin.R > log
R --vanilla --quiet --args /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/${celltype}_count_mtx.txt /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/${celltype}_offset.txt  /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/covariates/${celltype}_cvrts_counts.txt < 7_convert_binary/7_txt2bin.R > log
done