#BSUB -G compute-ccdg
#BSUB -q ccdg
#BSUB -a "docker(ndatth/rasqual:v0.0.0)"
#BSUB -J "l.paullee RASQUAL"
#BSUB -e error/RASQUAL_top_perm_NUMSTART_CELLTYPE_COVNO.err
#BSUB -o out/RASQUAL_top_perm_NUMSTART_CELLTYPE_COVNO.out
#BSUB -n 1
#BSUB -g /l.paullee/rasqual
#BSUB -M 50GB
#BSUB -R 'select[mem>50GB] rusage[mem=50GB] span[hosts=1]'

for i in $(seq NUMSTART NUMEND)
do
Y="/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/CELLTYPE_count_mtx.bin"
K="/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/CELLTYPE_offset.bin"
#X="/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/QTL_covariates/atac_covariates/cvrts_no_counts.bin"
X="/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/covariates/CELLTYPE_cvrts_counts.bin"
vcf_file="/ccdg/Active/analysis/l.paullee/multiome_qtl/genotyping/ASVCFs/CELLTYPE_combined_filtered_no_chr.vcf.gz"
#vcf_file="/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/data/VCF/meta_atac_ASVCFs/${peaks}/CELLTYPE.all.dose.filter.correct.reformat.vcf.gz"
#out_dir="/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/results/QTL_Meta_8_10_24/CELLTYPE/QTL_Perm_covCOVNO"
out_dir="/ccdg/Active/analysis/l.paullee/multiome_qtl/benchmark_qtl/pseudobulk/rasqual_qtl/results/CELLTYPE"
peaks="/ccdg/Active/analysis/l.paullee/multiome_qtl/common_files/peaks_metadata/peaks.saf"
#Param file for 10k
#param_file="/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/2b_Meta_ATAC_QTL/params/${peaks}/atac_param_10kb.txt"
param_file="/ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/params/atac_param_10kb.txt"
line_num=$i

param=$(cat $param_file | sed "${line_num}q;d")
gene_id=$(echo $param | cut -d' ' -f 2)
gene_name=$(echo $param | cut -d' ' -f 2)
region=$(echo $param | cut -d' ' -f 9)
#Feature SNP = SNP in Peaks
#Refernece SNPs = SNPs being tested
n_rsnp=$(echo $param | cut -d' ' -f 6)
n_fsnp=$(echo $param | cut -d' ' -f 5)
exon_start_positions=$(echo $param | cut -d' ' -f 3)
exon_end_positions=$(echo $param | cut -d' ' -f 4)
feat_id=$(grep $gene_id -n $peaks | cut -d":" -f1,1)
window_size=20000
n_sample=$(cat /ccdg/Active/analysis/l.paullee/multiome_qtl/qtl/data/celltype_count_mtx/CELLTYPE_sample_list.tsv | wc -l)
echo id: $gene_id 
echo name: $gene_name 
echo region: $region
echo reference snps: $n_rsnp
echo feature snps: $n_fsnp
echo feature id: $feat_id

#if grep -q $gene_name $out_dir/caQTL_CELLTYPE_perm_cov_full.txt
#then 
#echo "Found"
#else
#echo "Not Found"
for i in 1 2 3 4
do
tabix $vcf_file $region | \
rasqual \
-y $Y -k $K -x $X \
-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
-s $exon_start_positions -e $exon_end_positions \
--cis-window-size $window_size \
-f $gene_name --n_threads 2 \
--random-permutation -t -p COVNO \
--force -v | awk -v i=$i -F'\t' '{print $0 "\t" i}' >> $out_dir/caQTL_CELLTYPE_perm_cov_full.txt
done
done

#Permutation Run x 5
#for i in 1 2 3 4 5
#do
#tabix $vcf_file $region | \
#rasqual \
#-y $Y -k $K -x $X \
#-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
#-s $exon_start_positions -e $exon_end_positions \
#--cis-window-size $window_size \
#-f $gene_name --n_threads 2 \
#--random-permutation -t -p 10 \
#--force -v >> $out_dir/caQTL_CELLTYPE_perm.txt
#done
