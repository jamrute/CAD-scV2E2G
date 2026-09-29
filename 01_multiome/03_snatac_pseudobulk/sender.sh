#Job Submission Script: Loops over all bam files.
for peaks in 1k_peaks
do
ls $HOME/ccdg/analysis/l.paullee/coronary_multiome_Final/1c_snATAC_Preprocess/pseudobulk_bams/CB_annotations/$peaks | while read i
do
id=$(echo $i | cut -d. -f1)
echo $id
mkdir -p $HOME/ccdg/analysis/l.paullee/coronary_multiome_Final/1c_snATAC_Preprocess/pseudobulk_bams/out_atac/$peaks/$id
ls $HOME/ccdg/analysis/l.paullee/Miller_SC/cellranger_outputs/$id/*/*.bam | cut -d'/' -f 5- | while read f
do
echo $f
sed -e "s|filelocation|$f|g" -e "s|sampleid|$id|g" -e "s|type|atac|g" -e "s|peaks|$peaks|g" < 2_subset_bam.sh | cat
sed -e "s|filelocation|$f|g" -e "s|sampleid|$id|g" -e "s|type|atac|g" -e "s|peaks|$peaks|g" < 2_subset_bam.sh | bsub
done
echo " "
done
done
