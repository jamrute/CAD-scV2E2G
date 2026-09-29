#BSUB -G compute-ccdg
#BSUB -q ccdg
#BSUB -n 1
#BSUB -a "docker(pegi3s/feature-counts:latest)"
#BSUB -J "l.paullee mergecounts"  
#BSUB -e error/mergecounts_SMC.err 
#BSUB -o out/mergecounts_SMC.out
#BSUB -M 100GB
#BSUB -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]'

cell="SMC"
touch /scratch1/fs1/ccdg/l.paullee/SMC_counts/NUMSTART_${cell}_total_count.txt
count_path="/ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/counts"
#cat ${count_path}/${cell}/D1127_AAACCAACACCGTTCC-1.txt | tail -n +3 | cut -f 1 > ${count_path}/${cell}_total_count.txt
cat /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/Dynamic_QTL/CB_list.txt | sed -n NUMSTART,NUMENDp | while read CB
do
cat ${count_path}/${cell}/${CB}.txt | tail -n +3 | cut -f 7 > ${count_path}/${cell}/${CB}.col7
paste -d'\t' /scratch1/fs1/ccdg/l.paullee/SMC_counts/NUMSTART_${cell}_total_count.txt ${count_path}/${cell}/${CB}.col7 > NUMSTART_tmpout && mv NUMSTART_tmpout /scratch1/fs1/ccdg/l.paullee/SMC_counts/NUMSTART_${cell}_total_count.txt
rm ${count_path}/${cell}/${CB}.col7
done


#CD into directory
for file in *_SMC_total_count.txt; do
    awk '{sub(/^\t/, ""); sub(/\t$/, ""); print}' "$file" > "${file}_processed"
done

paste $(ls -1 *processed| sort -n -k1,1) > SMC_total_count.txt

rm *_processed