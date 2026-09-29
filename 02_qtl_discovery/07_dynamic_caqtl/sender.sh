#Send Featurecounts
#numbers=($(seq 1 1000 60754))
#for numstart in "${numbers[@]}"
#do
#numend=$(($numstart + 999))
#if [ "$numend" -gt 60754 ]; then
#    numend=60754
#fi
#echo "Run"
#echo $numstart
#echo $numend
#sed -e "s/NUMSTART/$numstart/g" -e "s/NUMEND/$numend/g" < 1_FeatureCounts/FeatureCountsCBSub.sh | bsub
#done
#done

#Send Mergecounts
numbers=($(seq 1 1000 60754))
for numstart in "${numbers[@]}"
do
numend=$(($numstart + 999))
if [ "$numend" -gt 60754 ]; then
    numend=60754
fi
echo "Run"
echo $numstart
echo $numend
sed -e "s/NUMSTART/$numstart/g" -e "s/NUMEND/$numend/g" < 2_merge_counts/mergecounts.sh | bsub
done
