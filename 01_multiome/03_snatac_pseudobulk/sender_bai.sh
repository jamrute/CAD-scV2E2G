cat full_bamlist.txt | while read bamfile
do 
	sed -e "s|bamfile|$bamfile|g" < 2_AddBaiSub.sh | bsub
done

