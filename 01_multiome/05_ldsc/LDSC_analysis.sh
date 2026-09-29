#Partitioned heritability tutorial
#https://hgen471.hakyimlab.org/post/2022/02/22/partition-heritability/


cd /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/LDSC

#Generate Cell Specific Bed
R Generate_Cell_Specific_Bed.R

#Munge Sumstats
python /ldsc/munge_sumstats.py \
--sumstats GCST90132314_b38.assoc \
--merge-alleles w_hm3.snplist \
--chunksize 50000 \
--out CAD_Aragam \
--a1-inc

#Make annotation
#--bed-file peaks/${ct}.peaks.bed \
for ct in SMC Myeloid Endothelium Fibroblast Adipocyte
do
for i in {1..22}
do
python make_annot.py \
--bed-file peaks/${ct}.peaks.FDR001.bed \
--bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim \
--annot-file annots/${ct}/${ct}.FDR001.chr${i}.annot.gz
done
done

#Make cell specific LD scores
for ct in SMC Myeloid Endothelium Fibroblast Adipocyte
do
for i in {1..22}
do
python /ldsc/ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} --ld-wind-cm 1 --annot annots/${ct}/${ct}.FDR001.chr${i}.annot.gz --thin-annot --out annots/${ct}/${ct}.FDR001.chr${i} --print-snps snplist.snp
done
done

#Partitioned LD score for celltypes
python /ldsc/ldsc.py \
    --h2-cts sumstats/CAD_Aragam.sumstats.gz \
    --ref-ld-chr annots/baseline_v1.2/baseline. \
    --out CAD_Multiome_CellType_LDSC_FDR001_All \
    --ref-ld-chr-cts cad_multiome_FDR001.ldcts \
    --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.