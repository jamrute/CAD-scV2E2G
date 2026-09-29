#BSUB -G compute-ccdg
#BSUB -q ccdg 
#BSUB -n 5
#BSUB -a "docker(pl92297/sinto:latest)"
#BSUB -J "l.paullee cell specific bams"  
#BSUB -e cell_specific_sampleid_type.err 
#BSUB -o cell_specific_sampleid_type.out 
#BSUB -M 100GB
#BSUB -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]'

#Script to generate cell type specific sample per sample using Sinto


sinto filterbarcodes -b /ccdg/Active/filelocation -p 8 -c /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/1c_snATAC_Preprocess/pseudobulk_bams/CB_annotations/peaks/sampleid.txt --barcodetag "CB" --outdir /ccdg/Active/analysis/l.paullee/coronary_multiome_Final/1c_snATAC_Preprocess/pseudobulk_bams/out_atac/peaks/sampleid/
