#BSUB -G compute-ccdg
#BSUB -q ccdg
#BSUB -a "docker(bhuvic/singlecell:v10)"
#BSUB -J "l.paullee generate_cvrt"
#BSUB -e error/generate_cvrt.err
#BSUB -o out/generate_cvrt.out
#BSUB -M 100GB
#BSUB -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]'

export RASQUALDIR="/rasqual"
echo $RASQUALDIR

R < 6_generate_covariates/6_generate_cvrt_matrix.R