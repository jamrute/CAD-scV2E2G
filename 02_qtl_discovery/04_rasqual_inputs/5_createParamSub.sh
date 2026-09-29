#BSUB -G compute-ccdg
#BSUB -q ccdg
#BSUB -n 5
#BSUB -a "docker(ndatth/rasqual:v0.0.0)"
#BSUB -J "l.paullee generate_param"
#BSUB -e error/generate_param.err
#BSUB -o out/generate_param.out
#BSUB -R 'select[mem>200GB] rusage[mem=200GB] span[hosts=1]'
#BSUB -M 200GB

export RASQUALDIR="/rasqual"
echo $RASQUALDIR

R --vanilla --quiet < 5_create_params/5_createParam.R > log
