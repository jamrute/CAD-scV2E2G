#!/bin/bash
#SBATCH --job-name=fithic_merge
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48                 # parallel chrom workers within ONE sample/res task
#SBATCH --mem-per-cpu=2G                    # per srun); adjust upward if needed
#SBATCH --time=12:00:00
#SBATCH --array=0-65%13             # 66 total tasks; add "%13" after to run at most 13 at once (important because they also spawn per-chr so 24 more tasks, only have 320 total nodes, so 320/24=13 max concurrent)
#SBATCH --output=<path>/FULL_SAMPLES/FITHIC_2026/logs/%x.%A_%a.out
#SBATCH --error=<path>/FULL_SAMPLES/FITHIC_2026/logs/%x.%A_%a.err

set -euo pipefail
IFS=$'\n\t'

################################################################################
# USER CONFIG
################################################################################
SAMPLE_ROOT="<path>/FULL_SAMPLES"
FINAL_ROOT="<path>/FULL_SAMPLES/FITHIC_2026/LOOPS"
IMAGE="<path>/fithic.sif"
UTILITYFOLDER="/opt/fithic-utils/"
FAST_TMP_BASE="/scratch/$USER/fithic_merge_tmp"
FDR=0.05

# Your sample list (explicit, stable ordering)
SAMPLES=(SAMPLE_01 SAMPLE_02 SAMPLE_03 SAMPLE_04 SAMPLE_05 SAMPLE_06 SAMPLE_07 SAMPLE_08 SAMPLE_09 SAMPLE_10 SAMPLE_11 SAMPLE_12 SAMPLE_13 SAMPLE_14 SAMPLE_15 SAMPLE_16 SAMPLE_17 SAMPLE_18 SAMPLE_19 SAMPLE_20 SAMPLE_21 SAMPLE_22)

# Resolutions to process
RESOLUTIONS=(10000 5000 2000)

# Within-task parallelism across chromosomes
MAX_PARALLEL=${SLURM_NTASKS:-48}

################################################################################
# Derive (sample,res) for this array task
################################################################################
NS=${#SAMPLES[@]}
NR=${#RESOLUTIONS[@]}
NT=$((NS * NR))

task=${SLURM_ARRAY_TASK_ID}

if (( task < 0 || task >= NT )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${task} out of range [0,$((NT-1))]"
  exit 2
fi

sample_index=$(( task / NR ))
res_index=$(( task % NR ))
sample="${SAMPLES[$sample_index]}"
res="${RESOLUTIONS[$res_index]}"

echo "Task ${task}/${NT}: sample=${sample} res=${res}"
echo "MAX_PARALLEL=${MAX_PARALLEL} IMAGE=${IMAGE}"
echo

mkdir -p "${FINAL_ROOT}"
mkdir -p "<path>/FULL_SAMPLES/FITHIC_2026/logs"

################################################################################
# Find input significances file
################################################################################
RUN_BASE_GLOB="${SAMPLE_ROOT}/${sample}/fithic/${res}/run_*"
rundir="$(ls -1td ${RUN_BASE_GLOB} 2>/dev/null | head -n1 || true)"

if [ -z "${rundir}" ]; then
  echo "SKIP: no run dir for ${sample} ${res} at ${SAMPLE_ROOT}/${sample}/fithic/${res}"
  exit 0
fi

SIGNIF="${rundir}/FitHiC.spline_pass2.res${res}.significances.txt.gz"
if [ ! -s "${SIGNIF}" ]; then
  echo "SKIP: significance file not found: ${SIGNIF}"
  exit 0
fi

################################################################################
# Use fast tmp for intermediates
################################################################################
FAST_TMP_CURRENT="${FAST_TMP_BASE}/${sample}/${res}_$(date +%s)_${SLURM_JOB_ID}"
mkdir -p "${FAST_TMP_CURRENT}"
echo "Using fast tmp: ${FAST_TMP_CURRENT}"

cleanup() {
  rm -rf "${FAST_TMP_CURRENT}"
}
trap cleanup EXIT

cp -v "${SIGNIF}" "${FAST_TMP_CURRENT}/signif.gz"
INPUT_NAME="signif.gz"
OUTDIR_FAST="${FAST_TMP_CURRENT}/merged"
mkdir -p "${OUTDIR_FAST}"

################################################################################
# Create per-chrom jobs using merge-filter-parallelized.sh
################################################################################
echo "Running merge-filter-parallelized.sh..."
apptainer exec -B "${FAST_TMP_CURRENT}:/data" -B "${rundir}:/rundir" "${IMAGE}" \
  bash -lc "bash ${UTILITYFOLDER}/merge-filter-parallelized.sh /data/${INPUT_NAME} ${res} /data/merged ${FDR} ${UTILITYFOLDER}"

JOBFILES=( $(find "${OUTDIR_FAST}" -maxdepth 2 -type f -name 'fithic_*.job' | sort) )
if [ ${#JOBFILES[@]} -eq 0 ]; then
  echo "WARNING: no jobfiles produced for ${sample} ${res}"
  exit 1
fi
echo "Found ${#JOBFILES[@]} per-chrom jobs; launching up to ${MAX_PARALLEL} concurrently"

################################################################################
# Run per-chrom jobs in parallel on this allocation
################################################################################
idx=0
for job in "${JOBFILES[@]}"; do
  cmd="$(sed -n '1,200p' "${job}")"

  srun --exclusive -N1 -n1 bash -lc \
    "apptainer exec -B '${FAST_TMP_CURRENT}:/data' -B '${rundir}:/rundir' '${IMAGE}' bash -lc \"${cmd}\"" &

  idx=$((idx+1))
  if (( idx % MAX_PARALLEL == 0 )); then
    wait
  fi
done
wait
echo "All per-chrom jobs done for ${sample} ${res}"

################################################################################
# Aggregate post-merged outputs and write final file
################################################################################
AGG_OUT_FAST="${OUTDIR_FAST}/postmerged_all.res${res}.gz"
rm -f "${AGG_OUT_FAST}"

first=1
for f in "${OUTDIR_FAST}"/*/postmerged_fithic_*.gz; do
  [ -f "$f" ] || continue
  if (( first )); then
    zcat "$f"
    first=0
  else
    zcat "$f" | awk 'NR>1'
  fi
done | gzip -c > "${AGG_OUT_FAST}"

FINAL_NAME="${FINAL_ROOT}/${sample}.fithic.${res}.loops.gz"
cp -v "${AGG_OUT_FAST}" "${FINAL_NAME}"
echo "Wrote: ${FINAL_NAME}"
ls -lh "${FINAL_NAME}"

echo "DONE sample=${sample} res=${res}"
