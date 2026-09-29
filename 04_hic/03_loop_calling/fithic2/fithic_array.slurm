#!/bin/bash
#SBATCH --job-name=fithic_array
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-65
# ^ 22 samples * 3 resolutions = 66 tasks (0..65)
#   %12 limits concurrency; tune to taste / quota. (I eliminated this, "%12" came after 0-62 above, but I am happy to have all run at once)

set -euo pipefail

# ---------- config ----------
BASE_ROOT="<path>/FULL_SAMPLES"
IMAGE="<path>/fithic.sif"

# samples in the order you listed
SAMPLES=(SAMPLE_01 SAMPLE_02 SAMPLE_03 SAMPLE_04 SAMPLE_05 SAMPLE_06 SAMPLE_07 SAMPLE_08 SAMPLE_09 SAMPLE_10 SAMPLE_11 SAMPLE_12 SAMPLE_13 SAMPLE_14 SAMPLE_15 SAMPLE_16 SAMPLE_17 SAMPLE_18 SAMPLE_19 SAMPLE_20 SAMPLE_21 SAMPLE_22)
RESOLUTIONS=(10000 5000 2000)

# FitHiC params
PASSES=2
LOWER_MULT=2          # lower bound = 2 * RES
UPPER=2000000         # 2Mb; consider changing per res later if you want

# chromosome filtering rules
STD_CHR_REGEX='^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$|M$)'
DROP_REGEX='(_alt|_random|_fix|Un)'

# environment hygiene
export MPLCONFIGDIR="/tmp/matplotlib"
export TMPDIR="/tmp"
export APPTAINERENV_MPLCONFIGDIR="/tmp/matplotlib"
export APPTAINERENV_TMPDIR="/tmp"

# ---------- map array index -> (sample, res) ----------
n_samples=${#SAMPLES[@]}
n_res=${#RESOLUTIONS[@]}
task=${SLURM_ARRAY_TASK_ID}

sample_idx=$(( task / n_res ))
res_idx=$(( task % n_res ))

if (( sample_idx >= n_samples )); then
  echo "ERROR: sample_idx out of range (task=${task}, sample_idx=${sample_idx}, n_samples=${n_samples})"
  exit 2
fi

SAMPLE="${SAMPLES[$sample_idx]}"
RES="${RESOLUTIONS[$res_idx]}"

BASE_DIR="${BASE_ROOT}/${SAMPLE}/fithic/${RES}"

MATRIX="${BASE_DIR}/${SAMPLE}.${RES}.matrix"
BED="${BASE_DIR}/${SAMPLE}.${RES}_abs.bed"
BIAS="${BASE_DIR}/${SAMPLE}.${RES}_iced.matrix.biases"

# Output dirs
LOG_DIR="${BASE_DIR}/logs"
CONV_DIR="${BASE_DIR}/converted"
FILTER_DIR="${BASE_DIR}/filtered"
OUT_DIR="${BASE_DIR}/run_$(date +%Y%m%d)_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "${LOG_DIR}" "${CONV_DIR}" "${FILTER_DIR}" "${OUT_DIR}"

# Redirect stdout/err to per-task logs under the sample/res dir (optional but useful)
exec > >(tee -a "${LOG_DIR}/fithic.${SAMPLE}.${RES}.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}.out") \
     2> >(tee -a "${LOG_DIR}/fithic.${SAMPLE}.${RES}.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}.err" >&2)

echo "=== Job info ==="
echo "JobID:    ${SLURM_JOB_ID}"
echo "ArrayID:  ${SLURM_ARRAY_TASK_ID}"
echo "Node:     ${SLURM_NODELIST}"
echo "CPUs:     ${SLURM_CPUS_PER_TASK}"
echo "Start:    $(date)"
echo "Sample:   ${SAMPLE}"
echo "Res:      ${RES}"
echo "Base dir: ${BASE_DIR}"
echo

# ---------- input checks ----------
echo "=== Inputs ==="
for f in "${MATRIX}" "${BED}" "${BIAS}" "${IMAGE}"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing/empty required file: $f"
    exit 3
  fi
  ls -lh "$f"
done
echo

# Bind only what we need
BIND_OPTS="-B ${BASE_DIR}:/data"

# ---------- Step 1: Convert ----------
echo "=== Step 1: HiCPro2FitHiC conversion ==="
apptainer exec ${BIND_OPTS} "${IMAGE}" \
  HiCPro2FitHiC.py \
    -i "/data/$(basename "${MATRIX}")" \
    -b "/data/$(basename "${BED}")" \
    -s "/data/$(basename "${BIAS}")" \
    -o "/data/converted" \
    -r "${RES}"

echo
echo "=== Converted outputs (listing) ==="
ls -lah "${CONV_DIR}"
echo

FRAG_GZ="${CONV_DIR}/fithic.fragmentMappability.gz"
INT_GZ="${CONV_DIR}/fithic.interactionCounts.gz"
BIAS_GZ="${CONV_DIR}/fithic.biases.gz"

if [[ ! -s "${FRAG_GZ}" ]]; then FRAG_GZ="$(ls -1 ${CONV_DIR}/*frag*gz 2>/dev/null | head -n 1 || true)"; fi
if [[ ! -s "${INT_GZ}"  ]]; then INT_GZ="$(ls -1 ${CONV_DIR}/*inter*gz 2>/dev/null | head -n 1 || true)"; fi
if [[ ! -s "${BIAS_GZ}" ]]; then BIAS_GZ="$(ls -1 ${CONV_DIR}/*bias*gz 2>/dev/null | head -n 1 || true)"; fi

if [[ -z "${FRAG_GZ}" || -z "${INT_GZ}" || -z "${BIAS_GZ}" || ! -s "${FRAG_GZ}" || ! -s "${INT_GZ}" || ! -s "${BIAS_GZ}" ]]; then
  echo "ERROR: Could not find expected gz outputs in ${CONV_DIR}."
  echo "FRAG_GZ=${FRAG_GZ}"
  echo "INT_GZ=${INT_GZ}"
  echo "BIAS_GZ=${BIAS_GZ}"
  exit 4
fi

echo "Using fragment file:    ${FRAG_GZ}"
echo "Using interaction file: ${INT_GZ}"
echo "Using biases file:      ${BIAS_GZ}"
echo

# ---------- Step 1b: Prefix chr + filter ----------
echo "=== Step 1b: chr-prefix + standard chromosome filter ==="

FRAG_FILT="${FILTER_DIR}/fithic.fragmentMappability.filtered.gz"
INT_FILT="${FILTER_DIR}/fithic.interactionCounts.filtered.gz"
BIAS_FILT="${FILTER_DIR}/fithic.biases.filtered.gz"

gunzip -c "${FRAG_GZ}" \
  | awk 'BEGIN{OFS="\t"} { $1="chr"$1; print }' \
  | awk -v re="${STD_CHR_REGEX}" '$1 ~ re' \
  | grep -vE "${DROP_REGEX}" \
  | gzip -c > "${FRAG_FILT}"

gunzip -c "${BIAS_GZ}" \
  | awk 'BEGIN{OFS="\t"} { $1="chr"$1; print }' \
  | awk -v re="${STD_CHR_REGEX}" '$1 ~ re' \
  | grep -vE "${DROP_REGEX}" \
  | gzip -c > "${BIAS_FILT}"

gunzip -c "${INT_GZ}" \
  | awk 'BEGIN{OFS="\t"} { $1="chr"$1; $3="chr"$3; print }' \
  | awk -v re="${STD_CHR_REGEX}" '$1 ~ re && $3 ~ re' \
  | grep -vE "${DROP_REGEX}" \
  | gzip -c > "${INT_FILT}"

echo "Filtered fragment file:    ${FRAG_FILT}"
echo "Filtered interaction file: ${INT_FILT}"
echo "Filtered biases file:      ${BIAS_FILT}"
ls -lh "${FRAG_FILT}" "${INT_FILT}" "${BIAS_FILT}"
echo

echo "=== Sanity: frac of bins with bias in [0.5,2] (ignoring -1) ==="
gunzip -c "${BIAS_FILT}" | awk '
BEGIN{n=0;ok=0}
{ b=$3; if(b==-1) next; n++; if(b>=0.5 && b<=2) ok++ }
END{ printf "valid_n=%d ok=%d frac_ok=%.6f\n", n, ok, (n?ok/n:0) }'
echo

# ---------- Step 2: FitHiC ----------
LOWER=$(( LOWER_MULT * RES ))

echo "=== Step 2: Run FitHiC ==="
echo "Params: RES=${RES} LOWER=${LOWER} UPPER=${UPPER} PASSES=${PASSES}"
echo "OUT_DIR=${OUT_DIR}"
echo

apptainer exec ${BIND_OPTS} "${IMAGE}" \
  fithic \
    --visual \
    -f "/data/filtered/$(basename "${FRAG_FILT}")" \
    -i "/data/filtered/$(basename "${INT_FILT}")" \
    -t "/data/filtered/$(basename "${BIAS_FILT}")" \
    -r "${RES}" \
    -o "/data/$(basename "${OUT_DIR}")" \
    -L "${LOWER}" \
    -U "${UPPER}" \
    -p "${PASSES}" \
    -tL 0.5 \
    -tU 2

echo
echo "=== FitHiC outputs (listing) ==="
ls -lah "${OUT_DIR}"
echo

echo "=== Count q<=0.05 (pass2) ==="
if [[ -s "${OUT_DIR}/FitHiC.spline_pass2.res${RES}.significances.txt.gz" ]]; then
  zcat "${OUT_DIR}/FitHiC.spline_pass2.res${RES}.significances.txt.gz" \
    | awk 'NR>1 && $7<=0.05 {c++} END{print c+0}'
else
  echo "WARNING: pass2 significances file missing"
fi

echo
echo "=== Done ==="
echo "End: $(date)"
