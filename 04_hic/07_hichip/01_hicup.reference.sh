#!/usr/bin/env bash
# =============================================================================
# 04_hic/hichip/01_hicup.sh — HiCUP v0.7.4 on HiChIP reads (Arima: GATC + GANTC)
# Status  : REFERENCE IMPLEMENTATION — confirm HiCUP config values
# =============================================================================
set -euo pipefail
FASTA=${FASTA:?hg38 fasta}; BT2_IDX=${BT2_IDX:?}; FASTQ_DIR=${FASTQ_DIR:?}; OUT=${OUT:?}
THREADS=${THREADS:-16}
mkdir -p "${OUT}/digest"
( cd "${OUT}/digest" && hicup_digester --genome hg38 --re1 '^GATC,DpnII:G^ANTC,Arima' "${FASTA}" )
DIGEST=$(ls "${OUT}"/digest/Digest_hg38_*.txt | head -1)

for r1 in "${FASTQ_DIR}"/*_R1.fastq.gz; do
  s=$(basename "${r1}" _R1.fastq.gz); r2=${r1/_R1/_R2}
  mkdir -p "${OUT}/${s}"
  hicup --bowtie2 "$(command -v bowtie2)" --index "${BT2_IDX}" --digest "${DIGEST}" \
        --outdir "${OUT}/${s}" --threads "${THREADS}" --zip "${r1}" "${r2}"
done
