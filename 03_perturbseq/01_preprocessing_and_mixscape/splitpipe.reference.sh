#!/usr/bin/env bash
# =============================================================================
# 03_perturbseq/00_splitpipe.sh
# Purpose : Parse Biosciences split-pipe v1.6.3 — build GRCh38 reference,
#           process 8 sublibraries, combine.
# Status  : SCAFFOLD — confirm chemistry flag and sample-well list
# =============================================================================
set -euo pipefail
FASTQ_DIR=${FASTQ_DIR:?}; REF_DIR=${REF_DIR:?}; OUT=${OUT:?}
SAMPLE_LIST=${SAMPLE_LIST:?split-pipe --samp_list file (sample name + wells)}

# split-pipe --mode mkref --genome_name hg38 --fasta ... --genes ... --output_dir "${REF_DIR}"
for i in $(seq 1 8); do
  split-pipe --mode all --chemistry v3 --genome_dir "${REF_DIR}" \
    --fq1 "${FASTQ_DIR}/SL${i}_R1.fastq.gz" --fq2 "${FASTQ_DIR}/SL${i}_R2.fastq.gz" \
    --samp_list "${SAMPLE_LIST}" --output_dir "${OUT}/SL${i}"
done
split-pipe --mode comb --sublibraries "${OUT}"/SL{1..8} --output_dir "${OUT}/combined"
