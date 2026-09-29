#!/usr/bin/env bash
# =============================================================================
# 04_hic/chipseq/01_h3k27ac_chipseq.sh
# Purpose : H3K27ac ChIP-seq (HCASMC, HCAEC): BWA-MEM -> Picard 3.4.0
#           MarkDuplicates (remove) -> DFilter peaks (-ks=60 -bs=100 -lpval=8)
#           -> depth-normalised bigWig
# Status  : REFERENCE IMPLEMENTATION — confirm DFilter invocation and bigWig step
# =============================================================================
set -euo pipefail
BWA_IDX=${BWA_IDX:?}; FASTQ_DIR=${FASTQ_DIR:?}; OUT=${OUT:?}; DFILTER=${DFILTER:?}
THREADS=${THREADS:-16}
mkdir -p "${OUT}"
for r1 in "${FASTQ_DIR}"/*_R1.fastq.gz; do
  s=$(basename "${r1}" _R1.fastq.gz); r2=${r1/_R1/_R2}
  bwa mem -t "${THREADS}" "${BWA_IDX}" "${r1}" "${r2}" | samtools sort -@ 4 -o "${OUT}/${s}.bam"
  picard MarkDuplicates I="${OUT}/${s}.bam" O="${OUT}/${s}.dedup.bam" \
         M="${OUT}/${s}.dup_metrics.txt" REMOVE_DUPLICATES=true
  samtools index "${OUT}/${s}.dedup.bam"
done
# DFilter (per sample, with matched input if available) — TODO confirm flags
# "${DFILTER}/run_tag.sh" -d="${OUT}/${s}.dedup.bam" -c=<input.bam> -f=bam \
#     -o="${OUT}/${s}.dfilter.bed" -ks=60 -bs=100 -lpval=8
