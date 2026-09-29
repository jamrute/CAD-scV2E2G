#!/usr/bin/env bash
# =============================================================================
# 03_perturbseq/crispri_validation/02_amotl2_rnaseq_align.sh
# Purpose : AMOTL2 KD vs NTC bulk RNA-seq (stressed HCASMCs; PE150, NovaSeq X Plus)
#           trim-galore (length 50) -> STAR 2.7.0 (hg38, GENCODE v38) -> HTSeq-count 2.1.2
# Data    : GSE329169
# Status  : REFERENCE IMPLEMENTATION — confirm trim-galore options and strandedness
# =============================================================================
set -euo pipefail
FASTQ=${FASTQ:?}; STAR_IDX=${STAR_IDX:?}; GTF=${GTF:?}; OUT=${OUT:?}; THREADS=${THREADS:-8}
mkdir -p "${OUT}"/{trim,star,counts}

for r1 in "${FASTQ}"/*_R1.fastq.gz; do
  s=$(basename "${r1}" _R1.fastq.gz); r2=${r1/_R1/_R2}
  trim_galore --paired --hardtrim5 50 -o "${OUT}/trim" "${r1}" "${r2}"   # TODO confirm: "trimmed to sequence length of 50"
  STAR --runThreadN "${THREADS}" --genomeDir "${STAR_IDX}" \
       --readFilesIn "${OUT}/trim/${s}_R1.50bp_5prime.fq.gz" "${OUT}/trim/${s}_R2.50bp_5prime.fq.gz" \
       --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix "${OUT}/star/${s}."
  htseq-count -f bam -r pos -s reverse -t exon -i gene_id \
       "${OUT}/star/${s}.Aligned.sortedByCoord.out.bam" "${GTF}" > "${OUT}/counts/${s}.txt"   # TODO -s
done
