#!/usr/bin/env python
"""
04_hic/ctcf_allelic/02_design_masked_primers.py
Purpose : SYBR qPCR primers (75-150 bp amplicons) around candidate SNPs with
          Primer3 (v2.6.1 via primer3-py). Any primer pair whose binding site
          overlaps a variant in the merged 4-line VCF is discarded, so annealing
          is identical across lines. The target SNP must lie inside the amplicon
          but outside both primers.
Status  : REFERENCE IMPLEMENTATION — confirm Primer3 global settings used
Requires: primer3-py, pysam
"""
import argparse
import csv
import primer3
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--snps", required=True, help="TSV with snp, chr, pos (1-based)")
ap.add_argument("--fasta", required=True)
ap.add_argument("--vcf", required=True, help="merged 4-line VCF (bgzipped + indexed)")
ap.add_argument("--flank", type=int, default=200)
ap.add_argument("--n-return", type=int, default=20)
ap.add_argument("--out", required=True)
a = ap.parse_args()

fa = pysam.FastaFile(a.fasta)
vcf = pysam.VariantFile(a.vcf)

GLOBAL = {
    "PRIMER_PRODUCT_SIZE_RANGE": [[75, 150]],
    "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 63.0,
    "PRIMER_NUM_RETURN": a.n_return,
}


def variant_positions(chrom, start0, end0):
    return {r.pos for r in vcf.fetch(chrom, start0, end0)}  # 1-based positions


with open(a.snps) as fh, open(a.out, "w", newline="") as oh:
    w = csv.writer(oh, delimiter="\t")
    w.writerow(["snp", "rank", "left", "right", "left_start", "right_end", "product_size", "left_tm", "right_tm"])
    for row in csv.DictReader(fh, delimiter="\t"):
        chrom, pos = row["chr"], int(row["pos"])
        s0 = pos - 1 - a.flank
        seq = fa.fetch(chrom, s0, pos + a.flank).upper()
        snp_idx = pos - 1 - s0
        variants = variant_positions(chrom, s0, pos + a.flank)
        res = primer3.bindings.design_primers(
            {"SEQUENCE_ID": row["snp"], "SEQUENCE_TEMPLATE": seq,
             "SEQUENCE_TARGET": [snp_idx, 1]}, GLOBAL)
        kept = 0
        for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            ls, ll = res[f"PRIMER_LEFT_{i}"]
            rs, rl = res[f"PRIMER_RIGHT_{i}"]            # rs = 3'-most base (0-based)
            left_bases = {s0 + ls + k + 1 for k in range(ll)}
            right_bases = {s0 + rs - k + 1 for k in range(rl)}
            if (left_bases | right_bases) & variants:
                continue                                 # masked: primer overlaps a variant
            kept += 1
            w.writerow([row["snp"], kept, res[f"PRIMER_LEFT_{i}_SEQUENCE"],
                        res[f"PRIMER_RIGHT_{i}_SEQUENCE"], s0 + ls + 1, s0 + rs + 1,
                        res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                        round(res[f"PRIMER_LEFT_{i}_TM"], 1), round(res[f"PRIMER_RIGHT_{i}_TM"], 1)])
        if kept == 0:
            print(f"[warn] {row['snp']}: no variant-free primer pair")
