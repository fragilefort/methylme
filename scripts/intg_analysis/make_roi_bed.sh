#!/usr/bin/env bash
set -euo pipefail

GENES="/vol/COMPEPIWS/pipelines/references/mm10_reduced_chr18_chr19_genes.bed"

ATAC_PEAKS="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/peaks/consensus_peaks.mLb.clN.bed"

CHIP_SEG_DIR="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation"

WGBS_SEG_DIR="/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/segmentation"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI"


mkdir -p "$OUTDIR"

###############################################################################
# 1. Promoters: 500 bp upstream and 1,500 bp downstream of the TSS
###############################################################################

awk 'BEGIN {OFS="\t"}
$6 == "+" {
    start = $2 - 500
    end = $2 + 1500

    if (start < 0) start = 0

    if (end > start) {
        print $1, start, end
    }
}
$6 == "-" {
    start = $3 - 1500
    end = $3 + 500

    if (start < 0) start = 0

    if (end > start) {
        print $1, start, end
    }
}' "$GENES" \
| LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
| uniq \
> "$OUTDIR/promoters_500up_1500down.bed"

###############################################################################
# 2. Chromatin-accessible regions: ATAC-seq consensus peaks
###############################################################################

cut -f1-3 "$ATAC_PEAKS" \
| LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
| uniq \
> "$OUTDIR/atac_consensus_peaks.bed"

###############################################################################
# 3. Enhancers: all ChromHMM states with "Enh" in their state label
###############################################################################

awk 'BEGIN {OFS="\t"} $4 ~ /Enh/ {print $1, $2, $3}' \
    "$CHIP_SEG_DIR/liver_15_segments.bed" \
| LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
| uniq \
> "$OUTDIR/liver_ChromHMM_enhancers.bed"

awk 'BEGIN {OFS="\t"} $4 ~ /Enh/ {print $1, $2, $3}' \
    "$CHIP_SEG_DIR/kidney_15_segments.bed" \
| LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
| uniq \
> "$OUTDIR/kidney_ChromHMM_enhancers.bed"

###############################################################################
# 4. PMDs: MethylSeekR segments labelled PMD
###############################################################################

awk 'BEGIN {OFS="\t"} $4 == "PMD" {print $1, $2, $3}' \
    "$WGBS_SEG_DIR/liver_final.bed" \
| LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
| uniq \
> "$OUTDIR/liver_PMDs.bed"

awk 'BEGIN {OFS="\t"} $4 == "PMD" {print $1, $2, $3}' \
    "$WGBS_SEG_DIR/kidney_final.bed" \
| LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
| uniq \
> "$OUTDIR/kidney_PMDs.bed"

###############################################################################
# 5. Checks
###############################################################################

echo "ROI BED files created:"
wc -l "$OUTDIR"/*.bed

echo
echo "Liver ChromHMM enhancer states:"
awk '$4 ~ /Enh/ {print $4}' \
    "$CHIP_SEG_DIR/liver_15_segments.bed" \
| sort \
| uniq -c

echo
echo "Kidney ChromHMM enhancer states:"
awk '$4 ~ /Enh/ {print $4}' \
    "$CHIP_SEG_DIR/kidney_15_segments.bed" \
| sort \
| uniq -c
