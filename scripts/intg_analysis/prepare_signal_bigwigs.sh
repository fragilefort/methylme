#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals"

CHROM_SIZES="/vol/COMPEPIWS/pipelines/references/mm10_chromsizes_reduced.txt"

KIDNEY_WGBS="/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/signal/kidney_merged_signal.bedgraph"
LIVER_WGBS="/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/signal/liver_merged_signal.bedgraph"

KIDNEY_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/signal/kidney_counts.bed"
LIVER_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/signal/liver_counts.bed"

mkdir -p "$OUTDIR"

###############################################################################
# Check that the reference chromosome-size file contains our chromosomes.
###############################################################################

echo "Chromosome sizes used:"
grep -E '^chr(18|19)[[:space:]]' "$CHROM_SIZES"

###############################################################################
# 1. Convert merged WGBS methylation BedGraphs to BigWigs.
#
# Input columns already are:
# chromosome, start, end, methylation beta value
###############################################################################

LC_ALL=C sort -k1,1 -k2,2n \
    "$KIDNEY_WGBS" \
    > "$OUTDIR/kidney_WGBS_merged_signal.sorted.bedgraph"

bedGraphToBigWig \
    "$OUTDIR/kidney_WGBS_merged_signal.sorted.bedgraph" \
    "$CHROM_SIZES" \
    "$OUTDIR/kidney_WGBS_merged_signal.bw"

LC_ALL=C sort -k1,1 -k2,2n \
    "$LIVER_WGBS" \
    > "$OUTDIR/liver_WGBS_merged_signal.sorted.bedgraph"

bedGraphToBigWig \
    "$OUTDIR/liver_WGBS_merged_signal.sorted.bedgraph" \
    "$CHROM_SIZES" \
    "$OUTDIR/liver_WGBS_merged_signal.bw"

###############################################################################
# 2. Convert aggregate ATAC BED tracks to BigWigs.
#
# ATAC BED input columns:
# 1 = chromosome
# 2 = start
# 3 = end
# 4 = .
# 5 = aggregate accessibility count/signal
# 6 = .
#
# Create four-column BedGraph:
# chromosome, start, end, ATAC signal from column 5
###############################################################################

awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5}' \
    "$KIDNEY_ATAC" \
| LC_ALL=C sort -k1,1 -k2,2n \
> "$OUTDIR/kidney_ATAC_aggregate_signal.sorted.bedgraph"

bedGraphToBigWig \
    "$OUTDIR/kidney_ATAC_aggregate_signal.sorted.bedgraph" \
    "$CHROM_SIZES" \
    "$OUTDIR/kidney_ATAC_aggregate_signal.bw"

awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5}' \
    "$LIVER_ATAC" \
| LC_ALL=C sort -k1,1 -k2,2n \
> "$OUTDIR/liver_ATAC_aggregate_signal.sorted.bedgraph"

bedGraphToBigWig \
    "$OUTDIR/liver_ATAC_aggregate_signal.sorted.bedgraph" \
    "$CHROM_SIZES" \
    "$OUTDIR/liver_ATAC_aggregate_signal.bw"

###############################################################################
# 3. Confirm final BigWig files
###############################################################################

echo
echo "Created BigWig signal files:"
ls -lh "$OUTDIR"/*.bw
