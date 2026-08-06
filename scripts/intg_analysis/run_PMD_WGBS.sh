#!/usr/bin/env bash
set -euo pipefail

LIVER_PMDS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_PMDs.bed"

KIDNEY_PMDS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_PMDs.bed"

LIVER_WGBS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/liver_WGBS_merged_signal.bw"

KIDNEY_WGBS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/kidney_WGBS_merged_signal.bw"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/PMDs_WGBS"

mkdir -p "$OUTDIR"

###############################################################################
# Compare WGBS methylation across broad, variable-length PMDs.
#
# Each PMD body is scaled to 5,000 bp for plotting.
# The 10-kb flanks are not scaled.
#
# No --missingDataAsZero:
# Missing methylation calls remain missing; they are not interpreted as beta = 0.
###############################################################################

computeMatrix scale-regions \
    -b 10000 \
    -a 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    -R \
        "$LIVER_PMDS" \
        "$KIDNEY_PMDS" \
    -S \
        "$LIVER_WGBS" \
        "$KIDNEY_WGBS" \
    --samplesLabel \
        "Liver merged WGBS methylation" \
        "Kidney merged WGBS methylation" \
    --outFileName "$OUTDIR/PMDs_WGBS_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/PMDs_WGBS_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/PMDs_WGBS_matrix.gz" \
    --outFileName "$OUTDIR/PMDs_WGBS_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/PMDs_WGBS_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --startLabel "PMD start" \
    --endLabel "PMD end" \
    --regionsLabel \
        "Liver PMDs" \
        "Kidney PMDs" \
    --plotTitle "DNA methylation across partially methylated domains" \
    --colorMap Blues \
    --zMin 0 0 \
    --zMax 1 1 \
    --missingDataColor white \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/PMDs_WGBS_matrix.gz" \
    --outFileName "$OUTDIR/PMDs_WGBS_profile.png" \
    --startLabel "PMD start" \
    --endLabel "PMD end" \
    --regionsLabel \
        "Liver PMDs" \
        "Kidney PMDs" \
    --plotTitle "Average DNA methylation across partially methylated domains" \
    --perGroup
