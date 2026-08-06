#!/usr/bin/env bash
set -euo pipefail

LIVER_PMDS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_PMDs.bed"

KIDNEY_PMDS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_PMDs.bed"

LIVER_ATAC="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/liver_ATAC_aggregate_signal.bw"

KIDNEY_ATAC="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/kidney_ATAC_aggregate_signal.bw"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/PMDs_ATAC"

mkdir -p "$OUTDIR"

computeMatrix scale-regions \
    -b 10000 \
    -a 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    -R \
        "$LIVER_PMDS" \
        "$KIDNEY_PMDS" \
    -S \
        "$LIVER_ATAC" \
        "$KIDNEY_ATAC" \
    --samplesLabel \
        "Liver aggregate ATAC" \
        "Kidney aggregate ATAC" \
    --missingDataAsZero \
    --outFileName "$OUTDIR/PMDs_ATAC_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/PMDs_ATAC_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/PMDs_ATAC_matrix.gz" \
    --outFileName "$OUTDIR/PMDs_ATAC_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/PMDs_ATAC_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --startLabel "PMD start" \
    --endLabel "PMD end" \
    --regionsLabel \
        "Liver PMDs" \
        "Kidney PMDs" \
    --plotTitle "ATAC accessibility across partially methylated domains" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/PMDs_ATAC_matrix.gz" \
    --outFileName "$OUTDIR/PMDs_ATAC_profile.png" \
    --startLabel "PMD start" \
    --endLabel "PMD end" \
    --regionsLabel \
        "Liver PMDs" \
        "Kidney PMDs" \
    --plotTitle "Average ATAC accessibility across partially methylated domains" \
    --perGroup
