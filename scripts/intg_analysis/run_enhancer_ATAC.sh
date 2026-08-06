#!/usr/bin/env bash
set -euo pipefail

LIVER_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_ChromHMM_enhancers.bed"

KIDNEY_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_ChromHMM_enhancers.bed"

LIVER_ATAC="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/liver_ATAC_aggregate_signal.bw"

KIDNEY_ATAC="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/kidney_ATAC_aggregate_signal.bw"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/enhancers_ATAC"

mkdir -p "$OUTDIR"

computeMatrix reference-point \
    --referencePoint center \
    -b 1000 \
    -a 1000 \
    --binSize 25 \
    -R \
        "$LIVER_ENHANCERS" \
        "$KIDNEY_ENHANCERS" \
    -S \
        "$LIVER_ATAC" \
        "$KIDNEY_ATAC" \
    --samplesLabel \
        "Liver aggregate ATAC" \
        "Kidney aggregate ATAC" \
    --missingDataAsZero \
    --outFileName "$OUTDIR/enhancers_ATAC_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/enhancers_ATAC_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/enhancers_ATAC_matrix.gz" \
    --outFileName "$OUTDIR/enhancers_ATAC_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/enhancers_ATAC_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "Enhancer centre" \
    --regionsLabel \
        "Liver ChromHMM enhancers" \
        "Kidney ChromHMM enhancers" \
    --plotTitle "ATAC accessibility around ChromHMM enhancers" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/enhancers_ATAC_matrix.gz" \
    --outFileName "$OUTDIR/enhancers_ATAC_profile.png" \
    --refPointLabel "Enhancer centre" \
    --regionsLabel \
        "Liver ChromHMM enhancers" \
        "Kidney ChromHMM enhancers" \
    --plotTitle "Average ATAC accessibility around ChromHMM enhancers" \
    --perGroup
