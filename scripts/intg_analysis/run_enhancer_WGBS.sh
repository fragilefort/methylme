#!/usr/bin/env bash
set -euo pipefail

LIVER_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_ChromHMM_enhancers.bed"

KIDNEY_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_ChromHMM_enhancers.bed"

LIVER_WGBS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/liver_WGBS_merged_signal.bw"

KIDNEY_WGBS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/kidney_WGBS_merged_signal.bw"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/enhancers_WGBS"

mkdir -p "$OUTDIR"

computeMatrix reference-point \
    --referencePoint center \
    -b 1000 \
    -a 1000 \
    --binSize 50 \
    -R \
        "$LIVER_ENHANCERS" \
        "$KIDNEY_ENHANCERS" \
    -S \
        "$LIVER_WGBS" \
        "$KIDNEY_WGBS" \
    --samplesLabel \
        "Liver merged WGBS methylation" \
        "Kidney merged WGBS methylation" \
    --outFileName "$OUTDIR/enhancers_WGBS_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/enhancers_WGBS_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/enhancers_WGBS_matrix.gz" \
    --outFileName "$OUTDIR/enhancers_WGBS_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/enhancers_WGBS_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "Enhancer centre" \
    --regionsLabel \
        "Liver ChromHMM enhancers" \
        "Kidney ChromHMM enhancers" \
    --plotTitle "DNA methylation around ChromHMM enhancers" \
    --colorMap Blues \
    --zMin 0 0 \
    --zMax 1 1 \
    --missingDataColor white \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/enhancers_WGBS_matrix.gz" \
    --outFileName "$OUTDIR/enhancers_WGBS_profile.png" \
    --refPointLabel "Enhancer centre" \
    --regionsLabel \
        "Liver ChromHMM enhancers" \
        "Kidney ChromHMM enhancers" \
    --plotTitle "Average DNA methylation around ChromHMM enhancers" \
    --perGroup
