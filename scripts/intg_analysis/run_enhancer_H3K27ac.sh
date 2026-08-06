#!/usr/bin/env bash
set -euo pipefail

LIVER_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_ChromHMM_enhancers.bed"

KIDNEY_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_ChromHMM_enhancers.bed"

LIVER_H3K27AC="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K27ac_REP1.mLb.clN.bigWig"

KIDNEY_H3K27AC="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K27ac_REP1.mLb.clN.bigWig"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/enhancers_H3K27ac"

mkdir -p "$OUTDIR"

###############################################################################
# Measure H3K27ac around the centre of liver and kidney ChromHMM enhancers.
###############################################################################

computeMatrix reference-point \
    --referencePoint center \
    -b 1000 \
    -a 1000 \
    --binSize 25 \
    -R \
        "$LIVER_ENHANCERS" \
        "$KIDNEY_ENHANCERS" \
    -S \
        "$LIVER_H3K27AC" \
        "$KIDNEY_H3K27AC" \
    --samplesLabel \
        "Liver E14.5 H3K27ac" \
        "Kidney E14.5 H3K27ac" \
    --missingDataAsZero \
    --outFileName "$OUTDIR/enhancers_H3K27ac_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/enhancers_H3K27ac_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/enhancers_H3K27ac_matrix.gz" \
    --outFileName "$OUTDIR/enhancers_H3K27ac_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/enhancers_H3K27ac_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "Enhancer centre" \
    --regionsLabel \
        "Liver ChromHMM enhancers" \
        "Kidney ChromHMM enhancers" \
    --plotTitle "H3K27ac signal around ChromHMM enhancers" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/enhancers_H3K27ac_matrix.gz" \
    --outFileName "$OUTDIR/enhancers_H3K27ac_profile.png" \
    --refPointLabel "Enhancer centre" \
    --regionsLabel \
        "Liver ChromHMM enhancers" \
        "Kidney ChromHMM enhancers" \
    --plotTitle "Average H3K27ac signal around ChromHMM enhancers" \
    --perGroup
