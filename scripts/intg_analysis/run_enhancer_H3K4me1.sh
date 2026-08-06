#!/usr/bin/env bash
set -euo pipefail

LIVER_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_ChromHMM_enhancers.bed"

KIDNEY_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_ChromHMM_enhancers.bed"

LIVER_H3K4ME1="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K4me1_REP1.mLb.clN.bigWig"

KIDNEY_H3K4ME1="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me1_REP1.mLb.clN.bigWig"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/enhancers_H3K4me1"

mkdir -p "$OUTDIR"

###############################################################################
# Measure H3K4me1 around the centre of liver and kidney ChromHMM enhancers.
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
        "$LIVER_H3K4ME1" \
        "$KIDNEY_H3K4ME1" \
    --samplesLabel \
        "Liver E14.5 H3K4me1" \
        "Kidney E14.5 H3K4me1" \
    --missingDataAsZero \
    --outFileName "$OUTDIR/enhancers_H3K4me1_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/enhancers_H3K4me1_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/enhancers_H3K4me1_matrix.gz" \
    --outFileName "$OUTDIR/enhancers_H3K4me1_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/enhancers_H3K4me1_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "Enhancer centre" \
    --regionsLabel \
        "Liver ChromHMM enhancers" \
        "Kidney ChromHMM enhancers" \
    --plotTitle "H3K4me1 signal around ChromHMM enhancers" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/enhancers_H3K4me1_matrix.gz" \
    --outFileName "$OUTDIR/enhancers_H3K4me1_profile.png" \
    --refPointLabel "Enhancer centre" \
    --regionsLabel \
        "Liver ChromHMM enhancers" \
        "Kidney ChromHMM enhancers" \
    --plotTitle "Average H3K4me1 signal around ChromHMM enhancers" \
    --perGroup
