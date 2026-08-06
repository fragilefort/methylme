#!/usr/bin/env bash
set -euo pipefail

LIVER_PMDS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_PMDs.bed"

KIDNEY_PMDS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_PMDs.bed"

LIVER_H3K9ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K9me3_REP1.mLb.clN.bigWig"

KIDNEY_H3K9ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K9me3_REP1.mLb.clN.bigWig"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/PMDs_H3K9me3"

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
        "$LIVER_H3K9ME3" \
        "$KIDNEY_H3K9ME3" \
    --samplesLabel \
        "Liver E14.5 H3K9me3" \
        "Kidney E14.5 H3K9me3" \
    --missingDataAsZero \
    --outFileName "$OUTDIR/PMDs_H3K9me3_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/PMDs_H3K9me3_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/PMDs_H3K9me3_matrix.gz" \
    --outFileName "$OUTDIR/PMDs_H3K9me3_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/PMDs_H3K9me3_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --startLabel "PMD start" \
    --endLabel "PMD end" \
    --regionsLabel \
        "Liver PMDs" \
        "Kidney PMDs" \
    --plotTitle "H3K9me3 signal across partially methylated domains" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/PMDs_H3K9me3_matrix.gz" \
    --outFileName "$OUTDIR/PMDs_H3K9me3_profile.png" \
    --startLabel "PMD start" \
    --endLabel "PMD end" \
    --regionsLabel \
        "Liver PMDs" \
        "Kidney PMDs" \
    --plotTitle "Average H3K9me3 signal across partially methylated domains" \
    --perGroup
