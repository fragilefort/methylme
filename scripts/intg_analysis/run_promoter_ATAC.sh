#!/usr/bin/env bash
set -euo pipefail

GENES="/vol/COMPEPIWS/pipelines/references/mm10_reduced_chr18_chr19_genes.bed"

LIVER_ATAC="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/liver_ATAC_aggregate_signal.bw"

KIDNEY_ATAC="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/kidney_ATAC_aggregate_signal.bw"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/promoters_ATAC"

mkdir -p "$OUTDIR"

computeMatrix reference-point \
    --referencePoint TSS \
    -b 500 \
    -a 1500 \
    --binSize 25 \
    -R "$GENES" \
    -S \
        "$LIVER_ATAC" \
        "$KIDNEY_ATAC" \
    --samplesLabel \
        "Liver aggregate ATAC" \
        "Kidney aggregate ATAC" \
    --missingDataAsZero \
    --outFileName "$OUTDIR/promoters_ATAC_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/promoters_ATAC_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/promoters_ATAC_matrix.gz" \
    --outFileName "$OUTDIR/promoters_ATAC_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/promoters_ATAC_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "TSS" \
    --plotTitle "Aggregate ATAC accessibility at promoters" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/promoters_ATAC_matrix.gz" \
    --outFileName "$OUTDIR/promoters_ATAC_profile.png" \
    --refPointLabel "TSS" \
    --plotTitle "Average ATAC accessibility around promoters" \
    --perGroup
