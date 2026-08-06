#!/usr/bin/env bash
set -euo pipefail

GENES="/vol/COMPEPIWS/pipelines/references/mm10_reduced_chr18_chr19_genes.bed"

LIVER_WGBS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/liver_WGBS_merged_signal.bw"

KIDNEY_WGBS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/signals/kidney_WGBS_merged_signal.bw"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/promoters_WGBS"

mkdir -p "$OUTDIR"

###############################################################################
# Aggregate WGBS methylation beta values around strand-aware promoter TSSs.
#
# No --missingDataAsZero:
# missing CpG methylation data are not interpreted as unmethylated DNA.
###############################################################################

computeMatrix reference-point \
    --referencePoint TSS \
    -b 500 \
    -a 1500 \
    --binSize 50 \
    -R "$GENES" \
    -S \
        "$LIVER_WGBS" \
        "$KIDNEY_WGBS" \
    --samplesLabel \
        "Liver merged WGBS methylation" \
        "Kidney merged WGBS methylation" \
    --outFileName "$OUTDIR/promoters_WGBS_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/promoters_WGBS_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$OUTDIR/promoters_WGBS_matrix.gz" \
    --outFileName "$OUTDIR/promoters_WGBS_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/promoters_WGBS_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "TSS" \
    --plotTitle "DNA methylation at promoters" \
    --colorMap Blues \
    --zMin 0 0 \
    --zMax 1 1 \
    --missingDataColor white \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$OUTDIR/promoters_WGBS_matrix.gz" \
    --outFileName "$OUTDIR/promoters_WGBS_profile.png" \
    --refPointLabel "TSS" \
    --plotTitle "Average DNA methylation around promoters" \
    --perGroup
