#!/usr/bin/env bash
set -euo pipefail

GENES="/vol/COMPEPIWS/pipelines/references/mm10_reduced_chr18_chr19_genes.bed"

KIDNEY_H3K4ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me3_REP1.mLb.clN.bigWig"

LIVER_H3K4ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K4me3_REP1.mLb.clN.bigWig"

OUTDIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/promoters_H3K4me3"

mkdir -p "$OUTDIR"

###############################################################################
# Calculate H3K4me3 signal from 500 bp upstream to 1500 bp downstream of TSS.
#
# Original BED12 annotation is used instead of BED3 so deepTools can use
# column 6 (strand) to orient plus- and minus-strand genes correctly.
###############################################################################

computeMatrix reference-point \
    --referencePoint TSS \
    -b 500 \
    -a 1500 \
    --binSize 25 \
    -R "$GENES" \
    -S \
        "$LIVER_H3K4ME3" \
        "$KIDNEY_H3K4ME3" \
    --samplesLabel \
        "Liver E14.5 H3K4me3" \
        "Kidney E14.5 H3K4me3" \
    --missingDataAsZero \
    --outFileName "$OUTDIR/promoters_H3K4me3_matrix.gz" \
    --outFileSortedRegions "$OUTDIR/promoters_H3K4me3_sorted_regions.bed"

###############################################################################
# Heatmap: one row per promoter, ordered by mean H3K4me3 signal.
###############################################################################

plotHeatmap \
    --matrixFile "$OUTDIR/promoters_H3K4me3_matrix.gz" \
    --outFileName "$OUTDIR/promoters_H3K4me3_heatmap.png" \
    --outFileNameMatrix "$OUTDIR/promoters_H3K4me3_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "TSS" \
    --plotTitle "H3K4me3 signal at promoters" \
    --whatToShow "heatmap and colorbar"

###############################################################################
# Average profile: mean H3K4me3 signal across all promoters.
###############################################################################

plotProfile \
    --matrixFile "$OUTDIR/promoters_H3K4me3_matrix.gz" \
    --outFileName "$OUTDIR/promoters_H3K4me3_profile.png" \
    --refPointLabel "TSS" \
    --plotTitle "Average H3K4me3 signal around promoters" \
    --perGroup
