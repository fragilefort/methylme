#!/usr/bin/env bash
set -euo pipefail

# Original ATAC-seq BigWigs: 14.5, replicate 1 only.
LIVER_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/liver_14.5_REP1.mLb.clN.bigWig"
KIDNEY_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/kidney_14.5_REP1.mLb.clN.bigWig"

# Existing region inputs on the cluster.
GENES="/vol/COMPEPIWS/pipelines/references/mm10_reduced_chr18_chr19_genes.bed"
LIVER_PMDS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_PMDs.bed"
KIDNEY_PMDS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_PMDs.bed"
LIVER_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/liver_ChromHMM_enhancers.bed"
KIDNEY_ENHANCERS="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/kidney_ChromHMM_enhancers.bed"

# New output directory: this does not overwrite the aggregate-track results.
OUTROOT="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis_deeptools_ROI/ATACseq2_14_5_REP1"
PROMOTER_OUTDIR="$OUTROOT/promoters_ATAC"
ENHANCER_OUTDIR="$OUTROOT/enhancers_ATAC"
PMD_OUTDIR="$OUTROOT/PMDs_ATAC"
mkdir -p "$PROMOTER_OUTDIR" "$ENHANCER_OUTDIR" "$PMD_OUTDIR"

for file in "$LIVER_ATAC" "$KIDNEY_ATAC" "$GENES" "$LIVER_PMDS" "$KIDNEY_PMDS" "$LIVER_ENHANCERS" "$KIDNEY_ENHANCERS"; do
    [[ -s "$file" ]] || { echo "Missing or empty input: $file" >&2; exit 1; }
done

# 1. Promoters: TSS -500 bp to +1,500 bp.
computeMatrix reference-point \
    --referencePoint TSS \
    -b 500 \
    -a 1500 \
    --binSize 25 \
    -R "$GENES" \
    -S "$LIVER_ATAC" "$KIDNEY_ATAC" \
    --samplesLabel "Liver 14.5 REP1 ATAC" "Kidney 14.5 REP1 ATAC" \
    --missingDataAsZero \
    --outFileName "$PROMOTER_OUTDIR/promoters_ATAC_matrix.gz" \
    --outFileSortedRegions "$PROMOTER_OUTDIR/promoters_ATAC_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$PROMOTER_OUTDIR/promoters_ATAC_matrix.gz" \
    --outFileName "$PROMOTER_OUTDIR/promoters_ATAC_heatmap.png" \
    --outFileNameMatrix "$PROMOTER_OUTDIR/promoters_ATAC_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "TSS" \
    --plotTitle "ATAC accessibility at promoters: 14.5 REP1" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$PROMOTER_OUTDIR/promoters_ATAC_matrix.gz" \
    --outFileName "$PROMOTER_OUTDIR/promoters_ATAC_profile.png" \
    --refPointLabel "TSS" \
    --plotTitle "Average ATAC accessibility around promoters: 14.5 REP1"

# 2. ChromHMM enhancers: 1 kb on either side of enhancer centre.
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 \
    -a 1000 \
    --binSize 25 \
    -R "$LIVER_ENHANCERS" "$KIDNEY_ENHANCERS" \
    -S "$LIVER_ATAC" "$KIDNEY_ATAC" \
    --samplesLabel "Liver 14.5 REP1 ATAC" "Kidney 14.5 REP1 ATAC" \
    --missingDataAsZero \
    --outFileName "$ENHANCER_OUTDIR/enhancers_ATAC_matrix.gz" \
    --outFileSortedRegions "$ENHANCER_OUTDIR/enhancers_ATAC_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$ENHANCER_OUTDIR/enhancers_ATAC_matrix.gz" \
    --outFileName "$ENHANCER_OUTDIR/enhancers_ATAC_heatmap.png" \
    --outFileNameMatrix "$ENHANCER_OUTDIR/enhancers_ATAC_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --refPointLabel "Enhancer centre" \
    --regionsLabel "Liver ChromHMM enhancers" "Kidney ChromHMM enhancers" \
    --plotTitle "ATAC accessibility around ChromHMM enhancers: 14.5 REP1" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$ENHANCER_OUTDIR/enhancers_ATAC_matrix.gz" \
    --outFileName "$ENHANCER_OUTDIR/enhancers_ATAC_profile.png" \
    --refPointLabel "Enhancer centre" \
    --regionsLabel "Liver ChromHMM enhancers" "Kidney ChromHMM enhancers" \
    --plotTitle "Average ATAC accessibility around ChromHMM enhancers: 14.5 REP1" \
    --perGroup

# 3. PMDs: 10 kb flanks and a region body scaled to 5 kb.
computeMatrix scale-regions \
    -b 10000 \
    -a 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    -R "$LIVER_PMDS" "$KIDNEY_PMDS" \
    -S "$LIVER_ATAC" "$KIDNEY_ATAC" \
    --samplesLabel "Liver 14.5 REP1 ATAC" "Kidney 14.5 REP1 ATAC" \
    --missingDataAsZero \
    --outFileName "$PMD_OUTDIR/PMDs_ATAC_matrix.gz" \
    --outFileSortedRegions "$PMD_OUTDIR/PMDs_ATAC_sorted_regions.bed"

plotHeatmap \
    --matrixFile "$PMD_OUTDIR/PMDs_ATAC_matrix.gz" \
    --outFileName "$PMD_OUTDIR/PMDs_ATAC_heatmap.png" \
    --outFileNameMatrix "$PMD_OUTDIR/PMDs_ATAC_values.tab" \
    --sortRegions descend \
    --sortUsing mean \
    --startLabel "PMD start" \
    --endLabel "PMD end" \
    --regionsLabel "Liver PMDs" "Kidney PMDs" \
    --plotTitle "ATAC accessibility across partially methylated domains: 14.5 REP1" \
    --whatToShow "heatmap and colorbar"

plotProfile \
    --matrixFile "$PMD_OUTDIR/PMDs_ATAC_matrix.gz" \
    --outFileName "$PMD_OUTDIR/PMDs_ATAC_profile.png" \
    --startLabel "PMD start" \
    --endLabel "PMD end" \
    --regionsLabel "Liver PMDs" "Kidney PMDs" \
    --plotTitle "Average ATAC accessibility across PMDs: 14.5 REP1" \
    --perGroup

echo "Finished. Results: $OUTROOT"
