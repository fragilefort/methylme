#!/usr/bin/env bash
#SBATCH --job-name=deeptools_replot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=deeptools_replot_%j.out
#SBATCH --error=deeptools_replot_%j.err

set -euo pipefail

ROI_DIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/integrative_analysis_deeptools_ROI"
OUTROOT="$ROI_DIR/unified_deeptools_option1_14_5_REP1"

command -v plotHeatmap >/dev/null || {
    echo "Missing command: plotHeatmap" >&2
    exit 1
}

# Re-render existing computeMatrix outputs only. This script does not rerun computeMatrix.
# The _improved files preserve the original figures and matrices.
plot_matrix() {
    local matrix="$1"
    local outdir="$2"
    local prefix="$3"
    local title="$4"
    local colour_map="$5"
    local y_label="$6"
    shift 6

    if [[ ! -s "$matrix" ]]; then
        echo "Skipping missing matrix: $matrix" >&2
        return 0
    fi

    plotHeatmap \
        --matrixFile "$matrix" \
        --outFileName "$outdir/${prefix}_improved.png" \
        --outFileNameMatrix "$outdir/${prefix}_improved_values.tab" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap "$colour_map" \
        --colorNumber 256 \
        --zMin auto \
        --zMax auto \
        --heatmapWidth 14 \
        --heatmapHeight 16 \
        --dpi 300 \
        --interpolationMethod nearest \
        --boxAroundHeatmaps no \
        --whatToShow "plot, heatmap and colorbar" \
        --legendLocation upper-right \
        --labelRotation 0 \
        --yAxisLabel "$y_label" \
        --plotTitle "$title" \
        "$@"

    plotHeatmap \
        --matrixFile "$matrix" \
        --outFileName "$outdir/${prefix}_improved.pdf" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap "$colour_map" \
        --colorNumber 256 \
        --zMin auto \
        --zMax auto \
        --heatmapWidth 14 \
        --heatmapHeight 16 \
        --dpi 300 \
        --interpolationMethod nearest \
        --boxAroundHeatmaps no \
        --whatToShow "plot, heatmap and colorbar" \
        --legendLocation upper-right \
        --labelRotation 0 \
        --yAxisLabel "$y_label" \
        --plotTitle "$title" \
        "$@"
}

# Promoters
plot_matrix \
    "$OUTROOT/promoters_ATAC/matrix.gz" \
    "$OUTROOT/promoters_ATAC" \
    "promoters_ATAC_heatmap_profile" \
    "ATAC accessibility at promoters" \
    "magma" \
    "Mean ATAC signal" \
    --refPointLabel "TSS" \
    --xAxisLabel "Distance from TSS (bp)"

plot_matrix \
    "$OUTROOT/promoters_WGBS/matrix.gz" \
    "$OUTROOT/promoters_WGBS" \
    "promoters_WGBS_heatmap_profile" \
    "DNA methylation at promoters" \
    "viridis" \
    "Mean methylation signal" \
    --refPointLabel "TSS" \
    --xAxisLabel "Distance from TSS (bp)"

plot_matrix \
    "$OUTROOT/promoters_H3K4me3/matrix.gz" \
    "$OUTROOT/promoters_H3K4me3" \
    "promoters_H3K4me3_heatmap_profile" \
    "H3K4me3 signal at promoters" \
    "Blues" \
    "Mean H3K4me3 signal" \
    --refPointLabel "TSS" \
    --xAxisLabel "Distance from TSS (bp)"

# Enhancers
plot_matrix \
    "$OUTROOT/enhancers_ATAC/matrix.gz" \
    "$OUTROOT/enhancers_ATAC" \
    "enhancers_ATAC_heatmap_profile" \
    "ATAC accessibility at enhancers" \
    "magma" \
    "Mean ATAC signal" \
    --refPointLabel "Centre" \
    --xAxisLabel "Distance from enhancer centre (bp)" \
    --regionsLabel "Liver enhancers" "Kidney enhancers"

plot_matrix \
    "$OUTROOT/enhancers_WGBS/matrix.gz" \
    "$OUTROOT/enhancers_WGBS" \
    "enhancers_WGBS_heatmap_profile" \
    "DNA methylation at enhancers" \
    "viridis" \
    "Mean methylation signal" \
    --refPointLabel "Centre" \
    --xAxisLabel "Distance from enhancer centre (bp)" \
    --regionsLabel "Liver enhancers" "Kidney enhancers"

plot_matrix \
    "$OUTROOT/enhancers_H3K4me1/matrix.gz" \
    "$OUTROOT/enhancers_H3K4me1" \
    "enhancers_H3K4me1_heatmap_profile" \
    "H3K4me1 signal at enhancers" \
    "YlGn" \
    "Mean H3K4me1 signal" \
    --refPointLabel "Centre" \
    --xAxisLabel "Distance from enhancer centre (bp)" \
    --regionsLabel "Liver enhancers" "Kidney enhancers"

plot_matrix \
    "$OUTROOT/enhancers_H3K27ac/matrix.gz" \
    "$OUTROOT/enhancers_H3K27ac" \
    "enhancers_H3K27ac_heatmap_profile" \
    "H3K27ac signal at enhancers" \
    "YlOrBr" \
    "Mean H3K27ac signal" \
    --refPointLabel "Centre" \
    --xAxisLabel "Distance from enhancer centre (bp)" \
    --regionsLabel "Liver enhancers" "Kidney enhancers"

# PMDs: one region group per plot, so omit --regionsLabel to prevent text crowding.
for tissue in liver kidney; do
    if [[ "$tissue" == "liver" ]]; then
        LABEL="Liver"
    else
        LABEL="Kidney"
    fi

    plot_matrix \
        "$OUTROOT/PMDs_ATAC/$tissue/matrix.gz" \
        "$OUTROOT/PMDs_ATAC/$tissue" \
        "${tissue}_PMDs_ATAC_heatmap_profile" \
        "$LABEL ATAC accessibility across PMDs" \
        "magma" \
        "Mean ATAC signal" \
        --startLabel "Start" \
        --endLabel "End" \
        --xAxisLabel "PMD body: 5 kb; flanks: +/-10 kb"

    plot_matrix \
        "$OUTROOT/PMDs_WGBS/$tissue/matrix.gz" \
        "$OUTROOT/PMDs_WGBS/$tissue" \
        "${tissue}_PMDs_WGBS_heatmap_profile" \
        "$LABEL DNA methylation across PMDs" \
        "viridis" \
        "Mean methylation signal" \
        --startLabel "Start" \
        --endLabel "End" \
        --xAxisLabel "PMD body: 5 kb; flanks: +/-10 kb"

    plot_matrix \
        "$OUTROOT/PMDs_H3K9me3/$tissue/matrix.gz" \
        "$OUTROOT/PMDs_H3K9me3/$tissue" \
        "${tissue}_PMDs_H3K9me3_heatmap_profile" \
        "$LABEL H3K9me3 signal across PMDs" \
        "Purples" \
        "Mean H3K9me3 signal" \
        --startLabel "Start" \
        --endLabel "End" \
        --xAxisLabel "PMD body: 5 kb; flanks: +/-10 kb"
done

echo "Replotting finished. Improved figures are saved beside the original figures under: $OUTROOT"
