#!/usr/bin/env bash
#SBATCH --job-name=tissue_matched_deeptools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=tissue_matched_deeptools_%j.out
#SBATCH --error=tissue_matched_deeptools_%j.err

set -euo pipefail

# One tissue per figure: Liver regions with Liver signal only;
# Kidney regions with Kidney signal only.
# Outputs are written to a NEW directory, preserving prior combined analyses.

export PATH="/vol/COMPEPIWS/groups/wgbs2/methylme/conda_deeptools/bin:$PATH"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
ROI_DIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/integrative_analysis_deeptools_ROI"
OUTROOT="$ROI_DIR/tissue_matched_deeptools_option1_14_5_REP1"

GENES="/vol/COMPEPIWS/pipelines/references/mm10_reduced_chr18_chr19_genes.bed"
LIVER_PMDS="$ROI_DIR/liver_PMDs.bed"
KIDNEY_PMDS="$ROI_DIR/kidney_PMDs.bed"
LIVER_ENHANCERS="$ROI_DIR/liver_ChromHMM_enhancers.bed"
KIDNEY_ENHANCERS="$ROI_DIR/kidney_ChromHMM_enhancers.bed"

LIVER_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/liver_14.5_REP1.mLb.clN.bigWig"
KIDNEY_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/kidney_14.5_REP1.mLb.clN.bigWig"
LIVER_WGBS="$ROI_DIR/signals/liver_WGBS_merged_signal.bw"
KIDNEY_WGBS="$ROI_DIR/signals/kidney_WGBS_merged_signal.bw"
LIVER_H3K4ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K4me3_REP1.mLb.clN.bigWig"
KIDNEY_H3K4ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me3_REP1.mLb.clN.bigWig"
LIVER_H3K4ME1="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K4me1_REP1.mLb.clN.bigWig"
KIDNEY_H3K4ME1="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me1_REP1.mLb.clN.bigWig"
LIVER_H3K27AC="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K27ac_REP1.mLb.clN.bigWig"
KIDNEY_H3K27AC="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K27ac_REP1.mLb.clN.bigWig"
LIVER_H3K9ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K9me3_REP1.mLb.clN.bigWig"
KIDNEY_H3K9ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K9me3_REP1.mLb.clN.bigWig"

for cmd in computeMatrix plotHeatmap; do
    command -v "$cmd" >/dev/null || { echo "Missing command: $cmd" >&2; exit 1; }
done

for file in \
    "$GENES" "$LIVER_PMDS" "$KIDNEY_PMDS" "$LIVER_ENHANCERS" "$KIDNEY_ENHANCERS" \
    "$LIVER_ATAC" "$KIDNEY_ATAC" "$LIVER_WGBS" "$KIDNEY_WGBS" \
    "$LIVER_H3K4ME3" "$KIDNEY_H3K4ME3" "$LIVER_H3K4ME1" "$KIDNEY_H3K4ME1" \
    "$LIVER_H3K27AC" "$KIDNEY_H3K27AC" "$LIVER_H3K9ME3" "$KIDNEY_H3K9ME3"
do
    [[ -s "$file" ]] || { echo "Missing or empty input: $file" >&2; exit 1; }
done

plot_single() {
    local matrix="$1"
    local outdir="$2"
    local prefix="$3"
    local title="$4"
    local colour_map="$5"
    local y_label="$6"
    shift 6

    plotHeatmap \
        --matrixFile "$matrix" \
        --outFileName "$outdir/${prefix}.png" \
        --outFileNameMatrix "$outdir/${prefix}_values.tab" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap "$colour_map" \
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
        --outFileName "$outdir/${prefix}.pdf" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap "$colour_map" \
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

run_reference() {
    local regions="$1"
    local signal="$2"
    local outdir="$3"
    local prefix="$4"
    local title="$5"
    local colour_map="$6"
    local y_label="$7"
    local ref_label="$8"
    local x_label="$9"

    mkdir -p "$outdir"
    computeMatrix reference-point --referencePoint center -b 1000 -a 1000 --binSize 25 -p "$THREADS" \
        -R "$regions" \
        -S "$signal" \
        --samplesLabel "$LABEL" \
        --missingDataAsZero \
        --outFileName "$outdir/matrix.gz" \
        --outFileSortedRegions "$outdir/sorted_regions.bed"
    plot_single "$outdir/matrix.gz" "$outdir" "$prefix" "$title" "$colour_map" "$y_label" \
        --refPointLabel "$ref_label" \
        --xAxisLabel "$x_label"
}

run_promoter() {
    local signal="$1"
    local outdir="$2"
    local prefix="$3"
    local title="$4"
    local colour_map="$5"
    local y_label="$6"

    mkdir -p "$outdir"
    computeMatrix reference-point --referencePoint TSS -b 500 -a 1500 --binSize 25 -p "$THREADS" \
        -R "$GENES" \
        -S "$signal" \
        --samplesLabel "$LABEL" \
        --missingDataAsZero \
        --outFileName "$outdir/matrix.gz" \
        --outFileSortedRegions "$outdir/sorted_regions.bed"
    plot_single "$outdir/matrix.gz" "$outdir" "$prefix" "$title" "$colour_map" "$y_label" \
        --refPointLabel "TSS" \
        --xAxisLabel "Distance from TSS (bp)"
}

run_pmd() {
    local signal="$1"
    local outdir="$2"
    local prefix="$3"
    local title="$4"
    local colour_map="$5"
    local y_label="$6"

    mkdir -p "$outdir"
    computeMatrix scale-regions -b 10000 -a 10000 --regionBodyLength 5000 --binSize 100 -p "$THREADS" \
        -R "$PMDS" \
        -S "$signal" \
        --samplesLabel "$LABEL" \
        --missingDataAsZero \
        --outFileName "$outdir/matrix.gz" \
        --outFileSortedRegions "$outdir/sorted_regions.bed"
    plot_single "$outdir/matrix.gz" "$outdir" "$prefix" "$title" "$colour_map" "$y_label" \
        --startLabel "Start" \
        --endLabel "End" \
        --xAxisLabel "PMD body: 5 kb; flanks: +/-10 kb"
}

for tissue in liver kidney; do
    if [[ "$tissue" == "liver" ]]; then
        LABEL="Liver"
        ENHANCERS="$LIVER_ENHANCERS"
        PMDS="$LIVER_PMDS"
        ATAC="$LIVER_ATAC"
        WGBS="$LIVER_WGBS"
        H3K4ME3="$LIVER_H3K4ME3"
        H3K4ME1="$LIVER_H3K4ME1"
        H3K27AC="$LIVER_H3K27AC"
        H3K9ME3="$LIVER_H3K9ME3"
    else
        LABEL="Kidney"
        ENHANCERS="$KIDNEY_ENHANCERS"
        PMDS="$KIDNEY_PMDS"
        ATAC="$KIDNEY_ATAC"
        WGBS="$KIDNEY_WGBS"
        H3K4ME3="$KIDNEY_H3K4ME3"
        H3K4ME1="$KIDNEY_H3K4ME1"
        H3K27AC="$KIDNEY_H3K27AC"
        H3K9ME3="$KIDNEY_H3K9ME3"
    fi

    run_promoter "$ATAC" "$OUTROOT/promoters_ATAC/$tissue" "${tissue}_promoters_ATAC" "$LABEL ATAC accessibility at promoters" "magma" "Mean ATAC signal"
    run_promoter "$WGBS" "$OUTROOT/promoters_WGBS/$tissue" "${tissue}_promoters_WGBS" "$LABEL DNA methylation at promoters" "viridis" "Mean methylation signal"
    run_promoter "$H3K4ME3" "$OUTROOT/promoters_H3K4me3/$tissue" "${tissue}_promoters_H3K4me3" "$LABEL H3K4me3 signal at promoters" "Blues" "Mean H3K4me3 signal"

    run_reference "$ENHANCERS" "$ATAC" "$OUTROOT/enhancers_ATAC/$tissue" "${tissue}_enhancers_ATAC" "$LABEL ATAC accessibility at enhancers" "magma" "Mean ATAC signal" "Centre" "Distance from enhancer centre (bp)"
    run_reference "$ENHANCERS" "$WGBS" "$OUTROOT/enhancers_WGBS/$tissue" "${tissue}_enhancers_WGBS" "$LABEL DNA methylation at enhancers" "viridis" "Mean methylation signal" "Centre" "Distance from enhancer centre (bp)"
    run_reference "$ENHANCERS" "$H3K4ME1" "$OUTROOT/enhancers_H3K4me1/$tissue" "${tissue}_enhancers_H3K4me1" "$LABEL H3K4me1 signal at enhancers" "YlGn" "Mean H3K4me1 signal" "Centre" "Distance from enhancer centre (bp)"
    run_reference "$ENHANCERS" "$H3K27AC" "$OUTROOT/enhancers_H3K27ac/$tissue" "${tissue}_enhancers_H3K27ac" "$LABEL H3K27ac signal at enhancers" "YlOrBr" "Mean H3K27ac signal" "Centre" "Distance from enhancer centre (bp)"

    run_pmd "$ATAC" "$OUTROOT/PMDs_ATAC/$tissue" "${tissue}_PMDs_ATAC" "$LABEL ATAC accessibility across PMDs" "magma" "Mean ATAC signal"
    run_pmd "$WGBS" "$OUTROOT/PMDs_WGBS/$tissue" "${tissue}_PMDs_WGBS" "$LABEL DNA methylation across PMDs" "viridis" "Mean methylation signal"
    run_pmd "$H3K9ME3" "$OUTROOT/PMDs_H3K9me3/$tissue" "${tissue}_PMDs_H3K9me3" "$LABEL H3K9me3 signal across PMDs" "Purples" "Mean H3K9me3 signal"
done

echo "Finished. Tissue-matched results are in: $OUTROOT"
