#!/usr/bin/env bash
#SBATCH --job-name=integrative_deeptools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=integrative_deeptools_%j.out
#SBATCH --error=integrative_deeptools_%j.err

set -euo pipefail

# Unified deepTools workflow
# Promoters/enhancers compare Liver and Kidney together.
# PMDs use matched tissue only: Liver PMDs with Liver signal, Kidney PMDs with Kidney signal.

THREADS="${SLURM_CPUS_PER_TASK:-1}"

ROI_DIR="/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/integrative_analysis_deeptools_ROI"
OUTROOT="$ROI_DIR/unified_deeptools_option1_14_5_REP1"

GENES="/vol/COMPEPIWS/pipelines/references/mm10_reduced_chr18_chr19_genes.bed"
LIVER_PMDS="$ROI_DIR/liver_PMDs.bed"
KIDNEY_PMDS="$ROI_DIR/kidney_PMDs.bed"
LIVER_ENHANCERS="$ROI_DIR/liver_ChromHMM_enhancers.bed"
KIDNEY_ENHANCERS="$ROI_DIR/kidney_ChromHMM_enhancers.bed"

# Original ATAC BigWigs: embryonic day 14.5, replicate 1.
LIVER_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/liver_14.5_REP1.mLb.clN.bigWig"
KIDNEY_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq2/kidney_14.5_REP1.mLb.clN.bigWig"

# Tissue-merged WGBS methylation signal tracks.
LIVER_WGBS="$ROI_DIR/signals/liver_WGBS_merged_signal.bw"
KIDNEY_WGBS="$ROI_DIR/signals/kidney_WGBS_merged_signal.bw"

# Original ChIP-seq BigWigs: embryonic day 14.5, replicate 1.
LIVER_H3K4ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K4me3_REP1.mLb.clN.bigWig"
KIDNEY_H3K4ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me3_REP1.mLb.clN.bigWig"
LIVER_H3K4ME1="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K4me1_REP1.mLb.clN.bigWig"
KIDNEY_H3K4ME1="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me1_REP1.mLb.clN.bigWig"
LIVER_H3K27AC="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K27ac_REP1.mLb.clN.bigWig"
KIDNEY_H3K27AC="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K27ac_REP1.mLb.clN.bigWig"
LIVER_H3K9ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K9me3_REP1.mLb.clN.bigWig"
KIDNEY_H3K9ME3="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K9me3_REP1.mLb.clN.bigWig"

for cmd in computeMatrix plotHeatmap; do
    command -v "$cmd" >/dev/null || {
        echo "Missing command: $cmd" >&2
        exit 1
    }
done

for file in \
    "$GENES" "$LIVER_PMDS" "$KIDNEY_PMDS" "$LIVER_ENHANCERS" "$KIDNEY_ENHANCERS" \
    "$LIVER_ATAC" "$KIDNEY_ATAC" "$LIVER_WGBS" "$KIDNEY_WGBS" \
    "$LIVER_H3K4ME3" "$KIDNEY_H3K4ME3" "$LIVER_H3K4ME1" "$KIDNEY_H3K4ME1" \
    "$LIVER_H3K27AC" "$KIDNEY_H3K27AC" "$LIVER_H3K9ME3" "$KIDNEY_H3K9ME3"
do
    [[ -s "$file" ]] || {
        echo "Missing or empty input: $file" >&2
        exit 1
    }
done

mkdir -p \
    "$OUTROOT/promoters_ATAC" \
    "$OUTROOT/promoters_WGBS" \
    "$OUTROOT/promoters_H3K4me3" \
    "$OUTROOT/enhancers_ATAC" \
    "$OUTROOT/enhancers_WGBS" \
    "$OUTROOT/enhancers_H3K4me1" \
    "$OUTROOT/enhancers_H3K27ac" \
    "$OUTROOT/PMDs_ATAC/liver" \
    "$OUTROOT/PMDs_ATAC/kidney" \
    "$OUTROOT/PMDs_WGBS/liver" \
    "$OUTROOT/PMDs_WGBS/kidney" \
    "$OUTROOT/PMDs_H3K9me3/liver" \
    "$OUTROOT/PMDs_H3K9me3/kidney"

# Produces a figure with the average profile above the heatmap and a colour bar.
render_heatmap() {
    local matrix="$1"
    local outdir="$2"
    local prefix="$3"
    local title="$4"
    local colour_map="$5"
    shift 5

    plotHeatmap \
        --matrixFile "$matrix" \
        --outFileName "$outdir/${prefix}_heatmap_profile.png" \
        --outFileNameMatrix "$outdir/${prefix}_values.tab" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap "$colour_map" \
        --whatToShow "plot, heatmap and colorbar" \
        --plotTitle "$title" \
        "$@"

    plotHeatmap \
        --matrixFile "$matrix" \
        --outFileName "$outdir/${prefix}_heatmap_profile.pdf" \
        --sortRegions descend \
        --sortUsing mean \
        --colorMap "$colour_map" \
        --whatToShow "plot, heatmap and colorbar" \
        --plotTitle "$title" \
        "$@"
}

# -------------------------------------------------------------------
# PROMOTERS: TSS -500 bp to +1,500 bp; Liver and Kidney together.
# -------------------------------------------------------------------
computeMatrix reference-point --referencePoint TSS -b 500 -a 1500 --binSize 25 -p "$THREADS" \
    -R "$GENES" \
    -S "$LIVER_ATAC" "$KIDNEY_ATAC" \
    --samplesLabel "Liver E14.5 REP1 ATAC" "Kidney E14.5 REP1 ATAC" \
    --missingDataAsZero \
    --outFileName "$OUTROOT/promoters_ATAC/matrix.gz" \
    --outFileSortedRegions "$OUTROOT/promoters_ATAC/sorted_regions.bed"
render_heatmap "$OUTROOT/promoters_ATAC/matrix.gz" "$OUTROOT/promoters_ATAC" "promoters_ATAC" "ATAC accessibility at promoters" "YlOrRd" --refPointLabel "TSS"

computeMatrix reference-point --referencePoint TSS -b 500 -a 1500 --binSize 25 -p "$THREADS" \
    -R "$GENES" \
    -S "$LIVER_WGBS" "$KIDNEY_WGBS" \
    --samplesLabel "Liver merged WGBS" "Kidney merged WGBS" \
    --missingDataAsZero \
    --outFileName "$OUTROOT/promoters_WGBS/matrix.gz" \
    --outFileSortedRegions "$OUTROOT/promoters_WGBS/sorted_regions.bed"
render_heatmap "$OUTROOT/promoters_WGBS/matrix.gz" "$OUTROOT/promoters_WGBS" "promoters_WGBS" "DNA methylation at promoters" "RdYlBu_r" --refPointLabel "TSS"

computeMatrix reference-point --referencePoint TSS -b 500 -a 1500 --binSize 25 -p "$THREADS" \
    -R "$GENES" \
    -S "$LIVER_H3K4ME3" "$KIDNEY_H3K4ME3" \
    --samplesLabel "Liver E14.5 H3K4me3" "Kidney E14.5 H3K4me3" \
    --missingDataAsZero \
    --outFileName "$OUTROOT/promoters_H3K4me3/matrix.gz" \
    --outFileSortedRegions "$OUTROOT/promoters_H3K4me3/sorted_regions.bed"
render_heatmap "$OUTROOT/promoters_H3K4me3/matrix.gz" "$OUTROOT/promoters_H3K4me3" "promoters_H3K4me3" "H3K4me3 signal at promoters" "Blues" --refPointLabel "TSS"

# -------------------------------------------------------------------
# ENHANCERS: centre +/-1 kb; Liver and Kidney regions and signals together.
# -------------------------------------------------------------------
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 --binSize 25 -p "$THREADS" \
    -R "$LIVER_ENHANCERS" "$KIDNEY_ENHANCERS" \
    -S "$LIVER_ATAC" "$KIDNEY_ATAC" \
    --samplesLabel "Liver E14.5 REP1 ATAC" "Kidney E14.5 REP1 ATAC" \
    --missingDataAsZero \
    --outFileName "$OUTROOT/enhancers_ATAC/matrix.gz" \
    --outFileSortedRegions "$OUTROOT/enhancers_ATAC/sorted_regions.bed"
render_heatmap "$OUTROOT/enhancers_ATAC/matrix.gz" "$OUTROOT/enhancers_ATAC" "enhancers_ATAC" "ATAC accessibility around ChromHMM enhancers" "YlOrRd" --refPointLabel "Enhancer centre" --regionsLabel "Liver enhancers" "Kidney enhancers"

computeMatrix reference-point --referencePoint center -b 1000 -a 1000 --binSize 25 -p "$THREADS" \
    -R "$LIVER_ENHANCERS" "$KIDNEY_ENHANCERS" \
    -S "$LIVER_WGBS" "$KIDNEY_WGBS" \
    --samplesLabel "Liver merged WGBS" "Kidney merged WGBS" \
    --missingDataAsZero \
    --outFileName "$OUTROOT/enhancers_WGBS/matrix.gz" \
    --outFileSortedRegions "$OUTROOT/enhancers_WGBS/sorted_regions.bed"
render_heatmap "$OUTROOT/enhancers_WGBS/matrix.gz" "$OUTROOT/enhancers_WGBS" "enhancers_WGBS" "DNA methylation around ChromHMM enhancers" "RdYlBu_r" --refPointLabel "Enhancer centre" --regionsLabel "Liver enhancers" "Kidney enhancers"

computeMatrix reference-point --referencePoint center -b 1000 -a 1000 --binSize 25 -p "$THREADS" \
    -R "$LIVER_ENHANCERS" "$KIDNEY_ENHANCERS" \
    -S "$LIVER_H3K4ME1" "$KIDNEY_H3K4ME1" \
    --samplesLabel "Liver E14.5 H3K4me1" "Kidney E14.5 H3K4me1" \
    --missingDataAsZero \
    --outFileName "$OUTROOT/enhancers_H3K4me1/matrix.gz" \
    --outFileSortedRegions "$OUTROOT/enhancers_H3K4me1/sorted_regions.bed"
render_heatmap "$OUTROOT/enhancers_H3K4me1/matrix.gz" "$OUTROOT/enhancers_H3K4me1" "enhancers_H3K4me1" "H3K4me1 signal around ChromHMM enhancers" "Greens" --refPointLabel "Enhancer centre" --regionsLabel "Liver enhancers" "Kidney enhancers"

computeMatrix reference-point --referencePoint center -b 1000 -a 1000 --binSize 25 -p "$THREADS" \
    -R "$LIVER_ENHANCERS" "$KIDNEY_ENHANCERS" \
    -S "$LIVER_H3K27AC" "$KIDNEY_H3K27AC" \
    --samplesLabel "Liver E14.5 H3K27ac" "Kidney E14.5 H3K27ac" \
    --missingDataAsZero \
    --outFileName "$OUTROOT/enhancers_H3K27ac/matrix.gz" \
    --outFileSortedRegions "$OUTROOT/enhancers_H3K27ac/sorted_regions.bed"
render_heatmap "$OUTROOT/enhancers_H3K27ac/matrix.gz" "$OUTROOT/enhancers_H3K27ac" "enhancers_H3K27ac" "H3K27ac signal around ChromHMM enhancers" "Oranges" --refPointLabel "Enhancer centre" --regionsLabel "Liver enhancers" "Kidney enhancers"

# -------------------------------------------------------------------
# PMDS: matched tissue only. No Liver signal on Kidney PMDs or vice versa.
# -------------------------------------------------------------------
for tissue in liver kidney; do
    if [[ "$tissue" == "liver" ]]; then
        PMDS="$LIVER_PMDS"
        ATAC="$LIVER_ATAC"
        WGBS="$LIVER_WGBS"
        H3K9ME3="$LIVER_H3K9ME3"
        LABEL="Liver"
    else
        PMDS="$KIDNEY_PMDS"
        ATAC="$KIDNEY_ATAC"
        WGBS="$KIDNEY_WGBS"
        H3K9ME3="$KIDNEY_H3K9ME3"
        LABEL="Kidney"
    fi

    computeMatrix scale-regions -b 10000 -a 10000 --regionBodyLength 5000 --binSize 100 -p "$THREADS" \
        -R "$PMDS" \
        -S "$ATAC" \
        --samplesLabel "$LABEL E14.5 REP1 ATAC" \
        --missingDataAsZero \
        --outFileName "$OUTROOT/PMDs_ATAC/$tissue/matrix.gz" \
        --outFileSortedRegions "$OUTROOT/PMDs_ATAC/$tissue/sorted_regions.bed"
    render_heatmap "$OUTROOT/PMDs_ATAC/$tissue/matrix.gz" "$OUTROOT/PMDs_ATAC/$tissue" "${tissue}_PMDs_ATAC" "$LABEL ATAC accessibility across $LABEL PMDs" "YlOrRd" --startLabel "PMD start" --endLabel "PMD end" --regionsLabel "$LABEL PMDs"

    computeMatrix scale-regions -b 10000 -a 10000 --regionBodyLength 5000 --binSize 100 -p "$THREADS" \
        -R "$PMDS" \
        -S "$WGBS" \
        --samplesLabel "$LABEL merged WGBS" \
        --missingDataAsZero 
