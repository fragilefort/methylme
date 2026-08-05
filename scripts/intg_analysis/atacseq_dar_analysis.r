#!/usr/bin/env Rscript

library(dplyr)
library(GenomicRanges)
library(rtracklayer)

wgbs_table <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_tiling.csv"
atac_master <- "/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/differential_peaks_annotated.tsv"

diff_df <- read.csv(wgbs_table)

dmr_df <- diff_df %>%
    filter(comb.p.adj.fdr < 0.05 & abs(mean.mean.diff) >= 0.2)

dmr_gr <- makeGRangesFromDataFrame(
    dmr_df,
    seqnames.field = "Chromosome",
    start.field = "Start",
    end.field = "End",
    keep.extra.columns = TRUE
)

cat("Total filtered DMRs:", length(dmr_gr), "\n\n")

cat("=== Filtering Master ATAC TSV for DARs ===\n")
atac_df <- read.delim(atac_master)

# Filter for statistically significant Differentially Accessible Regions (FDR < 0.05 & |log2FC| >= 1)
dars_df <- atac_df %>%
    filter(FDR < 0.05 & abs(log2FoldChange) >= 1)

# Separate by tissue directionality
kidney_dars_df <- dars_df %>% filter(log2FoldChange > 1)
liver_dars_df <- dars_df %>% filter(log2FoldChange < -1)

# Convert to GRanges objects
dars_gr <- makeGRangesFromDataFrame(
    dars_df,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
)

kidney_dars_gr <- makeGRangesFromDataFrame(
    kidney_dars_df,
    seqnames.field = "chr",
    start.field    = "start",
    end.field      = "end"
)

liver_dars_gr <- makeGRangesFromDataFrame(
    liver_dars_df,
    seqnames.field = "chr",
    start.field    = "start",
    end.field      = "end"
)

# Merge overlapping DAR intervals into a clean unified set
unified_dars_gr <- reduce(dars_gr)

cat("Kidney-accessible DARs:", length(kidney_dars_gr), "\n")
cat("Liver-accessible DARs: ", length(liver_dars_gr), "\n")
cat("Total filtered DARs:   ", length(dars_gr), "\n")
cat("Unified merged DARs:   ", length(unified_dars_gr), "\n\n")

# Compare DMRs with Filtered DARs
dmrs_in_kidney_dars <- subsetByOverlaps(dmr_gr, kidney_dars_gr)
dmrs_in_liver_dars <- subsetByOverlaps(dmr_gr, liver_dars_gr)
dmrs_in_all_dars <- subsetByOverlaps(dmr_gr, unified_dars_gr)

cat("DMRs overlapping Kidney-accessible peaks:", length(dmrs_in_kidney_dars), "\n")
cat("DMRs overlapping Liver-accessible peaks: ", length(dmrs_in_liver_dars), "\n")
cat("Total DMRs overlapping any filtered DAR: ", length(dmrs_in_all_dars), "\n")
