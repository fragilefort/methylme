#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

diff_meth <- read.csv("WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_tiling.csv")
atac_master <- read.delim("/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/differential_peaks_annotated.tsv")
out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/atac-wgbs"

# Filter DMRs (FDR < 0.05 & |mean.mean.diff| >= 0.2)
dmrs <- diff_meth %>%
    filter(comb.p.adj.fdr < 0.05 & abs(mean.mean.diff) >= 0.2)

# Convert to GRanges
dmr_gr <- makeGRangesFromDataFrame(
    dmrs,
    seqnames.field = "Chromosome",
    start.field = "Start",
    end.field = "End",
    keep.extra.columns = TRUE
)

atac_gr <- makeGRangesFromDataFrame(
    atac_master,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
)

# Find Overlaps between DMRs and ATAC Peaks
overlaps <- findOverlaps(dmr_gr, atac_gr)

# Build a joint dataframe for overlapping regions
matched_data <- data.frame(
    meth_diff   = dmr_gr$mean.mean.diff[queryHits(overlaps)], # Methylation difference (Kidney - Liver)
    atac_log2fc = atac_gr$log2FoldChange[subjectHits(overlaps)], # Accessibility fold-change
    atac_fdr    = atac_gr$FDR[subjectHits(overlaps)]
) %>%
    mutate(
        Accessibility_Status = case_when(
            atac_fdr < 0.05 & atac_log2fc > 1 ~ "Kidney Accessible",
            atac_fdr < 0.05 & atac_log2fc < -1 ~ "Liver Accessible",
            TRUE ~ "Not Differentially Accessible"
        )
    )

# Scatter Plot (Methylation Difference vs ATAC Fold Change)
scatter_plot <- ggplot(matched_data, aes(x = meth_diff, y = atac_log2fc)) +
    geom_point(aes(color = Accessibility_Status), alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c(
        "Kidney Accessible" = "#D95F02",
        "Liver Accessible" = "#7570B3",
        "Not Differentially Accessible" = "gray70"
    )) +
    theme_minimal(base_size = 14) +
    labs(
        title = "DNA Methylation vs. Chromatin Accessibility",
        subtitle = "Regions overlapping DMRs and ATAC peaks",
        x = "Methylation Difference (Kidney - Liver)",
        y = "ATAC-seq log2 Fold Change",
        color = "ATAC Status"
    )

ggsave(file.path(out_dir, "methylation_vs_accessibility_scatter.png"), plot = scatter_plot, width = 8, height = 6, dpi = 300)

box_plot <- ggplot(matched_data, aes(x = Accessibility_Status, y = meth_diff, fill = Accessibility_Status)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c(
        "Kidney Accessible" = "#D95F02",
        "Liver Accessible" = "#7570B3",
        "Not Differentially Accessible" = "gray70"
    )) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none") +
    labs(
        title = "Methylation Differences across Accessibility States",
        x = "Chromatin Accessibility Status",
        y = "Methylation Difference (Kidney - Liver)"
    )

ggsave(file.path(out_dir, "methylation_by_accessibility_boxplot.png"), plot = box_plot, width = 7, height = 5, dpi = 300)
