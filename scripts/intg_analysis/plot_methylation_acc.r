#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

wgbs_table <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_tiling.csv"
atac_master <- "/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/differential_peaks_annotated.tsv"
out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/atac-wgbs"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

diff_meth <- read.csv(wgbs_table)
atac_master <- read.delim(atac_master)

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

overlaps <- findOverlaps(dmr_gr, atac_gr)

matched_data <- data.frame(
    meth_diff   = dmr_gr$mean.mean.diff[queryHits(overlaps)], # Methylation difference (Kidney - Liver)
    atac_log2fc = atac_gr$log2FoldChange[subjectHits(overlaps)], # Accessibility fold-change
    atac_fdr    = atac_gr$padj[subjectHits(overlaps)]
) %>%
    mutate(
        Accessibility_Status = case_when(
            atac_fdr < 0.05 & atac_log2fc > 1 ~ "Kidney Accessible",
            atac_fdr < 0.05 & atac_log2fc < -1 ~ "Liver Accessible",
            TRUE ~ "Not Differentially Accessible"
        ),
        # 4 Quadrant biological categories (Kidney vs Liver active states)
        Regulatory_State = case_when(
            meth_diff < 0 & atac_log2fc > 0 ~ "Kidney Active (Hypomethylated & Open ATAC)",
            meth_diff > 0 & atac_log2fc < 0 ~ "Liver Active (Kidney Hypermethylated & Closed ATAC)",
            meth_diff > 0 & atac_log2fc > 0 ~ "Atypical (Hypermethylated & Open ATAC)",
            meth_diff < 0 & atac_log2fc < 0 ~ "Atypical (Hypomethylated & Closed ATAC)"
        )
    )

pearson_res <- cor.test(matched_data$meth_diff, matched_data$atac_log2fc, method = "pearson")
spearman_res <- cor.test(matched_data$meth_diff, matched_data$atac_log2fc, method = "spearman")

cat("\n========================================================\n")
cat(" CORRELATION STATISTICS (Overlapping DMR-ATAC Regions)\n")
cat("========================================================\n")
cat(sprintf("Total Overlapping Regions: %d\n", nrow(matched_data)))
cat(sprintf("Pearson  r: %6.3f  (p-value: %e)\n", pearson_res$estimate, pearson_res$p.value))
cat(sprintf("Spearman rho: %6.3f  (p-value: %e)\n", spearman_res$estimate, spearman_res$p.value))
cat("========================================================\n\n")

cat("=== REGION BREAKDOWN BY TISSUE STATE (QUADRANTS) ===\n")
quadrant_summary <- matched_data %>%
    group_by(Regulatory_State) %>%
    summarise(
        Count = n(),
        Percentage = round((n() / nrow(matched_data)) * 100, 2)
    )
print(as.data.frame(quadrant_summary))
cat("\n")

stats_text <- sprintf(
    "Pearson r = %.3f (p = %.2e)\nSpearman rho = %.3f (p = %.2e)",
    pearson_res$estimate, pearson_res$p.value,
    spearman_res$estimate, spearman_res$p.value
)

scatter_plot <- ggplot(matched_data, aes(x = meth_diff, y = atac_log2fc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = Regulatory_State), alpha = 0.65, size = 2) +
    geom_smooth(method = "lm", color = "black", linetype = "solid", se = TRUE) +
    annotate("label",
        x = min(matched_data$meth_diff) * 0.85, y = max(matched_data$atac_log2fc) * 0.85,
        label = stats_text, fill = "white", size = 3.8, hjust = 0
    ) +
    scale_color_manual(values = c(
        "Kidney Active (Hypomethylated & Open ATAC)" = "#D95F02", # Orange
        "Liver Active (Kidney Hypermethylated & Closed ATAC)" = "#7570B3", # Purple
        "Atypical (Hypermethylated & Open ATAC)" = "#E7298A", # Pink
        "Atypical (Hypomethylated & Closed ATAC)" = "#66A61E" # Green
    )) +
    theme_minimal(base_size = 13) +
    theme(
        legend.position = "bottom",
        legend.title = element_blank()
    ) +
    guides(color = guide_legend(ncol = 2)) +
    labs(
        title = "DNA Methylation Difference vs. Chromatin Accessibility",
        subtitle = "Regions overlapping DMRs and ATAC peaks (Kidney vs. Liver)",
        x = "Methylation Difference (Kidney - Liver)",
        y = "ATAC-seq log2 Fold Change (Kidney / Liver)"
    )

ggsave(file.path(out_dir, "methylation_vs_accessibility_scatter.png"), plot = scatter_plot, width = 8.5, height = 7, dpi = 300)

# Boxplot (Methylation Difference by Accessibility Status)
box_plot <- ggplot(matched_data, aes(x = Accessibility_Status, y = meth_diff, fill = Accessibility_Status)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c(
        "Kidney Accessible"             = "#D95F02",
        "Liver Accessible"              = "#7570B3",
        "Not Differentially Accessible" = "gray70"
    )) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none") +
    labs(
        title = "Methylation Differences across Accessibility States",
        x = "Chromatin Accessibility Status",
        y = "Methylation Difference (Kidney - Liver)"
    )

ggsave(file.path(out_dir, "methylation_by_accessibility_boxplot.png"), plot = box_plot, width = 7.5, height = 5.5, dpi = 300)
