#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(GenomicRanges)

wgbs_promoters <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"
atac_master    <- "/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/differential_peaks_annotated.tsv"
out_dir        <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/atac-wgbs"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

cat("=== Loading Promoter WGBS and Master ATAC Data ===\n")
diff_promoters <- read.csv(wgbs_promoters)
atac_master    <- read.delim(atac_master)

# Filter Differentially Methylated Promoters (FDR < 0.05 & |mean.mean.diff| >= 0.2)
promoter_dmrs <- diff_promoters %>%
    filter(comb.p.adj.fdr < 0.05 & abs(mean.mean.diff) >= 0.2)

# Convert to GRanges
promoter_gr <- makeGRangesFromDataFrame(
    promoter_dmrs,
    seqnames.field = "Chromosome",
    start.field    = "Start",
    end.field      = "End",
    keep.extra.columns = TRUE
)

atac_gr <- makeGRangesFromDataFrame(
    atac_master,
    seqnames.field = "chr",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE
)

overlaps <- findOverlaps(promoter_gr, atac_gr)

matched_data <- data.frame(
    meth_diff   = promoter_gr$mean.mean.diff[queryHits(overlaps)],  # Methylation difference (Kidney - Liver)
    atac_log2fc = atac_gr$log2FoldChange[subjectHits(overlaps)]    # Accessibility fold-change
) %>%
    mutate(
        Regulatory_State = case_when(
            meth_diff < 0 & atac_log2fc > 0 ~ "Kidney Active Promoter (Hypomethylated & Open ATAC)",
            meth_diff > 0 & atac_log2fc < 0 ~ "Liver Active Promoter (Kidney Hypermethylated & Closed ATAC)",
            meth_diff > 0 & atac_log2fc > 0 ~ "Atypical (Hypermethylated & Open ATAC)",
            meth_diff < 0 & atac_log2fc < 0 ~ "Atypical (Hypomethylated & Closed ATAC)"
        )
    )

cat("\n=== PROMOTER BREAKDOWN BY TISSUE STATE (QUADRANTS) ===\n")
quadrant_summary <- matched_data %>%
    group_by(Regulatory_State) %>%
    summarise(
        Count = n(),
        Percentage = round((n() / nrow(matched_data)) * 100, 2)
    )
print(as.data.frame(quadrant_summary))
cat("\n")

scatter_plot <- ggplot(matched_data, aes(x = meth_diff, y = atac_log2fc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = Regulatory_State), alpha = 0.7, size = 2.2) +
    scale_color_manual(values = c(
        "Kidney Active Promoter (Hypomethylated & Open ATAC)"            = "#D95F02", # Orange
        "Liver Active Promoter (Kidney Hypermethylated & Closed ATAC)"  = "#7570B3", # Purple
        "Atypical (Hypermethylated & Open ATAC)"                        = "#E7298A", # Pink
        "Atypical (Hypomethylated & Closed ATAC)"                       = "#66A61E"  # Green
    )) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title    = element_text(size = 11, face = "bold"), # Smaller title size
        plot.subtitle = element_text(size = 9.5),
        legend.position = "bottom",
        legend.title  = element_blank()
    ) +
    guides(color = guide_legend(ncol = 2)) +
    labs(
        title = "Promoter DNA Methylation Difference vs. Chromatin Accessibility",
        subtitle = "Overlapping Promoter DMRs and ATAC peaks (Kidney vs. Liver)",
        x = "Promoter Methylation Difference (Kidney - Liver)",
        y = "ATAC-seq log2 Fold Change (Kidney / Liver)"
    )

output_file <- file.path(out_dir, "promoter_methylation_vs_accessibility_scatter.png")
ggsave(output_file, plot = scatter_plot, width = 8.5, height = 7, dpi = 300)

cat(sprintf("Promoter scatter plot successfully saved to:\n%s\n", output_file))
