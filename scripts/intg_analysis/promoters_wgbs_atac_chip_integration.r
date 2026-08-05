#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

wgbs_promoters <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"
atac_master    <- "/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/differential_peaks_annotated.tsv"
kidney_seg_bed <- "/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation/kidney_15_segments.bed"
liver_seg_bed  <- "/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation/liver_15_segments.bed"
out_dir        <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/atac-wgbs"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

cat("=== 1. Loading WGBS Promoters & ATAC-seq Peaks ===\n")
diff_promoters <- read.csv(wgbs_promoters)
atac_master    <- read.delim(atac_master)

# Filter Differentially Methylated Promoters (FDR < 0.05 & |mean.mean.diff| >= 0.2)
promoter_dmrs <- diff_promoters %>%
    filter(comb.p.adj.fdr < 0.05 & abs(mean.mean.diff) >= 0.2)

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

kidney_seg <- import(kidney_seg_bed)
liver_seg  <- import(liver_seg_bed)

# Filter specifically for Heterochromatin states (Het_P and Het_S)
kidney_het_gr <- kidney_seg[grepl("Het", kidney_seg$name)]
liver_het_gr  <- liver_seg[grepl("Het", liver_seg$name)]

cat(sprintf("Kidney Heterochromatin regions: %d\n", length(kidney_het_gr)))
cat(sprintf("Liver Heterochromatin regions:  %d\n\n", length(liver_het_gr)))

# Find Overlaps & Layer Tissue-Specific Heterochromatin
cat("=== 3. Overlapping Promoters, ATAC, and Heterochromatin States ===\n")
overlaps <- findOverlaps(promoter_gr, atac_gr)

# Subset promoter GRanges for overlapping pairs
matched_promoters_gr <- promoter_gr[queryHits(overlaps)]

# Query overlap with heterochromatin in each tissue
matched_promoters_gr$is_kidney_het <- overlapsAny(matched_promoters_gr, kidney_het_gr)
matched_promoters_gr$is_liver_het  <- overlapsAny(matched_promoters_gr, liver_het_gr)

# Build unified dataframe
matched_df <- data.frame(
    meth_diff     = matched_promoters_gr$mean.mean.diff,  # Kidney - Liver Meth
    atac_log2fc   = atac_gr$log2FoldChange[subjectHits(overlaps)],
    is_kidney_het = matched_promoters_gr$is_kidney_het,
    is_liver_het  = matched_promoters_gr$is_liver_het
) %>%
    mutate(
        Regulatory_State = case_when(
            meth_diff < 0 & atac_log2fc > 0 ~ "Kidney Active (Hypo-meth & Open ATAC)",
            meth_diff > 0 & atac_log2fc < 0 ~ "Liver Active (Kidney Hyper & Closed ATAC)",
            meth_diff > 0 & atac_log2fc > 0 ~ "Atypical (Hyper-meth & Open ATAC)",
            meth_diff < 0 & atac_log2fc < 0 ~ "Atypical (Hypo-meth & Closed ATAC)"
        ),
        # Evaluate tissue-matching heterochromatin status
        In_Tissue_Heterochromatin = case_when(
            Regulatory_State == "Kidney Active (Hypo-meth & Open ATAC)"  ~ is_kidney_het,
            Regulatory_State == "Liver Active (Kidney Hyper & Closed ATAC)" ~ is_liver_het,
            TRUE                                                          ~ (is_kidney_het | is_liver_het)
        )
    )

cat("========================================================================\n")
cat(" PROMOTER HETEROCHROMATIN OVERLAP BREAKDOWN BY REGULATORY STATE\n")
cat("========================================================================\n")

summary_table <- matched_df %>%
    group_by(Regulatory_State) %>%
    summarise(
        Total_Promoters    = n(),
        In_Heterochromatin = sum(In_Tissue_Heterochromatin),
        Not_Heterochromatin= sum(!In_Tissue_Heterochromatin),
        `Pct_Heterochromatin (%)` = round((sum(In_Tissue_Heterochromatin) / n()) * 100, 2),
        `Pct_Not_Het (%)`         = round((sum(!In_Tissue_Heterochromatin) / n()) * 100, 2)
    )

print(as.data.frame(summary_table))

bar_plot <- ggplot(matched_df, aes(x = Regulatory_State, fill = In_Tissue_Heterochromatin)) +
    geom_bar(position = "fill", alpha = 0.85, width = 0.6) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(
        values = c("FALSE" = "#2B8CBE", "TRUE" = "#E41A1C"),
        labels = c("FALSE" = "Euchromatin (Not Het)", "TRUE" = "Heterochromatin (Het_P/Het_S)")
    ) +
    theme_minimal(base_size = 13) +
    theme(
        axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank()
    ) +
    labs(
        title = "Heterochromatin Proportion Across Promoter Regulatory States",
        subtitle = "Promoter DMRs overlapping ATAC peaks evaluated against ChIP-seq ChromHMM",
        x = "Promoter Regulatory Category",
        y = "Percentage of Promoters"
    )

output_bar_file <- file.path(out_dir, "promoter_heterochromatin_proportions_barplot.png")
ggsave(output_bar_file, plot = bar_plot, width = 8.5, height = 6, dpi = 300)

scatter_het <- ggplot(matched_df, aes(x = meth_diff, y = atac_log2fc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = In_Tissue_Heterochromatin), alpha = 0.7, size = 2) +
    scale_color_manual(
        values = c("FALSE" = "#2B8CBE", "TRUE" = "#E41A1C"),
        labels = c("FALSE" = "Euchromatin (Not Het)", "TRUE" = "Heterochromatin (Het_P/Het_S)")
    ) +
    theme_minimal(base_size = 13) +
    theme(
        legend.position = "bottom",
        legend.title = element_blank()
    ) +
    labs(
        title = "Promoter Methylation vs. Accessibility (ChIP-seq Heterochromatin Overlay)",
        x = "Promoter Methylation Difference (Kidney - Liver)",
        y = "ATAC-seq log2 Fold Change (Kidney / Liver)"
    )

output_scatter_file <- file.path(out_dir, "promoter_methylation_vs_accessibility_het_overlay.png")
ggsave(output_scatter_file, plot = scatter_het, width = 8, height = 6.5, dpi = 300)

cat(sprintf("Plots saved to:\n1. %s\n2. %s\n", output_bar_file, output_scatter_file))
