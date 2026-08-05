#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(GenomicRanges)

wgbs_mystery_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/diff_meth_mystery_regions_annotated.csv"
atac_master_file  <- "/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/differential_peaks_annotated.tsv"
out_dir           <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/atac-wgbs"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

diff_meth   <- read.csv(wgbs_mystery_file)
atac_master <- read.delim(atac_master_file)

cat(sprintf("Loaded Mystery Regions: %d rows\n", nrow(diff_meth)))
cat(sprintf("Loaded ATAC Peaks:     %d rows\n\n", nrow(atac_master)))

# Convert to GRanges using exact chromosome/start/end column names
dmr_gr <- makeGRangesFromDataFrame(
    diff_meth,
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

# Exact Coordinate Matching (type = "equal")
overlaps <- findOverlaps(dmr_gr, atac_gr, type = "equal")

cat(sprintf("Total exact matched mystery regions: %d\n\n", length(overlaps)))

if (length(overlaps) == 0) {
    stop("No exact matches found! Please verify chromosome notation between both tables.")
}

matched_dmr  <- diff_meth[queryHits(overlaps), ]
matched_atac <- atac_master[subjectHits(overlaps), ]

# Kidney: mean.mean.g1 (WGBS) vs. meanVstCountGrp1_kidney (ATAC)
kidney_df <- data.frame(
    Region_ID     = paste0(matched_dmr$Chromosome, ":", matched_dmr$Start, "-", matched_dmr$End),
    Tissue        = "Kidney",
    Methylation   = matched_dmr$mean.mean.g1,
    Accessibility = matched_atac$meanVstCountGrp1_kidney
)

# Liver: mean.mean.g2 (WGBS) vs. meanVstCountGrp2_liver (ATAC)
liver_df <- data.frame(
    Region_ID     = paste0(matched_dmr$Chromosome, ":", matched_dmr$Start, "-", matched_dmr$End),
    Tissue        = "Liver",
    Methylation   = matched_dmr$mean.mean.g2,
    Accessibility = matched_atac$meanVstCountGrp2_liver
)

plot_df <- rbind(kidney_df, liver_df)

stats_df <- plot_df %>%
    group_by(Tissue) %>%
    summarise(
        pearson_r = cor(Methylation, Accessibility, method = "pearson", use = "complete.obs"),
        p_val     = cor.test(Methylation, Accessibility, method = "pearson")$p.value,
        .groups   = "drop"
    ) %>%
    mutate(
        label = sprintf("Pearson r = %.3f\np = %.2e", pearson_r, p_val)
    )

cat("========================================================\n")
cat(" PEARSON CORRELATION (Mystery Regions)\n")
cat("========================================================\n")
print(as.data.frame(stats_df))
cat("========================================================\n\n")

p <- ggplot(plot_df, aes(x = Methylation, y = Accessibility)) +
    geom_point(aes(color = Tissue), alpha = 0.6, size = 1.8) +
    geom_smooth(method = "lm", color = "black", linetype = "solid", se = TRUE) +
    geom_text(
        data = stats_df,
        aes(x = -Inf, y = Inf, label = label),
        hjust = -0.1, vjust = 1.2,
        size = 3.8,
        inherit.aes = FALSE
    ) +
    facet_wrap(~ Tissue, scales = "free_y") +
    scale_color_manual(values = c("Kidney" = "#D95F02", "Liver" = "#7570B3")) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title       = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "gray92", color = NA),
        strip.text       = element_text(size = 10, face = "bold"),
        panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        legend.position  = "none"
    ) +
    labs(
        title = "Correspondence of DNA Methylation and Chromatin Accessibility",
        subtitle = "Mystery Regions matched to ATAC Peaks (findOverlaps type = 'equal')",
        x = "DNA Methylation Level",
        y = "Normalized ATAC Accessibility Score (VST Count)"
    )

output_file <- file.path(out_dir, "mystery_regions_methylation_vs_accessibility_two_panel.png")
ggsave(output_file, plot = p, width = 9, height = 4.5, dpi = 300)
