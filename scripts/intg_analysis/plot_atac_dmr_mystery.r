#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(GenomicRanges)

wgbs_mystery_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/diff_meth_mystery_regions_annotated.csv"
atac_master_file <- "/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/differential_peaks_annotated.tsv"
out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/atac-wgbs"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

diff_meth <- read.csv(wgbs_mystery_file, stringsAsFactors = FALSE)
atac_master <- read.delim(atac_master_file, stringsAsFactors = FALSE)

cat(sprintf("Loaded Mystery Regions: %d rows\n", nrow(diff_meth)))
cat(sprintf("Loaded ATAC Peaks:     %d rows\n\n", nrow(atac_master)))

# Clean strings and format numeric coordinates
diff_meth$Chromosome <- trimws(gsub('"', "", as.character(diff_meth$Chromosome)))
diff_meth$Start <- as.numeric(diff_meth$Start)
diff_meth$End <- as.numeric(diff_meth$End)

atac_master$chr <- trimws(gsub('"', "", as.character(atac_master$chr)))
atac_master$start <- as.numeric(atac_master$start)
atac_master$end <- as.numeric(atac_master$end)

# Ensure "chr" prefix consistency
diff_meth$Chromosome <- ifelse(startsWith(diff_meth$Chromosome, "chr"), diff_meth$Chromosome, paste0("chr", diff_meth$Chromosome))
atac_master$chr <- ifelse(startsWith(atac_master$chr, "chr"), atac_master$chr, paste0("chr", atac_master$chr))

# Convert to GRanges
dmr_gr <- makeGRangesFromDataFrame(
    diff_meth,
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

# First try strict equal match
overlaps <- findOverlaps(dmr_gr, atac_gr, type = "equal", ignore.strand = TRUE)

# Fallback to maxgap = 1 if 0-based BED vs 1-based GRanges offset exists
if (length(overlaps) == 0) {
    cat("Notice: Strict type = 'equal' returned 0 matches. Retrying with maxgap = 1 (0-based/1-based offset tolerance)...\n")
    overlaps <- findOverlaps(dmr_gr, atac_gr, maxgap = 1, ignore.strand = TRUE)
}

if (length(overlaps) == 0) {
    cat("Notice: Retrying with range overlap (type = 'any')...\n")
    overlaps <- findOverlaps(dmr_gr, atac_gr, type = "any", ignore.strand = TRUE)
}

cat(sprintf("Successfully matched mystery regions: %d\n\n", length(overlaps)))

if (length(overlaps) == 0) {
    stop("Execution halted: No matching coordinates found between WGBS and ATAC tables.")
}

matched_dmr <- diff_meth[queryHits(overlaps), ]
matched_atac <- atac_master[subjectHits(overlaps), ]

kidney_df <- data.frame(
    Region_ID     = paste0(matched_dmr$Chromosome, ":", matched_dmr$Start, "-", matched_dmr$End),
    Tissue        = "Kidney",
    Methylation   = matched_dmr$mean.mean.g1,
    Accessibility = matched_atac$meanVstCountGrp1_kidney
)

liver_df <- data.frame(
    Region_ID     = paste0(matched_dmr$Chromosome, ":", matched_dmr$Start, "-", matched_dmr$End),
    Tissue        = "Liver",
    Methylation   = matched_dmr$mean.mean.g2,
    Accessibility = matched_atac$meanVstCountGrp2_liver
)

plot_df <- rbind(kidney_df, liver_df)

# Format p-value to prevent "p = 0e+00"
stats_df <- plot_df %>%
    group_by(Tissue) %>%
    summarise(
        pearson_r = cor(Methylation, Accessibility, method = "pearson", use = "complete.obs"),
        p_val     = cor.test(Methylation, Accessibility, method = "pearson")$p.value,
        .groups   = "drop"
    ) %>%
    mutate(
        p_formatted = ifelse(p_val < 2.2e-16, "p < 2.2e-16", sprintf("p = %.2e", p_val)),
        label       = sprintf("Pearson r = %.3f\n%s", pearson_r, p_formatted)
    )

cat(" PEARSON CORRELATION RESULTS (Mystery Regions)\n")
print(as.data.frame(stats_df))

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
    facet_wrap(~Tissue, scales = "free_y") +
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
        subtitle = "Mystery Regions matched to ATAC Peaks",
        x = "DNA Methylation Level",
        y = "Normalized ATAC Accessibility Score (VST Count)"
    )

output_file <- file.path(out_dir, "mystery_regions_methylation_vs_accessibility_two_panel.png")
ggsave(output_file, plot = p, width = 9, height = 4.5, dpi = 300)
