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

map_chromhmm_category <- function(state) {
    case_when(
        state == "Pr_A"                       ~ "Active Promoter (Pr_A)",
        state == "Pr_B"                       ~ "Bivalent Promoter (Pr_B)",
        state %in% c("Pr_F", "Pr_W")          ~ "Flanking/Weak Promoter (Pr_F/W)",
        state %in% c("Enh_A", "Enh_P", 
                     "Enh_W", "Enh_ds")       ~ "Enhancer (Enh_*)",
        state %in% c("Tx_I", "Tx_S", "Tx_W")  ~ "Transcription (Tx_*)",
        state %in% c("Het_P", "Het_S")        ~ "Heterochromatin (Het_P/S)",
        state %in% c("Mix", "NS")             ~ "Mixed/Low Signal (Mix/NS)",
        TRUE                                  ~ "Unannotated"
    )
}

cat("=== 2. Overlapping Promoters with ATAC and Mapping ChromHMM States ===\n")
overlaps <- findOverlaps(promoter_gr, atac_gr)

matched_promoters_gr <- promoter_gr[queryHits(overlaps)]

# Map overlapping ChromHMM state names for each tissue (first overlap match)
get_tissue_state <- function(query_gr, target_seg) {
    olaps <- findOverlaps(query_gr, target_seg)
    states <- rep("None", length(query_gr))
    first_hits <- !duplicated(queryHits(olaps))
    states[queryHits(olaps)[first_hits]] <- target_seg$name[subjectHits(olaps)[first_hits]]
    return(states)
}

matched_promoters_gr$kidney_state <- get_tissue_state(matched_promoters_gr, kidney_seg)
matched_promoters_gr$liver_state  <- get_tissue_state(matched_promoters_gr, liver_seg)

# Build unified dataframe
matched_df <- data.frame(
    meth_diff    = matched_promoters_gr$mean.mean.diff,  # Kidney - Liver Meth
    atac_log2fc  = atac_gr$log2FoldChange[subjectHits(overlaps)],
    kidney_state = matched_promoters_gr$kidney_state,
    liver_state  = matched_promoters_gr$liver_state
) %>%
    mutate(
        Regulatory_State = case_when(
            meth_diff < 0 & atac_log2fc > 0 ~ "Kidney Active (Hypo-meth & Open ATAC)",
            meth_diff > 0 & atac_log2fc < 0 ~ "Liver Active (Kidney Hyper & Closed ATAC)",
            meth_diff > 0 & atac_log2fc > 0 ~ "Atypical (Hyper-meth & Open ATAC)",
            meth_diff < 0 & atac_log2fc < 0 ~ "Atypical (Hypo-meth & Closed ATAC)"
        ),
        # Extract tissue-relevant ChromHMM state based on regulatory direction
        Relevant_Tissue_State = case_when(
            Regulatory_State == "Kidney Active (Hypo-meth & Open ATAC)"  ~ kidney_state,
            Regulatory_State == "Liver Active (Kidney Hyper & Closed ATAC)" ~ liver_state,
            TRUE                                                          ~ kidney_state
        ),
        ChromHMM_Category = map_chromhmm_category(Relevant_Tissue_State)
    )

facet_scatter_plot <- ggplot(matched_df, aes(x = meth_diff, y = atac_log2fc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = Regulatory_State), alpha = 0.8, size = 2.2) +
    facet_wrap(~ ChromHMM_Category, ncol = 3, scales = "fixed") +
    scale_color_manual(values = c(
        "Kidney Active (Hypo-meth & Open ATAC)"           = "#D95F02", # Orange
        "Liver Active (Kidney Hyper & Closed ATAC)" = "#7570B3", # Purple
        "Atypical (Hyper-meth & Open ATAC)"              = "#E7298A", # Pink
        "Atypical (Hypo-meth & Closed ATAC)"             = "#66A61E"  # Green
    )) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title       = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "gray92", color = NA),
        strip.text       = element_text(size = 9.5, face = "bold"),
        panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        legend.position  = "bottom",
        legend.title     = element_blank()
    ) +
    guides(color = guide_legend(ncol = 2)) +
    labs(
        title = "Promoter Methylation vs. Accessibility across ChromHMM Categories",
        x = "Promoter Methylation Difference (Kidney - Liver)",
        y = "ATAC-seq log2 Fold Change (Kidney / Liver)"
    )

output_grid_file <- file.path(out_dir, "promoter_methylation_vs_accessibility_chromhmm_grid.png")
ggsave(output_grid_file, plot = facet_scatter_plot, width = 10, height = 8, dpi = 300)

cat(sprintf("\nFaceted grid plot successfully saved to:\n%s\n", output_grid_file))
