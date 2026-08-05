#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(GenomicRanges)

wgbs_mystery_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/diff_meth_mystery_regions_annotated.csv"
gene_ref_file <- "/vol/COMPEPIWS/pipelines/references/mm10_reduced_chr18_chr19_genes.bed"
deg_bed_file <- "/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq2/DEGs/DEGs_complete.bed"
out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/rna-wgbs"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

diff_meth <- read.csv(wgbs_mystery_file, stringsAsFactors = FALSE)
gene_ref <- read.delim(gene_ref_file, header = FALSE, stringsAsFactors = FALSE)
deg_df <- read.delim(deg_bed_file,
    header = FALSE, stringsAsFactors = FALSE,
    col.names = c("chr", "start", "end", "gene_name", "baseMean", "log2FC_rna", "expr_kidney", "expr_liver", "gene_id")
)

diff_meth$Chromosome <- trimws(gsub('"', "", as.character(diff_meth$Chromosome)))
diff_meth$Start <- as.numeric(diff_meth$Start)
diff_meth$End <- as.numeric(diff_meth$End)
diff_meth$Chromosome <- ifelse(startsWith(diff_meth$Chromosome, "chr"), diff_meth$Chromosome, paste0("chr", diff_meth$Chromosome))

colnames(gene_ref)[1:4] <- c("chr", "start", "end", "gene_id")
gene_ref$chr <- trimws(gsub('"', "", as.character(gene_ref$chr)))
gene_ref$start <- as.numeric(gene_ref$start)
gene_ref$end <- as.numeric(gene_ref$end)
gene_ref$chr <- ifelse(startsWith(gene_ref$chr, "chr"), gene_ref$chr, paste0("chr", gene_ref$chr))

cat(sprintf("Loaded Mystery Regions: %d rows\n", nrow(diff_meth)))
cat(sprintf("Loaded Reference Genes: %d rows\n", nrow(gene_ref)))
cat(sprintf("Loaded RNA-seq DEGs:    %d rows\n\n", nrow(deg_df)))

# Convert to GRanges
dmr_gr <- makeGRangesFromDataFrame(
    diff_meth,
    seqnames.field = "Chromosome",
    start.field = "Start",
    end.field = "End",
    keep.extra.columns = TRUE
)

gene_ref_gr <- makeGRangesFromDataFrame(
    gene_ref,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
)

nearest_hits <- distanceToNearest(dmr_gr, gene_ref_gr, ignore.strand = TRUE)

matched_dmr <- diff_meth[queryHits(nearest_hits), ]
matched_gene <- gene_ref[subjectHits(nearest_hits), ]
dist_val <- mcols(nearest_hits)$distance

annotated_df <- matched_dmr %>%
    mutate(
        DMR_ID            = paste0(Chromosome, ":", Start, "-", End),
        Annotated_Gene_ID = matched_gene$gene_id,
        Distance_to_Gene  = dist_val
    ) %>%
    left_join(deg_df, by = c("Annotated_Gene_ID" = "gene_id"))

# Save annotated table
write.csv(annotated_df, file.path(out_dir, "mystery_regions_annotated_reference_genes.csv"), row.names = FALSE)

total_dmrs <- length(unique(annotated_df$DMR_ID))
unique_genes <- length(unique(annotated_df$Annotated_Gene_ID))

cat(" DMR TO GENE MAPPING STATISTICS\n")
cat(sprintf("Total Unique DMRs analyzed : %d\n", total_dmrs))
cat(sprintf("Total Unique Genes mapped  : %d\n\n", unique_genes))

dmrs_per_gene <- annotated_df %>%
    group_by(Annotated_Gene_ID) %>%
    summarise(num_dmrs = n_distinct(DMR_ID), .groups = "drop")

freq_table <- dmrs_per_gene %>%
    count(num_dmrs, name = "gene_count") %>%
    mutate(percentage = (gene_count / sum(gene_count)) * 100)

print("Frequency of DMRs per Gene:")
print(as.data.frame(freq_table))

p <- ggplot(dmrs_per_gene, aes(x = factor(num_dmrs))) +
    geom_bar(fill = "#2B8CBE", color = "black", alpha = 0.8, width = 0.7) +
    geom_text(
        stat = "count",
        aes(label = after_stat(count)),
        vjust = -0.5,
        size = 3.8
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title       = element_text(size = 11, face = "bold"),
        panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        axis.text        = element_text(color = "black")
    ) +
    labs(
        title = "Distribution of Differentially Methylated Regions (DMRs) per Gene",
        subtitle = "Frequency of reference genes associated with single or multiple mystery regions",
        x = "Number of DMRs Assigned to a Gene",
        y = "Number of Genes"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

output_plot <- file.path(out_dir, "dmrs_per_gene_barplot.png")
ggsave(output_plot, plot = p, width = 8, height = 5, dpi = 300)
