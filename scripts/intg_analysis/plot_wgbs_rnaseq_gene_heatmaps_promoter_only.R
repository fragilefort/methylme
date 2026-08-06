suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# WGBS promoter differential-methylation table.
promoter_wgbs_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"

# RNA BED: col4 gene symbol; col5 -log10(FDR); col6 log2FC;
# col7 Kidney CPM; col8 Liver CPM; col9 Ensembl gene ID.
deg_bed_file <- "/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq2/DEGs/DEGs_complete.bed"

out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/rnaseq-wgbs/promoter_gene_heatmaps_clustered_by_rna"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wgbs_mean_diff_cutoff <- 0.20
wgbs_fdr_cutoff <- 0.05
rna_log2fc_cutoff <- 1
rna_fdr_cutoff <- 0.05
rna_minus_log10_fdr_cutoff <- -log10(rna_fdr_cutoff)

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

# WGBS promoter table already contains Ensembl IDs and gene symbols.
# mean.mean.kidney and mean.mean.liver are the promoter methylation values.
wgbs_promoters <- read_csv(promoter_wgbs_file, show_col_types = FALSE) %>%
  transmute(
    ensembl_gene_id = as.character(id),
    gene_symbol = as.character(symbol),
    chr = as.character(Chromosome),
    start = as.integer(Start),
    end = as.integer(End),
    kidney_methylation = as.numeric(mean.mean.kidney),
    liver_methylation = as.numeric(mean.mean.liver),
    mean_difference = as.numeric(mean.mean.diff),
    fdr = as.numeric(comb.p.adj.fdr),
    combined_rank = as.numeric(combinedRank)
  ) %>%
  filter(
    mean_difference > wgbs_mean_diff_cutoff,
    fdr <= wgbs_fdr_cutoff,
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  )

message("Promoter WGBS rows passing filters: ", nrow(wgbs_promoters))

# Median-collapse any duplicate promoter entries per Ensembl gene ID.
wgbs_gene <- wgbs_promoters %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    wgbs_gene_symbol = first_nonmissing(gene_symbol),
    kidney_methylation = median(kidney_methylation, na.rm = TRUE),
    liver_methylation = median(liver_methylation, na.rm = TRUE),
    median_mean_difference = median(mean_difference, na.rm = TRUE),
    median_wgbs_fdr = median(fdr, na.rm = TRUE),
    median_combined_rank = median(combined_rank, na.rm = TRUE),
    n_promoter_rows = n(),
    .groups = "drop"
  )

message("Unique promoter WGBS genes after median collapse: ", nrow(wgbs_gene))

# Read, filter, and median-collapse RNA data per Ensembl gene ID.
rna_gene <- read_tsv(deg_bed_file, col_names = FALSE, show_col_types = FALSE) %>%
  transmute(
    rna_gene_symbol = as.character(X4),
    minus_log10_fdr = as.numeric(X5),
    log2fc = as.numeric(X6),
    kidney_cpm = as.numeric(X7),
    liver_cpm = as.numeric(X8),
    ensembl_gene_id = as.character(X9)
  ) %>%
  filter(
    log2fc > rna_log2fc_cutoff,
    minus_log10_fdr >= rna_minus_log10_fdr_cutoff,
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    rna_gene_symbol = first_nonmissing(rna_gene_symbol),
    kidney_cpm = median(kidney_cpm, na.rm = TRUE),
    liver_cpm = median(liver_cpm, na.rm = TRUE),
    median_log2fc = median(log2fc, na.rm = TRUE),
    median_minus_log10_fdr = median(minus_log10_fdr, na.rm = TRUE),
    n_rna_rows = n(),
    .groups = "drop"
  )

message("RNA-seq genes passing filters: ", nrow(rna_gene))

# Use Ensembl IDs for the join. This is more reliable than joining symbols.
integrated_gene_data <- wgbs_gene %>%
  inner_join(rna_gene, by = "ensembl_gene_id") %>%
  mutate(
    gene_symbol = case_when(
      !is.na(wgbs_gene_symbol) & wgbs_gene_symbol != "" ~ wgbs_gene_symbol,
      !is.na(rna_gene_symbol) & rna_gene_symbol != "" ~ rna_gene_symbol,
      TRUE ~ ensembl_gene_id
    )
  ) %>%
  arrange(desc(median_mean_difference), desc(median_log2fc))

message("Final shared promoter-WGBS/RNA genes: ", nrow(integrated_gene_data))

if (nrow(integrated_gene_data) < 2) {
  stop("Fewer than two shared genes remain after promoter-WGBS and RNA filtering.")
}

write_csv(
  integrated_gene_data,
  file.path(out_dir, "integrated_promoter_WGBS_RNA_gene_table.csv")
)

wgbs_matrix <- integrated_gene_data %>%
  select(kidney_methylation, liver_methylation) %>%
  as.matrix()
colnames(wgbs_matrix) <- c("Kidney", "Liver")
rownames(wgbs_matrix) <- paste0(integrated_gene_data$gene_symbol, " | ", integrated_gene_data$ensembl_gene_id)

rna_matrix <- integrated_gene_data %>%
  select(kidney_cpm, liver_cpm) %>%
  mutate(across(everything(), ~ log2(.x + 1))) %>%
  as.matrix()
colnames(rna_matrix) <- c("Kidney", "Liver")
rownames(rna_matrix) <- rownames(wgbs_matrix)

# RNA-BASED clustering: cluster genes by log2(Kidney CPM + 1) and
# log2(Liver CPM + 1), using Euclidean distance and complete linkage.
rna_row_clustering <- hclust(dist(rna_matrix), method = "complete")
rna_row_order <- rna_row_clustering$order

show_gene_names <- nrow(rna_matrix) <= 100
plot_height_cm <- max(12, min(100, 0.30 * nrow(rna_matrix)))
plot_height_px <- max(1600, min(12000, 35 * nrow(rna_matrix)))

column_annotation <- HeatmapAnnotation(
  Tissue = c("Kidney", "Liver"),
  col = list(Tissue = c(Kidney = "#3B82F6", Liver = "#D97706")),
  annotation_name_side = "left"
)

methylation_colors <- colorRamp2(
  c(0, 0.5, 1),
  c("#2166AC", "white", "#B2182B")
)

rna_max <- max(rna_matrix, na.rm = TRUE)
rna_colors <- colorRamp2(
  c(0, rna_max / 2, rna_max),
  c("white", "#FDB863", "#B2182B")
)

# WGBS heatmap uses RNA-defined order without independently reclustering.
ht_wgbs <- Heatmap(
  wgbs_matrix,
  name = "Promoter methylation",
  col = methylation_colors,
  top_annotation = column_annotation,
  cluster_rows = FALSE,
  row_order = rna_row_order,
  cluster_columns = FALSE,
  show_row_names = show_gene_names,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = "Mean promoter methylation (ordered by RNA clusters)",
  heatmap_legend_param = list(title = "Methylation")
)

# RNA heatmap defines and displays the row clustering dendrogram.
ht_rna <- Heatmap(
  rna_matrix,
  name = "log2(CPM + 1)",
  col = rna_colors,
  top_annotation = column_annotation,
  cluster_rows = as.dendrogram(rna_row_clustering),
  cluster_columns = FALSE,
  show_row_names = show_gene_names,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = "Mean RNA-seq CPM (RNA expression clustering)",
  heatmap_legend_param = list(title = "log2(CPM + 1)")
)

pdf(file.path(out_dir, "promoter_WGBS_methylation_RNA_clustered_heatmap.pdf"),
    width = 8, height = plot_height_cm / 2.54)
draw(ht_wgbs)
dev.off()

pdf(file.path(out_dir, "RNA_CPM_RNA_clustered_heatmap.pdf"),
    width = 8, height = plot_height_cm / 2.54)
draw(ht_rna)
dev.off()

png(file.path(out_dir, "promoter_WGBS_methylation_RNA_clustered_heatmap.png"),
    width = 2400, height = plot_height_px, res = 300)
draw(ht_wgbs)
dev.off()

png(file.path(out_dir, "RNA_CPM_RNA_clustered_heatmap.png"),
    width = 2400, height = plot_height_px, res = 300)
draw(ht_rna)
dev.off()

pdf(file.path(out_dir, "promoter_WGBS_and_RNA_aligned_heatmaps_clustered_by_RNA.pdf"),
    width = 14, height = plot_height_cm / 2.54)
draw(ht_wgbs + ht_rna,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

message("Finished.")
message("Output directory: ", out_dir)
