suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------

wgbs_mystery_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/diff_meth_mystery_regions_annotated.csv"

closest_gene_file <- paste0(
  "/vol/COMPEPIWS/groups/wgbs2/methylme/",
  "integrative_analysis/rnaseq-wgbs/",
  "mystery_regions_closest_annotated_genes.tsv"
)

deg_bed_file <- "/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq2/DEGs/DEGs_complete.bed"

out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/rnaseq-wgbs/gene_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Thresholds
# ------------------------------------------------------------------

wgbs_mean_diff_cutoff <- 0.20
wgbs_fdr_cutoff <- 0.05

rna_log2fc_cutoff <- 1
rna_fdr_cutoff <- 0.05
rna_minus_log10_fdr_cutoff <- -log10(rna_fdr_cutoff)

# ------------------------------------------------------------------
# 1. Read and filter WGBS DMRs
# g1 = Kidney
# g2 = Liver
# ------------------------------------------------------------------

wgbs <- read_csv(wgbs_mystery_file, show_col_types = FALSE) %>%
  transmute(
    chr = Chromosome,
    start = as.integer(Start),
    end = as.integer(End),
    kidney_methylation = mean.mean.g1,
    liver_methylation = mean.mean.g2,
    mean_difference = mean.mean.diff,
    fdr = comb.p.adj.fdr
  ) %>%
  filter(
    mean_difference > wgbs_mean_diff_cutoff,
    fdr <= wgbs_fdr_cutoff
  )

message("WGBS DMRs passing filters: ", nrow(wgbs))

# ------------------------------------------------------------------
# 2. Read closest-transcript file
# Column 1-3: DMR coordinates
# Column 9: RefSeq NM_ transcript ID
# ------------------------------------------------------------------

closest_transcript <- read_tsv(
  closest_gene_file,
  col_names = FALSE,
  show_col_types = FALSE
) %>%
  transmute(
    chr = X1,
    start = as.integer(X2),
    end = as.integer(X3),
    refseq_nm = sub("\\..*$", "", X9)
  ) %>%
  filter(str_detect(refseq_nm, "^NM_")) %>%
  distinct(chr, start, end, refseq_nm)

# Coordinate join: filtered DMRs + closest NM_ transcript
wgbs_with_nm <- wgbs %>%
  inner_join(closest_transcript, by = c("chr", "start", "end"))

message("Filtered WGBS DMR-transcript rows after coordinate join: ", nrow(wgbs_with_nm))
message("Unique NM_ transcripts after coordinate join: ", n_distinct(wgbs_with_nm$refseq_nm))

# ------------------------------------------------------------------
# 3. Map RefSeq NM_ transcript IDs to Ensembl mouse gene IDs
# ------------------------------------------------------------------

refseq_to_ensembl <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(wgbs_with_nm$refseq_nm),
  keytype = "REFSEQ",
  columns = c("ENSEMBL", "SYMBOL")
) %>%
  transmute(
    refseq_nm = REFSEQ,
    ensembl_gene_id = ENSEMBL,
    gene_symbol = SYMBOL
  ) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  distinct()

wgbs_mapped <- wgbs_with_nm %>%
  inner_join(refseq_to_ensembl, by = "refseq_nm")

message("WGBS rows mapped to Ensembl gene IDs: ", nrow(wgbs_mapped))
message("Unique mapped WGBS Ensembl genes: ", n_distinct(wgbs_mapped$ensembl_gene_id))

# Median-collapse all DMRs/transcripts mapping to the same Ensembl gene.
wgbs_gene <- wgbs_mapped %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    gene_symbol = first(na.omit(gene_symbol), default = NA_character_),
    kidney_methylation = median(kidney_methylation, na.rm = TRUE),
    liver_methylation = median(liver_methylation, na.rm = TRUE),
    median_mean_difference = median(mean_difference, na.rm = TRUE),
    n_wgbs_rows = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------------
# 4. Read and filter RNA-seq DEG BED
#
# Input columns:
# 1 chr
# 2 start
# 3 end
# 4 gene symbol
# 5 -log10(FDR)
# 6 log2FC
# 7 Kidney mean CPM
# 8 Liver mean CPM
# 9 Ensembl gene ID
# ------------------------------------------------------------------

rna <- read_tsv(
  deg_bed_file,
  col_names = FALSE,
  show_col_types = FALSE
) %>%
  transmute(
    chr = X1,
    start = as.integer(X2),
    end = as.integer(X3),
    rna_gene_symbol = X4,
    minus_log10_fdr = as.numeric(X5),
    log2fc = as.numeric(X6),
    kidney_cpm = as.numeric(X7),
    liver_cpm = as.numeric(X8),
    ensembl_gene_id = X9
  ) %>%
  filter(
    log2fc > rna_log2fc_cutoff,
    minus_log10_fdr >= rna_minus_log10_fdr_cutoff,
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  )

message("RNA-seq genes passing filters: ", nrow(rna))

# Median-collapse duplicate Ensembl IDs.
rna_gene <- rna %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    rna_gene_symbol = first(na.omit(rna_gene_symbol), default = NA_character_),
    kidney_cpm = median(kidney_cpm, na.rm = TRUE),
    liver_cpm = median(liver_cpm, na.rm = TRUE),
    median_log2fc = median(log2fc, na.rm = TRUE),
    n_rna_rows = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------------
# 5. Keep only genes passing both WGBS and RNA filters
# ------------------------------------------------------------------

integrated_gene_data <- wgbs_gene %>%
  inner_join(rna_gene, by = "ensembl_gene_id") %>%
  mutate(
    gene_label = case_when(
      !is.na(gene_symbol) & gene_symbol != "" ~ paste0(gene_symbol, " | ", ensembl_gene_id),
      !is.na(rna_gene_symbol) & rna_gene_symbol != "" ~ paste0(rna_gene_symbol, " | ", ensembl_gene_id),
      TRUE ~ ensembl_gene_id
    )
  ) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  arrange(desc(median_mean_difference), desc(median_log2fc))

message("Final shared genes for both heatmaps: ", nrow(integrated_gene_data))

if (nrow(integrated_gene_data) < 2) {
  stop(
    "Fewer than two genes remain after all filters. ",
    "Check thresholds, coordinate matching, or NM_-to-Ensembl mapping."
  )
}

write_csv(
  integrated_gene_data,
  file.path(out_dir, "integrated_WGBS_RNA_gene_table.csv")
)

# ------------------------------------------------------------------
# 6. Create aligned matrices
# WGBS heatmap uses methylation fractions, range 0 to 1.
# RNA heatmap uses log2(CPM + 1).
# ------------------------------------------------------------------

wgbs_matrix <- integrated_gene_data %>%
  select(kidney_methylation, liver_methylation) %>%
  as.matrix()

colnames(wgbs_matrix) <- c("Kidney", "Liver")
rownames(wgbs_matrix) <- integrated_gene_data$gene_label

rna_matrix <- integrated_gene_data %>%
  select(kidney_cpm, liver_cpm) %>%
  mutate(across(everything(), ~ log2(.x + 1))) %>%
  as.matrix()

colnames(rna_matrix) <- c("Kidney", "Liver")
rownames(rna_matrix) <- integrated_gene_data$gene_label

# Cluster the rows once using WGBS methylation.
row_clustering <- hclust(dist(wgbs_matrix), method = "complete")
row_order <- row_clustering$order

show_gene_names <- nrow(wgbs_matrix) <= 100
plot_height_cm <- max(12, min(100, 0.30 * nrow(wgbs_matrix)))
plot_height_px <- max(1600, min(12000, 35 * nrow(wgbs_matrix)))

# ------------------------------------------------------------------
# 7. Heatmap settings
# ------------------------------------------------------------------

methylation_colors <- colorRamp2(
  c(0, 0.5, 1),
  c("#2166AC", "white", "#B2182B")
)

rna_max <- max(rna_matrix, na.rm = TRUE)

rna_colors <- colorRamp2(
  c(0, rna_max / 2, rna_max),
  c("white", "#FDB863", "#B2182B")
)

column_annotation <- HeatmapAnnotation(
  Tissue = c("Kidney", "Liver"),
  col = list(Tissue = c(Kidney = "#3B82F6", Liver = "#D97706")),
  annotation_name_side = "left"
)

ht_wgbs <- Heatmap(
  wgbs_matrix,
  name = "Mean methylation",
  col = methylation_colors,
  top_annotation = column_annotation,
  cluster_rows = as.dendrogram(row_clustering),
  cluster_columns = FALSE,
  show_row_names = show_gene_names,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = "Mean methylation of WGBS DMR-associated genes",
  heatmap_legend_param = list(title = "Methylation")
)

ht_rna <- Heatmap(
  rna_matrix,
  name = "log2(CPM + 1)",
  col = rna_colors,
  top_annotation = column_annotation,
  cluster_rows = FALSE,
  row_order = row_order,
  cluster_columns = FALSE,
  show_row_names = show_gene_names,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = "Mean RNA-seq CPM of the same genes",
  heatmap_legend_param = list(title = "log2(CPM + 1)")
)

# ------------------------------------------------------------------
# 8. Save separate PDFs
# ------------------------------------------------------------------

pdf(
  file.path(out_dir, "WGBS_mean_methylation_heatmap.pdf"),
  width = 8,
  height = plot_height_cm / 2.54
)
draw(ht_wgbs)
dev.off()

pdf(
  file.path(out_dir, "RNA_mean_CPM_heatmap.pdf"),
  width = 8,
  height = plot_height_cm / 2.54
)
draw(ht_rna)
dev.off()

# ------------------------------------------------------------------
# 9. Save separate PNGs
# ------------------------------------------------------------------

png(
  file.path(out_dir, "WGBS_mean_methylation_heatmap.png"),
  width = 2400,
  height = plot_height_px,
  res = 300
)
draw(ht_wgbs)
dev.off()

png(
  file.path(out_dir, "RNA_mean_CPM_heatmap.png"),
  width = 2400,
  height = plot_height_px,
  res = 300
)
draw(ht_rna)
dev.off()

# ------------------------------------------------------------------
# 10. Save both aligned heatmaps in one PDF
# ------------------------------------------------------------------

pdf(
  file.path(out_dir, "WGBS_and_RNA_aligned_gene_heatmaps.pdf"),
  width = 14,
  height = plot_height_cm / 2.54
)
draw(
  ht_wgbs + ht_rna,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()

message("Finished.")
message("Output directory: ", out_dir)
