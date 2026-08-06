suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Input files
wgbs_mystery_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/diff_meth_mystery_regions_annotated.csv"
closest_gene_file <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/rnaseq-wgbs/mystery_regions_closest_annotated_genes.tsv"
gene_gtf_file <- "/vol/COMPEPIWS/pipelines/references/mm10.reduced.refGene.gtf"
deg_bed_file <- "/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq2/DEGs/DEGs_complete.bed"

out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/rnaseq-wgbs/gene_symbol_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Thresholds
wgbs_mean_diff_cutoff <- 0.20
wgbs_fdr_cutoff <- 0.05
rna_log2fc_cutoff <- 1
rna_fdr_cutoff <- 0.05
rna_minus_log10_fdr_cutoff <- -log10(rna_fdr_cutoff)

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

# WGBS: mean.mean.g1 = Kidney; mean.mean.g2 = Liver
wgbs <- read_csv(wgbs_mystery_file, show_col_types = FALSE) %>%
  transmute(
    chr = as.character(Chromosome),
    start = as.integer(Start),
    end = as.integer(End),
    kidney_methylation = as.numeric(mean.mean.g1),
    liver_methylation = as.numeric(mean.mean.g2),
    mean_difference = as.numeric(mean.mean.diff),
    fdr = as.numeric(comb.p.adj.fdr)
  ) %>%
  filter(
    mean_difference > wgbs_mean_diff_cutoff,
    fdr <= wgbs_fdr_cutoff
  )

message("WGBS DMRs passing filters: ", nrow(wgbs))

# Headerless whitespace-delimited bedtools output:
# columns 1-3 = DMR coordinates; column 9 = closest NM_ transcript.
closest_raw <- read.table(
  closest_gene_file,
  header = FALSE,
  sep = "",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

if (ncol(closest_raw) < 9) {
  stop("Closest-gene file has fewer than 9 fields; expected RefSeq NM_ transcript in field 9.")
}

closest_transcript <- closest_raw %>%
  transmute(
    chr = trimws(as.character(V1)),
    start = as.integer(V2),
    end = as.integer(V3),
    refseq_nm = sub("\\..*$", "", trimws(as.character(V9)))
  ) %>%
  filter(str_detect(refseq_nm, "^NM_")) %>%
  distinct(chr, start, end, refseq_nm)

# Exact coordinates were confirmed to match between these two WGBS-derived files.
wgbs_with_nm <- wgbs %>%
  inner_join(closest_transcript, by = c("chr", "start", "end")) %>%
  distinct()

message("Filtered WGBS DMR-transcript rows after coordinate join: ", nrow(wgbs_with_nm))
message("Unique closest NM_ transcripts: ", n_distinct(wgbs_with_nm$refseq_nm))

if (nrow(wgbs_with_nm) == 0) {
  stop("The coordinate join returned zero rows. Check parsing of closest_gene_file.")
}

# Local GTF mapping: transcript_id NM_ -> gene_id/gene_name.
# In this refGene GTF, gene_id and gene_name are mouse gene symbols.
gtf_transcripts <- read_tsv(
  gene_gtf_file,
  col_names = FALSE,
  comment = "#",
  show_col_types = FALSE
) %>%
  filter(X3 == "transcript") %>%
  transmute(
    refseq_nm = str_match(X9, 'transcript_id "([^"]+)"')[, 2],
    gene_id_gtf = str_match(X9, 'gene_id "([^"]+)"')[, 2],
    gene_name_gtf = str_match(X9, 'gene_name "([^"]+)"')[, 2]
  ) %>%
  mutate(
    refseq_nm = sub("\\..*$", "", refseq_nm),
    gene_symbol = coalesce(gene_name_gtf, gene_id_gtf)
  ) %>%
  filter(!is.na(refseq_nm), !is.na(gene_symbol), gene_symbol != "") %>%
  select(refseq_nm, gene_symbol) %>%
  distinct()

wgbs_mapped <- wgbs_with_nm %>%
  inner_join(gtf_transcripts, by = "refseq_nm")

message("WGBS rows mapped from NM_ to gene symbols: ", nrow(wgbs_mapped))
message("Unique mapped WGBS gene symbols: ", n_distinct(wgbs_mapped$gene_symbol))

if (nrow(wgbs_mapped) == 0) {
  stop("No NM_ transcripts mapped in the local GTF. Check gene_gtf_file and transcript IDs.")
}

# Median-collapse all DMR/transcript entries belonging to the same gene symbol.
wgbs_gene <- wgbs_mapped %>%
  group_by(gene_symbol) %>%
  summarise(
    kidney_methylation = median(kidney_methylation, na.rm = TRUE),
    liver_methylation = median(liver_methylation, na.rm = TRUE),
    median_mean_difference = median(mean_difference, na.rm = TRUE),
    median_wgbs_fdr = median(fdr, na.rm = TRUE),
    n_wgbs_rows = n(),
    .groups = "drop"
  )

# RNA BED has no header:
# col 4 = gene symbol; col 5 = -log10(FDR); col 6 = log2FC;
# col 7 = Kidney mean CPM; col 8 = Liver mean CPM; col 9 = Ensembl ID.
# The Ensembl ID is retained in the output table but the merge uses gene symbols.
rna <- read_tsv(deg_bed_file, col_names = FALSE, show_col_types = FALSE) %>%
  transmute(
    gene_symbol = as.character(X4),
    minus_log10_fdr = as.numeric(X5),
    log2fc = as.numeric(X6),
    kidney_cpm = as.numeric(X7),
    liver_cpm = as.numeric(X8),
    ensembl_gene_id = as.character(X9)
  ) %>%
  filter(
    log2fc > rna_log2fc_cutoff,
    minus_log10_fdr >= rna_minus_log10_fdr_cutoff,
    !is.na(gene_symbol),
    gene_symbol != ""
  )

message("RNA-seq rows passing filters: ", nrow(rna))

# Median-collapse duplicate RNA rows per gene symbol.
rna_gene <- rna %>%
  group_by(gene_symbol) %>%
  summarise(
    ensembl_gene_id = first_nonmissing(ensembl_gene_id),
    kidney_cpm = median(kidney_cpm, na.rm = TRUE),
    liver_cpm = median(liver_cpm, na.rm = TRUE),
    median_log2fc = median(log2fc, na.rm = TRUE),
    median_minus_log10_fdr = median(minus_log10_fdr, na.rm = TRUE),
    n_rna_rows = n(),
    .groups = "drop"
  )

# Keep genes that pass both WGBS and RNA filters.
integrated_gene_data <- wgbs_gene %>%
  inner_join(rna_gene, by = "gene_symbol") %>%
  arrange(desc(median_mean_difference), desc(median_log2fc))

message("Final shared genes for both heatmaps: ", nrow(integrated_gene_data))

if (nrow(integrated_gene_data) < 2) {
  stop("Fewer than two shared genes remain after filtering and gene-symbol matching.")
}

write_csv(
  integrated_gene_data,
  file.path(out_dir, "integrated_WGBS_RNA_gene_symbol_table.csv")
)

# Matrices share exactly the same gene-symbol rows.
wgbs_matrix <- integrated_gene_data %>%
  select(kidney_methylation, liver_methylation) %>%
  as.matrix()
colnames(wgbs_matrix) <- c("Kidney", "Liver")
rownames(wgbs_matrix) <- integrated_gene_data$gene_symbol

rna_matrix <- integrated_gene_data %>%
  select(kidney_cpm, liver_cpm) %>%
  mutate(across(everything(), ~ log2(.x + 1))) %>%
  as.matrix()
colnames(rna_matrix) <- c("Kidney", "Liver")
rownames(rna_matrix) <- integrated_gene_data$gene_symbol

# Apply the same row order to both plots.
row_clustering <- hclust(dist(wgbs_matrix), method = "complete")
row_order <- row_clustering$order
show_gene_names <- nrow(wgbs_matrix) <= 100
plot_height_cm <- max(12, min(100, 0.30 * nrow(wgbs_matrix)))
plot_height_px <- max(1600, min(12000, 35 * nrow(wgbs_matrix)))

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
  column_title = "Mean methylation of filtered WGBS DMR-associated genes",
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
  column_title = "Mean RNA-seq CPM of the same filtered genes",
  heatmap_legend_param = list(title = "log2(CPM + 1)")
)

pdf(file.path(out_dir, "WGBS_mean_methylation_gene_symbol_heatmap.pdf"),
    width = 8, height = plot_height_cm / 2.54)
draw(ht_wgbs)
dev.off()

pdf(file.path(out_dir, "RNA_mean_CPM_gene_symbol_heatmap.pdf"),
    width = 8, height = plot_height_cm / 2.54)
draw(ht_rna)
dev.off()

png(file.path(out_dir, "WGBS_mean_methylation_gene_symbol_heatmap.png"),
    width = 2400, height = plot_height_px, res = 300)
draw(ht_wgbs)
dev.off()

png(file.path(out_dir, "RNA_mean_CPM_gene_symbol_heatmap.png"),
    width = 2400, height = plot_height_px, res = 300)
draw(ht_rna)
dev.off()

pdf(file.path(out_dir, "WGBS_and_RNA_aligned_gene_symbol_heatmaps.pdf"),
    width = 14, height = plot_height_cm / 2.54)
draw(ht_wgbs + ht_rna,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

message("Finished.")
message("Output directory: ", out_dir)
