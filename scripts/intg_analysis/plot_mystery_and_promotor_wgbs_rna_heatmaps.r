suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# -------------------------- Input files ---------------------------
base_out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/rnaseq-wgbs"
mystery_wgbs_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/diff_meth_mystery_regions_annotated.csv"
promoter_wgbs_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"
gene_gtf_file <- "/vol/COMPEPIWS/pipelines/references/mm10.reduced.refGene.gtf"
deg_bed_file <- "/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq2/DEGs/DEGs_complete.bed"

mystery_out_dir <- file.path(base_out_dir, "mystery_regions_RNA_clustered_green_black_red")
promoter_out_dir <- file.path(base_out_dir, "promoters_RNA_clustered_green_black_red")
dir.create(mystery_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(promoter_out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- Thresholds ----------------------------
wgbs_mean_diff_cutoff <- 0.20
wgbs_fdr_cutoff <- 0.05
rna_log2fc_cutoff <- 1
rna_fdr_cutoff <- 0.05
rna_minus_log10_fdr_cutoff <- -log10(rna_fdr_cutoff)

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

# ------------------------- Shared RNA data -------------------------
rna_gene <- read_tsv(deg_bed_file, col_names = FALSE, show_col_types = FALSE) %>%
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
    !is.na(ensembl_gene_id), ensembl_gene_id != "",
    !is.na(gene_symbol), gene_symbol != ""
  ) %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    gene_symbol = first_nonmissing(gene_symbol),
    kidney_cpm = median(kidney_cpm, na.rm = TRUE),
    liver_cpm = median(liver_cpm, na.rm = TRUE),
    median_log2fc = median(log2fc, na.rm = TRUE),
    median_minus_log10_fdr = median(minus_log10_fdr, na.rm = TRUE),
    n_rna_rows = n(),
    .groups = "drop"
  )

message("RNA-seq genes passing filters: ", nrow(rna_gene))

# ------------------------ Common plotting --------------------------
plot_heatmap_pair <- function(dat, out_dir, prefix, methylation_title) {
  if (nrow(dat) < 2) stop(prefix, ": fewer than two shared genes remain.")

  dat <- dat %>% arrange(desc(median_mean_difference), desc(median_log2fc))

  gene_label <- paste0(dat$gene_symbol, " | ", dat$ensembl_gene_id)

  wgbs_matrix <- dat %>%
    select(kidney_methylation, liver_methylation) %>%
    as.matrix()
  colnames(wgbs_matrix) <- c("Kidney", "Liver")
  rownames(wgbs_matrix) <- gene_label

  rna_matrix <- dat %>%
    select(kidney_cpm, liver_cpm) %>%
    mutate(across(everything(), ~ log2(.x + 1))) %>%
    as.matrix()
  colnames(rna_matrix) <- c("Kidney", "Liver")
  rownames(rna_matrix) <- gene_label

  # RNA expression defines the clustering and row order for both heatmaps.
  rna_cluster <- hclust(dist(rna_matrix), method = "complete")
  rna_order <- rna_cluster$order

  show_gene_names <- nrow(dat) <= 40
  height_cm <- max(12, min(100, 0.30 * nrow(dat)))
  height_px <- max(1600, min(12000, 35 * nrow(dat)))

  tissue_annotation <- HeatmapAnnotation(
    Tissue = c("Kidney", "Liver"),
    col = list(Tissue = c(Kidney = "#1B9E77", Liver = "#D73027")),
    annotation_name_side = "left"
  )

  # Same low-middle-high palette for both assays: green -> black -> red.
  methylation_colours <- colorRamp2(c(0, 0.5, 1), c("#1B9E77", "black", "#D73027"))
  rna_max <- max(rna_matrix, na.rm = TRUE)
  rna_colours <- colorRamp2(c(0, rna_max / 2, rna_max), c("#1B9E77", "black", "#D73027"))

  ht_wgbs <- Heatmap(
    wgbs_matrix,
    name = "Methylation",
    col = methylation_colours,
    top_annotation = tissue_annotation,
    cluster_rows = FALSE,
    row_order = rna_order,
    cluster_columns = FALSE,
    show_row_names = show_gene_names,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title = paste0(methylation_title, " (ordered by RNA clusters)"),
    heatmap_legend_param = list(title = "Methylation")
  )

  ht_rna <- Heatmap(
    rna_matrix,
    name = "log2(CPM + 1)",
    col = rna_colours,
    top_annotation = tissue_annotation,
    cluster_rows = as.dendrogram(rna_cluster),
    cluster_columns = FALSE,
    show_row_names = show_gene_names,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title = "Mean RNA-seq CPM (RNA expression clustering)",
    heatmap_legend_param = list(title = "log2(CPM + 1)")
  )

  write_csv(dat, file.path(out_dir, paste0(prefix, "_integrated_gene_table.csv")))

  pdf(file.path(out_dir, paste0(prefix, "_methylation_heatmap.pdf")), width = 8, height = height_cm / 2.54)
  draw(ht_wgbs, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  png(file.path(out_dir, paste0(prefix, "_methylation_heatmap.png")), width = 2400, height = height_px, res = 300)
  draw(ht_wgbs, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  pdf(file.path(out_dir, paste0(prefix, "_RNA_CPM_heatmap.pdf")), width = 8, height = height_cm / 2.54)
  draw(ht_rna, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  png(file.path(out_dir, paste0(prefix, "_RNA_CPM_heatmap.png")), width = 2400, height = height_px, res = 300)
  draw(ht_rna, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  pdf(file.path(out_dir, paste0(prefix, "_aligned_heatmaps.pdf")), width = 15, height = height_cm / 2.54)
  draw(ht_wgbs + ht_rna, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  png(file.path(out_dir, paste0(prefix, "_aligned_heatmaps.png")), width = 4200, height = height_px, res = 300)
  draw(ht_wgbs + ht_rna, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  message(prefix, ": finished; shared genes = ", nrow(dat))
}

# ------------------- 1. Mystery-region WGBS data ------------------
if (Sys.which("bedtools") == "") stop("bedtools is not available in PATH.")

mystery_dmrs <- read_csv(mystery_wgbs_file, show_col_types = FALSE) %>%
  transmute(
    dmr_id = paste(Chromosome, Start, End, sep = ":"),
    chr = as.character(Chromosome),
    start = as.integer(Start),
    end = as.integer(End),
    kidney_methylation = as.numeric(mean.mean.g1),
    liver_methylation = as.numeric(mean.mean.g2),
    mean_difference = as.numeric(mean.mean.diff),
    fdr = as.numeric(comb.p.adj.fdr)
  ) %>%
  filter(mean_difference > wgbs_mean_diff_cutoff, fdr <= wgbs_fdr_cutoff) %>%
  arrange(chr, start, end)

message("Mystery-region WGBS DMRs passing filters: ", nrow(mystery_dmrs))

gtf_transcripts <- read_tsv(gene_gtf_file, col_names = FALSE, comment = "#", show_col_types = FALSE) %>%
  filter(X3 == "transcript") %>%
  transmute(
    chr = as.character(X1),
    start_bed = as.integer(X4) - 1L,
    end = as.integer(X5),
    strand = as.character(X7),
    refseq_nm = str_match(X9, 'transcript_id "([^"]+)"')[, 2],
    gene_id_gtf = str_match(X9, 'gene_id "([^"]+)"')[, 2],
    gene_name_gtf = str_match(X9, 'gene_name "([^"]+)"')[, 2]
  ) %>%
  mutate(gene_symbol = coalesce(gene_name_gtf, gene_id_gtf)) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  select(chr, start_bed, end, gene_symbol, refseq_nm, strand) %>%
  distinct() %>%
  arrange(chr, start_bed, end)

dmr_bed <- tempfile(fileext = ".bed")
gtf_bed <- tempfile(fileext = ".bed")
closest_out <- tempfile(fileext = ".tsv")
write_tsv(mystery_dmrs %>% transmute(chr, start = start - 1L, end, dmr_id), dmr_bed, col_names = FALSE)
write_tsv(gtf_transcripts, gtf_bed, col_names = FALSE)

status <- system2("bedtools", c("closest", "-a", dmr_bed, "-b", gtf_bed, "-D", "a", "-t", "first"), stdout = closest_out)
if (status != 0) stop("bedtools closest failed for mystery-region DMRs.")

closest_gene <- read_tsv(closest_out, col_names = FALSE, show_col_types = FALSE) %>%
  transmute(dmr_id = as.character(X4), gene_symbol = as.character(X8)) %>%
  filter(!is.na(gene_symbol), gene_symbol != ".") %>%
  distinct(dmr_id, .keep_all = TRUE)
unlink(c(dmr_bed, gtf_bed, closest_out))

mystery_gene <- mystery_dmrs %>%
  inner_join(closest_gene, by = "dmr_id") %>%
  group_by(gene_symbol) %>%
  summarise(
    kidney_methylation = median(kidney_methylation, na.rm = TRUE),
    liver_methylation = median(liver_methylation, na.rm = TRUE),
    median_mean_difference = median(mean_difference, na.rm = TRUE),
    median_wgbs_fdr = median(fdr, na.rm = TRUE),
    n_wgbs_rows = n(),
    .groups = "drop"
  ) %>%
  inner_join(rna_gene, by = "gene_symbol")

plot_heatmap_pair(mystery_gene, mystery_out_dir, "mystery_regions", "Mean mystery-region methylation")

# ---------------------- 2. Promoter WGBS data ---------------------
promoter_gene <- read_csv(promoter_wgbs_file, show_col_types = FALSE) %>%
  transmute(
    ensembl_gene_id = as.character(id),
    wgbs_gene_symbol = as.character(symbol),
    kidney_methylation = as.numeric(mean.mean.kidney),
    liver_methylation = as.numeric(mean.mean.liver),
    mean_difference = as.numeric(mean.mean.diff),
    fdr = as.numeric(comb.p.adj.fdr),
    combined_rank = as.numeric(combinedRank)
  ) %>%
  filter(
    mean_difference > wgbs_mean_diff_cutoff,
    fdr <= wgbs_fdr_cutoff,
    !is.na(ensembl_gene_id), ensembl_gene_id != ""
  ) %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    wgbs_gene_symbol = first_nonmissing(wgbs_gene_symbol),
    kidney_methylation = median(kidney_methylation, na.rm = TRUE),
    liver_methylation = median(liver_methylation, na.rm = TRUE),
    median_mean_difference = median(mean_difference, na.rm = TRUE),
    median_wgbs_fdr = median(fdr, na.rm = TRUE),
    median_combined_rank = median(combined_rank, na.rm = TRUE),
    n_wgbs_rows = n(),
    .groups = "drop"
  ) %>%
  inner_join(rna_gene, by = "ensembl_gene_id") %>%
  mutate(gene_symbol = coalesce(wgbs_gene_symbol, gene_symbol)) %>%
  select(-wgbs_gene_symbol)

plot_heatmap_pair(promoter_gene, promoter_out_dir, "promoters", "Mean promoter methylation")

message("All mystery-region and promoter heatmaps finished.")
message("Mystery results: ", mystery_out_dir)
message("Promoter results: ", promoter_out_dir)
