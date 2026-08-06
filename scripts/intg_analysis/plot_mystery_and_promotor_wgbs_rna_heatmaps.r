suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# -------------------------- Input files ---------------------------
out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/rnaseq-wgbs"
mystery_wgbs_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/diff_meth_mystery_regions_annotated.csv"
promoter_wgbs_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv"
gene_gtf_file <- "/vol/COMPEPIWS/pipelines/references/mm10.reduced.refGene.gtf"
deg_bed_file <- "/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq2/DEGs/DEGs_complete.bed"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

# Ensembl IDs are sometimes stored with a version suffix (e.g. ".3") in one
# file but not the other. Strip versions before joining so genes aren't
# silently dropped by a mismatched ID format.
strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)
rna_gene <- rna_gene %>% mutate(ensembl_gene_id = strip_ensembl_version(ensembl_gene_id))

# ---------------------- Publication colour scheme -------------------
# Colour-blind-safe diverging palette (Blue - White - Red), smoothly
# interpolated so gradients look continuous rather than banded.
methylation_colours <- colorRamp2(
  seq(0, 1, length.out = 7),
  c("#053061", "#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B", "#67001F")
)

tissue_palette <- c(Kidney = "#1B9E77", Liver = "#D95F02")

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
  height_px <- max(1800, min(11000, 32 * nrow(dat)))

  tissue_annotation <- HeatmapAnnotation(
    Tissue = c("Kidney", "Liver"),
    col = list(Tissue = tissue_palette),
    annotation_name_side = "left",
    simple_anno_size = unit(0.3, "cm"),
    annotation_legend_param = list(Tissue = list(title = "Tissue"))
  )

  rna_max <- max(rna_matrix, na.rm = TRUE)
  rna_colours <- colorRamp2(
    seq(0, rna_max, length.out = 7),
    c("#053061", "#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B", "#67001F")
  )

  font_family <- "sans"

  ht_wgbs <- Heatmap(
    wgbs_matrix,
    name = "Methylation",
    col = methylation_colours,
    top_annotation = tissue_annotation,
    cluster_rows = FALSE,
    row_order = rna_order,
    cluster_columns = FALSE,
    show_row_names = show_gene_names,
    row_names_gp = gpar(fontsize = 7, fontfamily = font_family),
    column_names_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = font_family),
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = methylation_title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = font_family),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    heatmap_legend_param = list(
      title = "Mean\nmethylation",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      legend_height = unit(3, "cm")
    )
  )

  ht_rna <- Heatmap(
    rna_matrix,
    name = "log2(CPM + 1)",
    col = rna_colours,
    top_annotation = tissue_annotation,
    cluster_rows = as.dendrogram(rna_cluster),
    cluster_columns = FALSE,
    show_row_names = show_gene_names,
    row_names_gp = gpar(fontsize = 7, fontfamily = font_family),
    column_names_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = font_family),
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = "RNA-seq expression",
    column_title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = font_family),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    heatmap_legend_param = list(
      title = "log2(CPM+1)",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      legend_height = unit(3, "cm")
    )
  )

  write_csv(dat, file.path(out_dir, paste0(prefix, "_integrated_gene_table.csv")))

  png(
    file.path(out_dir, paste0(prefix, "_heatmap.png")),
    width = 4200, height = height_px, res = 300, bg = "white"
  )
  draw(
    ht_wgbs + ht_rna,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    ht_gap = unit(6, "mm"),
    merge_legend = TRUE,
    padding = unit(c(4, 4, 4, 4), "mm")
  )
  dev.off()

  message(prefix, ": finished; shared genes = ", nrow(dat))
}

# Run one analysis end-to-end without letting a failure here block the
# other analysis. Prints a clear message on success or failure.
run_analysis <- function(label, expr) {
  message("\n=== Starting analysis: ", label, " ===")
  result <- tryCatch(
    { force(expr); message("=== ", label, " completed successfully ===") ; TRUE },
    error = function(e) {
      message("!!! ", label, " FAILED: ", conditionMessage(e))
      FALSE
    }
  )
  invisible(result)
}

# ------------------- 1. Mystery-region WGBS data ------------------
run_analysis("mystery regions", {

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

message("Mystery-region genes after joining to RNA-seq: ", nrow(mystery_gene))

plot_heatmap_pair(mystery_gene, out_dir, "mystery_regions", "Mean mystery-region methylation")

})

# ---------------------- 2. Promoter WGBS data ---------------------
run_analysis("promoters", {

promoter_wgbs_raw <- read_csv(promoter_wgbs_file, show_col_types = FALSE) %>%
  transmute(
    ensembl_gene_id = strip_ensembl_version(as.character(id)),
    wgbs_gene_symbol = as.character(symbol),
    kidney_methylation = as.numeric(mean.mean.kidney),
    liver_methylation = as.numeric(mean.mean.liver),
    mean_difference = as.numeric(mean.mean.diff),
    fdr = as.numeric(comb.p.adj.fdr),
    combined_rank = as.numeric(combinedRank)
  )
message("Promoter WGBS rows read: ", nrow(promoter_wgbs_raw))

promoter_wgbs_filtered <- promoter_wgbs_raw %>%
  filter(
    mean_difference > wgbs_mean_diff_cutoff,
    fdr <= wgbs_fdr_cutoff,
    !is.na(ensembl_gene_id), ensembl_gene_id != ""
  )
message("Promoter WGBS rows passing mean-diff/FDR filters: ", nrow(promoter_wgbs_filtered))
if (nrow(promoter_wgbs_filtered) == 0) {
  stop("No promoter regions passed the WGBS filters (mean_difference > ",
       wgbs_mean_diff_cutoff, ", fdr <= ", wgbs_fdr_cutoff, "). Check the ",
       "input file and thresholds.")
}

promoter_gene <- promoter_wgbs_filtered %>%
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
  )
message("Unique promoter genes after collapsing duplicate rows: ", nrow(promoter_gene))

promoter_gene <- promoter_gene %>%
  inner_join(rna_gene, by = "ensembl_gene_id") %>%
  mutate(gene_symbol = coalesce(wgbs_gene_symbol, gene_symbol)) %>%
  select(-wgbs_gene_symbol)
message("Promoter genes after joining to RNA-seq: ", nrow(promoter_gene))

if (nrow(promoter_gene) < 2) {
  stop("Fewer than two promoter genes overlap with the RNA-seq DEG list after ",
       "joining on ensembl_gene_id. This is usually an ID-format mismatch ",
       "(e.g. gene symbol vs Ensembl ID, or stale annotation build) rather ",
       "than a true lack of overlap -- inspect promoter_wgbs_raw$ensembl_gene_id ",
       "and rna_gene$ensembl_gene_id directly.")
}

plot_heatmap_pair(promoter_gene, out_dir, "promoters", "Mean promoter methylation")

})

message("\nAll requested analyses attempted.")
message("Results written to: ", out_dir)
