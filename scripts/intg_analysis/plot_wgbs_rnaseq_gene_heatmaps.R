suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

wgbs_mystery_file <- "/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/diff_meth_mystery_regions_annotated.csv"
gene_gtf_file <- "/vol/COMPEPIWS/pipelines/references/mm10.reduced.refGene.gtf"
deg_bed_file <- "/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq2/DEGs/DEGs_complete.bed"

out_dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/integrative_analysis/rnaseq-wgbs/gene_symbol_heatmaps_from_gtf"
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

if (Sys.which("bedtools") == "") {
  stop("bedtools is not available in PATH. Load or activate an environment containing bedtools.")
}

# WGBS: g1 = Kidney, g2 = Liver.
wgbs <- read_csv(wgbs_mystery_file, show_col_types = FALSE) %>%
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
  filter(mean_difference > wgbs_mean_diff_cutoff, fdr <= wgbs_fdr_cutoff)

message("WGBS DMRs passing filters: ", nrow(wgbs))

# Local refGene GTF: retain transcript records only, then extract NM_ and gene symbols.
gtf_transcripts <- read_tsv(
  gene_gtf_file,
  col_names = FALSE,
  comment = "#",
  show_col_types = FALSE
) %>%
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
  mutate(
    refseq_nm = sub("\\..*$", "", refseq_nm),
    gene_symbol = coalesce(gene_name_gtf, gene_id_gtf)
  ) %>%
  filter(!is.na(refseq_nm), !is.na(gene_symbol), gene_symbol != "") %>%
  select(chr, start_bed, end, gene_symbol, refseq_nm, strand) %>%
  distinct()

message("Transcript records available in local GTF: ", nrow(gtf_transcripts))

# BED coordinates are 0-based half-open, whereas the WGBS CSV starts are 1-based.
dmr_bed <- tempfile(fileext = ".bed")
gtf_bed <- tempfile(fileext = ".bed")
closest_out <- tempfile(fileext = ".tsv")

write_tsv(
  wgbs %>% transmute(chr, start = start - 1L, end, dmr_id),
  dmr_bed,
  col_names = FALSE
)

write_tsv(gtf_transcripts, gtf_bed, col_names = FALSE)

# Assign each filtered DMR to its closest transcript using genomic coordinates.
# -t first makes ties deterministic; -D a also records signed distance.
status <- system2(
  "bedtools",
  args = c("closest", "-a", dmr_bed, "-b", gtf_bed, "-D", "a", "-t", "first"),
  stdout = closest_out
)

if (status != 0) {
  stop("bedtools closest failed.")
}

# Output fields: A has 4 columns; B has 6 columns; final column is distance.
closest_gene <- read_tsv(closest_out, col_names = FALSE, show_col_types = FALSE) %>%
  transmute(
    dmr_id = as.character(X4),
    closest_gene_symbol = as.character(X8),
    closest_refseq_nm = as.character(X9),
    closest_gene_distance_bp = as.numeric(X11)
  ) %>%
  filter(!is.na(closest_gene_symbol), closest_gene_symbol != ".") %>%
  distinct(dmr_id, .keep_all = TRUE)

unlink(c(dmr_bed, gtf_bed, closest_out))

wgbs_mapped <- wgbs %>%
  inner_join(closest_gene, by = "dmr_id")

message("WGBS DMRs assigned to closest GTF genes: ", nrow(wgbs_mapped))
message("Unique WGBS gene symbols: ", n_distinct(wgbs_mapped$closest_gene_symbol))

if (nrow(wgbs_mapped) == 0) {
  stop("No DMRs were assigned to genes by bedtools closest.")
}

# Median-collapse DMRs assigned to the same gene.
wgbs_gene <- wgbs_mapped %>%
  group_by(gene_symbol = closest_gene_symbol) %>%
  summarise(
    kidney_methylation = median(kidney_methylation, na.rm = TRUE),
    liver_methylation = median(liver_methylation, na.rm = TRUE),
    median_mean_difference = median(mean_difference, na.rm = TRUE),
    median_wgbs_fdr = median(fdr, na.rm = TRUE),
    median_distance_to_gene_bp = median(abs(closest_gene_distance_bp), na.rm = TRUE),
    n_wgbs_dmrs = n(),
    .groups = "drop"
  )

# RNA BED columns: 4 gene symbol; 5 -log10(FDR); 6 log2FC;
# 7 Kidney CPM; 8 Liver CPM; 9 Ensembl ID.
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

# Median-collapse duplicate RNA entries per gene symbol.
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

# Only genes passing both WGBS and RNA filters appear in both heatmaps.
integrated_gene_data <- wgbs_gene %>%
  inner_join(rna_gene, by = "gene_symbol") %>%
  arrange(desc(median_mean_difference), desc(median_log2fc))

message("Final shared genes for both heatmaps: ", nrow(integrated_gene_data))

if (nrow(integrated_gene_data) < 2) {
  stop("Fewer than two shared genes remain after WGBS/RNA filtering and gene-symbol matching.")
}

write_csv(integrated_gene_data, file.path(out_dir, "integrated_WGBS_RNA_gene_symbol_table.csv"))

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

methylation_colors <- colorRamp2(c(0, 0.5, 1), c("#2166AC", "white", "#B2182B"))
rna_max <- max(rna_matrix, na.rm = TRUE)
rna_colors <- colorRamp2(c(0, rna_max / 2, rna_max), c("white", "#FDB863", "#B2182B"))

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
  column_title = "Mean methylation of filtered DMR-associated genes",
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

pdf(file.path(out_dir, "WGBS_mean_methylation_heatmap.pdf"), width = 8, height = plot_height_cm / 2.54)
draw(ht_wgbs)
dev.off()

pdf(file.path(out_dir, "RNA_mean_CPM_heatmap.pdf"), width = 8, height = plot_height_cm / 2.54)
draw(ht_rna)
dev.off()

png(file.path(out_dir, "WGBS_mean_methylation_heatmap.png"), width = 2400, height = plot_height_px, res = 300)
draw(ht_wgbs)
dev.off()

png(file.path(out_dir, "RNA_mean_CPM_heatmap.png"), width = 2400, height = plot_height_px, res = 300)
draw(ht_rna)
dev.off()

pdf(file.path(out_dir, "WGBS_and_RNA_aligned_gene_heatmaps.pdf"), width = 14, height = plot_height_cm / 2.54)
draw(ht_wgbs + ht_rna, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

message("Finished.")
message("Output directory: ", out_dir)
