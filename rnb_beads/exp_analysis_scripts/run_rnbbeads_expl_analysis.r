###############################################################################
# Manual exploratory DNA methylation analysis after RnBeads preprocessing
# Dataset: mouse kidney and liver at E14.5 and E15.5
#
# This script DOES NOT rerun rnb.run.analysis().
# It loads the existing rnbSet_preprocessed object made by the RnBeads pipeline.
#
# IMPORTANT:
# - Density plots and boxplots use ALL complete CpG sites: no random sampling.
# - PCA and heatmaps intentionally select the most variable sites/regions.
# - Run this on a worker node inside a screen session because full-data plotting
#   can require substantial RAM.
###############################################################################

###############################################################################
# 0. Load packages
###############################################################################

library(RnBeads)
library(RnBeads.mm10)

library(ggplot2)
library(pheatmap)
library(matrixStats)

###############################################################################
# 1. Define paths
###############################################################################

report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"

plot.dir <- file.path(report.dir, "manual_plots")

dir.create(
    plot.dir,
    recursive = TRUE,
    showWarnings = FALSE
)

rnb.set.file <- file.path(
    report.dir,
    "rnbSet_preprocessed"
)

###############################################################################
# 2. RnBeads session settings
###############################################################################

rnb.options(
    assembly = "mm10",
    identifiers.column = "sampleId"
)

###############################################################################
# 3. Load the preprocessed RnBiseqSet
###############################################################################

rnb.set <- load.rnb.set(rnb.set.file)

# Sample-level metadata: one row per sample
pheno.tab <- pheno(rnb.set)

# Site-level methylation beta matrix:
# rows    = CpG sites
# columns = samples
# values  = methylation beta values between 0 and 1
beta <- meth(
    rnb.set,
    type = "sites",
    row.names = TRUE
)

print(pheno.tab)
print(colnames(pheno.tab))
print(dim(beta))
print(summarized.regions(rnb.set))

###############################################################################
# 4. Match sample annotation to methylation-matrix columns
###############################################################################

cat("\nMethylation matrix column names:\n")
print(colnames(beta))

cat("\nPhenotype-table row names:\n")
print(rownames(pheno.tab))

cat("\nSample IDs from sample sheet:\n")
print(pheno.tab$sampleId)

if (nrow(pheno.tab) != ncol(beta)) {
    stop(
        "Number of phenotype rows (", nrow(pheno.tab),
        ") does not equal number of methylation samples (", ncol(beta),
        "). Do not continue."
    )
}

# Preserve original RnBeads IDs as a metadata column.
pheno.tab$rnbeads.internal.id <- rownames(pheno.tab)

# Match metadata rows to methylation-matrix columns by position.
rownames(pheno.tab) <- colnames(beta)

stopifnot(identical(rownames(pheno.tab), colnames(beta)))

required.columns <- c(
    "sampleId",
    "tissue",
    "timepoint",
    "bedFile"
)

if (!all(required.columns %in% colnames(pheno.tab))) {
    stop(
        "Expected sample annotation columns are missing.\nFound columns: ",
        paste(colnames(pheno.tab), collapse = ", ")
    )
}

# Set biologically meaningful category order.
pheno.tab$tissue <- factor(
    pheno.tab$tissue,
    levels = c("kidney", "liver")
)

pheno.tab$timepoint <- factor(
    pheno.tab$timepoint,
    levels = c("E14.5", "E15.5")
)

# Combined tissue-timepoint grouping.
pheno.tab$group <- factor(
    paste(pheno.tab$tissue, pheno.tab$timepoint, sep = "_"),
    levels = c(
        "kidney_E14.5",
        "kidney_E15.5",
        "liver_E14.5",
        "liver_E15.5"
    )
)

# Name metadata vectors with the corresponding beta-matrix columns.
plot.sample.id <- pheno.tab$sampleId
names(plot.sample.id) <- colnames(beta)

tissue <- pheno.tab$tissue
names(tissue) <- colnames(beta)

timepoint <- pheno.tab$timepoint
names(timepoint) <- colnames(beta)

group <- pheno.tab$group
names(group) <- colnames(beta)

print(pheno.tab[, required.columns])
print(table(pheno.tab$tissue, pheno.tab$timepoint))

###############################################################################
# 5. Define colours
###############################################################################

tissue.cols <- c(
    kidney = "#1B9E77",
    liver = "#D95F02"
)

timepoint.cols <- c(
    "E14.5" = "#7570B3",
    "E15.5" = "#E7298A"
)

group.cols <- c(
    kidney_E14.5 = "#66C2A5",
    kidney_E15.5 = "#1B9E77",
    liver_E14.5 = "#FC8D62",
    liver_E15.5 = "#D95F02"
)

###############################################################################
# Create one complete-case CpG matrix for all later full-data density/boxplots.
# This retains every CpG with a valid beta value in all samples.
# No random sampling is performed anywhere in this script.
###############################################################################

complete.site.indices <- which(rowSums(is.na(beta)) == 0)

beta.complete <- beta[
    complete.site.indices,
    ,
    drop = FALSE
]

n.complete.sites <- nrow(beta.complete)

message(
    "Complete CpG sites retained for full-data plots: ",
    format(n.complete.sites, big.mark = ",")
)

###############################################################################
# TASK 2. Manual PCA and scree plot
###############################################################################

# Calculate variation across samples for every complete CpG.
site.sd <- matrixStats::rowSds(beta.complete)

# Select the 10,000 most variable CpGs for PCA.
n.pca.sites <- min(10000, nrow(beta.complete))

top.pca.indices <- order(
    site.sd,
    decreasing = TRUE
)[seq_len(n.pca.sites)]

beta.pca <- beta.complete[
    top.pca.indices,
    ,
    drop = FALSE
]

# PCA expects samples in rows and CpGs/features in columns.
pca <- prcomp(
    t(beta.pca),
    center = TRUE,
    scale. = FALSE
)

# Percentage variance explained by each PC.
pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca.df <- data.frame(
    sampleId = plot.sample.id[rownames(pca$x)],
    tissue = tissue[rownames(pca$x)],
    timepoint = timepoint[rownames(pca$x)],
    group = group[rownames(pca$x)],
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    row.names = rownames(pca$x)
)

p.pca <- ggplot(
    pca.df,
    aes(
        x = PC1,
        y = PC2,
        colour = tissue,
        shape = timepoint,
        label = sampleId
    )
) +
    geom_point(size = 4, alpha = 0.90) +
    geom_text(
        vjust = -0.8,
        size = 3,
        show.legend = FALSE
    ) +
    scale_colour_manual(values = tissue.cols) +
    scale_shape_manual(values = c("E14.5" = 16, "E15.5" = 17)) +
    labs(
        title = "PCA of the 10,000 most variable CpG sites",
        x = paste0("PC1 (", round(pve[1], 1), "% variance)"),
        y = paste0("PC2 (", round(pve[2], 1), "% variance)"),
        colour = "Tissue",
        shape = "Timepoint"
    ) +
    theme_classic(base_size = 14)

ggsave(
    filename = file.path(
        plot.dir,
        "PCA_top_10000_variable_CpGs.png"
    ),
    plot = p.pca,
    width = 10,
    height = 7,
    dpi = 300
)

# Scree plot: up to eight PCs because there are eight samples.
n.scree <- min(8, length(pve))

scree.df <- data.frame(
    PC = factor(seq_len(n.scree)),
    VarianceExplained = pve[seq_len(n.scree)]
)

p.scree <- ggplot(
    scree.df,
    aes(x = PC, y = VarianceExplained)
) +
    geom_col(fill = "#2C7FB8") +
    geom_text(
        aes(label = paste0(round(VarianceExplained, 1), "%")),
        vjust = -0.4,
        size = 4
    ) +
    labs(
        title = "Scree plot: variance explained by principal components",
        x = "Principal component",
        y = "Variance explained (%)"
    ) +
    coord_cartesian(
        ylim = c(0, max(scree.df$VarianceExplained) * 1.15)
    ) +
    theme_classic(base_size = 14)

ggsave(
    filename = file.path(
        plot.dir,
        "PCA_scree_plot.png"
    ),
    plot = p.scree,
    width = 9,
    height = 6,
    dpi = 300
)

###############################################################################
# TASK 3. Density plot by tissue: ALL complete CpG sites
###############################################################################

# Convert every complete CpG methylation value in every sample to long format.
# as.vector() reads by columns, so each metadata label is repeated once per CpG.
density.df <- data.frame(
    beta = as.vector(beta.complete),
    tissue = rep(tissue, each = n.complete.sites),
    timepoint = rep(timepoint, each = n.complete.sites)
)

p.density.tissue <- ggplot(
    density.df,
    aes(x = beta, colour = tissue, fill = tissue)
) +
    geom_density(alpha = 0.20, linewidth = 0.9) +
    scale_colour_manual(values = tissue.cols) +
    scale_fill_manual(values = tissue.cols) +
    labs(
        title = "Genome-wide methylation distributions by tissue",
        subtitle = paste(
            "All", format(n.complete.sites, big.mark = ","),
            "complete CpG sites per sample"
        ),
        x = "Methylation beta value",
        y = "Density",
        colour = "Tissue",
        fill = "Tissue"
    ) +
    theme_classic(base_size = 14)

ggsave(
    filename = file.path(
        plot.dir,
        "Density_methylation_by_tissue_FULL_DATA.png"
    ),
    plot = p.density.tissue,
    width = 10,
    height = 7,
    dpi = 300
)

# Optional: same full data, split into one panel per developmental timepoint.
p.density.tissue.time <- ggplot(
    density.df,
    aes(x = beta, colour = tissue, fill = tissue)
) +
    geom_density(alpha = 0.20, linewidth = 0.9) +
    facet_wrap(~timepoint) +
    scale_colour_manual(values = tissue.cols) +
    scale_fill_manual(values = tissue.cols) +
    labs(
        title = "Genome-wide methylation by tissue and timepoint",
        subtitle = paste(
            "All", format(n.complete.sites, big.mark = ","),
            "complete CpG sites per sample"
        ),
        x = "Methylation beta value",
        y = "Density",
        colour = "Tissue",
        fill = "Tissue"
    ) +
    theme_classic(base_size = 14)

ggsave(
    filename = file.path(
        plot.dir,
        "Density_methylation_tissue_by_timepoint_FULL_DATA.png"
    ),
    plot = p.density.tissue.time,
    width = 11,
    height = 7,
    dpi = 300
)

###############################################################################
# TASK 4. Density plot by CpG-island context: ALL eligible CpG sites
###############################################################################

site.anno <- annotation(
    rnb.set,
    type = "sites",
    add.names = TRUE
)

# Reorder annotation rows to exactly match beta rows.
site.anno <- site.anno[
    rownames(beta),
    ,
    drop = FALSE
]

print(colnames(site.anno))

# Find the annotation column containing Island/Shore/Shelf/OpenSea.
context.matches <- sapply(site.anno, function(x) {
    x <- as.character(x)
    sum(
        x %in% c("Island", "Shore", "Shelf", "OpenSea"),
        na.rm = TRUE
    )
})

if (max(context.matches) == 0) {
    stop(
        "Could not automatically detect the CpG-context annotation column.\n",
        "Inspect colnames(site.anno) and context.matches."
    )
}

context.column <- names(which.max(context.matches))

cat("Detected CpG-context column:", context.column, "\n")

cpg.context <- as.character(site.anno[[context.column]])

valid.context <- c(
    "Island",
    "Shore",
    "Shelf",
    "OpenSea"
)

# Retain all context-annotated CpGs with valid beta values in all samples.
# No random sampling.
context.site.indices <- which(
    cpg.context %in% valid.context &
    rowSums(is.na(beta)) == 0
)

beta.context <- beta[
    context.site.indices,
    ,
    drop = FALSE
]

n.context.sites <- nrow(beta.context)

message(
    "Complete CpG sites used for full context-density plot: ",
    format(n.context.sites, big.mark = ",")
)

print(table(
    cpg.context[context.site.indices],
    useNA = "ifany"
))

# Repeat every site context for every sample column.
context.df <- data.frame(
    beta = as.vector(beta.context),
    context = rep(
        cpg.context[context.site.indices],
        times = ncol(beta.context)
    )
)

context.df$context <- factor(
    context.df$context,
    levels = c("Island", "Shore", "Shelf", "OpenSea")
)

p.density.context <- ggplot(
    context.df,
    aes(x = beta, colour = context, fill = context)
) +
    geom_density(alpha = 0.18, linewidth = 0.8) +
    labs(
        title = "DNA methylation distribution by CpG context",
        subtitle = paste(
            "All", format(n.context.sites, big.mark = ","),
            "complete annotated CpG sites across all samples"
        ),
        x = "Methylation beta value",
        y = "Density",
        colour = "CpG context",
        fill = "CpG context"
    ) +
    theme_classic(base_size = 14)

ggsave(
    filename = file.path(
        plot.dir,
        "Density_methylation_by_CpG_context_FULL_DATA.png"
    ),
    plot = p.density.context,
    width = 10,
    height = 7,
    dpi = 300
)

###############################################################################
# TASK 5. Heatmap of the 1,000 most variable CpG sites
###############################################################################

n.top.sites <- min(1000, nrow(beta.complete))

top.1000.indices <- order(
    site.sd,
    decreasing = TRUE
)[seq_len(n.top.sites)]

heat.sites <- beta.complete[
    top.1000.indices,
    ,
    drop = FALSE
]

# Z-score each CpG across samples to show relative high/low methylation.
heat.sites.z <- t(scale(t(heat.sites)))
heat.sites.z[is.na(heat.sites.z)] <- 0

sample.annotation <- data.frame(
    Tissue = tissue,
    Timepoint = timepoint,
    Group = group,
    row.names = colnames(beta)
)

pheatmap(
    heat.sites.z,
    annotation_col = sample.annotation,
    annotation_colors = list(
        Tissue = tissue.cols,
        Timepoint = timepoint.cols,
        Group = group.cols
    ),
    show_rownames = FALSE,
    show_colnames = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = "1,000 most variable CpG sites",
    filename = file.path(
        plot.dir,
        "Heatmap_top_1000_variable_CpGs.png"
    ),
    width = 10,
    height = 12
)

###############################################################################
# TASK 6. Load, summarize, and save mystery regions
###############################################################################

mystery.annotation.file <- paste0(
    "/vol/COMPEPIWS/data/annotation/",
    "annotation_mm10_mystery.RData"
)

# Register the supplied mystery-region annotation under the name "mystery".
rnb.load.annotation(
    mystery.annotation.file,
    "mystery"
)

print(rnb.region.types("mm10"))

if (!("mystery" %in% rnb.region.types("mm10"))) {
    stop("The custom mystery annotation was not loaded successfully.")
}

# Calculate a coverage-weighted beta methylation value for each mystery region.
rnb.set <- summarize.regions(
    object = rnb.set,
    region.type = "mystery",
    aggregation = "coverage.weighted",
    overwrite = TRUE
)

print(summarized.regions(rnb.set))

# Required by save.rnb.set(): set ZIP command BEFORE saving.
Sys.setenv(
    R_ZIPCMD = "/vol/COMPEPIWS/conda/miniconda3/envs/core/bin/zip"
)

message("ZIP command: ", Sys.getenv("R_ZIPCMD"))
message("ZIP exists: ", file.exists(Sys.getenv("R_ZIPCMD")))
message("ZIP executable: ", file.access(Sys.getenv("R_ZIPCMD"), 1) == 0)

if (!file.exists(Sys.getenv("R_ZIPCMD"))) {
    stop("ZIP executable not found. Check the R_ZIPCMD path.")
}

rnb.set.mystery.file <- file.path(
    report.dir,
    "rnbSet_preprocessed_with_mystery"
)

# Creates rnbSet_preprocessed_with_mystery.zip.
# The original rnbSet_preprocessed object is not overwritten.
save.rnb.set(
    rnb.set,
    rnb.set.mystery.file
)

###############################################################################
# TASK 7. Heatmap of the 100 most variable mystery regions
###############################################################################

mystery.beta <- meth(
    rnb.set,
    type = "mystery",
    row.names = TRUE
)

# Retain mystery regions with beta values in all samples.
mystery.complete <- mystery.beta[
    rowSums(is.na(mystery.beta)) == 0,
    ,
    drop = FALSE
]

mystery.sd <- matrixStats::rowSds(mystery.complete)

n.top.mystery <- min(
    100,
    nrow(mystery.complete)
)

top.100.mystery.indices <- order(
    mystery.sd,
    decreasing = TRUE
)[seq_len(n.top.mystery)]

heat.mystery <- mystery.complete[
    top.100.mystery.indices,
    ,
    drop = FALSE
]

# Z-score each mystery region across samples.
heat.mystery.z <- t(scale(t(heat.mystery)))
heat.mystery.z[is.na(heat.mystery.z)] <- 0

pheatmap(
    heat.mystery.z,
    annotation_col = sample.annotation,
    annotation_colors = list(
        Tissue = tissue.cols,
        Timepoint = timepoint.cols,
        Group = group.cols
    ),
    show_rownames = FALSE,
    show_colnames = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = "100 most variable mystery regions",
    filename = file.path(
        plot.dir,
        "Heatmap_top_100_variable_mystery_regions.png"
    ),
    width = 10,
    height = 10
)

###############################################################################
# TASK 8. Boxplots comparing samples: ALL complete CpG sites
###############################################################################

# Convert every complete CpG beta value across all samples to long format.
# No random sampling.
box.df <- data.frame(
    beta = as.vector(beta.complete),
    sampleId = rep(plot.sample.id, each = n.complete.sites),
    tissue = rep(tissue, each = n.complete.sites),
    timepoint = rep(timepoint, each = n.complete.sites),
    group = rep(group, each = n.complete.sites)
)

# Preserve sample-sheet order on x-axis.
box.df$sampleId <- factor(
    box.df$sampleId,
    levels = plot.sample.id
)

box.df$group <- factor(
    box.df$group,
    levels = c(
        "kidney_E14.5",
        "kidney_E15.5",
        "liver_E14.5",
        "liver_E15.5"
    )
)

p.boxplot.sample <- ggplot(
    box.df,
    aes(x = sampleId, y = beta, fill = tissue)
) +
    geom_boxplot(
        outlier.size = 0.10,
        outlier.alpha = 0.05,
        linewidth = 0.3
    ) +
    scale_fill_manual(values = tissue.cols) +
    labs(
        title = "Genome-wide DNA methylation by sample",
        subtitle = paste(
            "All", format(n.complete.sites, big.mark = ","),
            "complete CpG sites per sample"
        ),
        x = "Sample",
        y = "Methylation beta value",
        fill = "Tissue"
    ) +
    theme_classic(base_size = 14) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        legend.position = "top"
    )

ggsave(
    filename = file.path(
        plot.dir,
        "Boxplot_global_methylation_by_sample_FULL_DATA.png"
    ),
    plot = p.boxplot.sample,
    width = 14,
    height = 7,
    dpi = 300
)

# Optional: same complete full data, grouped by tissue and developmental stage.
p.boxplot.group <- ggplot(
    box.df,
    aes(x = group, y = beta, fill = tissue)
) +
    geom_boxplot(
        outlier.size = 0.10,
        outlier.alpha = 0.05,
        linewidth = 0.3
    ) +
    scale_fill_manual(values = tissue.cols) +
    labs(
        title = "Genome-wide methylation by tissue and timepoint",
        subtitle = paste(
            "All", format(n.complete.sites, big.mark = ","),
            "complete CpG sites per sample"
        ),
        x = "Group",
        y = "Methylation beta value",
        fill = "Tissue"
    ) +
    theme_classic(base_size = 14) +
    theme(
        axis.text.x = element_text(
            angle = 30,
            hjust = 1
        ),
        legend.position = "top"
    )

ggsave(
    filename = file.path(
        plot.dir,
        "Boxplot_global_methylation_by_group_FULL_DATA.png"
    ),
    plot = p.boxplot.group,
    width = 10,
    height = 7,
    dpi = 300
)

###############################################################################
# 9. Confirm output files
###############################################################################

print(list.files(
    plot.dir,
    full.names = TRUE
))
