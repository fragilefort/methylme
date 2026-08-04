library(RnBeads)
library(ggplot2)
library(reshape)

rnb.options(assembly = "mm10")
rnb.load.annotation("/vol/COMPEPIWS/data/annotation/annotation_mm10_mystery.RData", "mystery")

report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
rnb.set.file <- file.path(report.dir, "rnbSet_preprocessed")
rnb.set <- load.rnb.set(rnb.set.file)

rnb.set <- summarize.regions(rnb.set, "mystery")

cpg_means <- rowMeans(meth(rnb.set, type = "sites"), na.rm = TRUE)
mystery_means <- rowMeans(meth(rnb.set, type = "mystery"), na.rm = TRUE)

meth_list <- list(
    "CpGs" = cpg_means,
    "Mystery Regions" = mystery_means
)
df_long <- melt(meth_list)
colnames(df_long) <- c("Methylation", "Feature")

ggplot(df_long, aes(x = Feature, y = Methylation, fill = Feature)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    labs(
        title = "Methylation Distribution: CpGs vs. Mystery Regions",
        x = "Genomic Feature",
        y = "Mean Methylation (β-value)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
