# Load libraries (same as your thesis style)
library(reshape2)
library(ggplot2)
library(tidyverse)

# Load the two tables (adjust file paths/names accordingly)
liver <- read.table("liver_14.5_mouse_1_2kbW_bed_counts.txt", header = TRUE)
kidney <- read.table("kidney_14.5_mouse_1_2kbW_bed_counts.txt", header = TRUE)

## Q4: Dimension of each table
dim(liver)
dim(kidney)

## Q5: Column names of each table
colnames(liver)
colnames(kidney)

## Q6: Genome length from each dataset
# Genome length = sum of bin widths (End - Start) per table
genome_length_liver <- sum(liver$End - liver$Start)
genome_length_kidney <- sum(kidney$End - kidney$Start)
genome_length_liver
genome_length_kidney

## Q7: Concatenate vertically + add cell_type column
liver$cell_type <- "liver"
kidney$cell_type <- "kidney"
combined_df <- rbind(liver, kidney)

## Q8: Dimension of the new data frame
dim(combined_df)

## Q9: Reshape wide to long (excluding chr, start, end), using cell_type as id
long_df <- melt(combined_df[, !(names(combined_df) %in% c("Chr", "Start", "End"))],
                 id.vars = "cell_type",
                 variable.name = "variable",
                 value.name = "value")
dim(long_df)

## Q10: Density plot per variable, faceted by cell_type, transparent colors, save as PDF
density_plot <- ggplot(long_df, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cell_type) +
  labs(title = "Density of histone marks by cell type",
       x = "Signal value", y = "Density") +
  theme_minimal()

ggsave("density_plot.pdf", plot = density_plot, width = 10, height = 6)

## Q11: Re-plot with x-axis limited to 100, save as PDF
density_plot_xlim <- density_plot + xlim(0, 100)

ggsave("density_plot_xlim100.pdf", plot = density_plot_xlim, width = 10, height = 6)

## Q12: Concatenate horizontally, keeping only H3K27me3, H3K36me3, H3K9me3
horiz_df <- cbind(
  liver[, c("H3K27me3", "H3K36me3", "H3K9me3")],
  kidney[, c("H3K27me3", "H3K36me3", "H3K9me3")]
)

## Q13: Rename columns
colnames(horiz_df) <- c("H3K27me3_liver", "H3K36me3_liver", "H3K9me3_liver",
                         "H3K27me3_kidney", "H3K36me3_kidney", "H3K9me3_kidney")

## Q14: PCA (same as your thesis: prcomp with center/scale)
pca_result <- prcomp(horiz_df, center = TRUE, scale. = TRUE)

## Q15: Variance explained by PC1 and PC2
pca_var <- pca_result$sdev^2
pca_var_proportion <- pca_var / sum(pca_var) * 100
pca_var_proportion[1]  # PC1
pca_var_proportion[2]  # PC2

## Q16: Screeplot and cumulative screeplot
pdf("screeplot.pdf", width = 8, height = 6)
screeplot(pca_result, type = "lines", main = "Scree Plot")
dev.off()

cum_var <- cumsum(pca_var_proportion)
pdf("cumulative_screeplot.pdf", width = 8, height = 6)
plot(cum_var, type = "b", xlab = "Principal Component", ylab = "Cumulative Variance Explained (%)",
     main = "Cumulative Scree Plot")
dev.off()

## Q17: Explore PC loadings - how many are there?
pca_result$rotation
dim(pca_result$rotation)  # rows = variables, columns = PCs (should be 6x6 here)

## Q18: Plot PC1 and PC2 loadings
loadings_df <- as.data.frame(pca_result$rotation)
loadings_df$variable <- rownames(loadings_df)

loadings_plot <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = variable)) +
  geom_point(size = 3, color = "darkred") +
  geom_text(vjust = -0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "PC1 vs PC2 Loadings") +
  theme_minimal()

ggsave("pc_loadings_plot.pdf", plot = loadings_plot, width = 8, height = 6)
