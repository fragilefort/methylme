# ===== Libraries =====
library(reshape2)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(GenomicRanges)

# ===== Load Data =====
liver <- read.table("/vol/COMPEPIWS/groups/shared/liver_14.5_mouse_1_2kbW_bed_counts.txt", header = TRUE)
kidney <- read.table("/vol/COMPEPIWS/groups/shared/kidney_14.5_mouse_1_2kbW_bed_counts.txt", header = TRUE)

## Q4: Dimension of each table
dim(liver)
dim(kidney)

## Q5: Column names of each table
colnames(liver)
colnames(kidney)

## Q6: Genome length from each dataset
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

## Q10: Density plot per variable, faceted by cell_type, transparent colors
density_plot <- ggplot(long_df, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cell_type) +
  labs(title = "Density of histone marks by cell type",
       x = "Signal value", y = "Density") +
  theme_minimal()

ggsave("density_plot.pdf", plot = density_plot, width = 10, height = 6)

## Q11: Re-plot with x-axis limited to 100
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

## Q14: PCA
pca_result <- prcomp(horiz_df, center = TRUE, scale. = TRUE)

## Q15: Variance explained by PC1 and PC2
pca_var <- pca_result$sdev^2
pca_var_proportion <- pca_var / sum(pca_var) * 100
pca_var_proportion[1]
pca_var_proportion[2]

## Q16: Screeplot and cumulative screeplot
pdf("screeplot.pdf", width = 8, height = 6)
screeplot(pca_result, type = "lines", main = "Scree Plot")
dev.off()

cum_var <- cumsum(pca_var_proportion)
pdf("cumulative_screeplot.pdf", width = 8, height = 6)
plot(cum_var, type = "b", xlab = "Principal Component", ylab = "Cumulative Variance Explained (%)",
     main = "Cumulative Scree Plot")
dev.off()

## Q17: Explore PC loadings
pca_result$rotation
dim(pca_result$rotation)

## Q18: Plot PC1 and PC2 loadings
loadings_df <- as.data.frame(pca_result$rotation)
loadings_df$variable <- rownames(loadings_df)
loadings_df$mark <- sub("_.*", "", loadings_df$variable)     
loadings_df$tissue <- sub(".*_", "", loadings_df$variable)  

loadings_plot <- ggplot(loadings_df, aes(x = PC1, y = PC2, color = mark, shape = tissue, label = variable)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 4, show.legend = FALSE, max.overlaps = 20, box.padding = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("H3K27me3" = "#E63946", "H3K36me3" = "#2A9D8F", "H3K9me3" = "#457B9D")) +
  labs(title = "PC1 vs PC2 Loadings", color = "Histone mark", shape = "Tissue") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("pc_loadings_plot.pdf", plot = loadings_plot, width = 9, height = 7)


## 1. Create two GRanges objects from liver and kidney data frames
## columns 4,5,6 (H3K27me3, H3K36me3, H3K9me3) as metadata

liver_gr <- GRanges(
  seqnames = liver$Chr,
  ranges = IRanges(start = liver$Start, end = liver$End),
  H3K27me3 = liver$H3K27me3,
  H3K36me3 = liver$H3K36me3,
  H3K9me3 = liver$H3K9me3
)

kidney_gr <- GRanges(
  seqnames = kidney$Chr,
  ranges = IRanges(start = kidney$Start, end = kidney$End),
  H3K27me3 = kidney$H3K27me3,
  H3K36me3 = kidney$H3K36me3,
  H3K9me3 = kidney$H3K9me3
)

liver_gr
kidney_gr

## 2. Total number of bases covered in each object
total_bases_liver <- sum(width(liver_gr))
total_bases_kidney <- sum(width(kidney_gr))
total_bases_liver
total_bases_kidney

# Alternative: bases covered by the UNION of overlapping ranges (merged/non-redundant)
total_bases_liver_reduced <- sum(width(reduce(liver_gr)))
total_bases_kidney_reduced <- sum(width(reduce(kidney_gr)))
total_bases_liver_reduced
total_bases_kidney_reduced

## 3. Subset one object to only chr2 lines
liver_chr2 <- liver_gr[seqnames(liver_gr) == "chr2"]
liver_chr2

## 4. Shift one object's ranges by 100bp upstream
# "Upstream" depends on strand; if no strand info, shift() moves in genomic coordinate direction
# shift() with negative value moves the range to lower coordinates (upstream on + strand)
liver_shifted <- shift(liver_gr, -100)
liver_shifted

## 5. Find overlapping regions between the two objects
overlaps <- findOverlaps(liver_gr, kidney_gr)
overlaps

# To extract the actual overlapping regions as a GRanges object:
overlapping_regions <- pintersect(liver_gr[queryHits(overlaps)], kidney_gr[subjectHits(overlaps)])
overlapping_regions

## 6. Make a GRangesList from the two GRanges objects
gr_list <- GRangesList(liver = liver_gr, kidney = kidney_gr)
gr_list
