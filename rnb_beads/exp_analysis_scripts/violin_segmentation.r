library(rtracklayer)
library(ggplot2)
library(purrr)
library(dplyr)

out.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/exploratory_analysis_3.3"
setwd(out.dir)

bed_files <- list.files(out.dir, pattern = "\\.bed$", full.names = TRUE)

df_sizes <- map_dfr(bed_files, function(bfile) {
  gr <- import(bfile, format = "BED")
  
  s_name <- if (grepl("kidney", bfile, ignore.case = TRUE)) "kidney" else "liver"
  
  data.frame(
    Sample = s_name,
    Category = as.character(mcols(gr)$name),
    Width_bp = width(gr)
  )
})

p_violin <- ggplot(df_sizes, aes(x = Category, y = Width_bp, fill = Category)) +
  geom_violin(scale = "width", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~ Sample) +
  labs(
    title = "Size Distribution of MethylSeekR Segmentation Categories",
    x = "Segmentation Category",
    y = "Region Width in Base Pairs (log10 scale)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

output_png <- file.path(out.dir, "segmentation_category_sizes_violin.png")
png(filename = output_png, width = 8, height = 6, units = "in", res = 300)
print(p_violin)
dev.off()
