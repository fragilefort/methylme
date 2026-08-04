library(RnBeads)
library(ggplot2)
library(purrr)
library(dplyr)

rnb.options(assembly = "mm10")

out.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/exploratory_analysis_3.3"
seg.dir <- file.path(out.dir, "segmentation")
setwd(out.dir)

merged_save_path <- file.path(out.dir, "rnbSet_merged")
rnb.set.merged <- load.rnb.set(merged_save_path)

sample_names <- samples(rnb.set.merged)
training_chr <- c("chr18", "chr19")

rnb.set.merged <- reduce(sample_names, function(rnb_obj, s) {
    message("Running segmentation for sample: ", s)

    rnb_obj <- rnb.execute.segmentation(
        rnb.set = rnb_obj,
        sample.name = s,
        chr.sel = training_chr,
        plot.path = seg.dir
    )

    rnb.bed.from.segmentation(
        rnb.set = rnb_obj,
        sample.name = s,
        type = "final",
        store.path = seg.dir
    )

    rnb_obj
}, .init = rnb.set.merged)

save.rnb.set(rnb.set.merged, path = merged_save_path, archive = FALSE)

categories <- c("UMRs", "LMRs", "HMDs", "PMDs")

df_all <- map_dfr(sample_names, function(s) {
    map_dfr(categories, function(cat_type) {
        rtype <- paste0(cat_type, "_", s)
        m_vals <- meth(rnb.set.merged, type = rtype)[, s]

        data.frame(
            Sample = s,
            Category = sub("s$", "", cat_type),
            Methylation = as.numeric(m_vals)
        )
    })
})

p_seg <- ggplot(df_all, aes(x = Category, y = Methylation, fill = Category)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    facet_wrap(~Sample) +
    labs(
        title = "Methylation Distribution Across Segmentation Categories",
        x = "Segmentation Category",
        y = "Mean Methylation (β-value)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

output_png <- file.path(seg.dir, "segmentation_categories_methylation.png")
png(filename = output_png, width = 8, height = 6, units = "in", res = 300)
print(p_seg)
dev.off()

message("Segmentation complete. Output saved in: ", seg.dir)
