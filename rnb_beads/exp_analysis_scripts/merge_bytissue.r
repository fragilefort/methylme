library(RnBeads)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(purrr)
library(dplyr)

rnb.options(assembly = "mm10")

report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
out.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/exploratory_analysis_3.3"
seg.dir <- file.path(out.dir, "segmentation")

if (!dir.exists(seg.dir)) {
    dir.create(seg.dir, showWarnings = FALSE, recursive = TRUE)
}

setwd(out.dir)

# Load preprocessed set and merge samples by tissue
message("Loading preprocessed dataset and merging samples by tissue...")
rnb.set.file <- file.path(report.dir, "rnbSet_preprocessed")
rnb.set <- load.rnb.set(rnb.set.file)
rnb.set.merged <- mergeSamples(rnb.set, grp.col = "tissue")

# Save merged dataset object
merged_save_path <- file.path(out.dir, "rnbSet_merged")
save.rnb.set(rnb.set.merged, path = merged_save_path, archive = FALSE)
message("Saved merged dataset to: ", merged_save_path)

# Export continuous signal bedGraph for each merged tissue
annot_sites <- annotation(rnb.set.merged, type = "sites")
meth_sites <- meth(rnb.set.merged, type = "sites")

for (s in samples(rnb.set.merged)) {
    message("Exporting continuous methylation signal for: ", s)

    df_signal <- data.frame(
        seqnames = annot_sites$Chromosome,
        start    = annot_sites$Start,
        end      = annot_sites$End,
        score    = meth_sites[, s]
    ) %>%
        filter(!is.na(score))

    gr_signal <- makeGRangesFromDataFrame(df_signal, keep.extra.columns = TRUE)

    out_bedgraph <- file.path(seg.dir, paste0(s, "_merged_signal.bedgraph"))
    export.bedGraph(gr_signal, con = out_bedgraph)

    message("Saved continuous signal: ", out_bedgraph)
}

sample_names <- samples(rnb.set.merged)
training_chr <- c("chr18", "chr19")

rnb.set.merged <- reduce(sample_names, function(rnb_obj, s) {
    message("Running segmentation for sample: ", s, " on ", paste(training_chr, collapse = ", "))

    rnb_obj <- rnb.execute.segmentation(
        rnb.set     = rnb_obj,
        sample.name = s,
        chr.sel     = training_chr,
        plot.path   = seg.dir
    )

    rnb.bed.from.segmentation(
        rnb.set     = rnb_obj,
        sample.name = s,
        type        = "final",
        store.path  = seg.dir
    )

    rnb_obj
}, .init = rnb.set.merged)

categories <- c("UMRs", "LMRs", "HMDs", "PMDs")

df_all <- map_dfr(sample_names, function(s) {
    map_dfr(categories, function(cat_type) {
        rtype <- paste0(cat_type, "_", s)
        m_vals <- meth(rnb.set.merged, type = rtype)[, s]

        data.frame(
            Sample      = s,
            Category    = sub("s$", "", cat_type), # Formats UMRs -> UMR
            Methylation = as.numeric(m_vals)
        )
    })
})

p_seg <- ggplot(df_all, aes(x = Category, y = Methylation, fill = Category)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    facet_wrap(~Sample) +
    labs(
        title    = "Methylation Distribution Across Segmentation Categories",
        subtitle = paste("Trained on:", paste(training_chr, collapse = ", ")),
        x        = "Segmentation Category",
        y        = "Mean Methylation (β-value)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

output_png <- file.path(seg.dir, "segmentation_categories_methylation.png")
png(filename = output_png, width = 8, height = 6, units = "in", res = 300)
print(p_seg)
dev.off()

message("Workflow complete! Output files generated in: ", seg.dir)
