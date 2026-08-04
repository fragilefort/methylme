library(RnBeads)
library(ggplot2)
library(GenomicRanges)
library(purrr)
library(dplyr)

rnb.options(assembly = "mm10")

report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
out.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/exploratory_analysis_3.3"
setwd(out.dir)

rnb.set.file <- file.path(report.dir, "rnbSet_preprocessed")
rnb.set <- load.rnb.set(rnb.set.file)
rnb.set.merged <- mergeSamples(rnb.set, grp.col = "tissue")

get_annotation_gr <- function(rnb_obj, type) {
    df_annot <- annotation(rnb_obj, type = type)
    makeGRangesFromDataFrame(
        df_annot,
        seqnames.field = "Chromosome",
        start.field = "Start",
        end.field = "End",
        strand.field = "Strand",
        keep.extra.columns = TRUE
    )
}

gr_sites <- get_annotation_gr(rnb.set.merged, type = "sites")

gr_promoters <- get_annotation_gr(rnb.set.merged, type = "promoters")
df_promoters <- data.frame(
    Region_Type = "Promoters",
    CpG_Count = countOverlaps(gr_promoters, gr_sites)
)

all_regions <- summarized.regions(rnb.set.merged)
mystery_types <- grep("mystery", all_regions, value = TRUE, ignore.case = TRUE)

df_counts_list <- list(df_promoters)

if (length(mystery_types) > 0) {
    gr_mystery <- get_annotation_gr(rnb.set.merged, type = mystery_types[1])
    df_mystery <- data.frame(
        Region_Type = "Mystery Regions",
        CpG_Count = countOverlaps(gr_mystery, gr_sites)
    )
    df_counts_list[[2]] <- df_mystery
}

df_site_freq <- bind_rows(df_counts_list)

p_hist <- ggplot(df_site_freq, aes(x = CpG_Count, fill = Region_Type)) +
    geom_histogram(binwidth = 2, color = "black", alpha = 0.7, position = "identity") +
    facet_wrap(~Region_Type, scales = "free") +
    labs(
        title = "Frequency of CpG Sites per Region",
        x = "Number of CpG Sites per Region",
        y = "Frequency (Number of Regions)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

output_png <- file.path(out.dir, "cpg_site_frequency_promoters.png")
png(filename = output_png, width = 8, height = 5, units = "in", res = 300)
print(p_hist)
dev.off()
