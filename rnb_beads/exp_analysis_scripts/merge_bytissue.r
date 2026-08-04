library(RnBeads)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

rnb.options(assembly = "mm10")

report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
out.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/exploratory_analysis_3.3"
seg.dir <- file.path(out.dir, "segmentation")

if (!dir.exists(seg.dir)) {
    dir.create(seg.dir, showWarnings = FALSE, recursive = TRUE)
}

setwd(out.dir)

rnb.set.file <- file.path(report.dir, "rnbSet_preprocessed")
rnb.set <- load.rnb.set(rnb.set.file)
rnb.set.merged <- mergeSamples(rnb.set, grp.col = "tissue")

merged_save_path <- file.path(out.dir, "rnbSet_merged")
save.rnb.set(rnb.set.merged, path = merged_save_path, archive = FALSE)
message("Merged RnBSet saved to: ", merged_save_path)

annot_sites <- annotation(rnb.set.merged, type = "sites")
meth_sites <- meth(rnb.set.merged, type = "sites")

for (s in samples(rnb.set.merged)) {
    message("Exporting merged methylation signal for: ", s)

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
