library(RnBeads)
library(dplyr)

rnb.options(assembly = "mm10")

report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
out.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/exploratory_analysis_3.3"
setwd(out.dir)

rnb.set.file <- file.path(report.dir, "rnbSet_preprocessed")
rnb.set <- load.rnb.set(rnb.set.file)

mystery_file <- "/vol/COMPEPIWS/data/annotation/annotation_mm10_mystery.RData"
rnb.load.annotation(mystery_file, "mystery")

rnb.set <- summarize.regions(rnb.set, "mystery")

diffmeth_obj <- rnb.execute.computeDiffMeth(
    rnb.set,
    pheno.cols = "tissue",
    region.types = "mystery"
)

comparisons <- get.comparisons(diffmeth_obj)

diff_table <- get.table(
    diffmeth_obj,
    comparison = comparisons[1],
    region.type = "mystery",
    return.data.frame = TRUE
)

# Add genomic coordinates to the differential methylation table
annot <- annotation(rnb.set, type = "mystery")

# Combine coordinate columns with the differential methylation results
coord_cols <- intersect(c("Chromosome", "Start", "End", "Strand"), colnames(annot))
diff_table_annotated <- cbind(annot[, coord_cols], diff_table)

diff_table_sorted <- diff_table_annotated %>%
    arrange(combinedRank)

message("--- Top 10 Ranked Mystery Regions ---")
print(head(diff_table_sorted[, c(coord_cols, "combinedRank", "mean.mean.diff", "comb.p.adj.fdr")], 10))

output_csv <- file.path(out.dir, "diff_meth_mystery_regions_annotated.csv")
write.csv(diff_table_sorted, file = output_csv, row.names = FALSE)

message("Annotated and sorted table successfully saved to: ", output_csv)
