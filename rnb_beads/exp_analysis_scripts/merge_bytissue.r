library(RnBeads)

rnb.options(assembly = "mm10")

report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
out.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/exploratory_analysis_3.3"
bedgraph.dir <- file.path(out.dir, "bedGraph")

if (!dir.exists(bedgraph.dir)) {
    dir.create(bedgraph.dir, recursive = TRUE)
}

rnb.set.file <- file.path(report.dir, "rnbSet_preprocessed")
rnb.set <- load.rnb.set(rnb.set.file)

rnb.set.merged <- mergeSamples(rnb.set, grp.col = "tissue")

cat("Merged samples:", samples(rnb.set.merged), "\n")

rnb.RnBSet.to.bedGraph(
    rnb.set = rnb.set.merged,
    out.dir = bedgraph.dir,
    reg.type = "sites"
)
