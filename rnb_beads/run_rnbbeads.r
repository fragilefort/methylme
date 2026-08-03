library(RnBeads)

data.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/CASINO/methylation_coverage"
report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
sample.sheet <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/sample_annotation.csv"

if (!dir.exists(report.dir)) {
    dir.create(report.dir, recursive = TRUE, showWarnings = FALSE)
}
rnb.initialize.reports(report.dir)

rnb.options("assembly" = "mm10")
rnb.options("import.bed.style" = "bismarkCov")
rnb.options("differential.comparison.columns" = "tissue")
rnb.options("differential.report.sites" = FALSE)
rnb.options("filtering.coverage.threshold" = 10)
rnb.options("filtering.low.coverage.masking" = TRUE)
rnb.options("filtering.high.coverage.outliers" = TRUE)
rnb.options("filtering.missing.value.quantile" = 1)

ds <- list(
    csv.file = sample.sheet,
    bed.dir = data.dir
)

rnb.run.analysis(
    dir.reports = report.dir,
    data.source = ds,
    data.type = "bs.bed.dir",
    initialize.reports = FALSE
)
