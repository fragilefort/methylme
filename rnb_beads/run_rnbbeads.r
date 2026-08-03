library(RnBeads)

data.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/CASINO/methylation_coverage"
report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
sample.sheet <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/sample_annotation.csv"

if (dir.exists(report.dir)) {
    unlink(report.dir, recursive = TRUE, force = TRUE)
}

rnb.options(
    assembly = "mm10",
    import.bed.style = "bismarkCov",
    differential.comparison.columns = c("tissue", "timepoint"),
    differential.report.sites = FALSE,
    filtering.coverage.threshold = 10,
    filtering.low.coverage.masking = TRUE,
    filtering.high.coverage.outliers = TRUE,
    filtering.missing.value.quantile = 0
)

ds <- list(
    dir = data.dir,
    sample.sheet = sample.sheet
)

rnb.run.analysis(
    dir.reports = report.dir,
    data.source = ds,
    data.type = "bed.dir",
    initialize.reports = TRUE
)
