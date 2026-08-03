library(RnBeads)

data.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/CASINO/methylation_coverage"
report.dir <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/rnbeads_reports"
sample.sheet <- "/vol/COMPEPIWS/groups/wgbs2/methylme/rnb_beads/sample_annotation.csv"

rnb.options(
    assembly = "mm10",
    import.bed.style = "bismarkCov",
    differential.comparison.columns = c("tissue"),
    differential.report.sites = FALSE,
    filtering.coverage.threshold = 10,
    filtering.low.coverage.masking = TRUE,
    filtering.high.coverage.outliers = TRUE,
    filtering.missing.value.quantile = 1,
    num.cores = 12
)

rnb.run.analysis(
    dir.reports = report.dir,
    sample.sheet = sample.sheet,
    data.dir = data.dir,
    data.type = "bismark.cov"
)
