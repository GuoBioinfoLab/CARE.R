# CARE.R

R package wrapping the analysis module of the **CARE** (Cancer Aging
Research Explorer for drug discovery) project.

## Installation

```r
# from the project root
install.packages("CARE_R",
                 repos = NULL, type = "source")
# or
install.packages("CARE.R_0.1.0.tar.gz", repos = NULL, type = "source")
```

Requires R >= 4.3.0, `GSVA` >= 1.50.0 from Bioconductor, and `jsonlite` >= 1.8.0.

```r
install.packages("jsonlite")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GSVA")
```

## Usage

```r
library(CARE.R)

res <- analyze_tpm(
  control_path   = "itcm_control_tpm.csv",
  treatment_path = "itcm_treatment_tpm.csv"
)

res$pebm              # EBM combined p-value (only if direction is consistent)
res$direction_consistent
res$hallmarks         # per-hallmark delta, one-sided Wilcoxon p-value, ssGSEA scores
res$distribution      # scaled chi-square density for the EBM tail plot
res$warnings          # small-sample and non-fatal input diagnostics
```

`run_analysis(control, treatment)` is a thin wrapper that accepts Shiny
`fileInput` objects, file paths, or genes-by-samples matrix/data.frame
objects. Matrix row names are gene IDs and matrix column names are sample IDs.

## Batch analysis

`batch_analyze()` runs the same CARE analysis for several compounds. The
expression input must be a genes-by-samples matrix (or CSV/TSV path), and the
metadata must contain a `sample_id` column whose values exactly match the
expression matrix column names. It also requires `compound` and `group`
columns, where `group` is `control` or `treatment`.

```r
demo_expression <- system.file(
  "extdata", "itcm_top5", "itcm_top5_expression.csv", package = "CARE.R"
)
demo_metadata <- system.file(
  "extdata", "itcm_top5", "itcm_top5_metadata.csv", package = "CARE.R"
)

batch <- batch_analyze(
  demo_expression,
  demo_metadata,
  compound_names = c(
    "Hupehenine", "Rosmarinic acid", "Z-Ligustilide",
    "Costunlide", "Fangchinoline"
  )
)
batch$summary[, c("compound", "status", "pebm", "control_samples",
                  "treatment_samples")]
batch$results[["Hupehenine"]]
```

The package includes this five-compound ITCM PEBM example under
`inst/extdata/itcm_top5/`. Failed compounds are isolated in the summary with
`status = "error"`; a mixed run returns `status = "partial"`.

## Visualization

The package includes static plots that mirror the web Analysis page. They
use base R graphics and do not add a `ggplot2` dependency.

```r
plot_hallmarks(res)                         # three hallmark boxplots
plot_ebm(res)                               # EBM density and right tail
plot_analysis(res, file = "care-analysis.png")
```

The `file` argument supports `.png`, `.pdf`, and `.svg`. With `file = NULL`,
the plot is drawn on the current graphics device. Each function returns the
analysis result invisibly.

## What the analysis does

1. Reads two genes-by-samples TPM matrices (CSV or TSV, first column =
   gene id, >= 2 samples), averages duplicated genes, validates
   non-negative finite values.
2. Intersects gene sets.
3. Runs ssGSEA with the current GSVA parameter-object API
   (`GSVA::ssgseaParam()` and `GSVA::gsva()`) on the three
   cancer-aging hallmarks shared by cancer and aging.
4. Computes per-hallmark median change (treatment − control), one-sided
   Wilcoxon `alternative = "less"` p-values, and direction consistency
   (all three medians must decrease).
5. When direction is consistent, combines the three p-values with
   Empirical Brown's Method (`empirical_browns_method()`, with Fisher
   fallback) into the single evidence score `pebm`, and returns the
   scaled chi-square distribution used for visualization.

Missing files, malformed gene-set JSON, duplicate sample names, non-numeric
TPM values, negative values, and insufficient hallmark coverage produce
actionable errors. Groups with fewer than four samples are analyzed but add
a warning to `res$warnings`; the result should be interpreted cautiously.

## Package layout

```
CARE_R/
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── R/
│   ├── analysis.R      # analyze_tpm(), run_analysis()
│   ├── EBM_metaP.R     # empirical_browns_method()
│   ├── batch.R          # batch_analyze()
│   └── visualization.R # plot_hallmarks(), plot_ebm(), plot_analysis()
└── inst/
    └── extdata/
        ├── cancer_aging_hallmark.json
        └── itcm_top5/
            ├── itcm_top5_expression.csv
            └── itcm_top5_metadata.csv
```
