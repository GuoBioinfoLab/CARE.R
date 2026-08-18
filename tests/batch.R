library(CARE.R)

control_path <- file.path("..", "public", "examples", "itcm_control_tpm.csv")
treatment_path <- file.path("..", "public", "examples", "itcm_treatment_tpm.csv")
gene_sets <- system.file("extdata", "cancer_aging_hallmark.json", package = "CARE.R")
if (!file.exists(control_path) || !file.exists(treatment_path)) {
  cat("BATCH TESTS SKIPPED: demo matrices not available\n")
  quit(status = 0)
}

control <- read.csv(control_path, check.names = FALSE, row.names = 1)
treatment <- read.csv(treatment_path, check.names = FALSE, row.names = 1)
expression <- cbind(
  setNames(control, paste0("CTRL_", seq_len(ncol(control)))),
  setNames(treatment, paste0("TRT_", seq_len(ncol(treatment))))
)
metadata <- data.frame(
  sample_id = colnames(expression),
  compound = "Demo compound",
  group = c(rep("control", ncol(control)), rep("treatment", ncol(treatment))),
  replicate = rep(seq_len(ncol(control)), 2),
  stringsAsFactors = FALSE
)

batch <- batch_analyze(expression, metadata, gene_set_path = gene_sets)
stopifnot(
  identical(batch$status, "ok"),
  nrow(batch$summary) == 1,
  identical(names(batch$results), "Demo compound"),
  batch$summary$control_samples == 3,
  batch$summary$treatment_samples == 3
)

bad_metadata <- metadata
bad_metadata$sample_id[1] <- "missing_sample"
bad_error <- tryCatch(
  batch_analyze(expression, bad_metadata, gene_set_path = gene_sets),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("sample_id", bad_error, fixed = TRUE))

demo_expression_path <- system.file(
  "extdata", "itcm_top5", "itcm_top5_expression.csv", package = "CARE.R"
)
demo_metadata_path <- system.file(
  "extdata", "itcm_top5", "itcm_top5_metadata.csv", package = "CARE.R"
)
if (file.exists(demo_expression_path) && file.exists(demo_metadata_path)) {
  top5 <- c(
    "Hupehenine", "Rosmarinic acid", "Z-Ligustilide",
    "Costunlide", "Fangchinoline"
  )
  demo_batch <- batch_analyze(
    demo_expression_path,
    demo_metadata_path,
    compound_names = top5,
    gene_set_path = gene_sets
  )
  stopifnot(
    identical(demo_batch$status, "ok"),
    identical(demo_batch$summary$compound, top5),
    all(demo_batch$summary$control_samples == 3),
    all(demo_batch$summary$treatment_samples == 3),
    all(vapply(demo_batch$results, Negate(is.null), logical(1)))
  )
}

cat("BATCH TESTS PASSED\n")
