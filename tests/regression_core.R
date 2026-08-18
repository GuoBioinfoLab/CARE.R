library(CARE.R)

control <- file.path("..", "public", "examples", "itcm_control_tpm.csv")
treatment <- file.path("..", "public", "examples", "itcm_treatment_tpm.csv")
gene_sets <- system.file("extdata", "cancer_aging_hallmark.json", package = "CARE.R")
if (!file.exists(control) || !file.exists(treatment)) {
  cat("CORE REGRESSION TESTS SKIPPED: demo matrices not available\n")
  quit(status = 0)
}

result <- analyze_tpm(control, treatment, gene_sets)
stopifnot(identical(result$status, "ok"))
stopifnot(isTRUE(result$direction_consistent))
stopifnot(result$input_summary$common_genes == 20030)
stopifnot(result$input_summary$control_samples == 3)
stopifnot(result$input_summary$treatment_samples == 3)
stopifnot(is.character(result$warnings))
stopifnot(any(grepl("constant", result$warnings, ignore.case = TRUE)))

control_matrix <- as.matrix(read.csv(control, check.names = FALSE, row.names = 1))
treatment_matrix <- as.matrix(read.csv(treatment, check.names = FALSE, row.names = 1))
matrix_result <- analyze_tpm(control_matrix, treatment_matrix, gene_sets)
stopifnot(isTRUE(all.equal(
  result$pebm,
  matrix_result$pebm,
  tolerance = 1e-8
)))
wrapper_result <- run_analysis(control_matrix, treatment_matrix, gene_set_path = gene_sets)
stopifnot(isTRUE(all.equal(
  result$pebm,
  wrapper_result$pebm,
  tolerance = 1e-8
)))

baseline <- vapply(result$hallmarks, function(x) x$delta, numeric(1))
raw_sets <- jsonlite::fromJSON(gene_sets, simplifyVector = FALSE)
reordered_path <- tempfile(fileext = ".json")
jsonlite::write_json(raw_sets[c(3, 1, 2, 5, 4, 6)], reordered_path, auto_unbox = TRUE)
reordered <- analyze_tpm(control, treatment, reordered_path)
stopifnot(isTRUE(all.equal(
  baseline,
  vapply(reordered$hallmarks, function(x) x$delta, numeric(1)),
  tolerance = 1e-8
)))

missing_error <- tryCatch(
  analyze_tpm("missing-control.csv", treatment, gene_sets),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("Control", missing_error, fixed = TRUE))

control_frame <- read.csv(control, nrows = 30, check.names = FALSE)
control_frame[1, 2] <- "not-a-number"
bad_control <- tempfile(fileext = ".csv")
write.csv(control_frame, bad_control, row.names = FALSE, quote = TRUE)
bad_error <- tryCatch(
  analyze_tpm(bad_control, treatment, gene_sets),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("Control", bad_error, fixed = TRUE))

small_control <- tempfile(fileext = ".csv")
small_treatment <- tempfile(fileext = ".csv")
write.csv(read.csv(control, check.names = FALSE)[, 1:3], small_control, row.names = FALSE)
write.csv(read.csv(treatment, check.names = FALSE)[, 1:3], small_treatment, row.names = FALSE)
small <- analyze_tpm(small_control, small_treatment, gene_sets)
stopifnot(any(grepl("small sample", small$warnings, ignore.case = TRUE)))

unlink(c(reordered_path, bad_control, small_control, small_treatment))
cat("CORE REGRESSION TESTS PASSED\n")
