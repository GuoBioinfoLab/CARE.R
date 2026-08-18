library(CARE.R)

control <- file.path("..", "public", "examples", "itcm_control_tpm.csv")
treatment <- file.path("..", "public", "examples", "itcm_treatment_tpm.csv")
gene_sets <- system.file("extdata", "cancer_aging_hallmark.json", package = "CARE.R")
if (!file.exists(control) || !file.exists(treatment)) {
  cat("VISUALIZATION TESTS SKIPPED: demo matrices not available\n")
  quit(status = 0)
}
result <- analyze_tpm(control, treatment, gene_sets)

stopifnot(identical(invisible(plot_hallmarks(result)), result))
stopifnot(identical(invisible(plot_ebm(result)), result))
stopifnot(identical(invisible(plot_analysis(result)), result))

outputs <- file.path(tempdir(), c("hallmarks.png", "ebm.pdf", "analysis.svg"))
plot_hallmarks(result, outputs[1])
plot_ebm(result, outputs[2])
plot_analysis(result, outputs[3])
stopifnot(all(file.exists(outputs), file.info(outputs)$size > 0))

swapped <- analyze_tpm(treatment, control, gene_sets)
plot_ebm(swapped)
cat("VISUALIZATION TESTS PASSED\n")
