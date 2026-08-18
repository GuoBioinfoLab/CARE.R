# Smoke test for CARE.R analysis (no external test framework required).
# Run: Rscript tests/smoke.R
library(CARE.R)

demo_dir <- file.path("..", "public", "examples")
control <- file.path(demo_dir, "itcm_control_tpm.csv")
treatment <- file.path(demo_dir, "itcm_treatment_tpm.csv")
if (!file.exists(control) || !file.exists(treatment)) {
  cat("SMOKE TEST SKIPPED: demo matrices not available\n")
  quit(status = 0)
}

res <- analyze_tpm(control, treatment)

stopifnot(
  is.list(res),
  res$status %in% c("ok", "direction_not_consistent"),
  is.logical(res$direction_consistent),
  length(res$hallmarks) == 3,
  all(vapply(res$hallmarks, function(h) {
    is.numeric(h$delta) && is.numeric(h$p_value) && is.integer(h$gene_coverage)
  }, logical(1)))
)

# If direction is consistent, EBM fields must be present and valid
if (res$direction_consistent) {
  stopifnot(
    is.numeric(res$pebm), res$pebm >= 0, res$pebm <= 1,
    is.numeric(res$fisher_p),
    is.numeric(res$fisher_statistic),
    is.numeric(res$brown_df), res$brown_df > 0,
    is.numeric(res$brown_scale), res$brown_scale > 0,
    length(res$distribution) == 180
  )
}

cat("status:", res$status, "\n")
cat("direction_consistent:", res$direction_consistent, "\n")
if (res$direction_consistent) {
  cat("pebm:", res$pebm, "\n")
  cat("fisher_p:", res$fisher_p, "\n")
  cat("brown_df:", res$brown_df, " brown_scale:", res$brown_scale, "\n")
}
cat("input_summary:", paste(names(res$input_summary), unlist(res$input_summary), sep = "=", collapse = ", "), "\n")
for (h in res$hallmarks) {
  cat(sprintf("  %-22s delta=%.6f p=%.4g coverage=%d\n", h$label, h$delta, h$p_value, h$gene_coverage))
}
cat("SMOKE TEST PASSED\n")
