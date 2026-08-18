# Compare CARE.R package output against the original backend analyze_tpm.R JSON.
# Usage: Rscript tests/compare_backend.R <original.json>
# Run from CARE_R/ with ../CARE_R_lib on .libPaths().
library(CARE.R)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("BACKEND COMPARISON SKIPPED: expected original JSON path\n")
  quit(status = 0)
}
original <- jsonlite::fromJSON(args[1], simplifyVector = FALSE)

res <- analyze_tpm(
  file.path("..", "public", "examples", "itcm_control_tpm.csv"),
  file.path("..", "public", "examples", "itcm_treatment_tpm.csv")
)

fields <- c("status", "direction_consistent", "pebm", "fisher_p",
            "fisher_statistic", "brown_df", "brown_scale", "fallback")
for (f in fields) {
  a <- original[[f]]; b <- res[[f]]
  if (is.null(a) || is.null(b)) { cat(sprintf("%-18s present-in-both=%s\n", f, !is.null(a) && !is.null(b))); next }
  ok <- if (is.numeric(a)) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-12)) else identical(a, b)
  cat(sprintf("%-18s equal=%-5s original=%s package=%s\n", f, ok,
              if (is.numeric(a)) format(as.numeric(a), digits = 12) else paste(a, collapse = ","),
              if (is.numeric(b)) format(as.numeric(b), digits = 12) else paste(b, collapse = ",")))
  if (!ok) stop("Mismatch in field: ", f)
}

# Per-hallmark comparison
for (i in seq_along(res$hallmarks)) {
  oh <- original$hallmarks[[i]]; ph <- res$hallmarks[[i]]
  stopifnot(identical(oh$name, ph$name),
            isTRUE(all.equal(oh$delta, ph$delta, tolerance = 1e-12)),
            isTRUE(all.equal(oh$p_value, ph$p_value, tolerance = 1e-12)),
            oh$gene_coverage == ph$gene_coverage)
  oscores_c <- vapply(oh$control, function(s) s$score, numeric(1))
  pscores_c <- vapply(ph$control, function(s) s$score, numeric(1))
  oscores_t <- vapply(oh$treatment, function(s) s$score, numeric(1))
  pscores_t <- vapply(ph$treatment, function(s) s$score, numeric(1))
  stopifnot(isTRUE(all.equal(oscores_c, pscores_c, tolerance = 1e-10)),
            isTRUE(all.equal(oscores_t, pscores_t, tolerance = 1e-10)))
  cat(sprintf("hallmark %-22s delta/p/coverage identical\n", oh$name))
}

# Input summary
stopifnot(original$input_summary$common_genes == res$input_summary$common_genes,
          original$input_summary$control_samples == res$input_summary$control_samples,
          original$input_summary$treatment_samples == res$input_summary$treatment_samples)

# Distribution curve (density points)
if (res$direction_consistent) {
  od <- vapply(original$distribution, function(p) p$density, numeric(1))
  pd <- vapply(res$distribution, function(p) p$density, numeric(1))
  stopifnot(length(od) == length(pd), isTRUE(all.equal(od, pd, tolerance = 1e-12)))
  cat("distribution curve identical (", length(pd), " points )\n", sep = "")
}

cat("BACKEND COMPARISON PASSED\n")
