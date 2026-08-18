#' Plot hallmark ssGSEA boxplots
#'
#' @param result A result returned by [analyze_tpm()].
#' @param file Optional `.png`, `.pdf`, or `.svg` output path.
#' @param width Output width in inches.
#' @param height Output height in inches.
#' @return `result`, invisibly.
#' @export
plot_hallmarks <- function(result, file = NULL, width = 11, height = 4.5) {
  validate_analysis_result(result)
  with_plot_file(file, width, height, function() draw_hallmarks(result))
  invisible(result)
}

#' Plot the combined EBM evidence distribution
#'
#' @param result A result returned by [analyze_tpm()].
#' @param file Optional `.png`, `.pdf`, or `.svg` output path.
#' @param width Output width in inches.
#' @param height Output height in inches.
#' @return `result`, invisibly.
#' @export
plot_ebm <- function(result, file = NULL, width = 7, height = 4.5) {
  validate_analysis_result(result)
  with_plot_file(file, width, height, function() draw_ebm(result))
  invisible(result)
}

#' Plot the complete CARE.R analysis
#'
#' @param result A result returned by [analyze_tpm()].
#' @param file Optional `.png`, `.pdf`, or `.svg` output path.
#' @param width Output width in inches.
#' @param height Output height in inches.
#' @return `result`, invisibly.
#' @export
plot_analysis <- function(result, file = NULL, width = 11, height = 8.5) {
  validate_analysis_result(result)
  with_plot_file(file, width, height, function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::layout(
      matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
      heights = c(1.2, 1)
    )
    draw_hallmarks(result, set_layout = FALSE)
    draw_ebm(result)
  })
  invisible(result)
}

validate_analysis_result <- function(result) {
  if (!is.list(result) || !is.list(result$hallmarks)) {
    stop("result must be the list returned by analyze_tpm()")
  }
  if (length(result$hallmarks) == 0) {
    stop("result contains no hallmark results")
  }
  invisible(result)
}

with_plot_file <- function(file, width, height, plotter) {
  if (is.null(file)) {
    plotter()
    return(invisible(NULL))
  }
  if (length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("Plot output path must be a non-empty string")
  }
  extension <- tolower(tools::file_ext(file))
  if (extension == "png") {
    grDevices::png(file, width = width, height = height, units = "in", res = 150)
  } else if (extension == "pdf") {
    grDevices::pdf(file, width = width, height = height)
  } else if (extension == "svg") {
    grDevices::svg(file, width = width, height = height)
  } else {
    stop("Unsupported plot extension; use .png, .pdf, or .svg")
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  plotter()
  invisible(NULL)
}

format_plot_value <- function(value) {
    if (length(value) == 0 || !is.finite(value)) return("-")
  if (value != 0 && abs(value) < 1e-3) {
    return(formatC(value, format = "e", digits = 2))
  }
  formatC(value, format = "f", digits = 2, big.mark = ",")
}

format_plot_p <- function(value) {
  if (length(value) == 0 || !is.finite(value)) return("P = -")
  paste0("P = ", format_plot_value(value))
}

draw_hallmarks <- function(result, set_layout = TRUE) {
  if (set_layout) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(1, length(result$hallmarks)))
  }
  graphics::par(mar = c(5, 4, 4, 1) + 0.1)
  control_color <- "#3459E6"
  treatment_color <- "#25A18E"

  for (index in seq_along(result$hallmarks)) {
    hallmark <- result$hallmarks[[index]]
    control <- vapply(hallmark$control, function(x) x$score, numeric(1))
    treatment <- vapply(hallmark$treatment, function(x) x$score, numeric(1))
    values <- c(control, treatment)
    finite_values <- values[is.finite(values)]
    if (length(finite_values) == 0) {
      graphics::plot.new()
      graphics::title(main = hallmark$label)
      graphics::text(0.5, 0.5, "No finite scores")
      next
    }
    value_range <- range(finite_values)
    padding <- if (diff(value_range) > 0) diff(value_range) * 0.12 else max(abs(value_range[1]) * 0.05, 1)
    y_limits <- value_range + c(-padding, padding)
    graphics::boxplot(
      list(control, treatment),
      names = c("Control", "Treatment"),
      col = grDevices::adjustcolor(c(control_color, treatment_color), alpha.f = 0.22),
      border = c(control_color, treatment_color),
      outline = FALSE,
      ylim = y_limits,
      ylab = "ssGSEA score",
      main = hallmark$label,
      las = 1
    )
    graphics::points(
      jitter(rep(1, length(control)), amount = 0.06),
      control,
      pch = 21,
      bg = control_color,
      col = "white",
      cex = 1.1
    )
    graphics::points(
      jitter(rep(2, length(treatment)), amount = 0.06),
      treatment,
      pch = 21,
      bg = treatment_color,
      col = "white",
      cex = 1.1
    )
    graphics::mtext(
      paste0("Delta ", format_plot_value(hallmark$delta), "   ", format_plot_p(hallmark$p_value)),
      side = 3,
      line = 0.2,
      cex = 0.78
    )
    graphics::mtext(
      paste0("Gene coverage: ", hallmark$gene_coverage %||% "-"),
      side = 1,
      line = 3.8,
      cex = 0.72,
      col = "#5F6B7A"
    )
  }
}

draw_ebm <- function(result) {
  distribution <- result$distribution
  if (!isTRUE(result$direction_consistent) || length(distribution) == 0) {
    graphics::plot.new()
    graphics::text(
      0.5,
      0.58,
      "EBM was not tested",
      cex = 1.1,
      col = "#3F4A59"
    )
    graphics::text(
      0.5,
      0.42,
      "The three hallmark changes are not all in the expected direction.",
      cex = 0.85,
      col = "#5F6B7A"
    )
    return(invisible(NULL))
  }

  x <- vapply(distribution, function(point) point$x, numeric(1))
  density <- vapply(distribution, function(point) point$density, numeric(1))
  observed <- result$fisher_statistic / result$brown_scale
  tail <- x >= observed
  y_max <- max(density, na.rm = TRUE)
  graphics::plot(
    x,
    density,
    type = "l",
    lwd = 2,
    col = "#3459E6",
    xlab = "Scaled Fisher statistic",
    ylab = "Density",
    main = "Combined EBM evidence",
    ylim = c(0, y_max * 1.12)
  )
  if (any(tail)) {
    graphics::polygon(
      c(x[tail], max(x[tail])),
      c(density[tail], 0),
      col = grDevices::adjustcolor("#25A18E", alpha.f = 0.28),
      border = NA
    )
  }
  graphics::abline(v = observed, col = "#25A18E", lwd = 2, lty = 2)
  graphics::legend(
    "topright",
    legend = c(
      paste0("Observed = ", format_plot_value(observed)),
      paste0("P EBM = ", format_plot_value(result$pebm))
    ),
    bty = "n",
    text.col = "#3F4A59"
  )
}

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0 || is.na(left)) right else left
}
