#' Analyze control vs treatment TPM matrices with cancer/aging hallmarks
#'
#' This is the packaged version of the CARE backend `analyze_tpm.R` script.
#' It performs, in order:
#' \enumerate{
#'   \item Reads and validates two genes-by-samples TPM matrices (CSV or TSV,
#'         first column = gene identifier, at least two sample columns).
#'   \item Restricts both matrices to their common genes (>= 20 required).
#'   \item Runs single-sample GSEA (ssGSEA, \code{GSVA::gsva}) on three
#'         cancer-aging hallmarks: genomic instability, epigenetic alterations
#'         and chronic inflammation.
#'   \item Computes per-hallmark median change (treatment - control), a
#'         one-sided Wilcoxon test (\code{alternative = "less"}), and the
#'         direction consistency of the changes.
#'   \item If all three hallmarks decrease under treatment, combines the three
#'         p-values with Empirical Brown's Method (\code{\link{empirical_browns_method}})
#'         into the EBM evidence score \code{pebm}, together with the scaled
#'         chi-square distribution used for plotting.
#' }
#'
#' The return value mirrors the JSON structure produced by the backend
#' `/api/analysis` endpoint.
#'
#' @param control_path Path to the control TPM matrix file.
#' @param treatment_path Path to the treatment TPM matrix file.
#' @param gene_set_path Path to the hallmark gene-set JSON file. Defaults to the
#'   gene sets shipped with the package
#'   (\code{system.file("extdata", "cancer_aging_hallmark.json", package = "CARE.R")}).
#' @return A named list with components \code{status},
#'   \code{direction_consistent}, \code{pebm}, \code{fisher_p},
#'   \code{fisher_statistic}, \code{brown_df}, \code{brown_scale},
#'   \code{fallback}, \code{input_summary}, \code{hallmarks} and
#'   \code{distribution}, mirroring the backend JSON output.
#' @export
#' @examples
#' \dontrun{
#' res <- analyze_tpm("control.csv", "treatment.csv")
#' res$pebm
#' }
analyze_tpm <- function(
    control_path,
    treatment_path,
    gene_set_path = system.file("extdata", "cancer_aging_hallmark.json", package = "CARE.R")
) {

  hallmarks <- c(
    "genomic_instability",
    "epigenetic_alterations",
    "chronic_inflammation"
  )

  read_tpm <- function(input, role) {
    if (is.list(input) && !is.null(input$datapath)) {
      input <- input$datapath
    }

    if (is.matrix(input) || is.data.frame(input)) {
      if (is.data.frame(input) && ncol(input) >= 3 &&
          (is.character(input[[1]]) || is.factor(input[[1]]))) {
        genes <- trimws(as.character(input[[1]]))
        sample_names <- names(input)[-1]
        sample_values <- input[, -1, drop = FALSE]
      } else {
        genes <- rownames(input)
        if (is.null(genes) || any(!nzchar(genes))) {
          stop(role, " matrix must have non-empty gene IDs in row names")
        }
        sample_names <- colnames(input)
        sample_values <- as.data.frame(input, check.names = FALSE)
      }
    } else {
      if (length(input) != 1L || is.na(input) || !nzchar(input)) {
        stop(role, " TPM path is empty")
      }
      if (!file.exists(input)) {
        stop(role, " TPM file does not exist: ", input)
      }
      if (file.access(input, 4) != 0) {
        stop(role, " TPM file is not readable: ", input)
      }
      first_line <- readLines(input, n = 1, warn = FALSE)
      separator <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
      frame <- read.table(
        input,
        header = TRUE,
        sep = separator,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "\"",
        comment.char = ""
      )
      if (ncol(frame) < 3) {
        stop("Each TPM matrix requires a gene column and at least two samples")
      }
      genes <- trimws(as.character(frame[[1]]))
      sample_names <- names(frame)[-1]
      sample_values <- frame[, -1, drop = FALSE]
    }

    if (length(sample_names) < 2) {
      stop(role, " TPM requires at least two sample columns")
    }
    if (any(is.na(sample_names) | !nzchar(sample_names))) {
      stop(role, " TPM contains an empty sample name")
    }
    if (anyDuplicated(sample_names) > 0) {
      stop(role, " TPM contains duplicate sample names")
    }

    parsed_columns <- lapply(seq_along(sample_names), function(index) {
      values <- suppressWarnings(as.numeric(as.character(sample_values[[index]])))
      if (anyNA(values) || any(!is.finite(values))) {
        stop(
          role,
          " TPM column '",
          sample_names[index],
          "' contains non-numeric or non-finite values"
        )
      }
      if (any(values < 0)) {
        stop(
          role,
          " TPM column '",
          sample_names[index],
          "' contains negative values"
        )
      }
      values
    })
    matrix <- do.call(cbind, parsed_columns)
    colnames(matrix) <- sample_names
    keep <- !is.na(genes) & nzchar(genes)
    genes <- genes[keep]
    matrix <- matrix[keep, , drop = FALSE]

    rownames(matrix) <- genes
    if (anyDuplicated(rownames(matrix)) > 0) {
      gene_n <- table(rownames(matrix))
      matrix <- rowsum(matrix, group = rownames(matrix), reorder = FALSE)
      matrix <- sweep(
        matrix,
        MARGIN = 1,
        STATS = as.numeric(gene_n[rownames(matrix)]),
        FUN = "/"
      )
    }
    matrix
  }

  control_tpm <- read_tpm(control_path, "Control")
  treatment_tpm <- read_tpm(treatment_path, "Treatment")
  diagnostics <- character(0)
  if (ncol(control_tpm) < 4) {
    diagnostics <- c(
      diagnostics,
      paste0(
        "Small sample size: Control has ",
        ncol(control_tpm),
        " samples; results may be unstable."
      )
    )
  }
  if (ncol(treatment_tpm) < 4) {
    diagnostics <- c(
      diagnostics,
      paste0(
        "Small sample size: Treatment has ",
        ncol(treatment_tpm),
        " samples; results may be unstable."
      )
    )
  }
  common_genes <- intersect(rownames(control_tpm), rownames(treatment_tpm))
  if (length(common_genes) < 20) {
    stop("Control and Treatment matrices share fewer than 20 genes")
  }

  control_tpm <- control_tpm[common_genes, , drop = FALSE]
  treatment_tpm <- treatment_tpm[common_genes, , drop = FALSE]
  combined_tpm <- cbind(treatment_tpm, control_tpm)

  if (length(gene_set_path) != 1L || is.na(gene_set_path) || !nzchar(gene_set_path)) {
    stop("Gene-set JSON path is empty")
  }
  if (!file.exists(gene_set_path)) {
    stop("Gene-set JSON file does not exist: ", gene_set_path)
  }
  raw_gene_sets <- tryCatch(
    jsonlite::fromJSON(gene_set_path, simplifyVector = FALSE),
    error = function(error) {
      stop("Gene-set JSON could not be parsed: ", conditionMessage(error))
    }
  )
  gene_set_names <- c(
    "genomic instability",
    "epigenetic alterations",
    "chronic inflammation"
  )
  if (!is.list(raw_gene_sets) || is.null(names(raw_gene_sets))) {
    stop("Gene-set JSON must be a named object")
  }
  missing_gene_sets <- setdiff(gene_set_names, names(raw_gene_sets))
  if (length(missing_gene_sets) > 0) {
    stop(
      "Gene-set JSON is missing: ",
      paste(missing_gene_sets, collapse = ", ")
    )
  }
  gene_sets <- raw_gene_sets[gene_set_names]
  gene_sets <- lapply(gene_sets, function(x) unique(as.character(x)))
  names(gene_sets) <- hallmarks
  coverage <- vapply(gene_sets, function(x) length(intersect(x, common_genes)), integer(1))
  if (any(coverage < 3)) {
    stop(
      "Insufficient hallmark coverage: ",
      paste(names(coverage), coverage, sep = "=", collapse = ", ")
    )
  }

  ssgsea_parameters <- GSVA::ssgseaParam(
    exprData = combined_tpm,
    geneSets = gene_sets,
    normalize = FALSE
  )
  gsva_warnings <- character(0)
  scores <- withCallingHandlers(
    GSVA::gsva(ssgsea_parameters, verbose = FALSE),
    warning = function(warning) {
      gsva_warnings <<- c(gsva_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  diagnostics <- unique(c(diagnostics, gsva_warnings))
  scores <- scores[hallmarks, , drop = FALSE]
  treatment_scores <- scores[, seq_len(ncol(treatment_tpm)), drop = FALSE]
  control_scores <- scores[
    ,
    ncol(treatment_tpm) + seq_len(ncol(control_tpm)),
    drop = FALSE
  ]

  changes <- apply(treatment_scores, 1, median) -
    apply(control_scores, 1, median)
  direction_consistent <- all(is.finite(changes) & changes < 0)
  p_values <- vapply(
    seq_along(hallmarks),
    function(index) {
      wilcox.test(
        treatment_scores[index, ],
        control_scores[index, ],
        alternative = "less",
        exact = TRUE
      )$p.value
    },
    numeric(1)
  )

  hallmark_results <- lapply(seq_along(hallmarks), function(index) {
    list(
      name = hallmarks[index],
      label = c(
        "Genomic instability",
        "Epigenetic alterations",
        "Chronic inflammation"
      )[index],
      gene_coverage = unname(coverage[index]),
      control = lapply(seq_len(ncol(control_scores)), function(sample_index) {
        list(
          sample_id = colnames(control_scores)[sample_index],
          score = control_scores[index, sample_index]
        )
      }),
      treatment = lapply(seq_len(ncol(treatment_scores)), function(sample_index) {
        list(
          sample_id = colnames(treatment_scores)[sample_index],
          score = treatment_scores[index, sample_index]
        )
      }),
      delta = unname(changes[index]),
      p_value = unname(p_values[index])
    )
  })

  if (direction_consistent) {
    brown <- empirical_browns_method(
      data_matrix = scores,
      p_values = p_values,
      min_reference_samples = 4,
      cap_at_fisher = TRUE
    )
    observed <- brown$Fisher_statistic / brown$Scale_Factor
    x_max <- max(qchisq(0.999, brown$DF_Brown), observed * 1.15)
    x <- seq(0, x_max, length.out = 180)
    distribution <- lapply(x, function(value) {
      list(
        x = value,
        density = dchisq(value, brown$DF_Brown),
        in_right_tail = value >= observed
      )
    })
    result <- list(
      status = "ok",
      direction_consistent = TRUE,
      pebm = brown$P_Brown,
      fisher_p = brown$P_Fisher,
      fisher_statistic = brown$Fisher_statistic,
      brown_df = brown$DF_Brown,
      brown_scale = brown$Scale_Factor,
      fallback = brown$Fallback,
      warnings = diagnostics,
      input_summary = list(
        common_genes = length(common_genes),
        control_samples = ncol(control_tpm),
        treatment_samples = ncol(treatment_tpm),
        hallmark_coverage = as.list(coverage)
      ),
      hallmarks = hallmark_results,
      distribution = distribution
    )
  } else {
    result <- list(
      status = "direction_not_consistent",
      direction_consistent = FALSE,
      input_summary = list(
        common_genes = length(common_genes),
        control_samples = ncol(control_tpm),
        treatment_samples = ncol(treatment_tpm),
        hallmark_coverage = as.list(coverage)
      ),
      hallmarks = hallmark_results,
      distribution = list(),
      fallback = "not_tested",
      warnings = diagnostics
    )
  }

  result
}

#' Run the analysis on uploaded file objects (e.g. from a Shiny app)
#'
#' Convenience wrapper around \code{\link{analyze_tpm}} that accepts file
#' objects as returned by Shiny's \code{fileInput} (\code{list(datapath = ...)})
#' or by R's \code{file.choose()}.
#'
#' @param control File object (or path) for the control matrix.
#' @param treatment File object (or path) for the treatment matrix.
#' @param ... Further arguments passed to \code{\link{analyze_tpm}}.
#' @return The same result list as \code{\link{analyze_tpm}}.
#' @export
run_analysis <- function(control, treatment, ...) {
  control_path <- if (is.list(control)) control$datapath else control
  treatment_path <- if (is.list(treatment)) treatment$datapath else treatment
  analyze_tpm(control_path, treatment_path, ...)
}
