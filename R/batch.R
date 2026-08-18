#' Run CARE.R analysis for multiple compounds
#'
#' @param expression Genes-by-samples expression matrix/data.frame or a CSV/TSV path.
#' @param metadata Long-format data.frame or CSV/TSV path with one row per
#'   expression sample and columns for sample ID, compound, and group.
#' @param sample_id_col Metadata column containing exact expression matrix
#'   column names.
#' @param compound_col Metadata column identifying compounds.
#' @param group_col Metadata column containing `control` or `treatment`.
#' @param compound_names Optional character vector selecting and ordering compounds.
#' @param gene_set_path Path to the hallmark gene-set JSON file.
#' @return A list with `status`, `summary`, `results`, and `warnings`.
#' @export
batch_analyze <- function(
    expression,
    metadata,
    sample_id_col = "sample_id",
    compound_col = "compound",
    group_col = "group",
    compound_names = NULL,
    gene_set_path = system.file(
      "extdata",
      "cancer_aging_hallmark.json",
      package = "CARE.R"
    )
) {
  expression_matrix <- read_batch_expression(expression)
  metadata_frame <- read_batch_metadata(metadata)

  required_columns <- c(sample_id_col, compound_col, group_col)
  missing_columns <- setdiff(required_columns, names(metadata_frame))
  if (length(missing_columns) > 0) {
    stop(
      "Metadata is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (is.null(colnames(expression_matrix)) ||
      any(is.na(colnames(expression_matrix))) ||
      any(!nzchar(colnames(expression_matrix)))) {
    stop("Expression matrix must have non-empty sample column names")
  }
  if (anyDuplicated(colnames(expression_matrix)) > 0) {
    stop("Expression matrix contains duplicate sample column names")
  }

  sample_ids <- trimws(as.character(metadata_frame[[sample_id_col]]))
  if (any(is.na(sample_ids) | !nzchar(sample_ids))) {
    stop("Metadata column '", sample_id_col, "' contains empty values")
  }
  unmatched <- setdiff(unique(sample_ids), colnames(expression_matrix))
  if (length(unmatched) > 0) {
    stop(
      "Metadata column '",
      sample_id_col,
      "' contains sample_id not found in expression matrix: ",
      paste(unmatched, collapse = ", ")
    )
  }

  groups <- tolower(trimws(as.character(metadata_frame[[group_col]])))
  invalid_groups <- setdiff(unique(groups), c("control", "treatment"))
  if (length(invalid_groups) > 0 || any(is.na(groups))) {
    stop(
      "Metadata column '",
      group_col,
      "' must contain only control or treatment"
    )
  }

  compounds <- as.character(metadata_frame[[compound_col]])
  if (any(is.na(compounds) | !nzchar(trimws(compounds)))) {
    stop("Metadata column '", compound_col, "' contains empty values")
  }
  available_compounds <- unique(compounds)
  if (is.null(compound_names)) {
    selected_compounds <- available_compounds
  } else {
    selected_compounds <- as.character(compound_names)
    missing_compounds <- setdiff(selected_compounds, available_compounds)
    if (length(missing_compounds) > 0) {
      stop(
        "Requested compound(s) are missing from metadata: ",
        paste(missing_compounds, collapse = ", ")
      )
    }
  }

  results <- setNames(vector("list", length(selected_compounds)), selected_compounds)
  summary_rows <- vector("list", length(selected_compounds))
  batch_warnings <- character(0)

  for (index in seq_along(selected_compounds)) {
    compound <- selected_compounds[index]
    rows <- which(compounds == compound)
    compound_metadata <- metadata_frame[rows, , drop = FALSE]
    compound_groups <- groups[rows]
    compound_sample_ids <- sample_ids[rows]
    control_ids <- compound_sample_ids[compound_groups == "control"]
    treatment_ids <- compound_sample_ids[compound_groups == "treatment"]

    if (anyDuplicated(control_ids) > 0 || anyDuplicated(treatment_ids) > 0) {
      error_message <- paste0(
        "Compound '", compound,
        "' has duplicate sample_id values within a group"
      )
      summary_rows[[index]] <- batch_error_summary(compound, error_message)
      results[[index]] <- NULL
      batch_warnings <- c(batch_warnings, error_message)
      next
    }
    if (length(control_ids) < 2L || length(treatment_ids) < 2L) {
      error_message <- paste0(
        "Compound '", compound,
        "' requires at least two control and two treatment samples"
      )
      summary_rows[[index]] <- batch_error_summary(compound, error_message)
      results[[index]] <- NULL
      batch_warnings <- c(batch_warnings, error_message)
      next
    }

    result <- tryCatch(
      analyze_tpm(
        expression_matrix[, control_ids, drop = FALSE],
        expression_matrix[, treatment_ids, drop = FALSE],
        gene_set_path = gene_set_path
      ),
      error = function(error) error
    )
    if (inherits(result, "error")) {
      error_message <- paste0("Compound '", compound, "': ", conditionMessage(result))
      summary_rows[[index]] <- batch_error_summary(compound, error_message)
      results[[index]] <- NULL
      batch_warnings <- c(batch_warnings, error_message)
      next
    }

    results[[index]] <- result
    summary_rows[[index]] <- batch_result_summary(compound, result)
    if (length(result$warnings) > 0) {
      batch_warnings <- c(
        batch_warnings,
        paste0(compound, ": ", result$warnings)
      )
    }
  }

  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  names(results) <- selected_compounds
  successful <- vapply(results, Negate(is.null), logical(1))
  list(
    status = if (all(successful)) "ok" else "partial",
    summary = summary,
    results = results,
    warnings = unique(batch_warnings)
  )
}

read_batch_expression <- function(input) {
  if (is.list(input) && !is.null(input$datapath)) {
    input <- input$datapath
  }
  if (is.matrix(input)) {
    return(input)
  }
  if (is.data.frame(input)) {
    if (ncol(input) < 3) {
      stop("Expression matrix requires a gene column and at least two samples")
    }
    if (is.character(input[[1]]) || is.factor(input[[1]])) {
      genes <- trimws(as.character(input[[1]]))
      values <- input[, -1, drop = FALSE]
      result <- as.matrix(values)
      rownames(result) <- genes
      return(result)
    }
    return(as.matrix(input))
  }
  if (length(input) != 1L || is.na(input) || !nzchar(input)) {
    stop("Expression matrix path is empty")
  }
  if (!file.exists(input)) {
    stop("Expression matrix file does not exist: ", input)
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
    stop("Expression matrix requires a gene column and at least two samples")
  }
  result <- as.matrix(frame[, -1, drop = FALSE])
  rownames(result) <- trimws(as.character(frame[[1]]))
  result
}

read_batch_metadata <- function(input) {
  if (is.data.frame(input)) {
    return(input)
  }
  if (length(input) != 1L || is.na(input) || !nzchar(input)) {
    stop("Metadata path is empty")
  }
  if (!file.exists(input)) {
    stop("Metadata file does not exist: ", input)
  }
  first_line <- readLines(input, n = 1, warn = FALSE)
  separator <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
  read.table(
    input,
    header = TRUE,
    sep = separator,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
}

batch_result_summary <- function(compound, result) {
  values <- setNames(
    vapply(result$hallmarks, function(x) x$delta, numeric(1)),
    vapply(result$hallmarks, function(x) x$name, character(1))
  )
  p_values <- setNames(
    vapply(result$hallmarks, function(x) x$p_value, numeric(1)),
    vapply(result$hallmarks, function(x) x$name, character(1))
  )
  data.frame(
    compound = compound,
    status = result$status,
    direction_consistent = result$direction_consistent,
    pebm = result$pebm %||% NA_real_,
    fisher_p = result$fisher_p %||% NA_real_,
    control_samples = result$input_summary$control_samples,
    treatment_samples = result$input_summary$treatment_samples,
    genomic_instability_delta = values[["genomic_instability"]] %||% NA_real_,
    epigenetic_alterations_delta = values[["epigenetic_alterations"]] %||% NA_real_,
    chronic_inflammation_delta = values[["chronic_inflammation"]] %||% NA_real_,
    genomic_instability_p = p_values[["genomic_instability"]] %||% NA_real_,
    epigenetic_alterations_p = p_values[["epigenetic_alterations"]] %||% NA_real_,
    chronic_inflammation_p = p_values[["chronic_inflammation"]] %||% NA_real_,
    warnings = paste(result$warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

batch_error_summary <- function(compound, message) {
  data.frame(
    compound = compound,
    status = "error",
    direction_consistent = NA,
    pebm = NA_real_,
    fisher_p = NA_real_,
    control_samples = NA_integer_,
    treatment_samples = NA_integer_,
    genomic_instability_delta = NA_real_,
    epigenetic_alterations_delta = NA_real_,
    chronic_inflammation_delta = NA_real_,
    genomic_instability_p = NA_real_,
    epigenetic_alterations_p = NA_real_,
    chronic_inflammation_p = NA_real_,
    warnings = message,
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0 || is.na(left)) right else left
}
