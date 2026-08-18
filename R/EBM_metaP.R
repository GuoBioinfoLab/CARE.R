#' Empirical Brown's Method for combining dependent p-values
#'
#' Combines p-values from multiple hallmarks while accounting for their
#' correlation estimated empirically from a reference data matrix, following
#' Poole et al. (2016). If the reference data is unsuitable (too few samples,
#' zero variance, non-finite covariance, ...), it falls back to Fisher's method.
#'
#' @param data_matrix A numeric matrix with one row per p-value (e.g. per
#'   hallmark) and one column per reference sample (e.g. ssGSEA scores across
#'   control + treatment samples).
#' @param p_values A numeric vector of p-values, one per row of
#'   `data_matrix`, in the same order.
#' @param min_reference_samples Minimum number of complete reference columns
#'   required before the empirical covariance estimate is attempted.
#' @param cap_at_fisher If `TRUE`, the Brown degrees of freedom are capped at
#'   Fisher's df (equivalent to scale factor 1) when the estimate would exceed
#'   it.
#' @return A list with components `P_Brown`, `P_Fisher`, `Fisher_statistic`,
#'   `Scale_Factor`, `DF_Brown`, `Brown_Variance`, `Covariance_Sum`,
#'   `Covariance_Matrix`, `Transformed_Data`, `Reference_N` and `Fallback`
#'   (character reason when the Fisher fallback was used, otherwise `NA`).
#' @export
empirical_browns_method <- function(
    data_matrix,
    p_values,
    min_reference_samples = 4,
    cap_at_fisher = TRUE
) {

  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "double"

  p_values <- as.numeric(p_values)

  if (nrow(data_matrix) != length(p_values)) {
    stop(
      "The number of rows in data_matrix must equal ",
      "the number of p-values."
    )
  }

  if (length(p_values) == 0) {
    stop("No p-values were supplied.")
  }

  if (any(!is.finite(p_values))) {
    stop("p_values contain NA, NaN or Inf.")
  }

  if (any(p_values < 0 | p_values > 1)) {
    stop("All p-values must be between 0 and 1.")
  }

  # Only keep samples complete across all hallmarks
  complete_columns <- apply(
    data_matrix,
    2,
    function(x) all(is.finite(x))
  )

  data_matrix <- data_matrix[
    ,
    complete_columns,
    drop = FALSE
  ]

  k <- length(p_values)

  p_values <- pmax(
    pmin(p_values, 1),
    .Machine$double.xmin
  )

  fisher_statistic <- -2 * sum(log(p_values))
  fisher_df <- 2 * k

  fisher_p <- pchisq(
    fisher_statistic,
    df = fisher_df,
    lower.tail = FALSE
  )

  fisher_result <- function(reason) {

    list(
      P_Brown = fisher_p,
      P_Fisher = fisher_p,
      Fisher_statistic = fisher_statistic,
      Scale_Factor = 1,
      DF_Brown = fisher_df,
      Brown_Variance = 4 * k,
      Covariance_Sum = 0,
      Covariance_Matrix = NULL,
      Transformed_Data = NULL,
      Reference_N = ncol(data_matrix),
      Fallback = reason
    )
  }

  if (k == 1) {
    return(fisher_result("Only one P-value"))
  }

  if (ncol(data_matrix) < min_reference_samples) {
    return(
      fisher_result(
        paste0(
          "Fewer than ",
          min_reference_samples,
          " complete reference samples"
        )
      )
    )
  }

  raw_variance <- apply(
    data_matrix,
    1,
    var,
    na.rm = TRUE
  )

  if (
    any(!is.finite(raw_variance)) ||
    any(raw_variance == 0)
  ) {
    return(
      fisher_result(
        "At least one hallmark has zero variance"
      )
    )
  }

  # EBM uses the right-tail empirical cumulative distribution
  empirical_transform <- function(x) {

    n <- length(x)

    right_tail_probability <- (
      n - rank(
        x,
        ties.method = "min",
        na.last = "keep"
      ) + 1
    ) / n

    right_tail_probability <- pmax(
      right_tail_probability,
      1 / n
    )

    -2 * log(right_tail_probability)
  }

  transformed_data <- t(
    apply(
      data_matrix,
      1,
      empirical_transform
    )
  )

  rownames(transformed_data) <- rownames(data_matrix)
  colnames(transformed_data) <- colnames(data_matrix)

  transformed_variance <- apply(
    transformed_data,
    1,
    var,
    na.rm = TRUE
  )

  if (
    any(!is.finite(transformed_variance)) ||
    any(transformed_variance == 0)
  ) {
    return(
      fisher_result(
        "Zero variance after empirical transformation"
      )
    )
  }

  covariance_matrix <- cov(
    t(transformed_data),
    use = "complete.obs"
  )

  covariance_matrix <- as.matrix(
    covariance_matrix
  )

  if (any(!is.finite(covariance_matrix))) {
    return(
      fisher_result(
        "Covariance matrix contains invalid values"
      )
    )
  }

  covariance_sum <- sum(
    covariance_matrix[
      upper.tri(covariance_matrix)
    ]
  )

  expected_fisher <- 2 * k

  brown_variance <- 4 * k +
    2 * covariance_sum

  if (
    !is.finite(brown_variance) ||
    brown_variance <= 0
  ) {
    return(
      fisher_result(
        "Invalid Brown variance"
      )
    )
  }

  scale_factor <- brown_variance /
    (2 * expected_fisher)

  brown_df <- 2 * expected_fisher^2 /
    brown_variance

  if (
    !is.finite(scale_factor) ||
    !is.finite(brown_df) ||
    scale_factor <= 0 ||
    brown_df <= 0
  ) {
    return(
      fisher_result(
        "Invalid Brown parameters"
      )
    )
  }

  if (
    cap_at_fisher &&
    brown_df > fisher_df
  ) {
    brown_df <- fisher_df
    scale_factor <- 1
    brown_variance <- 4 * k
  }

  brown_p <- pchisq(
    fisher_statistic / scale_factor,
    df = brown_df,
    lower.tail = FALSE
  )

  brown_p <- min(
    max(brown_p, 0),
    1
  )

  list(
    P_Brown = brown_p,
    P_Fisher = fisher_p,
    Fisher_statistic = fisher_statistic,
    Scale_Factor = scale_factor,
    DF_Brown = brown_df,
    Brown_Variance = brown_variance,
    Covariance_Sum = covariance_sum,
    Covariance_Matrix = covariance_matrix,
    Transformed_Data = transformed_data,
    Reference_N = ncol(data_matrix),
    Fallback = NA_character_
  )
}
