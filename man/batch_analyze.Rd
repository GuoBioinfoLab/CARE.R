\name{batch_analyze}
\alias{batch_analyze}
\title{Run CARE analysis for multiple compounds}
\usage{
batch_analyze(expression, metadata, sample_id_col = "sample_id",
  compound_col = "compound", group_col = "group", compound_names = NULL,
  gene_set_path = system.file("extdata", "cancer_aging_hallmark.json",
    package = "CARE.R"))
}
\arguments{
  \item{expression}{A genes-by-samples matrix or data frame, or a CSV/TSV
  path. Matrix row names identify genes and column names identify samples.}
  \item{metadata}{A data frame or CSV/TSV path. It must contain columns named
  by \code{sample_id_col}, \code{compound_col}, and \code{group_col}.}
  \item{sample_id_col}{Metadata column whose values exactly match expression
  matrix column names.}
  \item{compound_col}{Metadata column identifying compounds.}
  \item{group_col}{Metadata column containing \code{control} or
  \code{treatment}.}
  \item{compound_names}{Optional character vector selecting and ordering
  compounds.}
  \item{gene_set_path}{Path to the hallmark gene-set JSON file.}
}
\value{
A list with \code{status}, \code{summary}, \code{results}, and \code{warnings}.
Each successful entry in \code{results} is the same result object returned by
\code{analyze_tpm}. A compound-level failure is reported in \code{summary}
without stopping other compounds.
}
\description{
Runs the CARE ssGSEA and Empirical Brown's Method analysis independently for
each compound in a long-format metadata table.
}
\examples{
\dontrun{
expression <- matrix(runif(300), nrow = 100,
  dimnames = list(paste0("gene", 1:100), paste0("sample", 1:3)))
metadata <- data.frame(
  sample_id = colnames(expression),
  compound = "demo",
  group = c("control", "control", "treatment")
)
batch_analyze(expression, metadata)
}
}
