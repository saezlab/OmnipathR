% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keggrest.R
\name{kegg_ddi}
\alias{kegg_ddi}
\title{Find adverse drug-drug interactions in KEGG}
\usage{
kegg_ddi(...)
}
\arguments{
\item{...}{
  Arguments passed on to \code{\link[=kegg_query]{kegg_query}}
  \describe{
    \item{\code{operation}}{Character: one of the KEGG REST API operations.}
  }}
}
\value{
Data frame (tibble) of four columns with names "drug_a", "drug_b",
    "interaction" and "mechanism".
}
\description{
See \url{https://www.kegg.jp/kegg/rest/keggapi.html#ddi} for details.
}
\examples{
kegg_ddi(c("D00564", "D00100", "D00109"))

}
