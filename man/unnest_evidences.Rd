% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evidences.R
\name{unnest_evidences}
\alias{unnest_evidences}
\title{Separate evidences by direction and effect sign}
\usage{
unnest_evidences(data, longer = FALSE, .keep = FALSE)
}
\arguments{
\item{data}{An interaction data frame with "evidences" column.}

\item{longer}{Logical: If TRUE, the "evidences" column is split into rows.}

\item{.keep}{Logical: keep the "evidences" column. When unnesting to longer
data frame, the "evidences" column will contain the unnested evidences,
while the original column will be retained under the "all_evidences"
name (if `.keep = TRUE`).}
}
\value{
The data frame with new columns or new rows by direction and sign.
}
\description{
Separate evidences by direction and effect sign
}
\examples{
\dontrun{
op <- omnipath_interactions(fields = "evidences")
op <- unnest_evidences(op)
colnames(op)
}

}
\seealso{
\itemize{
    \item{\code{\link{only_from}}}
    \item{\code{\link{filter_evidences}}}
    \item{\code{\link{from_evidences}}}
}
}
