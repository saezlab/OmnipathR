% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intercell.R
\name{intercell_resources}
\alias{intercell_resources}
\alias{get_intercell_resources}
\title{Retrieves a list of intercellular communication resources available in
OmniPath}
\usage{
intercell_resources(dataset = NULL)
}
\arguments{
\item{dataset}{ignored at this query type}
}
\value{
character vector with the names of the databases
}
\description{
Retrieves a list of the databases from
\url{https://omnipathdb.org/intercell}.
}
\examples{
intercell_resources()

}
\seealso{
\itemize{
    \item{\code{\link{resources}}}
    \item{\code{\link{intercell}}}
    \item{\code{\link{filter_intercell}}}
    \item{\code{\link{intercell_categories}}}
    \item{\code{\link{intercell_generic_categories}}}
    \item{\code{\link{intercell_summary}}}
    \item{\code{\link{intercell_network}}}
}
}
