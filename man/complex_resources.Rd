% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexes.R
\name{complex_resources}
\alias{complex_resources}
\alias{get_complex_resources}
\title{Retrieve a list of complex resources available in Omnipath}
\usage{
complex_resources(dataset = NULL)
}
\arguments{
\item{dataset}{ignored for this query type}
}
\value{
character vector with the names of the databases
}
\description{
Get the names of the resources from \url{https://omnipathdb.org/complexes}
}
\examples{
complex_resources()

}
\seealso{
\itemize{
    \item{\code{\link{resources}}}
    \item{\code{\link{complexes}}}
}
}
