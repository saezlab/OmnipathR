% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.R
\name{enzsub_graph}
\alias{enzsub_graph}
\title{Enzyme-substrate graph}
\usage{
enzsub_graph(enzsub)
}
\arguments{
\item{enzsub}{Data frame created by \code{\link{enzyme_substrate}}}
}
\value{
An igraph directed graph object.
}
\description{
Transforms the a data frame with enzyme-substrate relationships
(obtained by \code{\link{enzyme_substrate}}) to an igraph graph object.
}
\examples{
enzsub <- enzyme_substrate(resources = c('PhosphoSite', 'SIGNOR'))
enzsub_g <- enzsub_graph(enzsub = enzsub)

}
\seealso{
\itemize{
    \item{\code{\link{enzyme_substrate}}}
    \item{\code{\link{giant_component}}}
    \item{\code{\link{find_all_paths}}}
}
}
