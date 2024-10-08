% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intercell.R
\name{simplify_intercell_network}
\alias{simplify_intercell_network}
\title{Simplify an intercell network}
\usage{
simplify_intercell_network(network, ...)
}
\arguments{
\item{network}{An intercell network data frame, as provided by
\code{\link{intercell_network}}.}

\item{...}{Optional, further columns to select.}
}
\value{
An intercell network data frame with some columns removed.
}
\description{
The intercellular communication network data frames, created by
\code{\link{intercell_network}}, are combinations of a network data
frame with two copies of the intercell annotation data frames, all of them
already having quite some columns. Here we keep only the names of the
interacting pair, their intercellular communication roles, and the minimal
information of the origin of both the interaction and the annotations.
Optionally further columns can be selected.
}
\examples{
icn <- intercell_network()
icn_s <- simplify_intercell_network(icn)

}
\seealso{
\itemize{
    \item{\code{\link{intercell_network}}}
    \item{\code{\link{filter_intercell_network}}}
    \item{\code{\link{unique_intercell_network}}}
    \item{\code{\link{intercell}}}
    \item{\code{\link{intercell_categories}}}
    \item{\code{\link{intercell_generic_categories}}}
    \item{\code{\link{intercell_summary}}}
}
}
