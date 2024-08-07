% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ramp.R
\name{ramp_table}
\alias{ramp_table}
\title{Return table from RaMP database}
\usage{
ramp_table(name, version = "2.5.4")
}
\arguments{
\item{name}{Character. The name of the RaMP table to fetch.}

\item{version}{Character. The version of RaMP to download.}
}
\value{
Character vector of table names in the RaMP SQLite database.
}
\description{
Return table from RaMP database
}
\examples{
ramp_table('source')

}
\seealso{
\itemize{
    \item{\code{\link{ramp_sqlite}}}
    \item{\code{\link{ramp_tables}}}
}
}
