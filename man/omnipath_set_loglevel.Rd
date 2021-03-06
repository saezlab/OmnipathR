% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log.R
\name{omnipath_set_loglevel}
\alias{omnipath_set_loglevel}
\title{Sets the log level for the package logger}
\usage{
omnipath_set_loglevel(level, target = "logfile")
}
\arguments{
\item{level}{Character or class `loglevel`. The desired log level.}

\item{target}{Character, either 'logfile' or 'console'}
}
\value{
Returns `NULL`.
}
\description{
Sets the log level for the package logger
}
\examples{
omnipath_set_loglevel(logger::FATAL, target = 'console')

}
