% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/options.R
\name{omnipath_get_config_path}
\alias{omnipath_get_config_path}
\title{Current config file path}
\usage{
omnipath_get_config_path(user = FALSE)
}
\arguments{
\item{user}{Logical: prioritize the user level config even if a config in
the current working directory is available.}
}
\value{
Character: path to the config file.
}
\description{
Current config file path
}
\examples{
omnipath_get_config_path()

}
