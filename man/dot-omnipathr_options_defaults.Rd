% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/options.R
\docType{data}
\name{.omnipathr_options_defaults}
\alias{.omnipathr_options_defaults}
\title{Default values for the package options}
\format{
An object of class \code{list} of length 35.
}
\usage{
.omnipathr_options_defaults
}
\value{
Nothing, this is not a function but a list.
}
\description{
These options describe the default settings for OmnipathR so you do not
need to pass these parameters at each function call.
Currently the only option useful for the public web service at
omnipathdb.org is ``omnipathr.license``. If you are a for-profit
user set it to ``'commercial'`` to make sure all the data you download
from OmniPath is legally allowed for commercial use. Otherwise just leave
it as it is: ``'academic'``.
If you don't use omnipathdb.org but within your organization you deployed
your own pypath server and want to share data whith a limited availability
to outside users, you may want to use a password. For this you can use
the ``omnipathr.password`` option.
Also if you want the R package to work from another pypath server instead
of omnipathdb.org, you can change the option ``omnipathr.url``.
}
\keyword{datasets}
