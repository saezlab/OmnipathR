% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nichenet.R
\name{nichenet_signaling_network_omnipath}
\alias{nichenet_signaling_network_omnipath}
\title{Builds signaling network for NicheNet using OmniPath}
\usage{
nichenet_signaling_network_omnipath(min_curation_effort = 0, ...)
}
\arguments{
\item{min_curation_effort}{Lower threshold for curation effort}

\item{...}{Passed to \code{\link{import_post_translational_interactions}}}
}
\value{
A network data frame (tibble) with signaling interactions
    suitable for use with NicheNet.
}
\description{
Retrieves network prior knowledge from OmniPath and provides it in
a format suitable for NicheNet.
This method never downloads the `ligrecextra` dataset because the
ligand-receptor interactions are supposed to come from \code{
\link{nichenet_lr_network_omnipath}}.
}
\examples{
# use interactions with at least 2 evidences (reference or database)
op_signaling_network <- nichenet_signaling_network_omnipath(
    min_curation_effort = 2
)

}
\seealso{
\itemize{
    \item{\code{\link{nichenet_signaling_network}}}
}
}
