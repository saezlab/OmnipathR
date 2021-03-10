% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nichenet.R
\name{nichenet_gr_network_omnipath}
\alias{nichenet_gr_network_omnipath}
\title{Builds gene regulatory network for NicheNet using OmniPath}
\usage{
nichenet_gr_network_omnipath(min_curation_effort = 0, ...)
}
\arguments{
\item{min_curation_effort}{Lower threshold for curation effort}

\item{...}{Passed to \code{\link{import_transcriptional_interactions}}}
}
\description{
Retrieves network prior knowledge from OmniPath and provides it in
a format suitable for NicheNet.
This method never downloads the `ligrecextra` dataset because the
ligand-receptor interactions are supposed to come from \code{
\link{nichenet_lr_network_omnipath}}.
}
\examples{
# use interactions up to confidence level "C" from DoRothEA:
op_gr_network <- nichenet_gr_network_omnipath(
    dorothea_levels = c('A', 'B', 'C')
)

}
\seealso{
\code{\link{nichenet_gr_network_evex},
    \link{nichenet_gr_network_harmonizome},
    \link{nichenet_gr_network_htridb},
    \link{nichenet_gr_network_omnipath},
    \link{nichenet_gr_network_pathwaycommons},
    \link{nichenet_gr_network_regnetwork},
    \link{nichenet_gr_network_remap},
    \link{nichenet_gr_network_trrust}}
}