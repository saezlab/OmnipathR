% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nichenet.R
\name{nichenet_signaling_network_vinayagam}
\alias{nichenet_signaling_network_vinayagam}
\title{NicheNet signaling network from Vinayagam}
\usage{
nichenet_signaling_network_vinayagam(...)
}
\arguments{
\item{...}{Ignored.}
}
\value{
A network data frame (tibble) with signaling interactions
    suitable for use with NicheNet.
}
\description{
Builds signaling network prior knowledge for NicheNet using Vinayagam 2011
Supplementary Table S6. Find out more at
\url{https://doi.org/10.1126/scisignal.2001699}.
}
\examples{
vi_signaling_network <- nichenet_signaling_network_vinayagam()

}
