% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactions.R
\name{import_small_molecule_protein_interactions}
\alias{import_small_molecule_protein_interactions}
\title{Interactions from the small molecule-protein dataset of OmniPath}
\usage{
import_small_molecule_protein_interactions(
  resources = NULL,
  organism = "human",
  fields = NULL,
  default_fields = TRUE,
  references_by_resource = TRUE,
  exclude = NULL,
  strict_evidences = FALSE,
  ...
)
}
\arguments{
\item{resources}{interactions not reported in these databases are
removed. See \code{\link{get_interaction_resources}} for more
information.}

\item{organism}{Character or integer: Name or NCBI Taxonomy ID of one or
organisms. The web service currently provides interactions for
human, mouse and rat.}

\item{fields}{Optional fields to be added.}

\item{default_fields}{whether to include the default fields (columns) for
the query type. If FALSE, only the fields defined by the user in the
`fields` argument will be added.}

\item{references_by_resource}{If \code{FALSE}, removes the resource name
prefixes from the references (PubMed IDs); this way the information
which reference comes from which resource will be lost and the PubMed
IDs will be unique.}

\item{exclude}{Character: datasets or resources to exclude.}

\item{strict_evidences}{Logical: restrict the evidences to the queried
datasets and resources. If set to FALSE, the directions and effect signs
and references might be based on other datasets and resources.}

\item{...}{optional additional arguments}
}
\value{
A dataframe of small molecule-protein interactions
}
\description{
Imports the dataset from:
\url{https://omnipathdb.org/interactions?datasets=small_molecule},
which contains small molecule-protein interactions. Small molecules
can be metabolites, intrinsic ligands or drug compounds.
}
\examples{
# What are the targets of aspirin?
interactions <-
    import_small_molecule_protein_interactions(
        sources = 'ASPIRIN'
    )
# The prostaglandin synthases:
interactions

}
\seealso{
\itemize{
    \item{\code{\link{get_interaction_resources}}}
    \item{\code{\link{import_all_interactions}}}
    \item{\code{\link{interaction_graph}}}
    \item{\code{\link{print_interactions}}}
}
}
