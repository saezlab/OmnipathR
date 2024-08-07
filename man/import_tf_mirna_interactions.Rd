% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactions.R
\name{import_tf_mirna_interactions}
\alias{import_tf_mirna_interactions}
\title{Imports interactions from the TF-miRNA dataset of OmniPath}
\usage{
import_tf_mirna_interactions(
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
removed. See \code{\link{get_interaction_resources}} for more information.}

\item{organism}{Character or integer: Name or NCBI Taxonomy ID of one or
organisms. The web service currently provides interactions for
human, mouse and rat.}

\item{fields}{The user can define here the fields to be added. If used, set
the next argument, `default_fields`, to FALSE.}

\item{default_fields}{whether to include the default fields (columns) for
the query type. If FALSE, only the fields defined by the user in the
`fields` argument will be added.}

\item{references_by_resource}{if FALSE, removes the resource name prefixes
from the references (PubMed IDs); this way the information which reference
comes from which resource will be lost and the PubMed IDs will be unique.}

\item{exclude}{Character: datasets or resources to exclude.}

\item{strict_evidences}{Logical: restrict the evidences to the queried
datasets and resources. If set to FALSE, the directions and effect signs
and references might be based on other datasets and resources.}

\item{...}{optional additional arguments}
}
\value{
A dataframe containing TF-miRNA interactions
}
\description{
Imports the dataset from:
\url{https://omnipathdb.org/interactions?datasets=tf_mirna},
which contains transcription factor-miRNA gene interactions
}
\examples{
interactions <-
    import_tf_mirna_interactions(
        resources = c('TransmiR')
    )

}
\seealso{
\itemize{
    \item{\code{\link{get_interaction_resources}}}
    \item{\code{\link{import_all_interactions}}}
    \item{\code{\link{interaction_graph}}}
    \item{\code{\link{print_interactions}}}
}
}
