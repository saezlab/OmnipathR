% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enz_sub.R
\name{import_omnipath_enzsub}
\alias{import_omnipath_enzsub}
\alias{import_Omnipath_PTMS}
\alias{import_OmniPath_PTMS}
\title{Imports enzyme-substrate relationships from OmniPath}
\usage{
import_omnipath_enzsub(
  resources = NULL,
  organism = 9606,
  fields = NULL,
  default_fields = TRUE,
  references_by_resource = TRUE,
  exclude = NULL,
  genesymbol_resource = NULL,
  ...
)
}
\arguments{
\item{resources}{PTMs not reported in these databases are
removed. See \code{\link{get_ptms_databases}} for more information.}

\item{organism}{Character or integer: Name or NCBI Taxonomy ID of one or
organisms. The web service currently provides interactions for
human, mouse and rat.}

\item{fields}{You can define here additional fields to be added to the
result. If used, set the next argument, \code{default_fields}, to
\code{FALSE}.}

\item{default_fields}{Whether to include the default fields (columns) for
the query type. If \code{FALSE}, only the fields defined by the user
in the \code{fields} argument will be added.}

\item{references_by_resource}{if FALSE, removes the resource name prefixes
from the references (PubMed IDs); this way the information which
reference comes from which resource will be lost and the PubMed IDs
will be unique.}

\item{exclude}{Character: datasets or resources to exclude.}

\item{genesymbol_resource}{Character: either "uniprot" or "ensembl". The
former leaves intact the gene symbols returned by the web service,
originally set from UniProt. The latter updates the gene symbols from
Ensembl, which uses a slightly different gene symbol standard. In this
case a few records will be duplicated, where Ensembl provides ambiguous
translation.}

\item{...}{Optional additional arguments.}
}
\value{
A data frame containing the information about ptms
}
\description{
Imports the enzyme-substrate (more exactly, enzyme-PTM) relationship
database from \url{https://omnipathdb.org/enzsub}
}
\examples{
enzsub <- import_omnipath_enzsub(
    resources = c('PhosphoSite', 'SIGNOR'),
    organism = 9606
)

}
\seealso{
\itemize{
    \item{\code{\link{get_enzsub_resources}}}
    \item{\code{\link{import_omnipath_interactions}}}
    \item{\code{\link{enzsub_graph}}}
    \item{\code{\link{print_interactions}}}
}
}
