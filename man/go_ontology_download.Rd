% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/go.R
\name{go_ontology_download}
\alias{go_ontology_download}
\title{The Gene Ontology tree}
\usage{
go_ontology_download(
  basic = TRUE,
  tables = TRUE,
  subset = NULL,
  relations = c("is_a", "part_of", "occurs_in", "regulates", "positively_regulates",
    "negatively_regulates")
)
}
\arguments{
\item{basic}{Logical: use the basic or the full version of GO. As written
on the GO home page: "the basic version of the GO is filtered such
that the graph is guaranteed to be acyclic and annotations can be
propagated up the graph. The relations included are is a, part of,
regulates, negatively regulates and positively regulates. This
version excludes relationships that cross the 3 GO hierarchies.
This version should be used with most GO-based annotation tools."}

\item{tables}{In the result return data frames or nested lists. These
later can be converted to each other if necessary. However converting
from table to list is faster.}

\item{subset}{Character: the GO subset (GO slim) name. GO slims are
subsets of the full GO which "give a broad overview of the ontology
content without the detail of the specific fine grained terms". This
option, if not \code{NULL}, overrides the \code{basic} parameter.
Available GO slims are: "agr" (Alliance for Genomics Resources),
"generic", "aspergillus", "candida", "drosophila", "chembl",
"metagenomic", "mouse", "plant", "pir" (Protein Information Resource),
"pombe" and "yeast".}

\item{relations}{Character vector: the relations to include in the
processed data.}
}
\value{
A list with the following elements: 1) "names" a list with
    terms as names and names as values; 2) "namespaces" a list with
    terms as names and namespaces as values; 3) "relations" a list with
    relations between terms: terms are keys, values are lists with
    relations as names and character vectors of related terms as
    values; 4) "subsets" a list with terms as keys and character
    vectors of subset names as values (or \code{NULL} if the term
    does not belong to any subset); 5) "obsolete" character vector
    with all the terms labeled as obsolete. If the \code{tables}
    parameter is \code{TRUE}, "names", "namespaces", "relations"
    and "subsets" will be data frames (tibbles).
}
\description{
The Gene Ontology tree
}
\examples{
# retrieve the generic GO slim, a small subset of the full ontology
go <- go_ontology_download(subset = 'generic')

}
