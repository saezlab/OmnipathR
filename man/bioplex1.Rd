% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioplex.R
\name{bioplex1}
\alias{bioplex1}
\title{Downloads the BioPlex version 1.0 interaction dataset}
\usage{
bioplex1()
}
\value{
Data frame (tibble) with interactions.
}
\description{
This dataset contains ~24,000 interactions detected in HEK293T cells
using 2,594 baits.
More details at \url{https://bioplex.hms.harvard.edu/interactions.php}.
}
\examples{
bioplex_interactions <- bioplex1()
nrow(bioplex_interactions)
# [1] 23744
colnames(bioplex_interactions)
# [1] "GeneA"         "GeneB"        "UniprotA"   "UniprotB"
# [5] "SymbolA"       "SymbolB"      "p_wrong"    "p_no_interaction"
# [9] "p_interaction"

}
\seealso{
\itemize{
    \item{\code{\link{bioplex2}}}
    \item{\code{\link{bioplex3}}}
    \item{\code{\link{bioplex_hct116_1}}}
    \item{\code{\link{bioplex_all}}}
}
}
