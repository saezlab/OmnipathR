% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioplex.R
\name{bioplex_hct116_1}
\alias{bioplex_hct116_1}
\title{Downloads the BioPlex HCT116 version 1.0 interaction dataset}
\usage{
bioplex_hct116_1()
}
\value{
Data frame (tibble) with interactions.
}
\description{
This dataset contains ~71,000 interactions detected in HCT116 cells
using 5,522 baits.
More details at \url{https://bioplex.hms.harvard.edu/interactions.php}.
}
\examples{
bioplex_interactions <- bioplex_hct116_1()
nrow(bioplex_interactions)
# [1] 70966
colnames(bioplex_interactions)
# [1] "GeneA"         "GeneB"        "UniprotA"   "UniprotB"
# [5] "SymbolA"       "SymbolB"      "p_wrong"    "p_no_interaction"
# [9] "p_interaction"

}
\seealso{
\itemize{
    \item{\code{\link{bioplex1}}}
    \item{\code{\link{bioplex2}}}
    \item{\code{\link{bioplex3}}}
    \item{\code{\link{bioplex_all}}}
}
}
