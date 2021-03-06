% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trrust.R
\name{trrust_download}
\alias{trrust_download}
\title{Downloads TF-target interactions from TRRUST}
\usage{
trrust_download(organism = "human")
}
\arguments{
\item{organism}{Character: either "human" or "mouse".}
}
\value{
A data frame of TF-target interactions.
}
\description{
TRRUST v2 (\url{https://www.grnpedia.org/trrust/}) is a database of
literature mined TF-target interactions for human and mouse.
}
\examples{
trrust_interactions <- trrust_download()
trrust_interactions
# # A tibble: 11,698 x 4
#    source_genesymbol target_genesymbol effect reference
#    <chr>             <chr>              <dbl> <chr>
#  1 AATF              BAX                   -1 22909821
#  2 AATF              CDKN1A                 0 17157788
#  3 AATF              KLK3                   0 23146908
#  4 AATF              MYC                    1 20549547
#  5 AATF              TP53                   0 17157788
#  6 ABL1              BAX                    1 11753601
#  7 ABL1              BCL2                  -1 11753601
# # . with 11,688 more rows

}
