% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harmonizome.R
\name{harmonizome_download}
\alias{harmonizome_download}
\title{Downloads a Harmonizome network dataset}
\usage{
harmonizome_download(dataset)
}
\arguments{
\item{dataset}{The dataset part of the URL. Please refer to the download
section of the Harmonizome webpage.}
}
\value{
Data frame (tibble) with interactions.
}
\description{
Downloads a single network dataset from Harmonizome
\url{https://maayanlab.cloud/Harmonizome}.
}
\examples{
harmonizome_data <- harmonizome_download('phosphositeplus')
harmonizome_data
# # A tibble: 6,013 x 7
#    source   source_desc source_id target target_desc target_id weight
#    <chr>    <chr>           <dbl> <chr>  <chr>           <dbl>  <dbl>
#  1 TP53     na               7157 STK17A na               9263      1
#  2 TP53     na               7157 TP53RK na             112858      1
#  3 TP53     na               7157 SMG1   na              23049      1
#  4 UPF1     na               5976 SMG1   na              23049      1
# # . with 6,003 more rows

}
