#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2023
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Alberto Valdeolivas
#                  Dénes Türei (turei.denes@gmail.com)
#                  Attila Gábor
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Downloads the list of transcription factors from TF census
#'
#' Vaquerizas et al. published in 2009 a list of transcription factors. This
#' function retrieves Supplementary Table 2 from the article
#' (\url{http://www.nature.com/nrg/journal/v10/n4/index.html}).
#'
#' @return A data frame (tibble) listing transcription factors.
#'
#' @examples
#' tfcensus <- tfcensus_download()
#' tfcensus
#' # # A tibble: 1,987 x 7
#' #    Class `Ensembl ID` `IPI ID` `Interpro DBD` `Interpro DNA-b.
#' #    <chr> <chr>        <chr>    <chr>          <chr>
#' #  1 a     ENSG0000000. IPI0021. NA             IPR001289
#' #  2 a     ENSG0000000. IPI0004. IPR000047;IPR. NA
#' #  3 a     ENSG0000000. IPI0001. IPR001356;IPR. NA
#' #  4 a     ENSG0000000. IPI0029. IPR000910;IPR. NA
#' #  5 a     ENSG0000000. IPI0001. IPR007087;IPR. IPR006794
#' # # . with 1,977 more rows, and 2 more variables: `HGNC symbol` <chr>,
#' # # `Tissue-specificity` <chr>
#'
#' @export
#' @importFrom magrittr %T>%
tfcensus_download <- function(){

    .slow_doctest()

    suppressWarnings(
        generic_downloader(
            url_key = 'tfcensus',
            reader_param = list(
                skip = 11
            ),
            resource = 'TF census'
        )
    ) %T>%
    load_success()

}
