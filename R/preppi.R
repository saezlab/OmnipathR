#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2021
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
#  Website: https://saezlab.github.io/omnipathr
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Interactions from PrePPI
#'
#' Retrieves predicted protein-protein interactions from the PrePPI
#' database (\url{http://honig.c2b2.columbia.edu/preppi}). The interactions
#' in this table are supposed to be correct with a > 0.5 probability.
#'
#' @param ... Minimum values for the scores. The available scores are:
#'     str, protpep, str_max, red, ort, phy, coexp, go, total, exp and final.
#'
#' @return A data frame (tibble) of interactions with scores, databases
#'     and literature references.
#'
#' @details
#' PrePPI is a combination of many prediction methods, each resulting a
#' score. For an explanation of the scores see
#' \url{https://honiglab.c2b2.columbia.edu/hfpd/help/Manual.html}.
#' The minimum and maximum values of the scores:
#'     | Score   | Minimum | Maximum           |
#'     | ------- | ------- | ----------------- |
#'     | str     |       0 |           6,495   |
#'     | protpep |       0 |          38,138   |
#'     | str_max |       0 |          38,138   |
#'     | red     |       0 |              24.4 |
#'     | ort     |       0 |           5,000   |
#'     | phy     |       0 |            2.42   |
#'     | coexp   |       0 |              45.3 |
#'     | go      |       0 |             181   |
#'     | total   |       0 | 106,197,000,000   |
#'     | exp     |       1 |           4,626   |
#'     | final   |     600 |           4.91e14 |
#'
#' @examples
#' preppi <- preppi_download()
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv cols
#' @importFrom dplyr filter
#' @importFrom rlang !!!
#' @export
#' @md
preppi_download <- function(...){

    'preppi' %>%
    archive_extractor(
        path = 'preppi_final600.txt',
        reader = read_tsv,
        reader_param = list(
            col_types = cols(),
            progress = FALSE
        ),
        resource = 'PrePPI'
    ) %T>%
    load_success()

}