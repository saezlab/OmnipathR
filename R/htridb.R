#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
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


#' Downloads TF-target interactions from HTRIdb
#'
#' HTRIdb (\url{https://www.lbbc.ibb.unesp.br/htri/}) is a database of
#' literature curated human TF-target interactions. As the database is
#' recently offline, the data is distributed by the OmniPath rescued data
#' repository (\url{https://rescued.omnipathdb.org/}).
#'
#' @return Data frame (tibble) with interactions.
#'
#' @examples
#' htridb_data <- htridb_download()
#' htridb_data
#' # # A tibble: 18,630 x 7
#' #      OID GENEID_TF SYMBOL_TF GENEID_TG SYMBOL_TG TECHNIQUE
#' #    <dbl>     <dbl> <chr>         <dbl> <chr>     <chr>
#' #  1 32399       142 PARP1           675 BRCA2     Electrophoretic Mobi.
#' #  2 32399       142 PARP1           675 BRCA2     Chromatin Immunoprec.
#' #  3 28907       196 AHR            1543 CYP1A1    Chromatin Immunoprec.
#' #  4 29466       196 AHR            1543 CYP1A1    Electrophoretic Mobi.
#' #  5 28911       196 AHR            1543 CYP1A1    Chromatin Immunoprec.
#' # # . with 18,620 more rows, and 1 more variable: PUBMED_ID <chr>
#'
#' @export
#' @importFrom readr cols col_skip
#' @importFrom magrittr %T>%
htridb_download <- function(){

    .slow_doctest()

    generic_downloader(
        url_key = 'htridb',
        reader = read_delim,
        reader_param = list(
            delim = ';',
            col_types = cols(
                X8 = col_skip()
            )
        ),
        resource = 'HTRIdb'
    ) %T>%
    load_success()

}
