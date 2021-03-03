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


#' Downloads TF-target interactions from HTRIdb
#'
#' HTRIdb (https://www.lbbc.ibb.unesp.br/htri/) is a database of literature
#' curated human TF-target interactions. As the database is recently offline,
#' the data is distributed by the OmniPath rescued data repository
#' (https://rescued.omnipathdb.org/).
#'
#' @export
#' @importFrom readr cols col_skip
#' @importFrom magrittr %T>%
htridb_download <- function(){

    suppressWarnings(
        generic_downloader(
            url_key = 'omnipath.htridb_url',
            reader = read_delim,
            reader_param = list(
                delim = ';',
                col_types = cols(
                    X8 = col_skip()
                )
            ),
            resource = 'HTRIdb'
        )
    ) %T>%
    load_success()

}