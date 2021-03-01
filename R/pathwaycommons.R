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


#' Interactions from PathwayCommons
#'
#' PathwayCommons (http://www.pathwaycommons.org/) provides molecular
#' interactions from a number of databases, in either BioPAX or SIF (simple
#' interaction format). This function retrieves all interactions in SIF
#' format. The data is limited to the interacting pair and the type of the
#' interaction.
#'
#' @importFrom readr cols
#' @export
pathwaycommons_download <- function(){

    generic_downloader(
        url_key = 'omnipath.pathwaycommons_url',
        reader_param = list(
            col_names = c('from', 'type', 'to'),
            col_types = cols()
        )
    )

}