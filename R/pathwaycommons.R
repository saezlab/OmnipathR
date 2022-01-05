#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2022
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
#' PathwayCommons (\url{http://www.pathwaycommons.org/}) provides molecular
#' interactions from a number of databases, in either BioPAX or SIF (simple
#' interaction format). This function retrieves all interactions in SIF
#' format. The data is limited to the interacting pair and the type of the
#' interaction.
#'
#' @return A data frame (tibble) with interactions.
#'
#' @examples
#' pc_interactions <- pathwaycommons_download()
#'  pc_interactions
#' # # A tibble: 1,884,849 x 3
#' #    from  type                        to
#' #    <chr> <chr>                       <chr>
#' #  1 A1BG  controls-expression-of      A2M
#' #  2 A1BG  interacts-with              ABCC6
#' #  3 A1BG  interacts-with              ACE2
#' #  4 A1BG  interacts-with              ADAM10
#' #  5 A1BG  interacts-with              ADAM17
#' # # . with 1,884,839 more rows
#'
#' @importFrom readr cols
#' @export
pathwaycommons_download <- function(){

    generic_downloader(
        url_key = 'pathwaycommons',
        reader_param = list(
            col_names = c('from', 'type', 'to'),
            col_types = cols()
        )
    )

}