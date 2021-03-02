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


#' Protein-protein interactions from Vinayagam 2011
#'
#' Retrieves the Supplementary Table S6 from Vinayagam et al. 2011.
#' Find out more at https://doi.org/10.1126/scisignal.2001699
#'
#' @importFrom magrittr %>%
#' @importFrom readxl read_xls
#' @importFrom utils unzip
#' @export
vinayagam_download <- function(){

    'omnipath.vinayagam_url' %>%
    archive_downloader() %>%
    `$`('path') %>%
    unzip(
        files = '2001699_Tables_S1_S2_S6.xls',
        exdir = tempdir()
    ) %>%
    `[`(1) %>%
    read_xls(
        sheet = 'S6',
        progress = FALSE
    )

}