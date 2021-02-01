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


#' Downloads a Harmonizome network dataset
#'
#' Downloads a single network dataset from Harmonizome
#' https://maayanlab.cloud/Harmonizome
#'
#' @param dataset The dataset part of the URL. Please refer to the download
#'     section of the Harmonizome webpage.
#' @importsFrom dplyr %>% mutate select
#' @importsFrom readr read_tsv read_lines
#' @export
harmonizome_download <- function(dataset){

    'omnipath.harmonizome_url' %>%
    options() %>%
    as.character() %>%
    sprintf(dataset) %>%
    url() %>%
    gzcon() %>%
    read_lines() %>%
    `[`(-2) %>%
    read_tsv(col_types = cols()) %>%
    select(from = source, to = target) %>%
    mutate(source = dataset)

}