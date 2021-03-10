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


#' Downloads ligand-receptor interactions from Ramilowski et al. 2015
#'
#' Curated ligand-receptor pairs from Supplementary Table 2 of the article
#' "A draft network of ligand-receptor mediated multicellular signaling in
#' human" (https://www.nature.com/articles/ncomms8866).
#'
#' @examples
#' rami_interactions <- ramilowski_download()
#'
#' @export
#' @importFrom magrittr %T>%
ramilowski_download <- function(){

    xls_downloader(
        url_key = 'omnipath.ramilowski_url',
        sheet = 'All.Pairs',
        resource = 'Ramilowski et al. 2015'
    ) %T>%
    load_success()

}