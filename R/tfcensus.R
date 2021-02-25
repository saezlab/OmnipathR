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


#' Downloads the list of transcription factors from TF census
#'
#' Vaquerizas et al. published in 2009 a list of transcription factors. This
#' function retrieves Supplementary Table 2 from the article
#' (http://www.nature.com/nrg/journal/v10/n4/index.html).
#'
#' @export
tfcensus_download <- function(){

    suppressWarnings(
        generic_downloader(
            url_key = 'omnipath.tfcensus_url',
            reader_param = list(
                skip = 11
            )
        )
    )

}