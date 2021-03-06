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


#' Downloads TF-target interactions from TRRUST
#'
#' TRRUST v2 (https://www.grnpedia.org/trrust/) is a database of literature
#' mined TF-target interactions for human and mouse.
#'
#' @param organism Character: either "human" or "mouse".
#'
#' @return A data frame of TF-target interactions.
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom dplyr mutate recode
#' @importFrom tidyr separate_rows
trrust_download <- function(organism = 'human'){

    generic_downloader(
        url_key = 'omnipath.trrust_url',
        url_param = list(organism),
        reader_param = list(
            col_names = c(
                'source_genesymbol',
                'target_genesymbol',
                'effect',
                'reference'
            )
        ),
        resource = 'TRRUST'
    ) %>%
    mutate(
        effect = recode(
            effect,
            Repression = -1,
            Unknown    =  0,
            Activation =  1
        )
    ) %>%
    separate_rows(
        reference,
        sep = ';'
    ) %T>%
    load_success()

}