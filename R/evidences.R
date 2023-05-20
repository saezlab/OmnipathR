#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2023
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


#' Converts the extra_attrs column from JSON encoded to list
#'
#' @param data A data frame from the OmniPath web service.
#'
#' @return The input data frame with the extra_attrs column converted
#'     to list.
#'
#' @importFrom magrittr %>%
#' @noRd
deserialize_evidences <- function(data){

    data %>%
    deserialize_json_col('evidences')

}


#' Tells if a data frame has a column "evidences"
#'
#' @importFrom magrittr %>%
#' @noRd
has_evidences <- function(data){

    data %>% has_column('evidences')

}
