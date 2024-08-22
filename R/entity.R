#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
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


#' Make sure the argument is a valid entity type
#'
#' @importFrom magrittr %>%
#' @noRd
ensure_entity_type <- function(entity_type = NULL) {

    valid_et <-
        entity_type %>%
        if_null('.default') %>%
        str_to_lower() %>%
        extract2(omnipathr.env$entity_types, .)

    if (is.null(valid_et)) {

        err <- sprintf('Unknown entity type: `%s`.', entity_type)
        log_error(err)
        stop(err)

    }

    return(valid_et)

}


#' The default entity type
#'
#' @noRd
default_entity_type <- function() {

    omnipathr.env$entity_types$.default

}


#' Load the entity type names and synonyms
#'
#' @importFrom yaml read_yaml
#' @noRd
.load_entity_types <- function(pkgname) {

    omnipathr.env$entity_types <-
        system.file(
            'internal',
            'entity_types.yaml',
            package = pkgname,
            mustWork = TRUE
        ) %>%
        read_yaml() %>%
        map2(
            names(.),
            ~set_names(
                rep(.y, length(.x) + 1L),
                c(.x, .y)
            )
        ) %>%
        unname %>%
        unlist

}
