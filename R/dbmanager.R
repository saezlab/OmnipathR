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


#' Creates an empty list to store the loaded databases
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom jsonlite fromJSON
#' @importFrom purrr map
#' @noRd
omnipath_init_db <- function(pkgname){

    omnipath.env$db <-
        system.file(
            'db',
            'db_def.json',
            package = pkgname,
            mustWork = TRUE
        ) %>%
        fromJSON() %>%
        map(
            function(dbdef){
                dbdef$lifetime %<>% if_empty(
                    getOption('omnipath.db_lifetime')
                )
                dbdef$package %<>% if_empty(pkgname)
                dbdef
            }
        )

}


#' Built in database definitions
#'
#' Databases are resources which might be costly to load but can be used many
#' times by functions which usually automatically load and retrieve them from
#' the database manager. Each database has a lifetime and will be unloaded
#' automatically upon expiry.
#'
#' @return A data frame with the build in databaase definitions.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest_wider
#' @export
omnipath_db_available <- function(){

    omnipath.env$db %>%
    tibble(db = .) %>%
    mutate(key = names(db)) %>%
    unnest_wider(db)

}