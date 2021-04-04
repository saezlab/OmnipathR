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


#' Populates the database register with the built in database definitions
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom jsonlite fromJSON
#' @importFrom purrr map
#' @importFrom later create_loop global_loop
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

    remove_tasks(pkgname)
    omnipath.env$.dbloop <- create_loop(parent = global_loop())

    db_lifetime_hook()

}


#' Removes the tasks scheduled in `later` by a certain package in order to
#' prevent their duplication. Necessary upon package reload.
#'
#' @param pkgname Character: name of the package. Only tasks with a callback
#'     belonging to this package's namespace will be canceled.
#'
#' @importFrom later destroy_loop
#' @importFrom magrittr %>%
#' @importFrom purrr walk
#' @importFrom rlang wref_key exec
#'
#' @noRd
remove_tasks <- function(pkgname){

    later%:::%.loops %>% ls %>% as.integer %>%
    walk(
        function(loop_id){
            envname <- ''
            for(task in exec(later%:::%list_queue_, loop_id)){
                envname <-
                    task$callback %>%
                    environment %>%
                    environmentName
                if(envname == pkgname){
                    exec(later%:::%cancel, as.character(task$id), loop_id)
                }
            }
            if(envname == pkgname){
                as.character(loop_id) %>%
                get(envir = later%:::%.loops) %>%
                wref_key %>%
                destroy_loop
            }
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


#' Remove unused databases
#'
#' In every 10 seconds unregisters the loaded databases if their lifetime
#' expired.
#'
#' @importFrom magrittr %<>%
#' @importFrom later later
#' @importFrom logger log_trace
#'
#' @noRd
db_lifetime_hook <- function(){

    log_trace('Cleaning up unused databases.')

    omnipath.env$db %<>%
        map(
            function(dbdef){
                if(
                    dbdef$loaded &&
                    as.numeric(
                        Sys.time() -
                        dbdef$last_used
                    ) > dbdef$lifetime
                ){
                    log_trace(
                        paste0(
                            'Removing database `%s` (not used in the ',
                            'past %d seconds)'
                        ),
                        dbdef$name,
                        dbdef$lifetime
                    )
                    dbdef$db <- list()
                    dbdef$loaded <- FALSE
                }
                dbdef
            }
        )

    later(db_lifetime_hook, delay = 10, loop = omnipath.env$.dbloop)

}