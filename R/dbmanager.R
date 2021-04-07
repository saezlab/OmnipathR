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
                dbdef$loaded <- FALSE
                dbdef$last_used <- NA
                dbdef$db <- NA
                dbdef$latest_param <- NA
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

    # NSE vs. R CMD check workaround
    .loops <- list_queue_ <- cancel <- NULL

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
#' @return A data frame with the built in databaase definitions.
#'
#' @examples
#' database_definitions <- omnipath_show_db()
#' database_definitions
#' # # A tibble: 14 x 10
#' #    name       last_used           lifetime package  loader    loader_p.
#' #    <chr>      <dttm>                 <dbl> <chr>    <chr>     <list>
#' #  1 Gene Onto. 2021-04-04 20:19:15      300 Omnipat. go_ontol. <named l.
#' #  2 Gene Onto. NA                       300 Omnipat. go_ontol. <named l.
#' #  3 Gene Onto. NA                       300 Omnipat. go_ontol. <named l.
#' #  4 Gene Onto. NA                       300 Omnipat. go_ontol. <named l.
#' #  5 Gene Onto. NA                       300 Omnipat. go_ontol. <named l.
#' # ... (truncated)
#' # # . with 4 more variables: latest_param <list>, loaded <lgl>, db <list>,
#' # #   key <chr>
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest_wider
#' @export
omnipath_show_db <- function(){

    # NSE vs. R CMD check workaround
    db <- NULL

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

    omnipath.env$db %<>%
        map(
            function(dbdef){
                if(
                    dbdef$loaded &&
                    as.numeric(
                        Sys.time() -
                        dbdef$last_used,
                        units = 'secs'
                    ) > dbdef$lifetime
                ){
                    log_trace(
                        paste0(
                            'Removing database `%s` (not ',
                            'used in the past %d seconds)'
                        ),
                        dbdef$name,
                        dbdef$lifetime
                    )
                    dbdef$db <- NA
                    dbdef$loaded <- FALSE
                }
                dbdef
            }
        )

    later(db_lifetime_hook, delay = 10, loop = omnipath.env$.dbloop)

}


#' Load a built in database
#'
#' @param key Character: the key of the database to load. For a list of
#'     available keys see \code{\link{omnipath_show_db}}.
#' @param param List: override the defaults or pass further parameters to
#'     the database loader function. See the loader functions and their
#'     default parameters in \code{\link{omnipath_show_db}}.
#'
#' @examples
#' load_db('go_slim')
#'
#' @importFrom magrittr %<>%
#' @importFrom logger log_fatal log_info
#' @importFrom rlang exec !!!
#' @export
#' @seealso  \code{\link{omnipath_show_db}}.
load_db <- function(key, param = list()){

    db_exists(key)

    dbdef <- omnipath.env$db[[key]]
    log_info('Loading database `%s`.', dbdef$name)
    loader <- get(dbdef$loader)
    param %<>% add_defaults(loader, dbdef$loader_param)
    db <- exec(loader, !!!param)
    omnipath.env$db[[key]]$db <- db
    omnipath.env$db[[key]]$latest_param <- param
    omnipath.env$db[[key]]$loaded <- TRUE
    omnipath.env$db[[key]]$last_used <- Sys.time()
    log_info('Loaded database `%s`.', dbdef$name)

}


#' Access a built in database
#'
#' Databases are resources which might be costly to load but can be used many
#' times by functions which usually automatically load and retrieve them from
#' the database manager. Each database has a lifetime and will be unloaded
#' automatically upon expiry.
#'
#' @param key Character: the key of the database to load. For a list of
#'     available keys see \code{\link{omnipath_show_db}}.
#' @param param List: override the defaults or pass further parameters to
#'     the database loader function. See the loader functions and their
#'     default parameters in \code{\link{omnipath_show_db}}. If the database
#'     is already loaded with different parameters it will be reloaded
#'     with the new parameters only if the \code{reload} option is
#'     \code{TRUE}.
#' @param reload Reload the database if \code{param} passed here is different
#'     from the parameters used the last time the database was loaded. If
#'     different functions with different parameters access the database
#'     repeatedly and request reload the frequent reloads might cost
#'     substantial time and resource use.
#'
#' @examples
#' goslim <- get_db('go_slim')
#'
#' @importFrom logger log_fatal log_info
#' @importFrom rlang exec !!!
#' @importFrom magrittr %<>%
#' @export
#' @seealso  \code{\link{omnipath_show_db}}.
get_db <- function(key, param = NULL, reload = FALSE){

    db_exists(key)

    omnipath.env$db[[key]]$last_used <- Sys.time()
    dbdef <- omnipath.env$db[[key]]

    if(!is.null(param) && reload && dbdef$loaded){

        loader <- get(dbdef$loader)
        param %<>% add_defaults(loader, dbdef$loader_param)
        reload <- !lists_identical(dbdef$latest_param, param)

    }

    if(!dbdef$loaded || reload){

        load_db(key, param = param)

    }

    return(omnipath.env$db[[key]]$db)

}


#' Raise an error if the database key does not exist.
#'
#' @param key Character: a database definition key.
#'
#' @noRd
db_exists <- function(key){

    if(!key %in% names(omnipath.env$db)){
        msg <- sprintf('No database defined with key `%s`.', key)
        log_fatal(msg)
        stop(msg)
    }

}