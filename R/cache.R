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


#' Sets up the cache once the module is loaded or after the cachedir option
#' has been changed
#'
#' @importFrom rappdirs user_cache_dir
omnipath_init_cache <- function(){

    cachedir <- omnipath_default_cachedir()
    omnipath_new_cachedir(cachedir)
    options(omnipath.cachedir = cachedir)
    logger::log_info('Initialized cache: `%s`.', cachedir)

    invisible(cachedir)

}


#' Returns the default cache directory path
#'
#' @importFrom magrittr %>%
omnipath_default_cachedir <- function(){

    user_cache_dir() %>%
    file.path('OmnipathR') %>%
    normalizePath()

}


#' Sets up a new cache directory
#'
#' @importFrom jsonlite toJSON
omnipath_new_cachedir <- function(path){

    if(!omnipath_is_cachedir(path)){

        logger::log_info('Setting up new cache directory `%s`.', path)
        dir.create(path, showWarnings = FALSE, recursive = TRUE)
        .omnipath_cache <<- list()
        write(
            jsonlite::toJSON(.omnipath_cache, pretty = TRUE),
            file.path(path, 'cache.json')
        )

    }

    invisible(path)

}


#' Tells if a directory looks like an OmnipathR cache directory
omnipath_is_cachedir <- function(path){

    dir.exists(path) && file.exists(file.path(cache, 'cache.json'))

}


#' Determines the path of the cache directory according to the settings
omnipath_get_cachedir <- function(){

    cachedir <- options('omnipath.cachedir')

    if(is.null(cachedir) || !omnipath_is_cachedir(cachedir)){

        cachedir <- omnipath_init_cache()

    }

    return(cachedir)

}


#' Creates a lock file in the cache directory in order to avoid simulatneous
#' reading and writing of the cache file database
omnipath_lock_cache_db <- function(){

    lockfile <- omnipath_cache_lock_path()

    for(i in 1:options('omnipath.cache_timeout')){

        if(file.exists(lockfile))
        Sys.sleep(1)

    }

}


#' Removes the lock file from the cache directory
#'
#' A lock file in the cache directory avoids simulatneous write and read.
#' It's supposed to be removed after each read and write operation. This
#' might not happen if the process crashes during such an operation. In
#' this case you can manually call this function.
#'
#' @importFrom magrittr %>%
#' @export
omnipath_unlock_cache_db <- function(){

    omnipath_cache_lock_path() %>%
    file.remove()

}


#' Tells the path to the cache lock file
#'
#' @importFrom magrittr %>%
omnipath_cache_lock_path <- function(){

    'omnipath.cachedir' %>%
    options() %>%
    file.path('cache.lock')

}


#' Throws a fatal error about locked cache.
omnipath_locked_cache_error <- function(){

    msg <- sprintf(
        paste0(
            'The cache directory `%s` is locked for longer than %d seconds. ',
            'If you keep experiencing this issue please call `omnipath_'
            'unlock_cache_db()`. '
        ),
        options('omnipath.cachedir'),
        options('omnipath.cache_timeout')
    )
    logger::log_fatal(msg)
    stop(msg)

}


#' Searches for cache items
omnipath_cache_search <- function(){



}


#' Removes contents from the cache directory
#'
#' According to the parameters, it can remove contents older than a certain
#' age, or contents having a more recent version, one specific item, or
#' wipe the entire cache.
omnipath_cache_remove <- function(){



}


#' Removes one item from the cache directory
omnipath_cache_remove_one <- function(){



}


#' Removes the items from the cache directory which are unknown by the cache
#' database
omnipath_cache_clean <- function(){



}


#' Retrieves one item from the cache directory
omnipath_cache_get <- function(){



}


#' Adds a new item to the cache
omnipath_cache_add <- function(){



}


#' Creates a record for the cache database, describing a cache item with its
#' metadata such as download date
omnipath_cache_record <- function(){



}


#' Generates a hash which identifies an element in the cache database
omnipath_cache_key <- function(){



}


#' Registers a new element in the cache database
omnipath_cache_db_add <- function(){



}


#' Deletes an element from the cache database
omnipath_cache_db_remove <- function(){



}


#' Modifies an element in the cache database
omnipath_cache_db_edit <- function(){



}


#' Reads the cache DB contents from the disk to the memory
#'
#' @importFrom jsonlite
omnipath_read_cache_db <- function(){



}


#' Writes the cache DB contents from the memory to the disk
#'
#' @importFrom jsonlite toJSON
omnipath_write_cache_db <- function(){

    omnipath_lock_cache_db()
    write(
        jsonlite::toJSON(.omnipath_cache, pretty = TRUE),
        file.path()
    )
    omnipath_unlock_cache_db()

}