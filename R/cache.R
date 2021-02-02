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
omnipath_init_cache <- function(){



}


#' Sets up a new cache directory
omnipath_new_cachedir <- function(){



}


#' Determines the path of the cache directory according to the settings
omnipath_get_cachedir <- function(){



}


#' Creates a lock file in the cache directory in order to avoid simulatneous
#' reading and writing of the cache file database
omnipath_lock_cache_db <- function(){



}


#' Removes the lock file from the cache directory
#'
#' A lock file in the cache directory avoids simulatneous write and read.
#' It's supposed to be removed after each read and write operation. This
#' might not happen if the process crashes during such an operation. In
#' this case you can manually call this function.
#'
#' @export
omnipath_unlock_cache_db <- function(){



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