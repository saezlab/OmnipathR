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
    omnipath_read_cache_db()
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

        if(file.exists(lockfile)){

            Sys.sleep(1)

        }else{

            file.create(lockfile)
            return(TRUE)

        }

    }

    omnipath_locked_cache_error()

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

    omnipath_get_cachedir() %>%
    file.path('cache.lock')

}


#' Tells the path to the cache database
#'
#' @importFrom magrittr %>%
omnipath_cache_db_path <- function(){

    omnipath_get_cachedir() %>%
    file.path('cache.json')

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
#'
#' @importFrom RCurl merge.list
omnipath_cache_add <- function(record){

    omnipath_lock_cache_db()
    omnipath_read_cache_db()

    record$versions[[1]]$dl_finished <- Sys.time()

    if(record$key %in% .omnipath_cache){

        .omnipath_cache[[record$key]]$versions <-
            RCurl::merge.list(
                .omnipath_cache[[record$key]]$versions,
                record$versions
            )

    }else{

        .omnipath_cache[[record$key]] <- record

    }

    omnipath_write_cache_db()
    omnipath_unlock_cache_db()

}


#' Creates a record for the cache database, describing a cache item with its
#' metadata such as download date
#'
#' @importFrom magrittr %>% %<>%
omnipath_cache_record <- function(
        key,
        url,
        version = 1,
        ext = NULL,
        post = NULL,
        payload = NULL,
        dl_start = Sys.time(),
        status = 'downloading'
    ){

    version %<>% as.character

    if(is.null(ext)) ext <- file_extension(url)

    path <-
        omnipath_get_cachedir() %>%
        file.path(
            sprintf('%s-%s', key, version)
        ) %>%
        file_add_extension(ext)


    list(
        key = key,
        url = url,
        post = post,
        payload = payload,
        ext = ext,
        versions =
            list(
                list(
                    number = version,
                    path = path,
                    dl_start = dl_start,
                    dl_finished = NULL,
                    status = status
                )
            ) %>%
            setNames(version)
    )

}


#' Generates a hash which identifies an element in the cache database
#'
#' @importFrom digest sha1_digest
#' @importFrom magrittr %>%
omnipath_cache_key <- function(url, post = NULL, payload = NULL){

    list(url, post, payload) %>%
    sha1_digest()

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

    .omnipath_cache <<-
        omnipath_cache_db_path() %>%
        jsonlite::fromJSON()

}


#' Writes the cache DB contents from the memory to the disk
#'
#' Never call this function directly.
#'
#' @importFrom jsonlite toJSON
#' @importFrom magrittr %>%
omnipath_write_cache_db <- function(){

    .omnipath_cache %>%
    jsonlite::toJSON(pretty = TRUE) %>%
    write(omnipath_cache_db_path())

}