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


.cache_status <- list(
    UNKNOWN <- 'unknown',
    STARTED <- 'started',
    READY <- 'ready',
    FAILED <- 'failed',
    DELETED <- 'deleted'
)

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
    file.remove(showWarnings = FALSE)

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
            'If you keep experiencing this issue please call `omnipath_',
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
#'
#' @param url URL pointing to the resource
#' @param post HTTP POST parameters as a list
#' @param payload HTTP data payload
#' @oaram create Create a new entry if doesn't exist yet
#' @return Cache record: an existing record if the entry already exists,
#' otherwise a newly created and inserted record
#' @export
omnipath_cache_get <- function(
    key = NULL,
    url = NULL,
    post = NULL,
    payload = NULL,
    create = TRUE,
    ...
){

    key <- omnipath_cache_ensure_key(
        key = key,
        url = url,
        post = post,
        payload = payload
    )

    if(!(key %in% names(.omnipath_cache)) && create){

        if(is.null(url)){

            msg <- 'Can not create cache record without an URL.'
            logger::log_fatal(msg)
            stop(msg)

        }

        record <- omnipath_cache_record(
            key = key,
            url = url,
            post = post,
            payload = payload,
            ...
        )
        omnipath_cache_add(record, new = TRUE)

    }

    record <- `if`(
        key %in% names(.omnipath_cache),
        .omnipath_cache[[key]],
        NULL
    )

    return(record)

}


#' Loads an R object from the cache
#'
#' Loads the object from RDS format.
#'
#' @param key Key of the cache item
#' @param url URL of the downloaded resource
#' @param post HTTP POST parameters as a list
#' @patam payload HTTP data payload
#' @param version Version of the cache item. If does not exist or NULL, the
#' latest version will be retrieved
#'
#' @importFrom logger log_info log_trace
#' @export
#' @seealso \code{\link{omnipath_cache_save}}
omnipath_cache_load <- function(
    key = NULL,
    version = NULL,
    url = NULL,
    post = NULL,
    payload = NULL
){

    record <- omnipath_cache_get(
        key = key,
        url = url,
        post = post,
        payload = payload,
        create = FALSE
    )

    if(is.null(record)){
        logger::log_info(
            'Cache record does not exist: '
        )
        return(NULL)
    }

    version <- `if`(
        (
            is.null(version) ||
            !(version %in% record$versions) ||
            !(record$versions[[version]]$status == .cache_status$READY)
        ),
        omnipath_cache_latest_version(record),
        version
    )

    if(!is.null(version)){

        path <- record$versions[[version]]$path
        data <- loadRDS(path)
        logger::log_trace('Loaded data from RDS `%s`.', path)
        return(data)

    }else{
        logger::log_info('No version is available for key `%s`.', key)
    }

}


#' Saves an R object to the cache
#'
#' Exports the object in RDS format, creates new cache record if necessary.
#'
#' @param data An object
#' @param key Key of the cache item
#' @param url URL of the downloaded resource
#' @param post HTTP POST parameters as a list
#' @patam payload HTTP data payload
#' @param version Version of the cache item. If does not exist a new version
#' item will be created
#'
#' @importFrom logger log_info
#' @export
#' @seealso \code{\link{omnipath_cache_move_in}}
omnipath_cache_save <- function(
    data,
    key = NULL,
    version = NULL,
    url = NULL,
    post = NULL,
    payload = NULL
){

    record <- omnipath_cache_get(
        key = key,
        url = url,
        post = post,
        payload = payload,
        ext = 'rds'
    )

    version <- omnipath_cache_update_status(
        key = record$key,
        version = version,
        status = .cache_status$STARTED
    )

    target_path <- .omnipath_cache[[key]]$versions[[version]]$path

    saveRDS(data, target_path)
    logger::log_trace('Exported RDS to `%s`.', target_path)

    omnipath_cache_download_ready(key, version)

}


#' Moves an existing file into the cache
#'
#' Either the key or the URL (with POST and payload) must be provided.
#'
#' @param path Path to the source file
#' @param key Key of the cache item
#' @param url URL of the downloaded resource
#' @param post HTTP POST parameters as a list
#' @patam payload HTTP data payload
#' @param version Version of the cache item. If does not exist a new version
#' item will be created
#' @param keep_original Whether to keep or remove the original file
#'
#' @importFrom logger log_info
#' @export
#' @seealso \code{\link{omnipath_cache_save}}
omnipath_cache_move_in <- function(
        path,
        key = NULL,
        version = NULL,
        url = NULL,
        post = NULL,
        payload = NULL,
        keep_original = FALSE
    ){

    version <- omnipath_cache_update_status(
        key = key,
        version = version,
        status = .cache_status$STARTED
    )

    target_path <- .omnipath_cache[[key]]$versions[[version]]$path

    file.copy(path, target_path)
    logger::log_trace('Copied `%s` to `%s`.', path, target_path)

    if(!keep_original){

        file.remove(path)

    }

    omnipath_cache_download_ready(key, version)

}


#' Sets the download status to ready for a cache item
#'
#' @param key Key of the cache item
#' @param version Version of the cache item. If does not exist a new version
#' item will be created
#' @export
omnipath_cache_download_ready <- function(key, version){

    omnipath_cache_update_status(
        key = key,
        version = version,
        status = 'ready',
        dl_finished = Sys.time()
    ) %>%
    invisible()

}


#' Updates the status of an existing cache record
#'
#' @param key Key of the cache item
#' @param version Version of the cache item. If does not exist a new version
#' item will be created
#' @param status The updated status value
#' @param dl_finished Timestamp for the time when download was finished,
#' if NULL the value remains unchanged
#'
#' @importFrom magrittr %<>%
#' @importFrom logger log_info log_warning
omnipath_cache_update_status <- cache_locked %@% function(
    key,
    version,
    status,
    dl_finished = NULL
){

    version %<>% as.character

    if(key %in% names(.omnipath_cache)){

        if(!(version %in% names(.omnipath_cache[[key]]))){

            version <- omnipath_cache_new_version(key, version = version)

        }

        old_status <- .omnipath_cache[[key]]$versions[[version]]$status
        .omnipath_cache[[key]]$versions[[version]]$status <- status

        if(!is.null(dl_finished)){
            .omnipath_cache[[key]]$versions[[version]]$dl_finished <-
                dl_finished
        }

        logger::log_info(
            'Cache item `%s` version %s: status changed from `%s` to `%s`.',
            key,
            version
            old_status,
            status
        )

    }else{

        logger::log_warning(
            'Failed to update cache: key `%s` does not exist.',
            key
        )

    }

    invisible(version)

}


#' Adds a new version item to an existing cache record
#'
#' @importFrom RCurl merge.list
omnipath_cache_new_version <- function(key, version = NULL){

    version <- `if`(
        is.null(version),
        omnipath_cache_new_version(key),
        version
    )

    record <- .omnipath_cache[[key]]

    new_version_record <- omnipath_cache_record(
        key,
        record$url,
        version = version,
        ext = record$ext,
        post = record$post,
        payload = record$payload,
        dl_started = NULL,
        status = .cache_status$UNKNOWN
    )

    .omnipath_cache[[key]]$versions <- merge.list(
        .omnipath_cache[[key]]$versions,
        new_version_record$versions
    )

    invisible(version)

}


#' Tells the next incremental version number for a cache record
#'
#' @return The next version number as character
#' @importFrom purrr map
#' @importFrom magrittr %>%
omnipath_cache_next_version <- function(key){

    .omnipath_cache[[key]]$versions %>%
    purrr::map(
        function(version){
            version$number
        }
    ) %>%
    unlist() %>%
    as.numeric() %>%
    max(0) %>%
    `+`(1) %>%
    as.character()

}


#' Finds the most recent version in a cache record
#'
#' @return The version ID with the most recent download finished time
#' @importFrom magrittr %<>%
omnipath_cache_latest_version <- function(record, max_age = NULL){

    versions <- names(record$versions)
    dates <- map(record$versions, function(v){v$dl_finished})
    versions <- versions[!sapply(dates, is.null)]
    dates %<>% unlist
    return(versions[which.max(dates)])
}


#' Adds a new item to the cache or updates an existing one
#'
#' @importFrom RCurl merge.list
omnipath_cache_add <- cache_locked %@% function(record, new = FALSE){

    if(record$key %in% .omnipath_cache && !new){

        .omnipath_cache[[record$key]]$versions <-
            RCurl::merge.list(
                .omnipath_cache[[record$key]]$versions,
                record$versions
            )

    }else{

        .omnipath_cache[[record$key]] <- record

    }

}


#' Returns the key if either the key or the URL is available
#'
#' @importFrom logger log_fatal
omnipath_cache_ensure_key <- function(
    key = NULL,
    url = NULL,
    post = NULL,
    payload = NULL
){

    if(is.null(key)){

        if(is.null(url)){

            msg <- 'Neither key or URL is available.'
            logger::log_fatal(msg)
            stop(msg)

        }

        key <- omnipath_cache_key(url = url, post = post, payload = payload)

    }

    return(key)

}


#' Executes a function with locking the cache database
cache_locked <- decorator %@% function(FUN){

    function(...){

        omnipath_lock_cache_db()
        omnipath_read_cache_db()

        result <- FUN(...)

        omnipath_write_cache_db()
        omnipath_unlock_cache_db()

        invisible(result)

    }

}


#' Creates a record for the cache database, describing a cache item with its
#' metadata such as download date
#'
#' @importFrom magrittr %>% %<>%
omnipath_cache_record <- function(
        key,
        url,
        version = '1',
        ext = NULL,
        post = NULL,
        payload = NULL,
        dl_started = Sys.time(),
        status = .cache_status$STARTED
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
                    dl_started = dl_started,
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


#' Reads the cache DB contents from the disk to the memory
#'
#' @importFrom jsonlite
omnipath_read_cache_db <- function(){

    .omnipath_cache <<-
        omnipath_cache_db_path() %>%
        jsonlite::fromJSON() %>%
        omnipath_cache_timestamps()

}


#' Converts the timestamps read from JSON to POSIXct objects
#'
#' @importFrom magrittr %<>%
omnipath_cache_timestamps <- function(cache_db){

    to_posixt <- function(value){

        `if`(is.null(value), NULL, as.POSIXct(value))

    }

    for(key in names(cache_db)){

        for(ver in names(cache_db[[key]]$versions)){

            cache_db[[key]]$versions[[ver]]$dl_started %<>% to_posixt()
            cache_db[[key]]$versions[[ver]]$dl_finished %<>% to_posixt()

        }

    }

    return(cache_db)

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