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


CACHE_STATUS <- list(
    UNKNOWN = 'unknown',
    STARTED = 'started',
    READY = 'ready',
    FAILED = 'failed',
    DELETED = 'deleted'
)

BEGINNING_OF_TIME <- as.POSIXct(0, origin = '1970-01-01')

PROTECTED_FILES <- c('cache.lock', 'cache.json')

#' Sets up the cache once the module is loaded or after the cachedir option
#' has been changed
#'
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
#' @importFrom rappdirs user_cache_dir
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

    dir.exists(path) && file.exists(file.path(path, 'cache.json'))

}


#' Determines the path of the cache directory according to the settings
omnipath_get_cachedir <- function(){

    cachedir <- options('omnipath.cachedir')[[1]]

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


#' Searches for cache items
#'
#' Searches the cache records by matching the URL against a string or regexp.
#'
#' @return List of cache records matching the pattern.
#'
#' @param pattern String or regular expression.
#' ... Passed to \code{\link{grep}}
#'
#' @importFrom purrr map_chr
#' @export
omnipath_cache_search <- function(pattern, ...){

    .omnipath_cache %>%
    map_chr('url') %>%
    grep(pattern = pattern, ...) %>%
    `[`(names(.omnipath_cache), .) %>%
    `[`(.omnipath_cache, .)


}


#' Removes contents from the cache directory
#'
#' According to the parameters, it can remove contents older than a certain
#' age, or contents having a more recent version, one specific item, or
#' wipe the entire cache.
#'
#' @param key The key of the cache record
#' @param url URL pointing to the resource
#' @param post HTTP POST parameters as a list
#' @param payload HTTP data payload
#' @param max_age Age of cache items in days. Remove everything that is older
#' than this age
#' @param min_age Age of cache items in days. Remove everything more recent
#' than this age
#' @param status Remove items having any of the states listed here
#' @param keep_latest Keep the latest
#' @param autoclean Remove the entries about failed downloads, the files in
#' the cache directory which are missing from the cache database, and the
#' entries without existing files in the cache directory
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom purrr map keep
#' @export
#' @seealso \code{\link{omnipath_cache_wipe}, \link{omnipath_cache_clean},
#' \link{omnipath_cache_autoclean}}
omnipath_cache_remove <- cache_locked %@% function(
    key = NULL,
    url = NULL,
    post = NULL,
    payload = NULL,
    max_age = NULL,
    min_age = NULL,
    status = NULL,
    keep_latest = FALSE,
    wipe = FALSE,
    autoclean = TRUE
){

    if(wipe){

        .omnipath_cache <<- list()
        omnipath_cache_clean()
        return(invisible(NULL))

    }

    key <- ifelse(
        is.null(key),
        omnipath_cache_key(url = url, post = post, payload = payload),
        key
    )

    .omnipath_cache %<>%
    {`if`(
        is.null(key),
        .,
        `if`(
            (
                is.null(max_age) &&
                is.null(min_age) &&
                is.null(status) &&
                !keep_latest
            ),
            `[`(., .omnipath_cache %>% names %>% setdiff(key)),
            `[`(., key)
        )
    )} %>%
    map(
        omnipath_cache_remove_versions,
        max_age = max_age,
        min_age = min_age,
        status = status,
        keep_latest = keep_latest
    ) %>%
    {`if`(
        is.null(key),
        .,
        map(
            .omnipath_cache,
            function(record){
                `if`(
                    record$key %in% names(.),
                    .[[record$key]],
                    record
                )
            }
        )
    )}

    .omnipath_cache <<- .omnipath_cache

    if(autoclean){

        omnipath_write_cache_db()
        omnipath_unlock_cache_db()

        omnipath_cache_clean_db()

        omnipath_lock_cache_db()

        omnipath_cache_autoclean()

    }

}


#' Removes version items from a cache record
#'
#' @return A cache record with the version items removed
#'
#' @importFrom magrittr %<>% %>%
omnipath_cache_remove_versions <- function(
    record,
    max_age = NULL,
    min_age = NULL,
    status = NULL,
    keep_latest = FALSE
){

    record$versions %<>%
        record$versions %>%
        `[`(
            omnipath_cache_filter_versions(
                record = record,
                latest = keep_latest,
                max_age = max_age,
                min_age = min_age,
                status = setdiff(CACHE_STATUS, status)
            )
        )

    invisible(record)

}


#' Removes the cache database entries without existing files
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map keep
omnipath_cache_clean_db <- cache_locked %@% function(){

    .omnipath_cache <<-
        .omnipath_cache %>%
        map(
            ~keep(
                .x$versions,
                function(v){file.exists(v$path)}
            )
        ) %>%
        keep(
            function(record){length(record$versions) > 0}
        )

}


#' Permanently removes all the cache contents
#'
#' After this operation the cache directory will be completely empty,
#' except an empty cache database file.
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \code{\link{omnipath_cache_remove}}
omnipath_cache_wipe <- cache_locked %@% function(){

    omnipath_get_cachedir() %>%
    list.files() %>%
    setdiff(PROTECTED_FILES) %>%
    file.path(omnipath_get_cachedir(), .) %>%
    file.remove()

    .omnipath_cache <<- list()
    return(invisible(NULL))

}


#' Removes the items from the cache directory which are unknown by the cache
#' database
#'
#' @importFrom purrr map_chr map
#' @importFrom magrittr %>%
omnipath_cache_clean <- function(){



    files_in_db <-
        .omnipath_cache %>%
        map(~map_chr(.x$versions, 'path')) %>%
        unlist() %>%
        basename()

    omnipath_get_cachedir() %>%
    list.files() %>%
    setdiff(files_in_db) %>%
    setdiff(PROTECTED_FILES) %>%
    file.path(omnipath_get_cachedir(), .) %>%
    file.remove()

    invisible(NULL)

}


#' Keeps only the latest versions of complete downloads
#'
#' Removes the old versions, the failed downloads and the files in the cache
#' directory which are missing from the database. For more flexible
#' operations use \link{\code{omnipath_cache_remove}} and
#' \link{\code{omnipath_cache_clean}}.
#'
#' @export
#'
#' Removes the items from the cache directory which are unknown by the cache
#' database
omnipath_cache_autoclean <- function(){

    omnipath_cache_remove(
        keep_most_recent = TRUE,
        status = CACHE_STATUS$READY
    )
    omnipath_cache_clean()

}


#' Retrieves one item from the cache directory
#'
#' @param key The key of the cache record
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
            !(record$versions[[version]]$status == CACHE_STATUS$READY)
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
        status = CACHE_STATUS$STARTED
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
        status = CACHE_STATUS$STARTED
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
            version,
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
        status = CACHE_STATUS$UNKNOWN
    )

    .omnipath_cache[[key]]$versions <- RCurl::merge.list(
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
#' @return Character: the version ID with the most recent download finished
#' time
#' @param record A cache record
#'
#' @export
omnipath_cache_latest_version <- function(record){

    omnipath_cache_filter_versions(record = record, latest = TRUE)

}

#' Filters the versions from one cache record
#'
#' Filters the versions based on multiple conditions: their age and status
#'
#' @return Character vector with version IDs, NA if no version satisfies the
#' conditions
#' @param record A cache record
#' @param latest Return the most recent version
#' @param max_age The maximum age in days (e.g. 5: 5 days old or more recent)
#' @param min_age The minimum age in days (e.g. 5: 5 days old or older)
#' @param status Character vector with status codes. By default only the
#' versions with `ready` (completed download) status are selected
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom purrr map map_chr
#' @export
omnipath_cache_filter_versions <- function(
    record,
    latest = FALSE,
    max_age = NULL,
    min_age = NULL,
    status = CACHE_STATUS$READY
){

    days_ago <- function(days){

        as.POSIXct(days * 24 * 3600 * -1, origin = Sys.time())

    }

    version_ids <- names(record$versions)
    versions <- record$versions
    statuses <- map_chr(versions, 'status')
    selection <- statuses %in% status

    if(!is.null(max_age)){
        max_age <- days_ago(max_age)
        selection %<>% `&`(which_dl_finished(versions, max_age))
    }

    if(!is.null(min_age)){
        min_age <- days_ago(min_age)
        selection %<>% `&`(which_dl_finished(versions, max_age, op = `<=`))
    }

    if(latest){
        t_latest <- versions %>%
            map(function(v){v$dl_finished}) %>%
            unlist()
        selection %<>% `&`(which_dl_finished, t_latest, op = `==`)
    }


    dates <- map(record$versions, function(v){v$dl_finished})
    versions <- versions[!sapply(dates, is.null)]
    dates %<>% unlist
    return(versions[which.max(dates)])

}


#' Selects the versions which fit a time threshold
#'
#' @param versions List of version entries
#' @param t Timestamp (POSIXct object)
#' @param op Operator, e.g. `>=` for larger than (more recent) or `<=` (older)
#'
#' @return Logical vector same length as `versions`
#'
#' @importFrom purrr map_lgl
which_dl_finished <- function(versions, t, op = `>=`){

    map_lgl(
        versions,
        function(v){
            !is.null(v$dl_finished) &&
            op(v$dl_finished, t)
        }
    )

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
        status = CACHE_STATUS$STARTED
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
#' @param url Character vector with URLs
#' @param post List with the HTTP POST parameters or a list of lists if
#' the url vector is longer than 1. NULL for queries without POST parameters.
#' @param payload HTTP data payload. List with multiple items if the url
#' vector is longer than 1. NULL for queries without data.
#'
#' @return Character vector of cache record keys
#'
#' @importFrom digest sha1_digest
#' @importFrom magrittr %>%
#' @export
omnipath_cache_key <- function(url, post = NULL, payload = NULL){

    url %>%
    omnipath_url_to_list(post, payload) %>%
    pmap_chr(sha1_digest)

}


#' Makes sure the post and the payload is the same length as url
#'
#' @param url Character vector with URLs
#' @param post List with the HTTP POST parameters or a list of lists if
#' the url vector is longer than 1. NULL for queries without POST parameters.
#' @param payload HTTP data payload. List with multiple items if the url
#' vector is longer than 1. NULL for queries without data.
#'
#' @return List of 3 lists or vectors: url, post and payload
#' @importFrom magrittr %<>% %>%
#' @importFrom purrr pmap_chr map map_lgl
omnipath_url_to_list <- function(url, post = NULL, payload = NULL){

    null_list <- function(value){

        # if the url is a character vector of length 1
        # and the post is a list of parameters
        if(
            length(url) == 1 &&
            is.list(value) &&
            value %>% map_lgl(is.list) %>% any %>% `!`
        ){
            value %<>% list
        }

        `if`(
            is.null(value),
            url %>% length %>% seq %>% map(function(x){NULL}),
            value
        )

    }

    if(length(url) != length(post) && !is.null(post)){
        post %<>% list
    }
    if(length(url) != length(payload) && !is.null(post)){
        post %<>% list
    }

    post %<>% null_list
    payload %<>% null_list

    list(url, post, payload)

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