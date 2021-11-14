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
#' @importFrom logger log_info
#'
#' @noRd
omnipath_init_cache <- function(){

    cachedir <- omnipath_default_cachedir()
    omnipath_new_cachedir(cachedir)
    options(omnipath.cachedir = cachedir)
    omnipath_read_cache_db()
    log_info('Initialized cache: `%s`.', cachedir)

    invisible(cachedir)

}


#' Returns the default cache directory path
#'
#' @importFrom magrittr %>%
#' @importFrom rappdirs user_cache_dir
#'
#' @noRd
omnipath_default_cachedir <- function(){

    user_cache_dir() %>%
    file.path('OmnipathR') %>%
    normalizePath(mustWork = FALSE) %>%
    suppressWarnings()

}


#' Sets up a new cache directory
#'
#' @importFrom jsonlite toJSON
#' @importFrom logger log_info
#'
#' @noRd
omnipath_new_cachedir <- function(path){

    if(!omnipath_is_cachedir(path)){

        logger::log_info('Setting up new cache directory `%s`.', path)
        dir.create(path, showWarnings = FALSE, recursive = TRUE)
        omnipath.env$cache <- list()
        write(
            jsonlite::toJSON(omnipath.env$cache, pretty = TRUE, null = 'null'),
            file.path(path, 'cache.json')
        )

    }

    invisible(path)

}


#' Change the cache directory
#'
#' @param path Character: path to the new cache directory. If don't exist,
#'     the directories will be created. If the path is an existing cache
#'     directory, the package's cache database for the current session will
#'     be loaded from the database in the directory. If \code{NULL}, the
#'     cache directory will be set to its default path.
#'
#' @return Returns \code{NULL}.
#'
#' @examples
#' tmp_cache <- tempdir()
#' omnipath_set_cachedir(tmp_cache)
#' # restore the default cache directory:
#' omnipath_set_cachedir()
#'
#' @importFrom magrittr %<>%
#' @export
omnipath_set_cachedir <- function(path = NULL){

    path %<>% if_null(omnipath_default_cachedir())
    options(omnipath.cachedir = path)
    omnipath_new_cachedir(path)
    omnipath_read_cache_db()
}


#' Tells if a directory looks like an OmnipathR cache directory
#'
#' @noRd
omnipath_is_cachedir <- function(path){

    dir.exists(path) && file.exists(file.path(path, 'cache.json'))

}


#' Determines the path of the cache directory according to the settings
#'
#' @noRd
omnipath_get_cachedir <- function(){

    cachedir <- options('omnipath.cachedir')[[1]]

    if(is.null(cachedir) || !omnipath_is_cachedir(cachedir)){

        cachedir <- omnipath_init_cache()

    }

    return(cachedir)

}


#' Creates a lock file in the cache directory in order to avoid simulatneous
#' reading and writing of the cache file database
#'
#' @noRd
omnipath_lock_cache_db <- function(){


    lockfile <- omnipath_cache_lock_path()

    for(i in seq(options('omnipath.cache_timeout')[[1]])){

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
#' @return Logical: returns TRUE if the cache was locked and now is
#'     unlocked; FALSE if it was not locked.
#'
#' @examples
#' omnipath_unlock_cache_db()
#'
#' @importFrom magrittr %>%
#' @export
omnipath_unlock_cache_db <- function(){

    omnipath_cache_lock_path() %>%
    {suppressWarnings(file.remove(.))}

}


#' Tells the path to the cache lock file
#'
#' @importFrom magrittr %>%
#'
#' @noRd
omnipath_cache_lock_path <- function(){

    omnipath_get_cachedir() %>%
    file.path('cache.lock')

}


#' Tells the path to the cache database
#'
#' @importFrom magrittr %>%
#'
#' @noRd
omnipath_cache_db_path <- function(){

    omnipath_get_cachedir() %>%
    file.path('cache.json')

}


#' Throws a fatal error about locked cache.
#'
#' @noRd
omnipath_locked_cache_error <- function(){

    msg <- sprintf(
        paste0(
            'The cache directory `%s` is locked for longer than %d seconds. ',
            'If you keep experiencing this issue please call `omnipath_',
            'unlock_cache_db()`. '
        ),
        omnipath_get_cachedir(),
        options('omnipath.cache_timeout')[[1]]
    )
    logger::log_fatal(msg)
    stop(msg)

}


#' Executes a function with locking the cache database
#'
#' @noRd
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
#' @param pattern String or regular expression.
#' @param ... Passed to \code{grep}
#'
#' @return List of cache records matching the pattern.
#'
#' @examples
#' # find all cache records from the BioPlex database
#' bioplex_cache_records <- omnipath_cache_search(
#'     'bioplex',
#'     ignore.case = TRUE
#' )
#'
#' @importFrom purrr map_chr
#' @export
omnipath_cache_search <- function(pattern, ...){

    omnipath.env$cache %>%
    map_chr('url') %>%
    grep(pattern = pattern, ...) %>%
    `[`(names(omnipath.env$cache), .) %>%
    `[`(omnipath.env$cache, .)


}


#' Removes contents from the cache directory
#'
#' According to the parameters, it can remove contents older than a certain
#' age, or contents having a more recent version, one specific item, or
#' wipe the entire cache.
#'
#' @usage omnipath_cache_remove(key = NULL, url = NULL, post = NULL,
#'     payload = NULL, max_age = NULL, min_age = NULL, status = NULL,
#'     only_latest = FALSE, wipe = FALSE, autoclean = TRUE)
#'
#' @param key The key of the cache record
#' @param url URL pointing to the resource
#' @param post HTTP POST parameters as a list
#' @param payload HTTP data payload
#' @param max_age Age of cache items in days. Remove everything that is older
#'     than this age
#' @param min_age Age of cache items in days. Remove everything more recent
#'     than this age
#' @param status Remove items having any of the states listed here
#' @param only_latest Keep only the latest version
#' @param autoclean Remove the entries about failed downloads, the files in
#'     the cache directory which are missing from the cache database, and the
#'     entries without existing files in the cache directory
#' @param wipe Logical: if TRUE, removes all files from the cache and the
#'     cache database. Same as calling \code{\link{omnipath_cache_wipe}}.
#'
#' @return Invisibly returns the cache database (list of cache records).
#'
#' @examples
#' \dontrun{
#' # remove all cache data from the BioPlex database
#' cache_records <- omnipath_cache_search(
#'     'bioplex',
#'     ignore.case = TRUE
#' )
#' omnipath_cache_remove(names(cache_records))
#'
#' # remove a record by its URL
#' regnetwork_url <- 'http://www.regnetworkweb.org/download/human.zip'
#' omnipath_cache_remove(url = regnetwork_url)
#'
#' # remove all records older than 30 days
#' omnipath_cache_remove(max_age = 30)
#'
#' # for each record, remove all versions except the latest
#' omnipath_cache_remove(only_latest = TRUE)
#' }
#'
#' bioc_url <- 'https://bioconductor.org/'
#' version <- omnipath_cache_latest_or_new(url = bioc_url)
#' httr::GET(bioc_url, httr::write_disk(version$path, overwrite = TRUE))
#' omnipath_cache_download_ready(version)
#' key <- omnipath_cache_key(bioc_url)
#' omnipath_cache_remove(key = key)
#'
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{omnipath_cache_wipe}}}
#'     \item{\code{\link{omnipath_cache_clean}}}
#'     \item{\code{\link{omnipath_cache_autoclean}}}
#' }
omnipath_cache_remove <- function(
    key = NULL,
    url = NULL,
    post = NULL,
    payload = NULL,
    max_age = NULL,
    min_age = NULL,
    status = NULL,
    only_latest = FALSE,
    wipe = FALSE,
    autoclean = TRUE
){

    do.call(.omnipath_cache_remove, as.list(environment()))

}


#' Workaround due to missing R CMD check option for ignoring codoc mismatches.
#'
#' We want to suppress warnings like "Codoc mismatches from documentation
#' object" because these are only a byproduct of the function is wrapped into
#' a decorator. Would be nicer if Roxygen and R CMD check had an option to
#' ignore certain arguments in the codoc mismatch check but unfortunately
#' they don't have.
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom purrr map keep map_lgl discard
#'
#' @noRd
.omnipath_cache_remove <- cache_locked %@% function(
    key = NULL,
    url = NULL,
    post = NULL,
    payload = NULL,
    max_age = NULL,
    min_age = NULL,
    status = NULL,
    only_latest = FALSE,
    wipe = FALSE,
    autoclean = FALSE
){

    if(wipe){

        omnipath.env$cache <- list()
        omnipath_cache_clean()
        return(invisible(NULL))

    }

    key %<>% ensure_list
    url %<>% ensure_list

    key <- ifelse(
        map_lgl(key, is.null) & !map_lgl(url, is.null),
        omnipath_cache_key(url = url, post = post, payload = payload),
        key
    )

    key %<>% list_null
    url %<>% list_null

    has_key <- !is.null(key)
    has_version_filtering <- (
        !is.null(max_age) ||
        !is.null(min_age) ||
        !is.null(status) ||
        only_latest
    )

    omnipath.env$cache %<>%
    {`if`(
        has_key,
        `if`(
            has_version_filtering,
            `[`(., key %>% unlist),
            `[`(., omnipath.env$cache %>% names %>% setdiff(key))
        ),
        .
    )} %>%
    map(
        omnipath_cache_keep_versions,
        max_age = max_age,
        min_age = min_age,
        status = status,
        only_latest = only_latest
    ) %>%
    {`if`(
        has_key && has_version_filtering,
        map(
            omnipath.env$cache,
            function(record){
                `if`(
                    record$key %in% names(.),
                    .[[record$key]],
                    record
                )
            }
        ),
        .
    )} %>%
    discard(
        function(record){
            length(record$versions) == 0L
        }
    )

    omnipath_cache_clean()

    if(autoclean){

        omnipath_write_cache_db()
        omnipath_unlock_cache_db()

        omnipath_cache_clean_db()
        omnipath_cache_autoclean()

        omnipath_lock_cache_db()

    }

    invisible(omnipath.env$cache)

}


#' Removes version items from a cache record
#'
#' @param record A cache record
#' @param max_age Age of cache items in days. Remove everything that is older
#' than this age
#' @param min_age Age of cache items in days. Remove everything more recent
#' than this age
#' @param status Keep items having any of the states listed here
#' @param only_latest Keep only the latest version
#'
#' @return A cache record with the version items removed.
#'
#' @importFrom magrittr %<>%
#'
#' @noRd
omnipath_cache_keep_versions <- function(
    record,
    max_age = NULL,
    min_age = NULL,
    status = NULL,
    only_latest = FALSE
){

    record$versions %<>%
        `[`(
            omnipath_cache_filter_versions(
                record = record,
                latest = only_latest,
                max_age = max_age,
                min_age = min_age,
                status = status
            )
        )

    invisible(record)

}


#' Removes the cache database entries without existing files
#'
#' @param ... Ignored.
#'
#' @return Returns `NULL`.
#'
#' @examples
#' omnipath_cache_clean_db()
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map keep discard
#' @importFrom logger log_trace
#'
#' @export
omnipath_cache_clean_db <- cache_locked %@% function(){

    omnipath.env$cache <-
        omnipath.env$cache %>%
        map(
            function(record){
                record$versions %<>%
                keep(
                    function(v){file.exists(v$path)}
                )
                record
            }
        ) %>%
        discard(
            function(record){
                no_versions <- length(record$versions) == 0L
                if(no_versions){
                    log_trace(
                        'Removing cache entry %s (%s): no files available.',
                        record$key,
                        record$url
                    )
                }
                return(no_versions)
            }
        )

}


#' Permanently removes all the cache contents
#'
#' After this operation the cache directory will be completely empty,
#' except an empty cache database file.
#'
#' @param ... Ignored.
#'
#' @return Returns `NULL`.
#'
#' @examples
#' \dontrun{
#' omnipath_cache_wipe()
#' # the cache is completely empty:
#' print(omnipath.env$cache)
#' # list()
#' list.files(omnipath_get_cachedir())
#' # [1] "cache.json"
#' }
#'
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom logger log_success
#' @seealso \code{\link{omnipath_cache_remove}}
omnipath_cache_wipe <- cache_locked %@% function(){

    omnipath_get_cachedir() %T>%
    {log_success('Removing all cache contents from `%s`.', .)} %>%
    list.files() %>%
    setdiff(PROTECTED_FILES) %>%
    file.path(omnipath_get_cachedir(), .) %>%
    file.remove()

    omnipath.env$cache <- list()
    invisible(NULL)

}


#' Removes the items from the cache directory which are unknown by the cache
#' database
#'
#' @return Returns `NULL`.
#'
#' @examples
#' omnipath_cache_clean()
#'
#' @importFrom purrr map_chr map walk
#' @importFrom magrittr %>%
#' @importFrom logger log_trace
#' @export
omnipath_cache_clean <- function(){

    files_in_db <-
        omnipath.env$cache %>%
        map(~map_chr(.x$versions, 'path')) %>%
        unlist() %>%
        empty_no_problem(basename)

    omnipath_get_cachedir() %>%
    list.files() %>%
    setdiff(files_in_db) %>%
    setdiff(PROTECTED_FILES) %>%
    file.path(omnipath_get_cachedir(), .) %>%
    walk(log_trace, fmt = 'Removing cache file `%s`.') %>%
    file.remove()

    invisible(NULL)

}


#' Keeps only the latest versions of complete downloads
#'
#' Removes the old versions, the failed downloads and the files in the cache
#' directory which are missing from the database. For more flexible
#' operations use \code{\link{omnipath_cache_remove}} and
#' \code{\link{omnipath_cache_clean}}.
#'
#' @return Invisibl returns the cache database (list of cache records).
#'
#' @examples
#' \dontrun{
#' omnipath_cache_autoclean()
#' }
#'
#' @export
omnipath_cache_autoclean <- function(){

    omnipath_cache_remove(
        only_latest = TRUE,
        status = CACHE_STATUS$READY,
        autoclean = FALSE
    )
    invisible(omnipath.env$cache)

}


#' Retrieves one item from the cache directory
#'
#' @param key The key of the cache record
#' @param url URL pointing to the resource
#' @param post HTTP POST parameters as a list
#' @param payload HTTP data payload
#' @param create Create a new entry if doesn't exist yet
#' @param ... Passed to \code{omnipath_cache_record} (internal function)
#'
#' @return Cache record: an existing record if the entry already exists,
#'     otherwise a newly created and inserted record
#'
#' @examples
#' # create an example cache record
#' bioc_url <- 'https://bioconductor.org/'
#' version <- omnipath_cache_latest_or_new(url = bioc_url)
#' omnipath_cache_remove(url = bioc_url) # cleaning up
#'
#' # retrieve the cache record
#' record <- omnipath_cache_get(url = bioc_url)
#' record$key
#' # [1] "41346a00fb20d2a9df03aa70cf4d50bf88ab154a"
#' record$url
#' # [1] "https://bioconductor.org/"
#'
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
    # omnipath.env$cache <- get('omnipath.env$cache', envir = .GlobalEnv)

    if(!(key %in% names(omnipath.env$cache)) && create){

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
        # omnipath.env$cache <- get('omnipath.env$cache', envir = .GlobalEnv)

    }

    record <- `if`(
        key %in% names(omnipath.env$cache),
        omnipath.env$cache[[key]],
        NULL
    )

    return(record)

}


#' The latest or a new version of a cache record
#'
#' Looks up a record in the cache and returns its latest valid version. If
#' the record doesn't exist or no valid version available, creates a new one.
#'
#' @param key The key of the cache record
#' @param url URL pointing to the resource
#' @param post HTTP POST parameters as a list
#' @param payload HTTP data payload
#' @param create Logical: whether to create and return a new version. If
#'     FALSE only the latest existing valid version is returned, if available.
#' @param ... Passed to \code{\link{omnipath_cache_get}}
#'
#' @return A cache version item.
#'
#' @examples
#' \dontrun{
#' # retrieve the latest version of the first cache record
#' # found by the search keyword "bioplex"
#' latest_bioplex <-
#'     omnipath_cache_latest_or_new(
#'         names(omnipath_cache_search('bioplex'))[1]
#'     )
#'
#' latest_bioplex$dl_finished
#' # [1] "2021-03-09 14:28:50 CET"
#' latest_bioplex$path
#' # [1] "/home/denes/.cache/OmnipathR/378e0def2ac97985f629-1.rds"
#' }
#'
#' # create an example cache record
#' bioc_url <- 'https://bioconductor.org/'
#' version <- omnipath_cache_latest_or_new(url = bioc_url)
#' omnipath_cache_remove(url = bioc_url) # cleaning up
#'
#' @importFrom rlang %||%
#' @importFrom logger log_info
#' @export
omnipath_cache_latest_or_new <- function(
    key = NULL,
    url = NULL,
    post = NULL,
    payload = NULL,
    create = TRUE,
    ...
){

    record <-
        omnipath_cache_get(
            key = key,
            url = url,
            post = post,
            payload = payload,
            ...
        )

    version <- omnipath_cache_latest_version(record)

    logger::log_info(
        'Looking up in cache `%s`: key=%s, %s.',
        url,
        record$key,
        `if`(
            is.null(version),
            'no version available',
            sprintf('latest version=%s', version)
        )
    )

    if(is.null(version) && create){

        version <-
            cache_locked(
                omnipath_cache_new_version
            )(record$key)

        logger::log_info(
            'Created new version for cache record %s: version %s.',
            record$key, version
        )

    }

    if(!is.null(version)){

        version <- omnipath.env$cache[[record$key]]$versions[[version]]

    }

    logger::log_trace('Cache file path: %s', version$path)

    return(version)

}


#' Loads an R object from the cache
#'
#' Loads the object from RDS format.
#'
#' @param key Key of the cache item
#' @param url URL of the downloaded resource
#' @param post HTTP POST parameters as a list
#' @param payload HTTP data payload
#' @param version Version of the cache item. If does not exist or NULL, the
#'     latest version will be retrieved
#'
#' @return Object loaded from the cache RDS file.
#'
#' @examples
#' url <- paste0(
#'     'https://omnipathdb.org/intercell?resources=Adhesome,Almen2009,',
#'     'Baccin2019,CSPA,CellChatDB&license=academic'
#' )
#' result <- read.delim(url, sep = '\t')
#' omnipath_cache_save(result, url = url)
#' # works only if you have already this item in the cache
#' intercell_data <- omnipath_cache_load(url = url)
#' class(intercell_data)
#' # [1] "data.frame"
#' nrow(intercell_data)
#' # [1] 16622
#' attr(intercell_data, 'origin')
#' # [1] "cache"
#'
#' # basic example of saving and loading to and from the cache:
#' bioc_url <- 'https://bioconductor.org/'
#' bioc_html <- readChar(url(bioc_url), nchars = 99999)
#' omnipath_cache_save(bioc_html, url = bioc_url)
#' bioc_html <- omnipath_cache_load(url = bioc_url)
#'
#' @importFrom logger log_info log_trace log_fatal
#' @importFrom rlang %||%
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
        log_info('Cache record does not exist: `%s`', key %||% url)
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

        if(!file.exists(path)){
            msg <- sprintf(
                paste0(
                    'Missing cache file: `%s`. Please run ',
                    '`omnipath_cache_clean_db()` and try again.'
                ),
                path
            )
            log_fatal(msg)
            stop(msg)
        }

        data <- readRDS(path)
        log_trace('Loaded data from RDS `%s`.', path)
        return(data %>% origin_cache)

    }else{
        log_info(
            'No version is available for key `%s` (%s).',
            record$key,
            url
        )
    }

}


#' Saves an R object to the cache
#'
#' Exports the object in RDS format, creates new cache record if necessary.
#'
#' @return Returns invisibly the data itself.
#'
#' @param data An object
#' @param key Key of the cache item
#' @param url URL of the downloaded resource
#' @param post HTTP POST parameters as a list
#' @param payload HTTP data payload
#' @param version Version of the cache item. If does not exist a new version
#' item will be created
#'
#' @return Invisibly returns the `data`.
#'
#' @examples
#' mydata <- data.frame(a = c(1, 2, 3), b = c('a', 'b', 'c'))
#' omnipath_cache_save(mydata, url = 'some_dummy_address')
#' from_cache <- omnipath_cache_load(url = 'some_dummy_address')
#' from_cache
#' #   a b
#' # 1 1 a
#' # 2 2 b
#' # 3 3 c
#' attr(from_cache, 'origin')
#' # [1] "cache"
#'
#' # basic example of saving and loading to and from the cache:
#' bioc_url <- 'https://bioconductor.org/'
#' bioc_html <- readChar(url(bioc_url), nchars = 99999)
#' omnipath_cache_save(bioc_html, url = bioc_url)
#' bioc_html <- omnipath_cache_load(url = bioc_url)
#'
#' @importFrom logger log_trace
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
    key <- record$key

    version <- omnipath_cache_update_status(
        key = key,
        version = version,
        status = CACHE_STATUS$STARTED
    )

    target_path <- omnipath.env$cache[[key]]$versions[[version]]$path

    log_trace('Exporting object to RDS: `%s`.', target_path)
    saveRDS(data, target_path)
    log_trace('Exported RDS to `%s`.', target_path)

    omnipath_cache_download_ready(version, key = key)

    invisible(data)

}


#' Moves an existing file into the cache
#'
#' Either the key or the URL (with POST and payload) must be provided.
#'
#' @param path Path to the source file
#' @param key Key of the cache item
#' @param url URL of the downloaded resource
#' @param post HTTP POST parameters as a list
#' @param payload HTTP data payload
#' @param version Version of the cache item. If does not exist a new version
#'     item will be created
#' @param keep_original Whether to keep or remove the original file
#'
#' @return Character: invisibly returns the version number of the cache
#'     version item.
#'
#' @examples
#' path <- tempfile()
#' saveRDS(rnorm(100), file = path)
#' omnipath_cache_move_in(path, url = 'the_download_address')
#'
#' # basic example of moving a file to the cache:
#'
#' bioc_url <- 'https://bioconductor.org/'
#' html_file <- tempfile(fileext = '.html')
#' httr::GET(bioc_url, httr::write_disk(html_file, overwrite = TRUE))
#' omnipath_cache_move_in(path = html_file, url = bioc_url)
#' omnipath_cache_remove(url = bioc_url) # cleaning up
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

    record <- omnipath_cache_get(
        key = key,
        url = url,
        post = post,
        payload = payload,
        ext = file_extension(path)
    )
    key <- record$key

    version <- omnipath_cache_update_status(
        key = key,
        version = version,
        status = CACHE_STATUS$STARTED
    )

    target_path <- omnipath.env$cache[[key]]$versions[[version]]$path

    file.copy(path, target_path)
    logger::log_trace('Copied `%s` to `%s`.', path, target_path)

    if(!keep_original){

        file.remove(path)

    }

    omnipath_cache_download_ready(version, key = key)

}


#' Sets the download status to ready for a cache item
#'
#' @param key Key of the cache item
#' @param version Version of the cache item. If does not exist a new version
#' item will be created
#'
#' @return Character: invisibly returns the version number of the cache
#'     version item.
#'
#' @examples
#' bioc_url <- 'https://bioconductor.org/'
#' # request a new version item (or retrieve the latest)
#' new_version <- omnipath_cache_latest_or_new(url = bioc_url)
#' # check if the version item is not a finished download
#' new_version$status
#' # [1] "unknown"
#' # download the file
#' httr::GET(bioc_url, httr::write_disk(new_version$path, overwrite = TRUE))
#' # report to the cache database that the download is ready
#' omnipath_cache_download_ready(new_version)
#' # now the status is ready:
#' version <- omnipath_cache_latest_or_new(url = bioc_url)
#' version$status
#' # "ready"
#' version$dl_finished
#' # [1] "2021-03-09 16:48:38 CET"
#' omnipath_cache_remove(url = bioc_url) # cleaning up
#'
#' @importFrom rlang %||%
#' @importFrom logger log_info
#' @export
omnipath_cache_download_ready <- function(version, key = NULL){

    key <- key %||% omnipath_cache_key_from_version(version)
    version <- `if`(is.list(version), version$number, version)

    logger::log_info('Download ready [key=%s, version=%s]', key, version)

    omnipath_cache_update_status(
        key = key,
        version = version,
        status = 'ready',
        dl_finished = Sys.time()
    ) %>%
    invisible()

}


#' From a cache version item extracts the cache record key
#'
#' @importFrom magrittr %>%
#'
#' @noRd
omnipath_cache_key_from_version <- function(version){

    version$path %>%
    basename() %>%
    substr(0, 40)


}


#' Updates the status of an existing cache record
#'
#' @usage omnipath_cache_update_status(key, version, status,
#'     dl_finished = NULL)
#'
#' @param key Key of the cache item
#' @param version Version of the cache item. If does not exist a new version
#'     item will be created
#' @param status The updated status value
#' @param dl_finished Timestamp for the time when download was finished,
#'     if `NULL` the value remains unchanged
#'
#' @return Character: invisibly returns the version number of the cache
#'     version item.
#'
#' @examples
#' bioc_url <- 'https://bioconductor.org/'
#' latest_version <- omnipath_cache_latest_or_new(url = bioc_url)
#' key <- omnipath_cache_key(bioc_url)
#' omnipath_cache_update_status(
#'     key = key,
#'     version = latest_version$number,
#'     status = 'ready',
#'     dl_finished = Sys.time()
#' )
#' omnipath_cache_remove(url = bioc_url) # cleaning up
#'
#' @export
omnipath_cache_update_status <- function(
    key,
    version,
    status,
    dl_finished = NULL
){

    do.call(.omnipath_cache_update_status, as.list(environment()))

}


#' Workaround due to missing R CMD check option for ignoring codoc mismatches.
#'
#' We want to suppress warnings like "Codoc mismatches from documentation
#' object" because these are only a byproduct of the function is wrapped into
#' a decorator. Would be nicer if Roxygen and R CMD check had an option to
#' ignore certain arguments in the codoc mismatch check but unfortunately
#' they don't have.
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom logger log_info log_warn
#'
#' @noRd
.omnipath_cache_update_status <- cache_locked %@% function(
    key,
    version,
    status,
    dl_finished = NULL
){

    version %<>%
        {`if`(
            is.null(.),
            .,
            as.character(.)
        )}

    if(key %in% names(omnipath.env$cache)){

        if(
            is.null(version) ||
            !(version %in% names(omnipath.env$cache[[key]]$versions))
        ){

            version <- omnipath_cache_new_version(key, version = version)

        }

        old_status <- omnipath.env$cache[[key]]$versions[[version]]$status
        omnipath.env$cache[[key]]$versions[[version]]$status <- status

        if(
            old_status != CACHE_STATUS$STARTED &&
            status == CACHE_STATUS$STARTED
        ){
            omnipath.env$cache[[key]]$versions[[version]]$dl_started <-
                Sys.time()
        }

        if(!is.null(dl_finished)){
            omnipath.env$cache[[key]]$versions[[version]]$dl_finished <-
                dl_finished
        }
        # omnipath.env$cache <- omnipath.env$cache

        logger::log_info(
            'Cache item `%s` version %s: status changed from `%s` to `%s`.',
            key,
            version,
            old_status,
            status
        )

    }else{

        logger::log_warn(
            'Failed to update cache: key `%s` does not exist.',
            key
        )

    }

    invisible(version)

}


#' Sets the file extension for a cache record
#'
#' @usage omnipath_cache_set_ext(key, ext)
#'
#' @param key Character: key for a cache item, alternatively a version entry.
#' @param ext Character: the file extension, e.g. "zip".
#'
#' @return Returns `NULL`.
#'
#' @examples
#' bioc_url <- 'https://bioconductor.org/'
#' version <- omnipath_cache_latest_or_new(url = bioc_url)
#' version$path
#' # [1] "/home/denes/.cache/OmnipathR/41346a00fb20d2a9df03-1"
#' httr::GET(bioc_url, httr::write_disk(version$path, overwrite = TRUE))
#' key <- omnipath_cache_key(url = bioc_url)
#' omnipath_cache_set_ext(key = key, ext = 'html')
#' version <- omnipath_cache_latest_or_new(url = bioc_url)
#' version$path
#' # [1] "/home/denes/.cache/OmnipathR/41346a00fb20d2a9df03-1.html"
#' record <- omnipath_cache_get(url = bioc_url)
#' record$ext
#' # [1] "html"
#' omnipath_cache_remove(url = bioc_url) # cleaning up
#'
#' @export
omnipath_cache_set_ext <- function(key, ext){

    do.call(.omnipath_cache_set_ext, as.list(environment()))

}


#' Workaround due to missing R CMD check option for ignoring codoc mismatches.
#'
#' We want to suppress warnings like "Codoc mismatches from documentation
#' object" because these are only a byproduct of the function is wrapped into
#' a decorator. Would be nicer if Roxygen and R CMD check had an option to
#' ignore certain arguments in the codoc mismatch check but unfortunately
#' they don't have.
#'
#' @return Returns `NULL`.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr map
#' @importFrom tools file_path_sans_ext
#'
#' @noRd
.omnipath_cache_set_ext <- cache_locked %@% function(key, ext){

    if(is.list(key)){

        key %<>% omnipath_cache_key_from_version

    }

    if((key %in% names(omnipath.env$cache))){

        omnipath.env$cache[[key]]$ext <- ext

        omnipath.env$cache[[key]]$versions %<>%
            map(
                function(version){

                    old_path <- version$path

                    version$path %<>%
                        file_path_sans_ext(compression = TRUE) %>%
                        paste(ext, sep = '.')

                    file.rename(old_path, version$path)

                    return(version)

                }
            )

    }

}


#' Adds a new version item to an existing cache record
#'
#' @return Character: key of the version item.
#'
#' @noRd
omnipath_cache_new_version <- function(key, version = NULL){

    version <- `if`(
        is.null(version),
        omnipath_cache_next_version(key),
        version
    )

    record <- omnipath.env$cache[[key]]

    new_version <- omnipath_cache_version(
        record,
        version = version,
        dl_started = NULL,
        status = CACHE_STATUS$UNKNOWN
    )

    omnipath.env$cache[[key]]$versions[[version]] <- new_version
    # omnipath.env$cache <- omnipath.env$cache

    invisible(version)

}


#' Tells the next incremental version number for a cache record
#'
#' @return The next version number as character
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @noRd
omnipath_cache_next_version <- function(key){

    omnipath.env$cache[[key]]$versions %>%
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
#' @param record A cache record
#'
#' @return Character: the version ID with the most recent download finished
#'     time
#'
#' @export
omnipath_cache_latest_version <- function(record){

    record %>%
    null_or_call(omnipath_cache_filter_versions, latest = TRUE) %>%
    empty_no_problem(identity)

}

#' Filters the versions from one cache record
#'
#' Filters the versions based on multiple conditions: their age and status
#'
#' @param record A cache record
#' @param latest Return the most recent version
#' @param max_age The maximum age in days (e.g. 5: 5 days old or more recent)
#' @param min_age The minimum age in days (e.g. 5: 5 days old or older)
#' @param status Character vector with status codes. By default only the
#' versions with `ready` (completed download) status are selected
#'
#' @return Character vector with version IDs, NA if no version satisfies the
#'     conditions.
#'
#' @examples
#' # creating an example cache record
#' bioc_url <- 'https://bioconductor.org/'
#' version <- omnipath_cache_latest_or_new(url = bioc_url)
#' httr::GET(bioc_url, httr::write_disk(version$path, overwrite = TRUE))
#' omnipath_cache_download_ready(version)
#' record <- dplyr::first(omnipath_cache_search('biocond'))
#'
#' # only the versions with status "ready"
#' version_numbers <- omnipath_cache_filter_versions(record, status = 'ready')
#' omnipath_cache_remove(url = bioc_url) # cleaning up
#'
#' @importFrom magrittr %<>% %>% and
#' @importFrom purrr map map_chr map_lgl
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

    versions <- record$versions
    statuses <- versions %>% map_chr('status')
    selection <- statuses %in% status

    if(!is.null(max_age)){
        max_age <- days_ago(max_age)
        selection %<>% and(which_dl_finished(versions, max_age))
    }

    if(!is.null(min_age)){
        min_age <- days_ago(min_age)
        selection %<>% and(which_dl_finished(versions, min_age, op = `<=`))
    }

    if(latest){
        t_latest <-
            versions %>%
            map('dl_finished') %>%
            unlist() %>%
            empty_no_problem(max)
        selection %<>% and(which_dl_finished(versions, t_latest, op = `==`))
    }

    versions %>%
    names %>%
    `[`(selection)

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
#'
#' @noRd
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
#' @noRd
omnipath_cache_add <- cache_locked %@% function(record, new = FALSE){

    if(record$key %in% omnipath.env$cache && !new){

        omnipath.env$cache[[record$key]]$versions <-
            merge_lists(
                omnipath.env$cache[[record$key]]$versions,
                record$versions
            )

    }else{

        omnipath.env$cache[[record$key]] <- record

    }

}


#' Returns the key if either the key or the URL is available
#'
#' @importFrom logger log_fatal
#'
#' @noRd
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
#'
#' @noRd
omnipath_cache_record <- function(
        key,
        url,
        ext = NULL,
        post = NULL,
        payload = NULL
    ){

    if(is.null(ext)) ext <- file_extension(url)

    post %<>% map(char_shorten)
    payload %<>% char_shorten(250)

    list(
        key = key,
        url = url,
        post = post,
        payload = payload,
        ext = ext,
        versions = list()
    )

}


#' In a character vector, shortens values longer than maxlen and adds a note
#' about the truncation and the original length
#'
#' @noRd
char_shorten <- function(value, maxlen = 50){

    ifelse(
        is.character(value) & nchar(value) > maxlen,
        sprintf(
            '%s... (truncated, full length: %d)',
            substr(value, 0, maxlen),
            nchar(value)
        ),
        value
    )

}


#' Creates a new version item for a cache record
#'
#' @return A version item
#'
#' @importFrom magrittr %>% %<>%
#'
#' @noRd
omnipath_cache_version <- function(
        record,
        version = '1',
        dl_started = Sys.time(),
        status = CACHE_STATUS$STARTED
    ){

    version %<>% as.character

    path <-
        omnipath_get_cachedir() %>%
        file.path(
            sprintf('%s-%s', record$key, version)
        ) %>%
        file_add_extension(record$ext)

    list(
        number = version,
        path = path,
        dl_started = dl_started,
        dl_finished = NULL,
        status = status
    )

}


#' Generates a hash which identifies an element in the cache database
#'
#' @param url Character vector with URLs
#' @param post List with the HTTP POST parameters or a list of lists if
#'     the url vector is longer than 1. NULL for queries without POST
#'     parameters.
#' @param payload HTTP data payload. List with multiple items if the url
#'     vector is longer than 1. NULL for queries without data.
#'
#' @return Character vector of cache record keys.
#'
#' @examples
#' bioc_url <- 'https://bioconductor.org/'
#' omnipath_cache_key(bioc_url)
#' # [1] "41346a00fb20d2a9df03aa70cf4d50bf88ab154a"
#'
#' @importFrom digest sha1_digest
#' @importFrom magrittr %>%
#' @importFrom purrr pmap
#' @export
omnipath_cache_key <- function(url, post = NULL, payload = NULL){

    url %>%
    omnipath_url_to_list(post, payload) %>%
    pmap(list) %>%
    map_chr(digest::sha1_digest)

}


#' Makes sure the post and the payload is the same length as url
#'
#' @param url Character vector with URLs
#' @param post List with the HTTP POST parameters or a list of lists if
#'     the url vector is longer than 1. NULL for queries without POST
#'     parameters.
#' @param payload HTTP data payload. List with multiple items if the url
#'     vector is longer than 1. NULL for queries without data.
#'
#' @return List of 3 lists or vectors: url, post and payload
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom purrr pmap_chr map map_lgl
#'
#' @noRd
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
    if(length(url) != length(payload) && !is.null(payload)){
        payload %<>% list
    }

    post %<>% null_list
    payload %<>% null_list

    list(url, post, payload)

}


#' Reads the cache DB contents from the disk to the memory
#'
#' @importFrom jsonlite fromJSON
#'
#' @noRd
omnipath_read_cache_db <- function(){

    omnipath.env$cache <-
        omnipath_cache_db_path() %>%
        jsonlite::fromJSON() %>%
        omnipath_cache_timestamps()

}


#' Converts the timestamps read from JSON to POSIXct objects
#'
#' @importFrom magrittr %<>%
#'
#' @noRd
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
#'
#' @noRd
omnipath_write_cache_db <- function(){

    omnipath.env$cache %>%
    jsonlite::toJSON(pretty = TRUE, null = 'null') %>%
    write(omnipath_cache_db_path())

}


#' Adds an attribute about cache origin
#'
#' Creates the attribute `origin` for a data frame (or other object) and
#' assigns the value "cache".
#'
#' @param data Data frame or any kind of object.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
origin_cache <- function(data, from_cache = TRUE){

    data %>%
    {`if`(from_cache, `attr<-`(., 'origin', 'cache'), .)}

}