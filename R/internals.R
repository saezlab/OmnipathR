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


READR_DEFAULTS = list(
    col_types = cols(),
    show_col_types = FALSE,
    progress = FALSE
)

CURL_DEBUG_TYPES = list(
    '* ',   # info
    '<=',  # hdr in
    '=>',  # hdr out
    '< ',   # data in
    ' >',   # data out
    '<*',  # tls data in
    '*>'   # tls data out
)


#' Prepend the current OmniPath server domain to an URL
#'
#' @param path_qs Character: part of the URL after the domain: the path and
#'     the query string.
#' @param notls Logical: use http instead of https (do not use TLS).
#'
#' @importFrom magrittr %>% or
#' @importFrom stringr str_replace
#' @noRd
omnipath_url <- function(path_qs, notls = FALSE){

    'omnipathr.notls_force' %>%
    getOption %>%
    or(notls) %>%
    `if`('notls_', '') %>%
    sprintf('omnipathr.%surl', .) %>%
    getOption %>%
    str_replace('/+$', '') %>%
    c(str_replace(path_qs, '^/+', '')) %>%
    paste(collapse = '/')

}


#' Retrieves an URL from the package's URL register
#'
#' @param key Character: name of the option containing the URL
#' @param key_param List: variables to insert into the `key`.
#'
#' @importFrom rlang exec !!!
#'
#' @noRd
get_url <- function(key, param = list()){

    key %>%
    exec(.fn = sprintf, !!!param) %>%
    `[[`(omnipathr.env$urls, .)

}


#' Retrieves an URL and inserts variable parameters
#'
#' @param url_key Character: name of the option containing the URL
#' @param url_key_param List: variables to insert into the `url_key`.
#' @param url_param List: variables to insert into the URL string (which is
#' returned from the options).
#'
#' @importFrom magrittr %>%
#' @importFrom utils URLencode
#'
#' @noRd
url_parser <- function(
    url_key,
    url_key_param = list(),
    url_param = list()
){

    url_key %>%
    get_url(url_key_param) %>%
    c(url_param) %>%
    do.call(what = sprintf) %>%
    URLencode()

}


#' Downloads an URL
#'
#' This function is convenient for appropriate resource retrieval. Following
#' http://bioconductor.org/developers/how-to/web-query/
#' It tries to retrieve the resource one or several times before failing.
#'
#' @param url Character: the URL to download.
#' @param fun The downloader function. Should be able to accept \code{url}
#'     as its first argument. Alternatively, it can be `NULL` or a function
#'     for reading the data, in the latter case `path` must be provided, or
#'     one of the `use_httr` or the `ignore_contents` parameters should be
#'     `TRUE`. In all these cases, `httr::GET` or `httr::POST` will be used
#'     to download the file.
#' @param post List with HTTP POST data. If \code{NULL}, the function
#'     \code{fun} will execute the download and potentially the reading from
#'     the retrieved data. If \code{post} is a list, \code{httr::POST} will
#'     be used to download the data and the contents will be channeled to
#'     \code{fun} to read it.
#' @param path Character: if not `NULL` the file will be downloaded
#'     to this path and the path will be returned.
#' @param http_headers List: a list of HTTP headers. Passed to
#'     `httr2::req_headers`, used only if the downloader function is set up
#'     here (see details at param `fun`).
#' @param init_url Character: retrieve first this URL, to obtain cookies
#'     or start a session.
#' @param init_headers List: HTTP headers for the initial request to
#'     `init_url`.
#' @param return_response Logical: return the response object from `httr`
#'     without any further processing. Used only if the downloader function
#'     is set up here (see details at param `fun`).
#' @param keep_headers Logical: add the response headers to the returned
#'     object as an attribute. If `ignore_contents` is `TRUE`, the returned
#'     object will be `FALSE`, but still might carry the headers in its
#'     `headers` attribute. Used only if the downloader method is set up here
#'     (see details at param `fun`).
#' @param ignore_contents Logical: do not extract the contents from the
#'     response object. The contents still might be saved to the disk if
#'     `path` is provided. The returned object will be the response object
#'     if `return_response` is `TRUE`, otherwise `FALSE` will be returned.
#' @param extract_headers Callable: a custom function which accepts the
#'     response object retrieved from `init_url` and returns a list of
#'     HTTP headers which can be used in the main request. Must be provided
#'     if `init_url` is not `NULL`.
#' @param use_httr Logical: use `httr::GET` to download the data even if no
#'     other argument or condition implies this.
#' @param ... Passed to \code{fun}.
#'
#' @return The output of the downloader function \code{fun}.
#'
#' @importFrom logger log_level log_warn log_trace
#' @importFrom httr2 request req_options req_body_form req_perform
#' @importFrom httr2 resp_body_string resp_status req_headers
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang !!! exec %||%
#' @importFrom curl new_handle
#'
#' @noRd
download_base <- function(
    url,
    fun = NULL,
    post = NULL,
    path = NULL,
    http_headers = list(),
    init_url = NULL,
    init_headers = NULL,
    return_response = FALSE,
    keep_headers = FALSE,
    ignore_contents = FALSE,
    extract_headers = NULL,
    use_httr = FALSE,
    ...
){

    on.exit(close_connection(url))

    if(!is.null(init_url)){

        init_response <- download_base(
            url = init_url,
            http_headers = init_headers,
            return_response = TRUE
        )

        if(!is.null(extract_headers)){

            http_headers <-
                init_response %>%
                extract_headers %>%
                {merge_lists(http_headers, .)}

        }

    }

    url_loglevel <- `if`(
        getOption('omnipathr.print_urls'),
        omnipath_console_loglevel(),
        logger::INFO
    )

    retries <- getOption('omnipathr.retry_downloads')

    log_level(level = url_loglevel, 'Retrieving URL: `%s`', url)

    if(
        !is.null(post) ||
        !is.null(path) ||
        ignore_contents ||
        is.null(fun) ||
        use_httr
){

        log_trace('Downloading by `httr2` in `download_base`.')

        reader <- fun %||% identity

        fun <- function(url, post = NULL, ...){

            log_trace('Preparing httr2 request.')

            http_headers %<>% list_null

            curlopt <- omnipath_new_handle(args_only = TRUE)

            req <-
                url %>%
                request() %>%
                req_headers(!!!http_headers) %>%
                req_options(!!!curlopt) %>%
                doif(!is.null(post), ~req_body_form(.x, !!!post))

            resp <-
                req %T>%
                {log_trace('Sending HTTP request.')} %>%
                req_perform(path = path)

            http_status <- resp %>% resp_status
            msg <- sprintf('HTTP %i', http_status)

            if(http_status != 200L){

                msg %<>% sprintf('%s: %s', resp %>% resp_status_desc)
                log_warn(msg)
                stop(msg)

            }

            log_trace(msg)

            if (!is.null(handle <- resp$request %>% attr('handle'))) {

                log_curl_stats(handle, url)

            } else {

                log_trace('No downlad stats available: no curl handle.')

            }

            result <- FALSE

            if (return_response) {

                result <- resp

            } else if (is.null(path) && !ignore_contents) {

                log_trace('Calling reader callback on response.')

                result <- resp %>% resp_body_string %>% reader(...)

            }

            if(keep_headers){

                attr(result, 'headers') <- headers(response)

            }

            return(result)

        }

    }

    args <- list(...)
    args$post <- post

    for(attempt in seq(retries)){

        log_trace('Attempt %d/%d: `%s`', attempt, retries, url)

        result <- tryCatch(
            exec(fun, url, !!!args),
            error = identity
        )

        if(inherits(result, 'error')){

            msg <-
                sprintf(
                    'Failed to download `%s` (attempt %d/%d); error: %s',
                    url,
                    attempt,
                    retries,
                    conditionMessage(result)
                )

            if(attempt == retries){

                log_error_with_info(msg)
                stop(result)

            }else{

                log_warn(msg)
                # to avoid too fast retries
                Sys.sleep(getOption('omnipathr.retry_downloads_in_seconds'))

            }

        }else{

            break

        }

    }

    return(result)

}


#' Downloads a file to the cache directory
#'
#' Retrieves a file by HTTP GET or POST and returns the path to the cache
#' file.
#'
#' @param url_key Character: name of the option containing the URL
#' @param url_key_param List: variables to insert into the `url_key`.
#' @param url_param List: variables to insert into the URL string (which is
#' returned from the options).
#' @param http_param List: further parameters for \code{httr::GET} or
#'     \code{httr::POST}.
#' @param ext Character: the file extension. If `NULL` a guess from the URL
#'     will be attempted.
#' @param post List with HTTP POST data. If \code{NULL}, \code{httr::GET}
#'     will be used to download the data, otherwise \code{httr::POST}.
#' @param http_headers List: a list of HTTP headers. Passed to
#'     `httr2::req_headers`, used only if the downloader function is set up
#'     here (see details at param `fun`).
#'
#' @noRd
download_to_cache <- function(
    url_key,
    url_key_param = list(),
    url_param = list(),
    ext = NULL,
    post = NULL,
    http_headers = NULL
){

    url <- url_parser(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param
    )

    version <- omnipath_cache_latest_or_new(url = url, post =post, ext = ext)

    from_cache <- version$status == CACHE_STATUS$READY

    if(!from_cache){

        download_base(
            url = url,
            fun = NULL,
            path = version$path,
            post = post,
            http_headers = http_headers
        )
        omnipath_cache_download_ready(version)

    }

    version$path %>%
    origin_cache(from_cache) %>%
    source_attrs(NULL, url)

}

#' Generic method to download a table
#'
#' Downloads a table which can be read by a function from the \code{readr}
#' package or other package.
#'
#' @param url_key Character: name of the option containing the URL
#' @param reader Function: the function to download and read the data.
#' @param url_key_param List: variables to insert into the `url_key`.
#' @param url_param List: variables to insert into the URL string (which is
#'     returned from the options).
#' @param reader_param List: options for the reader function.
#' @param resource Character: name of the resource.
#' @param post List with HTTP POST parameters.
#' @param use_httr Logical: force to use \code{httr::GET} instead of
#'     allowing \code{reader} to download the file.
#' @param ... Passed to \code{\link{download_base}}.
#'
#' @importFrom magrittr %>% %<>% equals
#' @importFrom readr cols
#' @importFrom rlang exec !!!
#' @importFrom logger log_trace
#'
#' @noRd
generic_downloader <- function(
    url_key,
    reader = curl_read_tsv,
    url_key_param = list(),
    url_param = list(),
    reader_param = list(),
    resource = NULL,
    post = NULL,
    use_httr = FALSE,
    ...
){

    log_trace('Downloading by `generic_downloader`.')

    reader_param %<>% doif(
        reader %>%
        environment %>%
        {isNamespace(.) && getNamespaceName(.) %>% equals('readr')},
        ~add_defaults(.x, reader, READR_DEFAULTS)
    )

    url <- url_parser(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param
    )

    log_trace('Looking up in cache: `%s`.', url)

    result <- omnipath_cache_load(url = url, post = post)

    if(is.null(result)){

        log_trace('Could not find in cache, initiating download: `%s`.', url)

        result <-
            exec(
                download_base,
                url,
                reader,
                post,
                !!!reader_param,
                use_httr = use_httr,
                ...
            ) %>%
            omnipath_cache_save(url = url, post = post)

    }

    result %>%
    source_attrs(resource, url)

}


#' Generic method to download an excel worksheet
#'
#' Downloads an xls or xlsx file and extracts a worksheet as a data frame.
#'
#' @param url_key Character: name of the option containing the URL
#' @param sheet Character or integer, passed to \code{read_excel}.
#' @param url_key_param List: variables to insert into the `url_key`.
#' @param url_param List: variables to insert into the URL string (which is
#' returned from the options).
#' @param ext Character: the file extension, either xls or xlsx.
#' @param http_param List: further parameters for \code{httr::GET} or
#'     \code{httr::POST}.
#' @param post List with HTTP POST data. If \code{NULL}, \code{httr::GET}
#'     will be used to download the data, otherwise \code{httr::POST}.
#'
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
#' @importFrom logger log_info
#'
#' @noRd
xls_downloader <- function(
    url_key,
    sheet = NULL,
    url_key_param = list(),
    url_param = list(),
    ext = 'xlsx',
    resource = NULL,
    http_param = list(),
    post = NULL
){

    path <- download_to_cache(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param,
        ext = ext,
        http_param = http_param,
        post = post
    )

    log_info('Reading XLS from `%s`', path)

    read_excel(path, sheet = sheet) %>%
    copy_source_attrs(path, resource = resource)

}


#' Generic method to download a zip archive
#'
#' Downloads a zip file or retrieves it from the cache. Returns the path
#' to the zip file and the list of paths in the archive.
#'
#' @param url_key Character: name of the option containing the URL
#' @param url_key_param List: variables to insert into the `url_key`.
#' @param url_param List: variables to insert into the URL string (which is
#'     returned from the options).
#' @param post List: passed to \code{curl::handle_setform}.
#' @param http_headers Named list with HTTP header keys and values.
#' @param cache_by_url Character: at the cache handling consider this URL
#'     and ignore the POST parameters or the data payload. This is useful if
#'     the download requires an access token which varies at each download
#'     but at reading from the cache no need for token.
#' @param ... Additional options for cURL. Passed to
#'     \code{\link{download_base}}.
#'
#' @importFrom utils unzip untar
#' @importFrom logger log_info log_warn log_trace
#' @importFrom magrittr %>%
#' @importFrom rlang exec !!!
#'
#' @noRd
archive_downloader <- function(
    url_key,
    url_key_param = list(),
    url_param = list(),
    post = NULL,
    http_headers = list(),
    cache_by_url = NULL,
    ...
){

    url <- url_parser(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param
    )

    cache_url <- `if`(is.null(cache_by_url), url, cache_by_url)
    cache_post <- `if`(is.null(cache_by_url), post, NULL)

    version <- omnipath_cache_latest_or_new(
        url = cache_url,
        post = cache_post
    )

    from_cache <- version$status == CACHE_STATUS$READY

    if(!from_cache){

        logger::log_info('Downloading `%s`', url)

        success <- download_base(
            url = url,
            path = version$path,
            post = post,
            http_headers = http_headers,
            ...
        )
        omnipath_cache_download_ready(version)
        key <- omnipath_cache_key_from_version(version)
        ext <- archive_type(version$path, url)
        if(is.null(ext)){
            log_warn('Could not state archive type: `%s`.', url)
        }else{
            log_trace('Archive type `%s`: `%s`', ext, url)
            omnipath_cache_set_ext(key, ext)
        }
        version <- omnipath_cache_latest_or_new(url = url)

    }

    key <- omnipath_cache_key_from_version(version)
    record <- omnipath_cache_get(key)
    extractor <- `if`(record$ext == 'zip', unzip, untar)

    list(
        path = version$path,
        url = cache_url,
        files = extractor(version$path, list = TRUE),
        ext = record$ext,
        from_cache = from_cache
    )

}


#' Generic method to download an archive and extract one file
#'
#' @param url_key Character: name of the option containing the URL
#' @param path Character: path to the file within the archive. If NULL, the
#'     first file will be extracted (not recommended).
#' @param url_key_param List: variables to insert into the `url_key`.
#' @param url_param List: variables to insert into the URL string (which is
#'     returned from the options).
#' @param post List: POST parameters. If NULL, a GET query is performed.
#' @param reader Optional, a function to read the connection.
#' @param reader_param List: arguments for the reader function.
#' @param cache_by_url Character: at the cache handling use this URL instead
#'     of the one provided in `url` and ignore the POST parameters or the
#'     data payload. This is useful if the download requires an access token
#'     which varies at each download but at reading from the cache no need
#'     for token.
#' @param resource Character: name of the resource.
#' @param extract_xls Logical: read worksheet from xls(x) file automatically.
#' @param to_tempdir Logical: if TRUE, the archive is extracted to a temporary
#'      directory from the zip file. Ignored for other archive types.
#' @param ... Passed to \code{\link{archive_downloader}}.
#'
#' @return A connection to the extracted file or a the contents read from
#'     the path inside the archive.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom logger log_fatal log_trace
#' @importFrom rlang %||% !!! exec
#' @importFrom readxl read_excel
#' @importFrom utils unzip
#'
#' @seealso \code{\link{archive_downloader}}
#'
#' @noRd
archive_extractor <- function(
    url_key,
    path = NULL,
    url_key_param = list(),
    url_param = list(),
    post = NULL,
    reader = NULL,
    reader_param = list(),
    cache_by_url = NULL,
    resource = NULL,
    extract_xls = TRUE,
    to_tempdir = FALSE,
    ...
){

    archive_data <- archive_downloader(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param,
        post = post,
        cache_by_url = cache_by_url,
        ...
    )

    # fallback to the first file
    path <- `if`(is.null(path), paths_in_archive(archive_data)[1], path)

    if(!path %in% paths_in_archive(archive_data)){

        log_trace(
            'First 8 bytes of `%s`: %s',
            archive_data$path,
            archive_data$path %>%
                readBin('raw', n = 8) %>%
                paste(collapse = ' ')
        )
        log_trace(
            'Archive type stated: `%s`; size: %d bytes.',
            archive_type(archive_data$path, archive_data$url),
            file.info(archive_data$path)$size
        )

        msg <- sprintf(
            'Path `%s` not found in archive `%s` (local file at `%s`)',
            path,
            archive_data$url,
            archive_data$path
        )
        logger::log_fatal(msg)
        stop(msg)

    }

    xls <- path %>% endsWith(c('xls', 'xlsx')) %>% any && extract_xls

    if(archive_data$ext == 'zip'){

        if(xls || to_tempdir){

            con <-
                archive_data$path %>%
                unzip(files = path, exdir = tempdir()) %>%
                extract(1L)

            if(xls) {
                reader <- read_excel
            }else{
                return(con)
            }

        }else{

            con <-
                archive_data$path %>%
                unz(path, open = 'rb')

        }

    }else{

        archive_data$path %>% untar(files = path, exdir = tempdir())
        ext_path <- tempdir() %>% file.path(path)

        if(xls){

            con <- ext_path
            reader <- read_excel

        }else{

            con <- ext_path %>% file(open = 'rb')

        }

    }

    if(is.null(reader)){

        return(con)

    }else{

        reader %>%
        exec(con, !!!reader_param) %>%
        origin_cache(archive_data$from_cache) %>%
        source_attrs(resource, archive_data$url) %T>%
        {if('connection' %in% class(con)) base::close(con)}

    }

}


#' Workaround for the different APIs of `unzip` and `untar`
#'
#' @return Character vector of paths in the archive.
#'
#' @param archive_data List: as returned by the
#' \code{\link{archive_downloader}} function.
#'
#' @seealso \code{\link{archive_downloader}}
#'
#' @noRd
paths_in_archive <- function(archive_data){

    `if`(
        'Name' %in% names(archive_data$files),
        archive_data$files$Name,
        archive_data$files
    )

}


#' Assigns attributes to a data frame about its sources
#'
#' @param data A data frame (or other object).
#' @param resource Character: the name of the resource.
#' @param url Character: the download URL.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang %||%
#'
#' @noRd
source_attrs <- function(data, resource, url){

    domain <-
        url %||% 'unknown domain' %>%
        domain_from_url

    source <- `if`(
        is.null(resource),
        domain,
        sprintf('%s (%s)', resource, domain)
    )

    data %>%
    `attr<-`('resource', resource) %>%
    `attr<-`('url', url) %>%
    `attr<-`('source', source)

}


#' Copies attributes about its sources from one object to another
#'
#' @param to The object to copy attributes to.
#' @param from The object to copy attributes from.
#' @param ... Passed to \code{\link{update_source_attrs}}.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
copy_source_attrs <- function(to, from, ...){

    to %>%
    copy_attrs(from, c('source', 'resource', 'url', 'origin')) %>%
    update_source_attrs(...)

}


#' Updates the source attributes of an object
#'
#' @importFrom rlang %||%
#' @importFrom magrittr %>%
#'
#' @noRd
update_source_attrs <- function(obj, ...){

    attrs <- list(...)
    origin <- attrs$origin %||% attr(obj, 'origin')
    from_cache <- (
        !is.null(origin) && origin == 'cache' ||
        attrs$from_cache %||% FALSE
    )
    resource <- attrs$resource %||% attr(obj, 'resource')
    url <- attrs$url %||% attr(obj, 'url')

    obj %>%
    origin_cache(from_cache) %>%
    source_attrs(resource, url)

}


#' Load the built-in magic byte database
#'
#' Called only in the package loading process.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest_longer unnest_wider
#' @importFrom purrr map
#'
#' @noRd
.load_magic_bytes <- function(pkgname){

    # NSE vs. R CMD check workaround
    magic <- NULL

    omnipathr.env$mb <-
        system.file(
            'internal',
            'magic_bytes.json',
            package = pkgname,
            mustWork = TRUE
        ) %>%
        safe_json(simplifyDataFrame = FALSE) %>%
        tibble(magic = .) %>%
        mutate(ext = names(magic)) %>%
        unnest_longer(magic) %>%
        unnest_wider(magic) %>%
        mutate(
            magic = map(magic, as.raw)
        )

}


#' Tells the type of an archive
#'
#' @param path Character: path to the archive.
#' @param url Character, optional: the URL the archive has been downloaded
#'     from.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr pmap detect
#'
#' @noRd
archive_type <- function(path, url = NULL){

    # NSE vs. R CMD check workaround
    ext <- NULL

    max_offset <- omnipathr.env$mb %>% {length(.$magic) + .$offset} %>% max

    header <- readBin(path, 'raw', n = max_offset)

    omnipathr.env$mb %>%
    pmap(list) %>%
    detect(
        function(row){
            all(
                row$magic ==
                header[row$offset:(row$offset + length(row$magic))]
            )
        }
    ) %>%
    `$`(ext) %>%
    {`if`(
        !is.null(.) && . %in% c('gz', 'bz', 'bz2', 'xz', 'zst'),
        `if`(
            any(grepl(sprintf('\\.tar\\.|\\.t%s', .), c(path, url))),
            sprintf('tar.%s', .),
            .
        ),
        .
    )}

}


#' Custom user-agent header from options
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr first
#' @noRd
user_agent <- function() {

    # NSE vs. R CMD check workaround
    'omnipathr.user_agent' %>%
    options %>%
    first %>%
    list('User-Agent' = .)

}


#' Download by curl, read by read_tsv
#'
#' @importFrom readr read_tsv
#' @noRd
curl_read_tsv <- function(url, curl_param = list(), ...) {

    args <- list(...) %>% add_defaults(read_tsv, READR_DEFAULTS)

    exec(omnipath_curl, url, curl_param, callback = read_tsv, !!!args)

}


#' Download by curl, read by stream_in
#'
# #' @importFrom jsonlite stream_in
#' @importFrom jsonlite fromJSON
#' @noRd
curl_read_json <- function(url, curl_param = list(), ...) {

    # this is a temporary replacement for jsonlite::stream_in
    # because the server returns JSON without any newlines
    # and then stream_in overflows the stack and crashes
    json_readlines <- function(con, ...) {

        jsonlite::fromJSON(readLines(con, warn = FALSE), ...)

    }

    omnipath_curl(url, curl_param, callback = json_readlines, ...)

}


#' Download by curl with OmnipathR specific options
#'
#' @param url Character: URL to download
#' @param curl_param List: parameters to pass to \code{curl::new_handle}
#' @param callback Optional, a function to call on the connection.
#' @param ... Passed to \code{callback}
#'
#' @importFrom readr read_tsv
#' @importFrom curl curl handle_data curl_fetch_disk
#' @importFrom fs path_ext
#' @importFrom logger log_trace
#' @importFrom rlang exec !!!
#' @importFrom magrittr %>% extract2
#' @noRd
omnipath_curl <- function(
        url,
        curl_param = list(),
        callback = NULL,
        compr = NULL,
        ...
    ) {

    handle <- exec(omnipath_new_handle, !!!curl_param)

    COMPR <- list(gz = gzfile, bz2 = bzfile, xz = xzfile)
    fname <- url %>% fname_from_url
    compr %<>% if_null(
        fname %>% path_ext %>% intersect(names(COMPR)) %>% if_null_len0(NULL)
    )

    if (!is.null(compr)) {

        path <- url %>% fname_from_url %>% file.path(tempdir(), .)
        log_trace(
            paste0(
               '`%s` compressed file, downloading to ',
               'temporary location: `%s`'
            ),
            compr,
            path
        )
        curl_fetch_disk(url, path, handle = handle)
        con <- COMPR %>% extract2(compr) %>% exec(path, open = 'rb')

    } else {

        con <- curl(url, open = 'rb', handle = handle)

    }

    if (!is.null(callback)) {

        result <-
            tryCatch(
                { callback(con, ...) },
                error = function(e) {
                    stop(e)
                },
                finally = {
                    close(con)
                }
            )

        log_curl_stats(handle, url)

        return(result)

    }

    return(con)

}


#' Log message with download statistics
#'
#' @importFrom magrittr %>% extract
#' @importFrom logger log_trace
#' @noRd
log_curl_stats <- function(handle, url) {

    handle_received <- `%:::%`('curl', 'handle_received')
    handle_speed <- `%:::%`('curl', 'handle_speed')
    domain <- url %>% domain_from_url
    stats <- handle %>% handle_data

    log_trace(
        paste(
            'Downloaded %s in %s from %s (%s/s); Redirect: %s, DNS look',
            'up: %s, Connection: %s, Pretransfer: %s, First byte at: %s'
        ),
        handle %>% handle_received %>% format_bytes,
        stats$times['total'] %>% format_period,
        domain,
        handle %>% handle_speed %>% extract(1L) %>% format_bytes,
        stats$times['redirect'] %>% format_period,
        stats$times['namelookup'] %>% format_period,
        stats$times['connect'] %>% format_period,
        stats$times['pretransfer'] %>% format_period,
        stats$times['starttransfer'] %>% format_period
    )

}


#' Curl debug callback
#'
#' @importFrom logger log_trace
#' @noRd
curl_debug <- function(type, data) {

    data %<>% {`if`(
        type == 2L || type == 3L || type == 5L || type == 6L,
        sprintf('%i bytes of data', length(.)),
        readBin(., 'character')
    )}

    log_trace('CURL DEBUG[%s]: %s', CURL_DEBUG_TYPES[[type + 1L]], data)

}


#' Create a new curl handle with OmnipathR specific options
#'
#' @param ... Curl options named by their curl/httr2 synonyms. See the output of
#'     \code{curl::curl_options}.
#' @param args_only Logical: if \code{TRUE} return only the arguments, otherwise
#'     return a curl handle.
#'
#' @return A list of arguments if \code{args_only = TRUE}, otherwise a curl handle.
#'
#' @importFrom rlang exec !!!
#' @importFrom magrittr %>% %<>% extract
#' @importFrom curl new_handle curl_options
#' @noRd
omnipath_new_handle <- function(args_only = FALSE, ...) {

    from_config <- list(
        connecttimeout =
            getOption('omnipathr.connect_timeout'),
        timeout = getOption('omnipathr.http_timeout'),
        debugfunction = curl_debug,
        verbose = getOption('omnipathr.curl_verbose'),
        tcp_keepalive = getOption('omnipathr.tcp_keepalive'),
        tcp_keepintvl = getOption('omnipathr.tcp_keepintvl'),
        tcp_keepidle = getOption('omnipathr.tcp_keepidle'),
        tcp_keepcnt = getOption('omnipathr.tcp_keepcnt'),
        upkeep_interval_ms = getOption('omnipathr.upkeep_interval_ms'),
        ssl_verifypeer = getOption('omnipathr.ssl_verifypeer'),
        ssl_verifyhost = getOption('omnipathr.ssl_verifyhost')
    )

    args <-
        list(...) %>%
        merge_lists(from_config)


    param_noavail <-
        args %>%
        names %>%
        setdiff(curl_options() %>% names)

    if (length(param_noavail) > 0L) {

        log_trace(
            'The following curl options are not available: %s',
            compact_repr(param_noavail)
        )

        args %<>% extract(names(.) %>% setdiff(param_noavail))

    }

    log_trace('Curl options: %s', compact_repr(args, limit = 99L))

    `if`(args_only, args, exec(new_handle, !!!args))

}
