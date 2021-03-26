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


#' Retrieves an URL from options and inserts variable parameters
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
    c(url_key_param) %>%
    do.call(what = sprintf) %>%
    options() %>%
    `[[`(1) %>%
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
#'     as its first argument.
#' @param ... Passed to \code{fun}.
#'
#' @return The output of the downloader function \code{fun}.
#'
#' @importFrom logger log_level log_error log_warn
#'
#' @noRd
download_base <- function(url, fun, ...){

    op <- options(timeout = 600)
    on.exit(options(op))

    url_loglevel <- `if`(
        getOption('omnipath.print_urls'),
        omnipath_console_loglevel(),
        logger::INFO
    )

    retries <- getOption('omnipath.retry_downloads')

    log_level(level = url_loglevel, 'Retrieving URL: `%s`', url)

    for(attempt in seq(retries)){

        result <- tryCatch(fun(url, ...), error = identity)

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

                log_error(msg)
                stop(msg)

            }else{

                log_warn(msg)

            }

        }else{

            break

        }

    }

    return(result)

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
#' returned from the options).
#' @param reader_param List: options for the reader function.
#' @param resource Character: name of the resource.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom readr read_tsv cols
#' @importFrom rlang exec !!!
#'
#' @noRd
generic_downloader <- function(
    url_key,
    reader = read_tsv,
    url_key_param = list(),
    url_param = list(),
    reader_param = list(col_types = cols()),
    resource = NULL
){

    reader_param %<>% add_defaults(
        reader,
        list(
            col_types = cols(),
            progress = FALSE
        )
    )

    url <- url_parser(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param
    )

    result <- omnipath_cache_load(url = url)

    if(is.null(result)){

        result <-
            exec(download_base, url, reader, !!!reader_param) %>%
            omnipath_cache_save(url = url)

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
#'
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
#' @importFrom utils download.file
#'
#' @noRd
xls_downloader <- function(
    url_key,
    sheet = NULL,
    url_key_param = list(),
    url_param = list(),
    ext = 'xlsx',
    resource = NULL
){

    url <- url_parser(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param
    )

    version <- omnipath_cache_latest_or_new(url = url, ext = ext)

    from_cache <- version$status == CACHE_STATUS$READY

    if(!from_cache){

        download.file(url = url, destfile = version$path, quiet = TRUE)
        omnipath_cache_download_ready(version)

    }

    read_excel(version$path, sheet = sheet) %>%
    origin_cache(from_cache) %>%
    source_attrs(resource, url)

}


#' Generic method to download a zip archive
#'
#' Downloads a zip file or retrieves it from the cache. Returns the path
#' to the zip file and the list of paths in the archive.
#'
#' @param url_key Character: name of the option containing the URL
#' @param url_key_param List: variables to insert into the `url_key`.
#' @param url_param List: variables to insert into the URL string (which is
#' returned from the options).
#' @param cache_by_url Logical: at the cache handling consider only the URL,
#'     and ignore the POST parameters or the data payload. This is useful if
#'     the download requires an access token which varies at each download
#'     but at reading from the cache no need for token.
#' @param ... Additional options for cURL. Passed to
#'     `curlSetOpt`.
#'
#' @importFrom utils unzip untar
#' @importFrom RCurl getCurlHandle CFILE curlSetOpt curlPerform close
#' @importFrom wand get_content_type
#' @importFrom logger log_info
#'
#' @noRd
archive_downloader <- function(
    url_key,
    url_key_param = list(),
    url_param = list(),
    post = NULL,
    curl_verbose = FALSE,
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
        # downloading the data
        curl_handle <- getCurlHandle()
        response <- CFILE(version$path, mode = 'wb')
        opt_set <- curlSetOpt(
            curl = curl_handle,
            url = url,
            verbose = curl_verbose,
            writedata = response@ref,
            ...
        )
        success <- curlPerform(curl = curl_handle)
        RCurl::close(response)
        omnipath_cache_download_ready(version)
        key <- omnipath_cache_key_from_version(version)
        content_type <- get_content_type(version$path)
        ext <- `if`('application/zip' %in% content_type, 'zip', 'tar.gz')
        omnipath_cache_set_ext(key, ext)
        version <- omnipath_cache_latest_or_new(url = url)

    }

    key <- omnipath_cache_key_from_version(version)
    record <- omnipath.env$cache[[key]]
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
#' @param ... Passed to \code{\link{archive_downloader}}.
#'
#' @return A connection to the extracted file or a the contents read from
#'     the path inside the archive.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom logger log_fatal
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

    if(!(path %in% paths_in_archive(archive_data))){

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

        if(xls){

            con <-
                archive_data$path %>%
                unzip(files = path, exdir = tempdir()) %>%
                `[`(1)

            reader <- read_excel

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
#' @importFrom httr parse_url
#'
#' @noRd
source_attrs <- function(data, resource, url){

    # NSE vs. R CMD check Workaround
    hostname <- NULL

    domain <-
        url %||% 'unknown domain' %>%
        parse_url %>%
        `$`(hostname) %||% 'unknown domain'

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
#'
#' @importFrom magrittr %>%
#'
#' @noRd
copy_source_attrs <- function(to, from){

    to %>%
    copy_attrs(from, c('source', 'resource', 'url', 'origin'))

}