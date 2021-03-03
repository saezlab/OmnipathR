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
#' @param resource Character: the name of the resource
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom readr read_tsv cols
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
            url %>%
            c(reader_param) %>%
            do.call(what = reader) %>%
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
#' @importFrom
xls_downloader <- function(
    url_key,
    sheet = NULL,
    url_key_param = list(),
    url_param = list(),
    ext = 'xlsx'
){

    url <- url_parser(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param
    )

    version <- omnipath_cache_latest_or_new(url = url, ext = ext)

    if(version$status != CACHE_STATUS$READY){

        download.file(url = url, destfile = version$path, quiet = TRUE)
        omnipath_cache_download_ready(version)

    }

    read_excel(version$path, sheet = sheet)

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

    if(version$status != CACHE_STATUS$READY){

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
    record <- .omnipath_cache[[key]]
    extractor <- `if`(record$ext == 'zip', unzip, untar)

    list(
        path = version$path,
        url = cache_url,
        files = extractor(version$path, list = TRUE),
        ext = record$ext
    )

}


#' Generic method to download a zip archive and extract one file
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
#'
#' @return A connection to the extracted file or a
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom logger log_fatal
#' @seealso \code{\link{archive_downloader}}
archive_extractor <- function(
    url_key,
    path = NULL,
    url_key_param = list(),
    url_param = list(),
    post = NULL,
    reader = NULL,
    reader_param = list(),
    cache_by_url = NULL,
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

    if(archive_data$ext == 'zip'){

        con <-
            archive_data$path %>%
            unz(path, open = 'rb')

    }else{

        archive_data$path %>% untar(files = path, exdir = tempdir())
        con <-
            tempdir() %>%
            file.path(path) %>%
            file(open = 'rb')

    }

    if(is.null(reader)){

        return(con)

    }else{

        reader_param %<>% c(list(con), .)
        result <- do.call(reader, reader_param)
        base::close(con)
        return(result)

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
#' @importFrom httr parse_url
source_attrs <- function(data, resource, url){

    domain <- parse_url(url)$hostname
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