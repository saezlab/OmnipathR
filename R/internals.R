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
#' @param post List with HTTP POST data. If \code{NULL}, the function
#'     \code{fun} will execute the download and potentially the reading from
#'     the retrieved data. If \code{post} is a list, \code{httr::POST} will
#'     be used to download the data and the contents will be channeled to
#'     \code{fun} to read it.
#' @param http_param List: further parameters for \code{httr::GET} or
#'     \code{httr::POST}. Used only if \code{path} is not \code{NULL} or
#'     if \code{post} is not \code{NULL}.
#' @param content_param List: further parameters for \code{httr::content}.
#'     Used only if \code{path} is \code{NULL} and \code{post} is not
#'     \code{NULL}.
#' @param path Character: if not `NULL` the file will be downloaded
#'     to this path and the path will be returned.
#' @param ... Passed to \code{fun}.
#'
#' @return The output of the downloader function \code{fun}.
#'
#' @importFrom logger log_level log_error log_warn log_trace
#' @importFrom httr POST GET content write_disk
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang !!! exec
#'
#' @noRd
download_base <- function(
    url,
    fun,
    post = NULL,
    http_param = list(),
    content_param = list(),
    path = NULL,
    ...
){

    op <- options(timeout = 600)
    on.exit(options(op))

    url_loglevel <- `if`(
        getOption('omnipath.print_urls'),
        omnipath_console_loglevel(),
        logger::INFO
    )

    retries <- getOption('omnipath.retry_downloads')

    log_level(level = url_loglevel, 'Retrieving URL: `%s`', url)

    if(!is.null(post) || !is.null(path)){

        reader <- fun

        fun <- function(url, post = NULL, ...){

            http_method <- `if`(is.null(post), GET, POST)
            http_param$body <- post

            if(!is.null(path)){
                http_param %<>% c(list(write_disk(path, overwrite = TRUE)))
            }

            exec(http_method, url = url, !!!http_param) %>%
            {`if`(
                is.null(path),
                exec(content, ., !!!content_param) %>%
                reader(...),
                path
            )}

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
#'
#' @noRd
download_to_cache <- function(
    url_key,
    url_key_param = list(),
    url_param = list(),
    http_param = list(),
    ext = NULL,
    post = NULL
){

    url <- url_parser(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param
    )

    version <- omnipath_cache_latest_or_new(url = url, ext = ext)

    from_cache <- version$status == CACHE_STATUS$READY

    if(!from_cache){

        download_base(
            url = url,
            fun = NULL,
            path = version$path,
            http_param = http_param,
            post = post
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
#' @importFrom magrittr %>% %<>%
#' @importFrom readr read_tsv cols
#' @importFrom rlang exec !!!
#' @importFrom logger log_trace
#'
#' @noRd
generic_downloader <- function(
    url_key,
    reader = read_tsv,
    url_key_param = list(),
    url_param = list(),
    reader_param = list(),
    resource = NULL,
    post = NULL,
    use_httr = FALSE,
    ...
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

    log_trace('Looking up in cache: `%s`.', url)

    result <- omnipath_cache_load(url = url, post = post)

    if(is.null(result)){

        log_trace('Could not find in cache, initiating download: `%s`.', url)

        result <-
            exec(download_base, url, reader, post, !!!reader_param, ...) %>%
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
#' @importFrom logger log_info log_warn log_trace
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
            ssl.cipher.list = 'HIGH:!ECDH',
            ...
        )
        success <- download_base(
            url = url,
            fun = function(url, ...){curlPerform(...)},
            curl = curl_handle
        )
        RCurl::close(response)
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
#' @importFrom jsonlite fromJSON
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest_longer unnest_wider
#' @importFrom purrr map
#'
#' @noRd
.load_magic_bytes <- function(pkgname){

    # NSE vs. R CMD check workaround
    magic <- NULL

    omnipath.env$mb <-
        system.file(
            'internal',
            'magic_bytes.json',
            package = pkgname,
            mustWork = TRUE
        ) %>%
        fromJSON(simplifyDataFrame = FALSE) %>%
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

    max_offset <- omnipath.env$mb %>% {length(.$magic) + .$offset} %>% max

    header <- readBin(path, 'raw', n = max_offset)

    omnipath.env$mb %>%
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