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
    do.call(what = sprintf)

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
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv cols
generic_downloader <- function(
    url_key,
    reader = read_tsv,
    url_key_param = list(),
    url_param = list(),
    reader_param = list(col_types = cols())
){

    if(!('col_types' %in% names(reader_param))){
        reader_param$col_types <- cols()
    }

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

    return(result)

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
#' @importFrom utils unzip
zip_downloader <- function(
    url_key,
    url_key_param = NULL
    url_param = NULL
){

    url <- url_parser(
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param
    )

    version <- omnipath_cache_latest_or_new(url = url, ext = 'zip')

    if(version$status != CACHE_STATUS$READY){

        download.file(url = url, destfile = version$path, quiet = TRUE)
        omnipath_cache_download_ready(version)

    }

    list(
        path = version$path,
        files = unzip(version$path, list = TRUE)
    )

}


zip_extractor <- function(

){

}