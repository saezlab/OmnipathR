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

    url <-
        url_key %>%
        do.call(sprintf, c(., url_key_param)) %>%
        options() %>%
        `[[`(1) %>%
        do.call(sprintf, url_param)

    result <- omnipath_cache_load(url = url)

    if(is.null(result)){

        result <- url %>% do.call(reader, reader_param)

        omnipath_cache_save(data = result, url = url)

    }

    return(result)

}