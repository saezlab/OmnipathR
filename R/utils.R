#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2023
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


#' Open one or more PubMed articles
#'
#' @param pmids Character or numberic vector of one or more PubMed IDs.
#' @param browser Character: name of the web browser executable. If `NULL`,
#'     the default web browser will be used.
#' @param sep Character: split the PubMed IDs by this separator.
#' @param max_pages Numeric: largest number of pages to open. This is to
#'     prevent opening hundreds or thousands of pages at once.
#'
#' @return Returns `NULL`.
#'
#' @examples
#' interactions <- import_omnipath_interactions()
#' pubmed_open(interactions$references[1])
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom logger log_error
#' @importFrom withr local_options
#' @importFrom stringr str_split str_replace_all
#' @importFrom purrr walk
#' @importFrom utils browseURL
#' @export
pubmed_open <- function(pmids, browser = NULL, sep = ';', max_pages = 25L){

    browser %<>% if_null(getOption('browser'))

    if(is_empty_2(browser)){

        logger::log_error(
            'To open pages in a web browser, set the `browser` option. ',
            'For example: options(browser = "firefox").'
        )
        return(invisible(NULL))

    }

    withr::local_options(browser = browser)

    logger::log_trace('Browser is: `%s`', browser)

    pmids %>%
    as.character %>%
    stringr::str_split(sep) %>%
    unlist %>%
    stringr::str_replace_all('[^:;]*+:(\\d+)', '\\1') %>%
    unique %>%
    {`if`(
        length(.) > max_pages,
        logger::log_error(
            paste0(
                'To open %i pages in a browser, increase the ',
                '`max_pages` parameter (currently %i).'
            ),
            length(.),
            max_pages
        ),
        `if`(
            length(.) > 1L,
            purrr::walk(., pubmed_open, browser = browser),
            sprintf('https://pubmed.ncbi.nlm.nih.gov/%s/', .) %>%
            utils::browseURL()
        )
    )}

}
