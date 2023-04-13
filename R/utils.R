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
    as.numeric %>%
    sort(decreasing = TRUE) %>%
    as.character %>%
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
        sprintf('https://pubmed.ncbi.nlm.nih.gov/%s/', .) %>%
        purrr::walk(., utils::browseURL),
    )}

}


#' Show evidences for an interaction
#'
#' @details
#' If the number of references is larger than `max_pages`, the most recent
#' ones will be opened. URLs are passed to the browser in order of decreasing
#' publication date, though browsers do not seem to respect the order at all.
#' In addition Firefox, if it's not open already, tends to randomly open empty
#' tab for the first or last URL, have no idea what to do about it.
#'
#' @param partner_a Identifier or name of one interacting partner. The order
#'     of the partners matter only if `directed` is `TRUE`. For both partners,
#'     vectors of more than one identifiers can be passed.
#' @param partner_b Identifier or name of the other interacting partner.
#' @param interactions An interaction data frame. If not provided, all
#'     interactions will be loaded within this function, but that takes
#'     noticeable time. If a `list` is provided, it will be used as
#'     parameters for \code{\link{import_omnipath_interactions}}. This way
#'     you can define the organism, datasets or the interaction type.
#' @param directed Logical: does the direction matter? If `TRUE`, only
#'     a → b interactions will be shown.
#' @param open Logical: open online articles in a web browser.
#' @param browser Character: override the web browser executable used
#'     to open online articles.
#' @param max_pages Numeric: largest number of pages to open. This is to
#'     prevent opening hundreds or thousands of pages at once.
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' evidences('CALM1', 'TRPC1', list(datasets = 'omnipath'))
#' }
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom rlang exec !!!
#' @importFrom logger log_success
#' @importFrom dplyr filter rename_with pull
#' @importFrom stringr str_split
#' @export
evidences <- function(
        partner_a,
        partner_b,
        interactions = NULL,
        directed = FALSE,
        open = TRUE,
        browser = NULL,
        max_pages = 25L
) {

    # R CMD check vs. NSE workaround
    source <- source_genesymbol <-
    target <- target_genesymbol <- references <- NULL

    source_side <- `if`(directed, partner_a, c(partner_a, partner_b))
    target_side <- `if`(directed, partner_b, c(partner_b, partner_a))

    interactions %>%
    {`if`(
        is.data.frame(.),
        .,
        exec(import_omnipath_interactions, !!!.)
    )} %>%
    rename_with(~sub('enzyme', 'source', .x)) %>%
    rename_with(~sub('substrate', 'target', .x)) %>%
    filter(
        (
            source %in% source_side |
            source_genesymbol %in% source_side
        ) & (
            target %in% target_side |
            target_genesymbol %in% target_side
        )
    ) %>%
    {`if`(
        nrow(.) > 0L,
        identity(.) %T>%
        {log_success(
            'Resources: %s.',
            .$sources %>%
            str_split(';') %>%
            unlist %>%
            sort %>%
            unique %>%
            paste(collapse = ', ')
        )} %>%
        strip_resource_labels %>%
        pull(references) %>%
        str_split(';') %>%
        unlist %>%
        unique %>%
        as.numeric %>%
        sort(decreasing = TRUE) %>%
        as.character %T>%
        {log_success('Found %i references.', length(.))} %>%
        head(n = max_pages) %>%
        pubmed_open,
        log_success(
            'No interaction between %s and %s.',
            partner_a,
            partner_b
        )
    )}

}
