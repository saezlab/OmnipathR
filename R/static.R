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


#' List the static tables available from OmniPath
#'
#' A few resources and datasets are available also as plain TSV files and
#' can be accessed without TLS. The purpose of these tables is to make the
#' most often used OmniPath data available on computers with configuration
#' issues. These tables are not the recommended way to access OmniPath
#' data, and a warning is issued each time they are accessed.
#'
#' @return A data frame listing the available tables.
#'
#' @examples
#' static_tables()
#'
#' @importFrom magrittr %>% extract
#' @importFrom tibble tibble
#' @importFrom tidyr separate_wider_regex separate_wider_delim
#' @importFrom stringr str_extract
#' @importFrom lubridate dmy hm
#' @export
static_tables <- function() {

    'omnipath_static' %>%
    download_to_cache %>%
    readLines %>%
    extract(8L:length(.) - 2L) %>%
    tibble(fname = .) %>%
    separate_wider_regex(
        fname,
        c(
           '<a [^>]+>',
           fname = '[A-z_\\d]+.tsv.gz',
           '</a>\\s+',
           date = '[^\\s]+',
           ' ',
           time = '[\\d:]+', '\\s+',
           size = '\\d+'
        )
    ) %>%
    separate_wider_delim(
        fname,
        '_',
        names = c('query', 'resource', 'organism'),
        cols_remove = FALSE
    ) %>%
    mutate(
        organism = str_extract(organism, '^\\d+'),
        date = dmy(date),
        time = hm(time),
        size = as.integer(size),
        url = paste0(url_parser('omnipath_static'), fname)
    )

}
