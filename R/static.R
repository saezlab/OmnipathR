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
#' @importFrom magrittr %>% extract extract2
#' @importFrom curl curl_fetch_memory
#' @importFrom stringi stri_split_lines1
#' @importFrom tibble tibble
#' @importFrom tidyr separate_wider_regex separate_wider_delim
#' @importFrom stringr str_extract
#' @importFrom lubridate dmy hm
#' @export
#' @seealso \code{\link{static_table}}
static_tables <- function() {

    # NSE vs. R CMD check workaround
    fname <- organism <- time <- size <- NULL

    'omnipath_static' %>%
    url_parser %>%
    curl_fetch_memory %>%
    extract2('content') %>%
    rawToChar %>%
    stri_split_lines1 %>%
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


#' Retrieve a static table from OmniPath
#'
#' A few resources and datasets are available also as plain TSV files and
#' can be accessed without TLS. The purpose of these tables is to make the
#' most often used OmniPath data available on computers with configuration
#' issues. These tables are not the recommended way to access OmniPath
#' data, and a warning is issued each time they are accessed.
#'
#' @param query Character: a query type such as "annotations" or
#'     "interactions".
#' @param resource Character: name of the resource or dataset, such as
#'     "CollecTRI" or "PROGENy".
#' @param organism Integer: NCBI Taxonomy of the organism: 9606 for human,
#'     10090 for mouse and 10116 for rat.
#' @param strict_evidences Logical: restrict the evidences to the queried
#'     datasets and resources. If set to FALSE, the directions and effect signs
#'     and references might be based on other datasets and resources.
#' @param wide Convert the annotation table to wide format, which
#'     corresponds more or less to the original resource. If the data comes
#'     from more than one resource a list of wide tables will be returned.
#'     See examples at \code{\link{pivot_annotations}}.
#' @param dorothea_levels Vector detailing the confidence levels of the
#'     interactions to be downloaded. In dorothea, every TF-target interaction
#'     has a confidence score ranging from A to E, being A the most reliable
#'     interactions.
#'     By default here we take A, B and C level interactions
#'     (\code{c("A", "B", "C")}).
#'     It is to note that E interactions are not available in OmnipathR.
#'
#' @return A data frame (tibble) with the requested resource.
#'
#' @examples
#' static_table("annotations", "PROGENy")
#'
#' @importFrom logger log_warn log_error
#' @importFrom magrittr %>% %<>% extract
#' @importFrom stringr str_to_lower str_detect
#' @importFrom dplyr pull filter
#' @importFrom readr read_tsv cols
#' @export
#' @seealso \code{\link{static_tables}}
static_table <- function(
        query,
        resource,
        organism = 9606L,
        strict_evidences = TRUE,
        wide = TRUE,
        dorothea_levels = c('A', 'B', 'C')
    ) {

    # NSE vs. R CMD check workaround
    .env <- NULL

    log_warn(
         sprintf('Accessing `%s` as a static table: this is not ', resource),
         'the recommended way to access OmniPath data; it is only a backup ',
         'plan for situations when our server or your computer is ',
         'experiencing issues.'
    )

    organism %<>% as.integer %>% as.character

    query_l <- str_to_lower(query)
    resource_l <- str_to_lower(resource)
    datasets <- `if`(
        resource_l %in% c('dorothea', 'collectri'),
        resource_l,
        NULL
    )
    resources <- `if`(is.null(datasets), resource_l, NULL)

    static_tables() %>%
    filter(
        query == query_l &
        str_to_lower(resource) == resource_l &
        organism == .env$organism
    ) %>%
    {`if`(
        nrow(.) < 1L,
        {
            msg <-
                paste0(
                    'Static table not available for query `%s`, resource ',
                    '`%s` and organism `%i`. For a list of available tables ',
                    'see `static_tables()`.'
                ) %>%
                sprintf(query, resource, organism)
            log_error(msg)
            stop(msg)
        },
        .
    )} %>%
    pull(url) %>%
    extract(1L) %>%
    omnipath_download(
        read_tsv,
        col_types = cols(),
        progress = FALSE,
        show_col_types = FALSE
    ) %>%
    omnipath_post_download(
         url = attr(., 'url'),
         strict_evidences = strict_evidences,
         param =
             list(
                 json_param = list(),
                 query_type = query_l,
                 datasets = datasets,
                 resources = resources
             ) %>%
             add_qt_message
    ) %>%
    {`if`(
        query_l == 'annotations' && wide,
        pivot_annotations(.),
        .
    )} %>%
    {`if`(
        resource_l == 'dorothea',
        filter(
            .,
            str_detect(
                dorothea_level,
                dorothea_levels %>% paste(collapse = '') %>% sprintf('[%s]', .)
            )
        ),
        .
    )}

}
