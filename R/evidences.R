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

EVIDENCES_KEYS <- c('positive', 'negative', 'directed', 'undirected')


#' Converts the evidences column from JSON encoded to list
#'
#' @param data A data frame from the OmniPath web service.
#' @param ... Passed to `jsonlite::fromJSON`.
#'
#' @return The input data frame with the "evidences" column converted
#'     to list.
#'
#' @importFrom magrittr %>%
#' @noRd
deserialize_evidences <- function(data, ...){

    data %>%
    deserialize_json_col('evidences', ...)

}


#' Tells if a data frame has a column "evidences"
#'
#' @importFrom magrittr %>% or
#' @noRd
has_evidences <- function(data, wide_ok = FALSE){

    data %>%
    has_column('evidences') %>%
    or(wide_ok && has_evidences_wide(data))

}


#' Tells if a data frame has evidence columns by direction and sign
#'
#' @importFrom magrittr %>% is_in
#' @noRd
has_evidences_wide <- function(data) {

    EVIDENCES_KEYS %>% is_in(colnames(data)) %>% all

}


#' Check for "evidences" column, throw an error if it is missing
#'
#' @importFrom magrittr %>% extract
#' @noRd
must_have_evidences <- function(data, wide_ok = FALSE, env = parent.frame()){

    if(!has_evidences(data, wide_ok = wide_ok)) {

        parent_call <- sys.call(-1L)

        m <-
            paste(
                '`%s` can be called only on data',
                'frames with `evidences` column.'
            ) %>%
            sprintf(parent_call %>% as.character %>% extract(1L))

        log_error(m)
        do.call(stop, list(m), envir = env)

    }

}


#' Separate evidences by direction and effect sign
#'
#' @param data An interaction data frame with "evidences" column.
#' @param longer Logical: If TRUE, the "evidences" column is split into rows.
#'
#' @return The data frame with new columns or new rows by direction and sign.
#'
#' @examples
#' omnipath <- import_omnipath_interactions(fields = 'evidences')
#' omnipath <- unnest_evidences(omnipath)
#' colnames(omnipath)
#'
#' @importFrom magrittr %>% extract
#' @importFrom dplyr mutate
#' @importFrom purrr map
#' @importFrom tidyr unnest_longer unnest_wider
#' @importFrom rlang exec !!!
#' @export
unnest_evidences <- function(data, longer = FALSE) {

    must_have_evidences(data)

    unnest_method <- `if`(longer, unnest_longer, unnest_wider)
    unnest_args <- `if`(longer, list(indices_to = 'direction'), list())

    data %>%
    mutate(evidences = map(evidences, extract, EVIDENCES_KEYS)) %>%
    exec(unnest_method, ., 'evidences', !!!unnest_args)

}


#' Filter evidences by dataset and resource
#'
#' @param data An interaction data frame with some columns containing evidences
#'      as nested lists.
#' @param ... The evidences columns to filter: tidyselect syntax is supported.
#'      By default the columns "evidences", "positive", "negative", "directed"
#'      and "undirected" are filtered, if present.
#' @param datasets A character vector of dataset names.
#' @param resources A character vector of resource names.
#'
#' @return The input data frame with the evidences in the selected columns
#'    filtered.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang expr
#' @importFrom tidyselect eval_select
#' @importFrom dplyr mutate across
#' @importFrom logger log_trace
#' @export
filter_evidences <- function(data, ..., datasets = NULL, resources = NULL) {

    columns <-
        expr(...) %>%
        eval_select(data) %>%
        names %>%
        if_null_len0(
            EVIDENCES_KEYS %>%
            c('evidences') %>%
            intersect(colnames(data))
        )

    log_trace(
        'Filtering evidence columns: %s; to datasets: %s; and resources: %s',
        paste(columns, collapse = ', '),
        paste(if_null(datasets, 'any'), collapse = ', '),
        paste(if_null(resources, 'any'), collapse = ', ')
    )

    data %>%
    mutate(across(columns, ~filter_evs_lst(.x, datasets, resources, columns)))

}


#' Filter evidences within a nested list
#'
#' @importFrom magrittr %>% is_in not
#' @importFrom purrr map2 keep map
#' @noRd
filter_evs_lst <- function(lst, datasets, resources, columns) {

    filter_evs <- function(x) {

        x %>%
        {`if`(
            names(.) %>% is.null %>% not,
            map2(
                 names(.),
                 .,
                 function(name, value) {
                     `if`(
                          name %>% is_in(EVIDENCES_KEYS),
                          value %>% filter_evs,
                          value
                     )
                 }
            ),
            keep(
                .,
                function(ev) {
                    (is.null(datasets) || ev$dataset %in% datasets) &&
                    (
                        is.null(resources) ||
                        ev$resource %in% resources ||
                        if_null(ev$via, '') %in% resources
                    )
                }
            )
        )}

    }


    lst %>%
    map(filter_evs)

}


#' Recreate interaction data frame based on certain datasets and resources
#'
#' @param data An interaction data frame from the OmniPath web service with
#'     evidences column.
#' @param datasets Character: a vector of dataset labels. Only evidences from
#'     these datasets will be used.
#' @param resources Character: a vector of resource labels. Only evidences
#'     from these resources will be used.
#'
#' @return A copy of the interaction data frame restricted to the given
#'     datasets and resources.
#'
#' @details
#' The OmniPath interactions database fully integrates all attributes from all
#' resources for each interaction. This comes with the advantage that
#' interaction data frames are ready for use in most of the applications;
#' however, it makes it impossible to know which of the resources and
#' references support the direction or effect sign of the interaction. This
#' information can be recovered from the "evidences" column. The "evidences"
#' column preserves all the details about interaction provenances. In cases
#' when you want to use a faithful copy of a certain resource or dataset, this
#' function will help you do so. Still, in most of the applications the best is
#' to use the interaction data as it is returned by the web service.
#'
#' @examples
#' ci <- collectri()
#' ci <- only_from(ci, datasets = 'collectri')
#'
#' @importFrom magrittr %>%
#' @importFrom rlang as_function
#' @importFrom purrr keep
#' @importFrom dplyr select
#' @importFrom tidyselect any_of
#' @export
only_from <- function(data, datasets = NULL, resources = NULL) {

    if(is.null(datasets) && is.null(resources)) {
        return(data)
    }

    must_have_evidences(data, wide_ok = TRUE)

    has_wide <- data %>% has_evidences_wide

    log_trace(
        'Restricting interaction records to datasets: %s; and resources %s',
        paste(if_null(datasets, 'any'), collapse = ', '),
        paste(if_null(resources, 'any'), collapse = ', ')
    )

    data %>%
    {`if`(has_wide, ., unnest_evidences(.))} %>%
    filter_evidences(datasets = datasets, resources = resources) %>%
    from_evidences() %>%
    {`if`(has_wide, ., select(., -any_of(EVIDENCES_KEYS)))}

}


#' Recreate interaction records from evidences columns
#'
#' @details
#' The OmniPath interaction data frames specify interactions primarily by
#' three columns: "is_directed", "is_stimulation" and "is_inhibition".
#' Besides these, there are the "sources" and "references" columns that are
#' always included in data frames created by OmnipathR and list the resources
#' and literature references for each interaction, respectively. The optional
#' "evidences" column is required to find out which of the resources and
#' references support the direction or effect sign of the interaction. To
#' properly recover information for arbitrary subsets of resources or
#' datasets, the evidences can be filtered first, and then the standard
#' data frame columns can be reconstructed from the selected evidences.
#' This function is able to do the latter. It expects either an "evidences"
#' column or evidences in their wide format 4 columns layout. It overwrites
#' the standard columns of interaction records based on data extracted from
#' the evidences.
#'
#' @importFrom stringr str_detect
#' @noRd
from_evidences <- function(data) {

    must_have_evidences(data, wide_ok = TRUE)

    prefix <- data$references[1L] %>% str_detect('^\\w+:')

    data %>%
    {`if`(has_evidences_wide(.), ., unnest_evidences(.))} %>%
    mutate(
        is_directed = as.integer(lgl_from(., positive, negative, directed)),
        is_stimulation = as.integer(map_lgl(positive, ~length(.x) > 0L)),
        is_inhibition = as.integer(map_lgl(negative, ~length(.x) > 0L)),
        sources = resources_from(., EVIDENCES_KEYS),
        references = references_from(., EVIDENCES_KEYS, prefix = prefix)
    )

}


#' Logical vector from a set of list columns
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom purrr pmap_lgl
#' @noRd
lgl_from <- function(data, ...) {

    data %>%
    select(...) %>%
    pmap_lgl(~length(c(...)) > 0L)

}


#' Extract resources from a set of evidence list columns
#'
#' @importFrom magrittr %>%
#' @noRd
resources_from <- function(data, ..., collapse = ';') {

    extract_resources <- function(ev) {
        paste0(ev$resource, `if`(is.null(ev$via), '', '_'), ev$via)
    }

    data %>%
    chr_from(..., fn = extract_resources, collapse = collapse)

}


#' Extract references from a set of evidence list columns
#'
#' @importFrom magrittr %>%
#' @noRd
references_from <- function(data, ..., prefix = TRUE, collapse = ';') {

    extract_references <- function(ev) {
        `if`(
            prefix,
            paste(ev$resource, ev$references, sep = ':'),
            ev$references, collapse = collapse
        ) %>%
        unique %>%
        sort %>%
        paste(collapse = collapse)
    }

    data %>%
    chr_from(..., fn = extract_references, collapse = collapse)

}


#' Character vector from a set of list columns
#'
#' @importFrom magrittr %>%
#' @importFrom purrr pmap_chr map_chr
#' @importFrom dplyr select
#' @noRd
chr_from <- function(data, ..., fn, collapse = ';') {

    map_fn <- function(...) {
        c(...) %>%
        map_chr(fn) %>%
        unique %>%
        sort %>%
        paste(collapse = collapse) %>%
        if_null_len0('')
    }

    data %>%
    select(...) %>%
    pmap_chr(map_fn)

}
