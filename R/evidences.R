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
#' @importFrom magrittr %>%
#' @noRd
has_evidences <- function(data){

    data %>% has_column('evidences')

}


#' Check for "evidences" column, throw an error if it is missing
#'
#' @importFrom magrittr %>% extract
#' @noRd
must_have_evidences <- function(data, env = parent.frame()){

    if(!has_evidences(data)) {

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

