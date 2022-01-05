#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2022
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


#' Retrieve the available resources for a given query type
#'
#' Collects the names of the resources available in OmniPath for a certain
#' query type and optionally for a dataset within that.
#'
#' @param query_type one of the query types `interactions`, `enz_sub`,
#'     `complexes`, `annotations` or `intercell`
#' @param datasets currently within the `interactions` query type only,
#'     multiple datasets are available: `omnipath`, `kinaseextra`,
#'     `pathwayextra`, `ligrecextra`, `dorothea`, `tf_target`, `tf_mirna`,
#'     `mirnatarget` and `lncrna_mrna`.
#' @param generic_categories for the `intercell` query type, restrict the
#'     search for some generic categories e.g. `ligand` or `receptor`.
#'
#' @return a character vector with resource names
#'
#' @export
#' @importFrom jsonlite fromJSON
#'
#' @examples
#' get_resources(query_type = 'interactions')
get_resources <- function(
    query_type,
    datasets = NULL,
    generic_categories = NULL
){

    null_or_matches <- function(
        res_data,
        values,
        key = deparse(substitute(values))
    ){

        qt_data <- res_data[['queries']][[query_type]]

        return(
            is.null(values) || (
                key %in% names(qt_data) &&
                length(intersect(qt_data[[key]], values)) > 0
            )
        )

    }

    query_type <- `if`(
        query_type %in% names(.omnipath_qt_synonyms),
        .omnipath_qt_synonyms[[query_type]],
        query_type
    )

    resources <- jsonlite::fromJSON(txt = 'https://omnipathdb.org/resources')

    return(
        sort(Filter(
            function(resource){
                query_type %in% names(resources[[resource]]$queries) &&
                null_or_matches(resources[[resource]], datasets) &&
                null_or_matches(resources[[resource]], generic_categories)
            },
            names(resources)
        ))
    )

}


########## ########## ########## ##########
########## RESOURCE FILTERING      ########
########## ########## ########## ##########
## Non exported functions (package internal functions) to filter PTMs,
## interactions, complexes and annotations according to the databases passed
## to the main functions

#' Filters OmniPath data by resources
#'
#' Keeps only those records which are supported by any of the resources of
#' interest.
#'
#' @param data A data frame downloaded from the OmniPath web service
#'     (interactions, enzyme-substrate or complexes).
#' @param resources Character vector with resource names to keep.
#'
#' @return The data frame filtered.
#'
#' @examples
#' interactions <- import_omnipath_interactions()
#' signor <- filter_by_resource(interactions, resources = 'SIGNOR')
#'
#' @importFrom logger log_success
#' @export
#'
#' @aliases filter_sources
filter_by_resource <- function(data, resources = NULL){

    if(!is.null(resources)){

        before <- nrow(data)

        field <- data %>% resources_colname

        data <- data[
            which(
                unlist(lapply(
                    strsplit(data[[field]], ';'),
                    function(res){
                        length(intersect(res, resources)) > 0
                    }
                ))
            ),
        ]

        after <- nrow(data)

        logger::log_success(
            'Filtering by resources: removed %d records.',
            before - after
        )

    }

    return(data)
}


#' Name of the column with the resources
#'
#' Unfortunately the column title is different across the various
#' query types in the OmniPath web service, so we need to guess.
#'
#' @param data A data frame downloaded by any \code{import_...} function
#'     in the current package.
#'
#' @return Character: the name of the column, if any of the column names
#'     matches.
#'
#' @examples
#' co <- import_omnipath_complexes()
#' resources_colname(co)
#' # [1] "sources"
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr first
#' @export
resources_colname <- function(data){

    intersect(
        c('sources', 'database', 'source'),
        colnames(data)
    ) %>%
    first

}


#' Alias for old function name.
#' @param ... Passed to \code{filter_by_resource}.
#'
#' @noRd
filter_sources <- function(...){
    .Deprecated("filter_by_resource")
    filter_by_resource(...)
}


#' Removes the resource labels from references (PubMed IDs) in the
#' interactions and enzyme-substrate data frames.
#'
#' @noRd
strip_resource_labels <- function(
    data,
    references_by_resource = FALSE,
    colname = 'references',
    inplace = TRUE,
    method = NULL
){

    if(!references_by_resource && colname %in% names(data)){

        result <- split_unique_join(
            gsub(
                '[-\\w]*:?(\\d+)',
                '\\1',
                data[[colname]],
                perl = TRUE
            ),
            method = method
        )

        if(inplace){
            data[[colname]] <- result
        }else{
            return(result)
        }

    }

    return(data)

}
