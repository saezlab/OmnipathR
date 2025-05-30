#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
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
#'
#' @examples
#' resources(query_type = "interactions")
#' @aliases get_resources
resources <- function(
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

    resources_url <- omnipath_url('resources')
    resources <- safe_json(path = resources_url)

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


# Aliases (old names) to be Deprecated
#' @rdname complexes
#' @param ... Passed to \code{resources}.
#' @export
#'
#' @noRd
get_resources <- function(...){
    .Deprecated('resources')
    resources(...)
}


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
#' interactions <- omnipath()
#' signor <- filter_by_resource(interactions, resources = "SIGNOR")
#'
#' @importFrom logger log_success
#' @export
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

        log_success(
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
#' co <- complexes()
#' resources_colname(co)
#' # [1] "sources"
#'
#' @importFrom magrittr %>% extract
#' @export
resources_colname <- function(data){

    intersect(
        c('sources', 'database', 'source'),
        colnames(data)
    ) %>%
    extract(1L)

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


#' OmniPath resource information
#'
#' The `resources` query type provides resource metadata in JSON format.
#' Here we retrieve this JSON and return it as a nested list structure.
#'
#' @return A nested list structure with resource metadata.
#'
#' @examples
#' resource_info()
#'
#' @export
resource_info <- function() {

    omnipath_query('resources', format = 'json')

}


#' OmniPath query parameters
#'
#' All parameter names and their possible values for a query type. Note:
#' parameters with `NULL` values have too many possible values to list
#' them.
#'
#' @param query_type Character: interactions, annotations, complexes, enz_sub
#'     or intercell.
#'
#' @return A named list with the parameter names and their possible values.
#'
#' @examples
#' ia_param <- query_info('interactions')
#' ia_param$datasets[1:5]
#' # [1] "dorothea"    "kinaseextra" "ligrecextra" "lncrna_mrna" "mirnatarget"
#'
#' @importFrom magrittr %>%
#' @export
query_info <- function(query_type) {

    query_type %>%
    sprintf('queries/%s', .) %>%
    omnipath_query(format = 'json')

}


#' Summary of the annotations and intercell database contents
#'
#' The `annotations_summary` and `intercell_summary` query types return
#' detailed information on the contents of these databases. It includes
#' all the available resources, fields and values in the database.
#'
#' @param query_type Character: either "annotations" or "intercell".
#' @param return_df Logical: return a data frame instead of list.
#'
#' @return Summary of the database contents: the available resources, fields,
#'     and their possible values. As a nested list if format is "json",
#'     otherwise a data frame.
#'
#' @examples
#' annotations_summary <- database_summary('annotations')
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map keep
#' @importFrom stringr str_split
#' @importFrom rlang set_names
#' @importFrom tidyr unnest_wider
#' @importFrom tibble tibble
#' @export
database_summary <- function(query_type, return_df = FALSE) {

    query_type %>%
    sprintf('%s_summary', .) %>%
    omnipath_query(format = 'json') %>%
    map(~modifyList(., set_names(str_split(.$value, '#'), 'value'))) %>%
    map(~modifyList(., list(value = keep(.$value, ~nchar(.x) > 0L)))) %>%
    `if`(return_df, tibble(.) %>% unnest_wider(1L), .)

}


#' Collect resource names from a data frame
#'
#' @param data A data frame from an OmniPath query.
#'
#' @return Character: resource names occuring in the data frame.
#'
#' @examples
#' pathways <- omnipath_interactions()
#' resources_in(pathways)
#'
#' @importFrom magrittr %>%
#' @export
resources_in <- function(data) {

    data %>%
    {`if`(
        has_evidences(.),
        resources_in_evidences(.),
        resources_in_simple(.)
    )}

}


#' Collect resource names from a data frame
#'
#' From the resources column of data frames returned by OmniPath queries.
#'
#' @param data A data frame from an OmniPath query.
#'
#' @return Character: resource names occuring in the data frame.
#'
#' @importFrom magrittr %>% extract
#' @noRd
resources_in_simple <- function(data) {

    col <- data %>% resources_colname

    data %>%
    pull(!!sym(col)) %>%
    str_split(';') %>%
    unique_sorted

}


#' Collect resource names from a data frame with evidences column
#'
#' @param data A data frame with "evidences" column, typically an interactions
#'     data frame.
#'
#' @return Character: resource names occuring in the data frame.
#'
#' @importFrom magrittr %>% extract
#' @importFrom dplyr pull
#' @importFrom purrr map keep
#' @noRd
resources_in_evidences <- function(data) {

    data %>%
    pull(evidences) %>%
    map(
        ~map(
            keep(.x, ~is.list(.x) && length(.x) > 0L),
            ~map(
                .x,
                ~paste(
                    unlist(extract(.x, c('resource', 'via'))),
                    collapse = '_'
                )
            )
        )
    ) %>%
    unique_sorted

}

#' Resources shared between the database and the query set
#'
#' @param database Character: resource names in the database (OmniPath); these
#'     names will be split by underscore to match the names of primary and
#'     secondary resources on their own; e.g. "SIGNOR_ProtMapper" will match
#'     both "SIGNOR", "ProtMapper" and "SIGNOR_ProtMapper".
#' @param query Character: resource names of interest; these will be left
#'     intact and matched against the first set.
#'
#' @return Character: resource names occuring both in database and query.
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @noRd
match_resources <- function(database, query) {

    database %>%
    c(str_split(., '_') %>% unlist) %>%
    intersect(query)

}
