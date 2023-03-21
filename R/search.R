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


#' Search index for OmniPath: load or build
#'
#' @param rebuild Logical: rebuild the index even if it could be loaded
#'     from the cache.
#'
#' @return Data frame of six columns: resource, dataset, key, value, where,
#'     database.
#'
#' @importFrom magrittr %>%
#' @noRd
build_index <- function(rebuild = FALSE) {

    `if`(rebuild, NULL, omnipath_cache_load('omnipath_search_index')) %>%
    if_null(build_index_impl())

}


#' Search index for OmniPath
#'
#' @return Data frame of six columns: resource, dataset, key, value, where,
#'     database.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom purrr keep map
#' @importFrom dplyr mutate bind_rows
#' @importFrom stringr str_to_lower str_extract
#' @importFrom logger log_success
#' @noRd
build_index_impl <- function() {

    log_success('Building search index, this may take a few minutes.')

    getNamespace('OmnipathR') %>%
    names %>%
    keep(endsWith, '_search_index') %>%
    map(
        function(index_builder){
            get(index_builder)() %>%
            mutate(database = str_extract(index_builder, '[^_]+'))
        }
    ) %>%
    bind_rows() %>%
    mutate(value_lower = str_to_lower(value)) %>%
    omnipath_cache_save(url = 'omnipath_search_index') %T>%
    {log_success('Finished building search index.')}

}


#' Search index for the OmniPath Annotations database
#'
#' @return Data frame of five columns: resource, dataset, key, value, where.
#'
#' @importFrom tidyr unnest
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows select distinct rename mutate
#' @importFrom logger log_info
#' @noRd
annotations_search_index <- function() {

    # NSE vs R CMD check workaround
    value <- label <- source <- resource <- NULL

    log_info('Building search index for the annotations database.')

    database_summary('annotations', return_df = TRUE) %>%
    unnest('value') %>%
    rename(resource = source) %>%
    mutate(dataset = NA) %>%
    bind_rows(

        select(., -value) %>%
        distinct %>%
        rename(value = label) %>%
        mutate(where = 'field_name', key = NA),

        rename(., key = label) %>%
        mutate(where = 'field_value')

    )

}


#' Search index for the OmniPath Intercell database
#'
#' @return Data frame of four columns: resource, dataset, value, where.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename bind_rows mutate distinct select
#' @importFrom logger log_info
#' @noRd
intercell_search_index <- function() {

    # NSE vs R CMD check workaround
    database <- category <- parent <- resource <- NULL

    log_info('Building search index for the intercell database.')

    database_summary('intercell', return_df = TRUE) %>%
    rename(resource = database) %>%
    mutate(dataset = NA) %>%
    bind_rows(

        rename(., key = parent, value = category) %>%
        distinct %>%
        mutate(where = 'category'),

        select(., -category) %>%
        rename(value = parent) %>%
        distinct %>%
        mutate(where = 'parent_category', key = NA)

    )

}

#' Search index for the OmniPath Interactions (Network) database
#'
#' @return Data frame of five columns: resource, dataset, key, value, where.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr unnest_longer separate_rows separate
#' @importFrom dplyr select rename bind_rows mutate distinct
#' @importFrom stringr str_extract str_replace
#' @importFrom logger log_info
#' @noRd
interactions_search_index <- function() {

    # NSE vs R CMD check workaround
    dataset <- datasets <- sources <- type <- resource <- NULL

    log_info('Building search index for the interactions database.')

    import_all_interactions(fields = c('extra_attrs', 'datasets')) %>%
    datasets_one_column(remove_logicals = TRUE) %>%
    select(resource = sources, type, extra_attrs, dataset = datasets) %>%
    unnest_longer(
        extra_attrs,
        indices_to = 'attr_key',
        transform = as.character
    ) %>%
    unnest_longer(extra_attrs) %>%
    unnest_longer(dataset) %>%
    separate_rows(resource, sep = ';') %>%
    bind_rows(

        # key is attr name, value is attr value
        select(., -resource, -type) %>%
        separate(
            attr_key,
            into = c('resource', 'key'),
            sep = '_',
            extra = 'merge'
        ) %>%
        rename(value = extra_attrs) %>%
        distinct %>%
        mutate(where = 'extra_attr'),

        # key is interaction type, value is extra attr name
        select(., dataset, key = type, value = attr_key) %>%
        mutate(
            resource = str_extract(value, '[^_]+'),
            where = 'extra_attr_name'
        ) %>%
        distinct,

        # key is NA, value is interaction type
        select(., resource, dataset, value = type) %>%
        distinct %>%
        mutate(key = NA, where = 'interaction_type'),

        # key is NA, value is secondary resource name
        select(., resource, dataset) %>%
        distinct %>%
        mutate(
            key = NA,
            value = str_replace(resource, '[^_]+_(.*)', '\\1'),
            resource = str_extract(resource, '[^_]+'),
            where = 'interaction_secondary_resource'
        ) %>%
        distinct

    )

}


#' Search index for the OmniPath Complexes database
#'
#' @return Data frame of five columns: resource, dataset, key, value, where.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr separate_rows
#' @importFrom dplyr select filter mutate distinct
#' @importFrom logger log_info
#' @noRd
complexes_search_index <- function() {

    log_info('Building search index for the complexes database.')

    import_omnipath_complexes() %>%
    select(resource = sources, value = name) %>%
    filter(!is.na(value)) %>%
    separate_rows(resource, sep = ';') %>%
    distinct %>%
    mutate(dataset = NA, key = NA, where = 'complex_name')

}


#' Search index for the OmniPath Enzyme-substrate database
#'
#' @return Data frame of five columns: resource, dataset, key, value, where.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr separate_rows
#' @importFrom dplyr select mutate distinct bind_rows
#' @importFrom stringr str_extract str_replace
#' @importFrom logger log_info
#' @noRd
enzsub_search_index <- function() {

    log_info('Building search index for the enzyme-substrate database.')

    import_omnipath_enzsub() %>%
    select(resource = sources, value = modification) %>%
    separate_rows(resource, sep = ';') %>%
    distinct %>%
    bind_rows(

        # key is NA, value is PTM type
        mutate(., dataset = NA, key = NA, where = 'ptm_type'),

        # key is NA, value is secondary resource name
        mutate(
            .,
            key = NA,
            value = str_replace(resource, '[^_]+_(.*)', '\\1'),
            resource = str_extract(resource, '[^_]+'),
            where = 'enzsub_secondary_resource'
        )

    )

}


#' Access the OmniPath search index
#'
#' @param rebuild Logical: rebuild the index even if it could be loaded
#'     from the cache.
#'
#' @return Data frame of six columns: resource, dataset, key, value, where,
#'     database.
#'
#' @noRd
get_index <- function(rebuild = FALSE) {

    get_db('search_index', param = list(rebuild = rebuild))

}


#' Search in OmniPath
#'
#' @param rebuild_index Logical: rebuild the index even if it could be loaded
#'     from the cache.
#' @param ... Ignored.
#'
#' Search for keywords in the resource and attribute names and metadata.
#'
search_omnipath <- function(..., rebuild_index = FALSE) {

    log_success('Search does not work yet.')
    index <- get_index(rebuild = rebuild_index)

}
