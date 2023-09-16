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

utils::globalVariables(
    c("category", "uniprot", "genesymbol", "annotations",
    "target", "database", "category_intercell_source", "target_genesymbol",
    "category_intercell_target", "parent_intercell_source",
    "database_intercell_target","parent_intercell_target","source_genesymbol",
    "is_stimulation", "is_inhibition", "consensus_direction",
    "consensus_stimulation","consensus_inhibition","dip_url","sources",
    "references", "curation_effort", "dorothea_level", "n_references",
    "n_resources", "ncbi_tax_id_target", "ncbi_tax_id_source")
)

.omnipath_qt_synonyms <- list(
    ptms = 'enzsub',
    enz_sub = 'enzsub',
    complex = 'complexes'
)

.omnipath_qt_messages <- list(
    interactions = 'interactions',
    enzsub = 'enzyme-substrate relationships',
    complexes = 'protein complexes',
    annotations = 'annotation records',
    intercell = 'intercellular communication role records'
)

.omnipath_qt_nolicense <- c(
    'annotations_summary',
    'intercell_summary'
)

.omnipath_arg_synonyms <- list(
    select_organism = 'organisms',
    filter_databases = 'resources',
    databases = 'resources',
    select_genes = 'proteins'
)

.omnipath_default_fields <- list(
    enzsub = c('sources', 'references', 'curation_effort'),
    interactions = c('sources', 'references', 'curation_effort')
)

.omnipath_querystring_param <- c(
    'genesymbols',
    'resources',
    'datasets',
    'organisms',
    'dorothea_levels',
    'dorothea_methods',
    'source_target',
    'fields',
    'format',
    'directed',
    'signed',
    'loops',
    'enzymes',
    'substrates',
    'partners',
    'proteins',
    'entity_types',
    'sources',
    'targets',
    'residues',
    'scope',
    'aspect',
    'source',
    'categories',
    'parent',
    'transmitter',
    'receiver',
    'secreted',
    'plasma_membrane_transmembrane',
    'plasma_membrane_peripheral',
    'topology',
    'causality',
    'license',
    'password',
    'types'
)


.omnipath_querystring_synonyms <- list(
    organism = 'organisms',
    resource = 'resources',
    databases = 'resources',
    database = 'resources',
    dorothea_level = 'dorothea_levels',
    tfregulons_levels = 'dorothea_levels',
    tfregulons_level = 'dorothea_levels',
    genesymbol = 'genesymbols',
    field = 'fields',
    dataset = 'datasets',
    directions = 'directed',
    entity_type = 'entity_types',
    signs = 'signed',
    modification = 'types',
    type = 'types'
)


.omnipath_param_misc_keys <- c(
    'query_type',
    'default_fields',
    'silent',
    'logicals',
    'download_args',
    'references_by_resource',
    'add_counts',
    'qt_message',
    'exclude',
    'extra_attrs',
    'evidences',
    'json_param',
    'strict_evidences'
)


#' Downloads data from the OmniPath web service
#'
#' Generic method for retrieval of a table and creating a data frame.
#' All methods specific for certain query types or datasets use this function
#' to manage the download.
#' Not exported.
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom tibble as_tibble
#' @importFrom readr read_tsv cols col_character
#' @importFrom utils modifyList
#' @importFrom rlang !!!
#'
#' @noRd
import_omnipath <- function(
    query_type,
    organism = 9606,
    resources = NULL,
    datasets = NULL,
    genesymbols = 'yes',
    fields = NULL,
    default_fields = TRUE,
    silent = FALSE,
    logicals = NULL,
    download_args = list(),
    references_by_resource = TRUE,
    add_counts = TRUE,
    license = NULL,
    password = NULL,
    exclude = NULL,
    json_param = list(),
    strict_evidences = FALSE,
    ...
){

    datasets %<>% setdiff(exclude)
    resources %<>% setdiff(exclude)

    param <- c(as.list(environment()), list(...))
    param <- omnipath_check_param(param)

    url <-
        param %>%
        omnipath_build_url %>%
        c(`if`(
            getOption('omnipath.notls_fallback') &&
            !getOption('omnipath.notls_force'),
            omnipath_build_url(param, notls = TRUE),
            NULL
        ))
    download_args_defaults <- list(
        url = url
    )
    dataframe_defaults <- list(
        fun = read_tsv,
        col_types = `if`(
            'dorothea_level' %in% param$fields,
            cols(dorothea_level = col_character()),
            cols()
        ),
        progress = FALSE,
        show_col_types = FALSE
    )
    json_defaults <- list(
        fun = safe_json,
        simplifyDataFrame = FALSE
    )
    download_args <- modifyList(
        `if`(
            !is.null(param$format) && param$format == 'json',
            json_defaults,
            dataframe_defaults
        ),
        download_args
    )
    download_args <- modifyList(
        download_args_defaults,
        download_args
    )

    result <-
        do.call(omnipath_download, download_args) %>%
        omnipath_post_download(
            url = url,
            logicals = logicals,
            references_by_resource = references_by_resource,
            strict_evidences = strict_evidences,
            exclude = exclude,
            param = param,
            add_counts = add_counts,
            silent = silent
        )

    return(result)

}


#' Post-processing of the data downloaded from OmniPath
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom tibble as_tibble
#' @importFrom rlang !!!
#' @noRd
omnipath_post_download <- function(
        result,
        url,
        param,
        logicals = NULL,
        references_by_resource = TRUE,
        strict_evidences = FALSE,
        exclude = NULL,
        add_counts = TRUE,
        silent = FALSE
    ) {

    omnipath_check_result(result, url)

    result %<>% cast_logicals(logicals)
    result %<>% strip_resource_labels(references_by_resource)
    result %<>% apply_exclude(exclude)
    result %<>% deserialize_extra_attrs(!!!param$json_param)
    result %<>% deserialize_evidences(!!!param$json_param)

    if(strict_evidences && param$query_type == 'interactions') {
        result %<>% only_from(
            datasets = param$datasets,
            resources = param$resources
        )
    }

    if(param$query_type %in% c('interactions', 'enzsub') && add_counts){
        result %<>% count_references
        result %<>% count_resources
    }

    if(is.data.frame(result)){
        result %<>% as_tibble
    }

    from_cache <- result %>% is_from_cache

    # reporting and returning result
    loglevel <- `if`(
        silent,
        logger::DEBUG,
        logger::SUCCESS
    )

    msg <- '%soaded %d %s%s.'

    logger::log_level(
        level = loglevel,
        msg,
        `if`(from_cache, 'L', 'Downl'),
        `if`(
            is.data.frame(result),
            nrow(result),
            length(result)
        ),
        param$qt_message,
        `if`(from_cache, ' from cache', '')
    )

    return(result)

}


#' Checks the arguments of \link{import_omnipath}, corrects some easy to
#' confuse or deprecated synonyms and selects the message printed by
#' the download function.
#' Not exported.
#'
#' @importFrom magrittr %<>%
#' @importFrom logger log_warn
#'
#' @noRd
omnipath_check_param <- function(param){

    # mapping query type synonyms
    param$query_type <- `if`(
        !is.null(param$query_type) &
        param$query_type %in% names(.omnipath_qt_synonyms),
        .omnipath_qt_synonyms[[param$query_type]],
        param$query_type
    )

    param %<>% add_qt_message

    # mapping the query string parameter synonyms
    for(name in names(param)){
        if(
            name %in% names(.omnipath_querystring_synonyms) &&
            !.omnipath_querystring_synonyms[[name]] %in% names(param)
        ){
            param[[.omnipath_querystring_synonyms[[name]]]] <- param[[name]]
        }
    }

    # checking DoRothEA confidence level values
    if(
        'dorothea_levels' %in% names(param) &&
        !all(param$dorothea_levels %in% c('A', 'B', 'C', 'D'))
    ){
        msg <- 'DoRothEA confidence levels available are A, B, C and D.'
        log_warn(msg)
        warning(msg)
    }

    # extra_attrs and evidences are accepted also as an argument
    for(name in c('extra_attrs', 'evidences')) {

        if(if_null(param[[name]], FALSE)) {

            param$fields %<>% union(name)
            param[[name]] <- NULL

        }

    }

    if(param$strict_evidences) {
        param$fields %<>% union('evidences')
    }

    # adding default fields if not disabled
    param$fields <- `if`(
        param$default_fields &&
        param$query_type %in% names(.omnipath_default_fields),
        unique(
            c(
                param$fields,
                .omnipath_default_fields[[param$query_type]],
                `if`(
                    'dorothea' %in% param$datasets,
                    'dorothea_level',
                    NULL
                )
            )
        ),
        param$fields
    )

    # removing some fields according to query type
    if(!param$query_type %in% c('interactions', 'enzsub')){
        param$genesymbols <- NULL
        param$organisms <- NULL
    }

    # checking for wrong resource names
    if(!is.null(param$resources)){

        all_resources <- get_resources(param$query_type)
        unknown_resources <- setdiff(param$resources, all_resources)

        if(length(unknown_resources) != 0){

            warning(
                sprintf(
                    paste(
                        'The following resources are not available: %s.',
                        'Check the resource names for spelling mistakes.'
                    ),
                    paste0(unknown_resources, collapse = ', ')
                )
            )

        }
    }

    # setting up generic defaults from options
    for(opt in c('license', 'password')){

        param[opt] <- `if`(
            (
                is.null(param[[opt]]) &&
                !param$query_type %in% .omnipath_qt_nolicense
            ),
            options(sprintf('omnipath.%s', opt)),
            param[[opt]]
        )

    }

    return(param)

}


#' Adds a message printed upon successful download
#'
#' @noRd
add_qt_message <- function(param) {

    # adding the message template which will be printed upon successful
    # download
    param$qt_message <- `if`(
        !is.null(param$query_type) &
        param$query_type %in% names(.omnipath_qt_messages),
        .omnipath_qt_messages[[param$query_type]],
        'records'
    )

    return(param)

}


#' Constructs the URL by creating a base URL according to the query type and
#' adding all user or package defined query string parameters.
#' Not exported.
#'
#' @importFrom magrittr %>%
#' @importFrom logger log_warn
#'
#' @noRd
omnipath_build_url <- function(param, notls = FALSE){

    baseurl <- omnipath_url(param$query_type, notls = notls)

    unknown_param <- setdiff(
        names(param),
        unique(c(
            .omnipath_querystring_param,
            names(.omnipath_querystring_synonyms),
            .omnipath_param_misc_keys
        ))
    )

    if(length(unknown_param) > 0L){

        log_warn(
            'Unknown %s: %s.',
            unknown_param %>% plural('parameter'),
            unknown_param %>% pretty_list
        )

    }

    url <- Reduce(
        function(url, key){
            omnipath_url_add_param(url, key, param[[key]])
        },
        .omnipath_querystring_param,
        init = baseurl
    )

    return(url)

}


#' Appends a query string parameter to the URL.
#' Not exported, used internally for assembling the URLs.
#'
#' @noRd
omnipath_url_add_param <- function(url, name, values = NULL){

    values <- `if`(
        is.null(values),
        NULL,
        `if`(
            identical(values, TRUE),
            'yes',
            `if`(
                identical(values, FALSE),
                'no',
                values
            )
        )
    )

    url <- `if`(
        all(is.null(values)) || all(is.na(values)),
        url,
        sprintf(
            '%s%s%s=%s',
            url,
            `if`(grepl('?', url, fixed = TRUE), '&', '?'),
            name,
            paste(values, collapse = ',')
        )
    )

    return(url)

}


#' Checks whether the response is real data or an error message.
#' In case of error stops the execution and prints the URL and the message
#' from the server.
#'
#' @noRd
omnipath_check_result <- function(result, url){

    if(length(result) == 1){
        server_msg <- paste(result[[1]], collapse = '\n')
        stop(
            sprintf(
                'Failed to download data from OmniPath:\nURL: %s\n%s\n',
                url,
                server_msg
            )
        )
    }

}

#' Makes sure the boolean variables, listed in argument `logicals`, are of
#' R logical type. Converts various string and numeric representations.
#' Checks only for TRUE values, whatever does not match remains FALSE.
#'
#' @noRd
cast_logicals <- function(data, logicals = NULL){

    true_values <- c('True', '1', 'TRUE', 'T', 'yes', 'YES', 'Y', 'y')

    for(name in logicals){
        data[[name]] <- (
            identical(data[[name]], TRUE) |
            data[[name]] %in% true_values |
            (is.numeric(data[[name]]) & data[[name]] > 0)
        )
    }

    return(data)

}


#' Removes records which are only from resources to be excluded
#'
#' @param data A data frame from the OmniPath web service.
#' @param exclude Character vector with the resource names to exclude.
#'
#' @return The input data frame with records removed according to the
#'     exclude list.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom rlang !! sym
#' @importFrom purrr map_lgl
#' @noRd
apply_exclude <- function(data, exclude){

    col <- data %>% resources_colname

    data %>%
    `if`(
        is.null(exclude),
        .,
        filter(
            .,
            map_lgl(
                str_split(!!sym(col), ';'),
                function(x){x %>% setdiff(exclude) %>% length %>% as.logical}
            )
        )
    )

}


#' For a character vector splits each element and re-joins sorted unique
#' values.
#'
#' @noRd
split_unique_join <- function(
    x,
    sep = ';',
    outsep = sep,
    method = NULL
){

    method <- `if`(
        is.null(method),
        function(values, outsep, ...){
            paste(sort(unique(values)), collapse = outsep)
        },
        method
    )

    return(
        split_apply(
            x = x,
            method = method,
            sep = sep,
            outsep
        )
    )

}


#' For a character vector splits each element and applies a method for
#' each sub vector.
#'
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @noRd
split_apply <- function(
    x,
    method,
    sep = ';',
    ...
){

    x %>%
    strsplit(sep) %>%
    map(method, ...) %>%
    unlist()

}


#' For an interactions or enzyme-substrate data frame adds a column
#' `n_resources` with the number of resources for each record.
#'
#' @noRd
count_resources <- function(data, only_primary = TRUE){

    data[['n_resources']] <- split_apply(
        data$sources,
        method = function(values, only_primary){
            if(only_primary){
                values <- values[!grepl('_', values)]
            }
            return(length(values))
        },
        sep = ';',
        only_primary = only_primary
    )

    return(data)

}


#' For an interactions or enzyme-substrate data frame adds a column
#' `n_references` with the number of references for each record.
#'
#' @importFrom dplyr n_distinct mutate
#' @importFrom magrittr %>%
#' @noRd
count_references <- function(data){

    # NSE vs. R CMD check workaround
    n_references <- references <- NULL

    data %>%
    mutate(
        n_references = ifelse(
            is.na(references),
            0,
            strip_resource_labels(
                .,
                inplace = FALSE,
                method = n_distinct
            )
        )
    )

}

#' For each undirected interaction adds a duplicate with the source and
#' target nodes swapped.
#'
#' @noRd
swap_undirected <- function(data){

    is_directed <- NULL

    data <- data %>%
        filter(is_directed == 0) %>%
        rename(
            source = target,
            target = source,
            source_genesymbol = target_genesymbol,
            target_genesymbol = source_genesymbol
        ) %>%
        {`if`(
            'ncbi_tax_id_source' %in% names(.),
            rename(
                .,
                ncbi_tax_id_source = ncbi_tax_id_target,
                ncbi_tax_id_target = ncbi_tax_id_source
            ),
            .
        )} %>%
        bind_rows(data)

    return(data)

}


#' Downloader dedicated to OmniPath web service URLs
#'
#' Just a thin wrapper around \code{download_base}.
#'
#' @param url Character: the URL to download. Elements after the first will
#'     be used as fallback URLs in case the first one fails.
#' @param fun The downloader function. Should be able to accept \code{url}
#'     as its first argument.
#' @param ... Passed to the internal function \code{download_base} and
#'     from there ultimately to \code{fun}.
#'
#' @importFrom logger log_trace log_info log_error
#' @noRd
omnipath_download <- function(url, fun, ...) {

    for(the_url in url) {

        log_trace('Looking up in cache: `%s`', the_url)
        from_cache <- omnipath_cache_load(url = the_url)

        if(!is.null(from_cache)){

            log_info('Loaded from cache: `%s`', the_url)
            return(from_cache)

        }

    }

    for(the_url in url){

        log_trace('Attempting `%s`', the_url)

        result <- tryCatch(
            download_base(the_url, fun, ...),
            error = function(e) {
                log_warn(
                    'Failed to download: `%s`; error: %s',
                    the_url,
                    conditionMessage(e)
                )
            }
        )

        if(!is.null(result)){

            log_info('Successfully retrieved: `%s`', the_url)
            omnipath_cache_save(data = result, url = the_url)
            return(result)

        }

    }

}
