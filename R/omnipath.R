#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2021
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


########## ########## ########## ##########
########## Generic non exported functions #
########## ########## ########## ##########

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
    'modification',
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
    'password'
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
    signs = 'signed'
)


.omnipath_param_misc_keys <- c(
    'query_type',
    'default_fields',
    'silent',
    'logicals',
    'download_args',
    'references_by_resource',
    'add_counts',
    'qt_message'
)


#' Downloads data from the OmniPath web service
#'
#' Generic method for retrieval of a table and creating a data frame.
#' All methods specific for certain query types or datasets use this function
#' to manage the download.
#' Not exported.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom readr read_tsv cols
#' @importFrom utils modifyList
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
    ...
){

    param <- c(as.list(environment()), list(...))
    param <- omnipath_check_param(param)

    url <- omnipath_url(param)
    download_args_defaults <- list(
        url = url
    )
    dataframe_defaults <- list(
        fun = read_tsv,
        col_types = cols()
    )
    json_defaults <- list(
        fun = jsonlite::fromJSON
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

    result <- do.call(omnipath_download, download_args)

    omnipath_check_result(result, url)

    msg <- '%soaded %d %s%s.'

    result <- cast_logicals(result, logicals)
    result <- strip_resource_labels(result, references_by_resource)
    if(param$query_type %in% c('interactions', 'enzsub') && add_counts){
        result <- count_references(result)
        result <- count_resources(result)
    }
    if(is.data.frame(result)){
        result %<>% as_tibble
    }
    from_cache <- result %>% is_from_cache

    loglevel <- `if`(
        silent,
        logger::DEBUG,
        logger::SUCCESS
    )

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

    # adding the message template which will be printed upon successful
    # download
    param$qt_message <- `if`(
        !is.null(param$query_type) &
        param$query_type %in% names(.omnipath_qt_messages),
        .omnipath_qt_messages[[param$query_type]],
        'records'
    )

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
            is.null(param[[opt]]),
            options(sprintf('omnipath.%s', opt)),
            param[[opt]]
        )

    }

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
omnipath_url <- function(param){

    baseurl <- options('omnipath.url')
    baseurl <- sprintf('%s%s', baseurl, param$query_type)

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
            outsep = outsep
        )
    )

}


#' For a character vector splits each element and applies a method for
#' each sub vector.
#'
#' @importFrom purrr map
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
#' @noRd
count_references <- function(data){

    data[['n_references']] <- strip_resource_labels(
        data,
        inplace = FALSE,
        method = function(refs, ...){
            length(unique(refs))
        }
    )

    return(data)

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

########## ########## ########## ##########
########## Enzyme-substrate      ##########
########## ########## ########## ##########

#' Imports enzyme-substrate relationships from OmniPath
#'
#' Imports the enzyme-substrate (more exactly, enzyme-PTM) relationship
#' database from \url{https://omnipathdb.org/enzsub}
#'
#' @return A data frame containing the information about ptms
#'
#' @param resources PTMs not reported in these databases are
#'     removed. See \code{\link{get_ptms_databases}} for more information.
#' @param organism PTMs are available for human, mouse and rat.
#'     Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields You can define here additional fields to be added to the
#'     result. If used, set the next argument, \code{default_fields}, to
#'     \code{FALSE}.
#' @param default_fields Whether to include the default fields (columns) for
#'     the query type. If \code{FALSE}, only the fields defined by the user
#'     in the \code{fields} argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#'     from the references (PubMed IDs); this way the information which
#'     reference comes from which resource will be lost and the PubMed IDs
#'     will be unique.
#' @param ... Optional additional arguments.
#'
#' @examples
#' enzsub <- import_omnipath_enzsub(
#'     resources = c('PhosphoSite', 'SIGNOR'),
#'     organism = 9606
#' )
#'
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_enzsub_resources}}}
#'     \item{\code{\link{import_omnipath_interactions}}}
#'     \item{\code{\link{enzsub_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_Omnipath_PTMS import_OmniPath_PTMS
import_omnipath_enzsub <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'enzsub',
        resources = resources,
        organism = organism,
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}

# Aliases (old names) to be Deprecated
#' @rdname import_omnipath_enzsub
#' @param ... Passed to \code{import_omnipath_enzsub}.
#' @export
#'
#' @noRd
import_Omnipath_PTMS <- function(...){
    .Deprecated("import_omnipath_enzsub")
    import_omnipath_enzsub(...)
}


# Aliases (old names) to be Deprecated
#' @rdname import_omnipath_enzsub
#' @param ... Passed to \code{import_omnipath_enzsub}.
#' @export
#'
#' @noRd
import_OmniPath_PTMS <- function(...){
    .Deprecated("import_omnipath_enzsub")
    import_omnipath_enzsub(...)
}


#' Retrieves a list of enzyme-substrate resources available in OmniPath
#'
#' Get the names of the enzyme-substrate relationship resources available
#' in \url{https://omnipath.org/enzsub}
#'
#' @param dataset ignored for this query type
#' @return character vector with the names of the enzyme-substrate resources
#' @export
#'
#' @examples
#' get_enzsub_resources()
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_omnipath_enzsub}}}
#' }
#'
#' @aliases get_ptms_databases
get_enzsub_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'enzsub', datasets = dataset))

}


# Aliases (old names) to be deprecated
#' @rdname get_enzsub_resources
#' @param ... Passed to \code{get_enzsub_resources}.
#' @export
#'
#' @noRd
get_ptms_databases <- function(...){
    .Deprecated("get_enzsub_resources")
    get_enzsub_resources(...)
}

########## ########## ########## ##########
########## INTERACTIONS          ##########
########## ########## ########## ##########

## Functions for importing interactions from OmniPath.
## The interactions database of OminPath consists of several datastes.

#' Imports interactions from the `omnipath` dataset of Omnipath
#'
#' Imports the database from \url{https://omnipathdb.org/interactions}, which
#' contains only interactions supported by literature references.
#' This part of the interaction database compiled a similar way as it has
#' been presented in the first paper describing OmniPath (Turei et al. 2016).
#'
#' @return A dataframe of protein-protein interactions
#' @export
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param datasets Names of the interaction datasts to download: omnipath 
#' (by default). Other possiblites are: pathwayextra, kinaseextra,
#' ligrecextra, dorothea,tf_target, mirnatarget, tf_mirna, lncrna_mrna.
#' The user can select multiple datasets as for example: c('omnipath',
#' 'pathwayextra', 'kinaseextra')
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments 
#'
#' @examples
#' interactions = import_omnipath_interactions(
#'     resources = c('SignaLink3'),
#'     organism = 9606
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_Omnipath_Interactions import_OmniPath_Interactions
import_omnipath_interactions <- function(
    resources = NULL,
    organism = 9606,
    datasets = 'omnipath',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = datasets,
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}

# Aliases (old names) to be deprecated
#' @rdname import_omnipath_interactions
#' @param ... Passed to \code{import_omnipath_interactions}.
#' @export
#'
#' @noRd
import_Omnipath_Interactions <- function(...){
    .Deprecated("import_omnipath_interactions")
    import_omnipath_interactions(...)
}


# Aliases (old names) to be deprecated
#' @rdname import_omnipath_interactions
#' @param ... Passed to \code{import_omnipath_interactions}.
#' @export
#'
#' @noRd
import_OmniPath_Interactions <- function(...){
    .Deprecated("import_omnipath_interactions")
    import_omnipath_interactions(...)
}


#' Imports interactions from the `pathway extra` dataset of Omnipath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=pathwayextra},
#' which contains activity flow interactions without literature reference.
#' The activity flow interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @return A dataframe containing activity flow interactions between proteins
#' without literature reference
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose one of those: 9606 human (default), 10116 rat or 10090 Mouse.
#' @param fields The user can define here the fields to be added. If used, set 
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments 
#'
#' @examples
#' interactions <-
#'     import_pathwayextra_interactions(
#'         resources = c('BioGRID', 'IntAct'),
#'         organism = 9606
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_PathwayExtra_Interactions
import_pathwayextra_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'pathwayextra',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


# Aliases (old names) to be deprecated
#'
#' @rdname import_pathwayextra_interactions
#' @param ... Passed to \code{import_pathwayextra_interactions}.
#' @export
#'
#' @noRd
import_PathwayExtra_Interactions <- function(...){
    .Deprecated("import_pathwayextra_interactions")
    import_pathwayextra_interactions(...)
}


#' Imports interactions from the `kinase extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=kinaseextra},
#' which contains enzyme-substrate interactions without literature reference.
#' The enzyme-substrate interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @return A dataframe containing enzyme-substrate interactions without
#' literature reference
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set 
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... Optional additional arguments.
#'
#' @examples
#' interactions <-
#'    import_kinaseextra_interactions(
#'        resources = c('PhosphoPoint', 'PhosphoSite'),
#'        organism = 9606
#'    )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_KinaseExtra_Interactions
import_kinaseextra_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL, 
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'kinaseextra',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_kinaseextra_interactions
#' @param ... Passed to \code{import_kinaseextra_interactions}.
#' @export
#'
#' @noRd
import_KinaseExtra_Interactions <- function(...){
    .Deprecated("import_kinaseextra_interactions")
    import_kinaseextra_interactions(...)
}


#' Imports interactions from the `ligrec extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=ligrecextra},
#' which contains ligand-receptor interactions without literature reference.
#' The ligand-receptor interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @return A dataframe containing ligand-receptor interactions including
#' the ones without literature references
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set 
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments 
#'
#' @examples
#' interactions <- import_ligrecextra_interactions(
#'     resources = c('HPRD', 'Guide2Pharma'),
#'     organism = 9606
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_LigrecExtra_Interactions
import_ligrecextra_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'ligrecextra',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_ligrecextra_interactions
#' @param ... Passed to \code{import_ligrecextra_interactions}.
#' @export
#'
#' @noRd
import_LigrecExtra_Interactions <- function(...){
    .Deprecated("import_ligrecextra_interactions")
    import_ligrecextra_interactions(...)
}


#' Imports all post-translational interactions from OmniPath
#'
#' Imports the dataset from all post-translational datasets of OmniPath.
#'
#' @return A dataframe containing post-translational interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param exclude datasets to exclude
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments
#'
#' @importFrom rlang %||% exec !!!
#' @importFrom magrittr %<>%
#' @examples
#' interactions <-
#'     import_post_translational_interactions(
#'         resources = c('BioGRID')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
import_post_translational_interactions <- function(
    resources = NULL,
    organism = 9606,
    exclude = NULL,
    references_by_resource = TRUE,
    ...
){

    args <- list(...)
    args$datasets <-
        args$datasets %||%
        c('omnipath', 'pathwayextra', 'kinaseextra', 'ligrecextra')
    args$datasets %<>% setdiff(exclude)

    args %<>% merge_lists(
        list(
            query_type = 'interactions',
            resources = resources,
            organism = organism,
            references_by_resource = references_by_resource
        )
    )


    exec(import_omnipath, !!!args)

}


#' From the OmniPath webservice imports interactions from the
#' DoRothEA dataset
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=dorothea}
#' which contains transcription factor (TF)-target interactions from DoRothEA
#' \url{https://github.com/saezlab/DoRothEA}
#'
#' @return A dataframe containing TF-target interactions from DoRothEA
#' @export
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param dorothea_levels Vector detailing the confidence levels of the
#' interactions to be downloaded. In dorothea, every TF-target interaction
#' has a confidence score ranging from A to E, being A the most reliable
#' interactions.
#' By default we take A and B level interactions (\code{c(A, B)}).
#' It is to note that E interactions are not available in OmnipathR.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments 
#'
#' @examples
#' interactions <- import_dorothea_interactions(
#'     resources = c('DoRothEA', 'ARACNe-GTEx_DoRothEA'),
#'     organism = 9606,
#'     dorothea_levels = c('A', 'B', 'C')
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_TFregulons_Interactions import_tfregulons_interactions
import_dorothea_interactions <- function(
    resources = NULL,
    organism = 9606,
    dorothea_levels = c('A', 'B'),
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        dorothea_levels = dorothea_levels,
        datasets = 'dorothea',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_dorothea_interactions
#' @param ... Passed to \code{import_dorothea_interactions}.
#' @export
#'
#' @noRd
import_TFregulons_Interactions <- function(...){
    .Deprecated("import_dorothea_interactions")
    import_dorothea_interactions(...)
}


#' @rdname import_dorothea_interactions
#' @param ... Passed to \code{import_dorothea_interactions}.
#' @export
#'
#' @noRd
import_tfregulons_interactions <- function(...){
    .Deprecated("import_dorothea_interactions")
    import_dorothea_interactions(...)
}


#' Imports interactions from the TF-target dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_target},
#' which contains transcription factor-target protein coding gene
#' interactions. Note: this is not the only TF-target dataset in OmniPath,
#' `dorothea` is the other one and the `tf_mirna` dataset provides
#' TF-miRNA gene interactions.
#'
#' @return A dataframe containing TF-target interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... Optional additional arguments
#'
#' @examples
#' interactions <-
#'     import_tf_target_interactions(
#'         resources = c('DoRothEA', 'SIGNOR')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
import_tf_target_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'tf_target',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


#' Imports all TF-target interactions from OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_target,dorothea},
#' which contains transcription factor-target protein coding gene
#' interactions.
#'
#' @return A dataframe containing TF-target interactions
#' @export
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param dorothea_levels Vector detailing the confidence levels of the
#' interactions to be downloaded. In dorothea, every TF-target interaction
#' has a confidence score ranging from A to E, being A the most reliable
#' interactions.
#' By default we take A and B level interactions (\code{c(A, B)}).
#' It is to note that E interactions are not available in OmnipathR.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments 
#'
#' @examples
#' interactions <-
#'     import_transcriptional_interactions(
#'         resources = c('PAZAR', 'ORegAnno', 'DoRothEA')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
import_transcriptional_interactions <- function(
    resources = NULL,
    organism = 9606,
    dorothea_levels = c('A', 'B'),
    references_by_resource = TRUE,
    ...
){

    is_directed <- NULL

    result <- rbind(
        import_dorothea_interactions(
            resources = resources,
            organism = organism,
            dorothea_levels = dorothea_levels,
            references_by_resource = references_by_resource,
            ...
        ),
        import_tf_target_interactions(
            resources = resources, 
            organism = organism, 
            references_by_resource = references_by_resource, 
            ...) %>% 
            mutate(dorothea_level = "") %>%
            select(source, target, source_genesymbol, target_genesymbol,
                is_directed, is_stimulation, is_inhibition, consensus_direction,
                consensus_stimulation, consensus_inhibition, dip_url, sources,
                references, curation_effort, dorothea_level,n_references, 
                n_resources)
    )

    return(result)

}


#' Imports interactions from the miRNA-target dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=mirnatarget},
#' which contains miRNA-mRNA interactions.
#'
#' @return A dataframe containing miRNA-mRNA interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set 
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments 
#'
#' @examples
#' interactions <-
#'     import_mirnatarget_interactions(
#'         resources = c('miRTarBase', 'miRecords')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_miRNAtarget_Interactions
import_mirnatarget_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'mirnatarget',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_mirnatarget_interactions
#' @param ... Passed to \code{import_mirnatarget_interactions}.
#' @export
#'
#' @noRd
import_miRNAtarget_Interactions <- function(...){
    .Deprecated("import_mirnatarget_interactions")
    import_mirnatarget_interactions(...)
}


#' Imports interactions from the TF-miRNA dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_mirna},
#' which contains transcription factor-miRNA gene interactions
#'
#' @return A dataframe containing TF-miRNA interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set 
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments 
#'
#' @examples
#' interactions <-
#'     import_tf_mirna_interactions(
#'         resources = c('TransmiR')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
import_tf_mirna_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL, 
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'tf_mirna',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


#' Imports interactions from the lncRNA-mRNA dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=lncrna_mrna},
#' which contains lncRNA-mRNA interactions
#'
#' @return A dataframe containing lncRNA-mRNA interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set 
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param ... optional additional arguments 
#'
#' @examples
#' interactions <-
#'     import_lncrna_mrna_interactions(
#'         resources = c('ncRDeathDB')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
import_lncrna_mrna_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'lncrna_mrna',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


#' Imports all interaction datasets available in OmniPath
#'
#' The interaction datasets currently available in OmniPath:
#'
#' omnipath: the OmniPath data as defined in the paper, an arbitrary optimum
#' between coverage and quality
#' pathwayextra: activity flow interactions without literature reference
#' kinaseextra: enzyme-substrate interactions without literature reference
#' ligrecextra: ligand-receptor interactions without literature reference
#' dorothea: transcription factor (TF)-target interactions from DoRothEA
#' tf_target: transcription factor (TF)-target interactions from other
#' resources
#' mirnatarget: miRNA-mRNA interactions
#' tf_mirna: TF-miRNA interactions
#' lncrna_mrna: lncRNA-mRNA interactions
#'
#' @return A dataframe containing all the datasets in the interactions query
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param dorothea_levels The confidence levels of the dorothea 
#' interactions (TF-target) which range from A to D. Set to A and B by
#' default.
#' @param exclude datasets to exclude
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE. 
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.

#' @param ... optional additional arguments 
#'
#' @examples
#' interactions <- import_all_interactions(
#'     resources = c('HPRD', 'BioGRID'),
#'     organism = 9606
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_AllInteractions
import_all_interactions <- function(
    resources = NULL,
    organism = 9606,
    dorothea_levels = c('A', 'B'),
    exclude = NULL,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    url <- paste0(options('omnipath.url'), 'queries/interactions?format=json')
    all_datasets <- jsonlite::fromJSON(txt = url)$datasets

    all_datasets <- setdiff(all_datasets, exclude)

    # it does not make sense without the type field
    fields <- unique(c(fields), 'type')

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        dorothea_levels = dorothea_levels,
        exclude = exclude,
        datasets = all_datasets,
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )
    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_all_interactions
#' @param ... Passed to \code{import_all_interactions}.
#' @export
import_AllInteractions <- function(...){
    .Deprecated("import_all_interactions")
    import_all_interactions(...)
}

#' Retrieve a list of interaction resources available in Omnipath
#'
#' Gets the names of the resources from
#' \url{https://omnipath.org/interactions}.
#'
#' @param dataset a dataset within the interactions query type. Currently
#' available datasets are `omnipath`, `kinaseextra`, `pathwayextra`,
#' `ligrecextra`, `dorothea`, `tf_target`, `tf_mirna`, `mirnatarget` and
#' `lncrna_mrna`
#'
#' @return character vector with the names of the interaction databases
#' @export
#'
#' @examples
#' get_interaction_resources()
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{import_omnipath_interactions}}}
#'     \item{\code{\link{import_pathwayextra_interactions}}}
#'     \item{\code{\link{import_kinaseextra_interactions}}}
#'     \item{\code{\link{import_ligrecextra_interactions}}}
#'     \item{\code{\link{import_mirnatarget_interactions}}}
#'     \item{\code{\link{import_dorothea_interactions}}}
#' }
#'
#' @aliases get_interaction_databases
get_interaction_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'interactions', datasets = dataset))

}

# Aliases (old names) to be deprecated
#' @rdname get_interaction_resources
#' @param ... Passed to \code{get_interaction_resources}.
#' @export
#'
#' @noRd
get_interaction_databases <- function(...){
    .Deprecated("get_interaction_resources")
    get_interaction_resources(...)
}


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
########## Complexes             ##########
########## ########## ########## ##########

#' Imports protein complexes from OmniPath
#'
#' Imports the complexes stored in Omnipath database from
#' \url{https://omnipathdb.org/complexes}.
#'
#' @return A dataframe containing information about complexes
#' @export
#'
#' @param resources complexes not reported in these databases are
#' removed. See \code{\link{get_complexes_databases}} for more information.
#' @param ... optional additional arguments 
#'
#' @examples
#' complexes = import_omnipath_complexes(
#'     resources = c('CORUM', 'hu.MAP')
#' )
#'
#' @seealso \itemize{\item{\code{\link{get_complexes_databases}}}}
#'
#' @aliases import_Omnipath_complexes import_OmniPath_complexes
import_omnipath_complexes <- function(
    resources = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'complexes',
        resources = resources,
        ...
    )

    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_omnipath_complexes
#' @param ... Passed to \code{import_omnipath_complexes}.
#' @export
#'
#' @noRd
import_Omnipath_complexes <- function(...){
    .Deprecated("import_omnipath_complexes")
    import_omnipath_complexes(...)
}


#' @rdname import_omnipath_complexes
#' @param ... Passed to \code{import_omnipath_complexes}.
#' @export
#'
#' @noRd
import_OmniPath_complexes <- function(...){
    .Deprecated("import_omnipath_complexes")
    import_omnipath_complexes(...)    
}


#' Retrieve a list of complex resources available in Omnipath
#'
#' Get the names of the resources from \url{https://omnipath.org/complexes}
#'
#' @param dataset ignored for this query type
#' @return character vector with the names of the databases
#'
#' @examples
#' get_complex_resources()
#'
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_omnipath_complexes}}}
#' }
#'
#' @aliases get_complexes_databases
get_complex_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'complexes', datasets = dataset))

}


# Aliases (old names) to be deprecated
#' @rdname get_complex_resources
#' @param ... Passed to \code{import_omnipath_enzsub}.
#' @export
#'
#' @noRd
get_complexes_databases <- function(...){
    .Deprecated("get_complex_resources")
    get_complex_resources(...)
}

########## ########## ########## ##########
########## Annotations           ##########
########## ########## ########## ##########

#' Imports annotations from OmniPath
#'
#' Imports protein annotations about function, localization, expression,
#' structure and other properties of proteins from OmniPath
#' \url{https://omnipathdb.org/annotations}.
#' Note: there might be also a few miRNAs annotated; a vast majority of
#' protein complex annotations are inferred from the annotations of the
#' members: if all members carry the same annotation the complex inherits.
#'
#' @details Downloading the full annotations
#' dataset is disabled by default because the size of this data is
#' around 1GB. We recommend to retrieve the annotations for a set of proteins
#' or only from a few resources, depending on your interest. You can always
#' download the full database from \url{
#' https://archive.omnipathdb.org/omnipath_webservice_annotations__recent.tsv}
#' using any standard R or \code{readr} method.
#'
#' @return A data frame containing different gene and complex annotations.
#' @export
#'
#' @param proteins Vector containing the genes or proteins for whom
#'     annotations will be retrieved (UniProt IDs or HGNC Gene Symbols or
#'     miRBase IDs). It is also possible to donwload annotations for protein
#'     complexes. To do so, write 'COMPLEX:' right before the genesymbols of
#'     the genes integrating the complex. Check the vignette for examples.
#' @param resources Load the annotations only from these databases.
#'     See \code{\link{get_annotation_resources}} for possible values.
#' @param wide Convert the annotation table to wide format, which
#'     corresponds more or less to the original resource. If the data comes
#'     from more than one resource a list of wide tables will be returned.
#'     See examples at \code{\link{pivot_annotations}}.
#' @param ... Additional arguments.
#'
#' @examples
#' annotations <- import_omnipath_annotations(
#'     proteins = c('TP53', 'LMNA'),
#'     resources = c('HPA_subcellular')
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_annotation_databases}}}
#'     \item{\code{\link{pivot_annotations}}}
#' }
#'
#' @aliases import_Omnipath_annotations import_OmniPath_annotations
import_omnipath_annotations <- function(
    proteins = NULL,
    resources = NULL,
    wide = FALSE,
    ...
){

    # account for the old argument name
    proteins <- c(proteins, list(...)$select_genes)

    if(
        is.null(proteins) &&
        is.null(resources)
    ){

        stop(
            paste(
                'Downloading the entire annotations database is not allowed',
                'by default because of its huge size (>1GB). If you really',
                'want to do that, you find static files at',
                'https://archive.omnipathdb.org/.',
                'However we recommend to query a set of proteins or a few',
                'resources, depending on your interest.'
            )
        )

    }

    if(length(proteins) < 600){

        result <- import_omnipath(
            query_type = 'annotations',
            proteins = proteins,
            resources = resources,
            ...
        )

        if(!is.null(proteins)){
            result <- result[
                which(
                    result$uniprot %in% proteins |
                    result$genesymbol %in% proteins
                ),
            ]
        }

    }else{

        parts <- list()

        proteins_chunks <-
            split(
                proteins,
                cut(
                    seq_along(proteins),
                    breaks = length(proteins) / 500,
                    labels = FALSE
                )
            )

        cat('Downloading ')

        for(proteins_chunk in proteins_chunks){

            cat('.')

            parts <- append(
                parts,
                list(
                    import_omnipath(
                        query_type = 'annotations',
                        proteins = proteins_chunk,
                        resources = resources,
                        recursive_call = TRUE,
                        silent = TRUE,
                        ...
                    )
                )
            )

        }

        cat(' ready.\n')

        result <- do.call(rbind, parts)

        logger::log_success(
            'Downloaded %d annotation records.',
            nrow(result)
        )

    }

    if(wide){

        result <- pivot_annotations(result)

    }

    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_omnipath_annotations
#' @param ... Passed to \code{import_omnipath_annotations}.
#' @export
#'
#' @noRd
import_Omnipath_annotations <- function(...){
    .Deprecated("import_omnipath_annotations")
    import_omnipath_annotations(...)
}


#' @rdname import_omnipath_annotations
#' @param ... Passed to \code{import_omnipath_annotations}.
#' @export
#'
#' @noRd
import_OmniPath_annotations <- function(...){
    .Deprecated("import_omnipath_annotations")
    import_omnipath_annotations(...)
}


#' Retrieves a list of available resources in the annotations database
#' of OmniPath
#'
#' Get the names of the resources from \url{https://omnipath.org/annotations}.
#'
#' @return character vector with the names of the annotation resources
#' @export
#' @param dataset ignored for this query type
#' @param ... optional additional arguments 
#'
#' @examples
#' get_annotation_resources()
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_omnipath_annotations}}}
#' }
#'
#' @aliases get_annotation_databases
get_annotation_resources <- function(dataset = NULL, ...){

    return(get_resources(query_type = 'annotations', datasets = dataset))

}


# Aliases (old names) to be deprecated
#' @rdname get_annotation_resources
#' @param ... Passed to \code{get_annotation_resources}.
#' @export
#'
#' @noRd
get_annotation_databases <- function(...){
    .Deprecated("get_annotation_resources")
    get_annotation_resources(...)
}


#' Annotation categories and resources
#'
#' A full list of annotation resources, keys and values.
#'
#' @return A data frame with resource names, annotation key labels and
#'     for each key all possible values.
#'
#' @examples
#' annot_cat <- annotation_categories()
#' annot_cat
#' # # A tibble: 46,307 x 3
#' #    source           label    value
#' #    <chr>            <chr>    <chr>
#' #  1 connectomeDB2020 role     ligand
#' #  2 connectomeDB2020 role     receptor
#' #  3 connectomeDB2020 location ECM
#' #  4 connectomeDB2020 location plasma membrane
#' #  5 connectomeDB2020 location secreted
#' #  6 KEGG-PC          pathway  Alanine, aspartate and glutamate metabolism
#' #  7 KEGG-PC          pathway  Amino sugar and nucleotide sugar metabolism
#' #  8 KEGG-PC          pathway  Aminoacyl-tRNA biosynthesis
#' #  9 KEGG-PC          pathway  Arachidonic acid metabolism
#' # 10 KEGG-PC          pathway  Arginine and proline metabolism
# . with 46,297 more rows
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr separate_rows
#' @export
annotation_categories <- function(){

    # NSE vs. R CMD check workaround
    value <- NULL

    import_omnipath('annotations_summary', license = NA) %>%
    separate_rows(value, sep = '#')

}


#' Converts annotation tables to a wide format
#'
#' Use this method to reconstitute the annotation tables into the format of
#' the original resources. With the `wide=TRUE` option
#' \code{\link{import_omnipath_annotations}} applies this function to the
#' downloaded data.
#'
#' @param annotations A data frame of annotations downloaded from the
#'    OmniPath web service by \code{\link{import_omnipath_annotations}}.
#'
#' @return A wide format data frame (tibble) if the provided data contains
#' annotations from one resource, otherwise a list of wide format tibbles.
#'
#' @examples
#' # single resource: the result is a data frame
#' disgenet <- import_omnipath_annotations(resources = 'DisGeNet')
#' disgenet <- pivot_annotations(disgenet)
#' disgenet
#' # # A tibble: 119,551 x 10
#' #    uniprot genesymbol entity_type disease score dsi   dpi   nof_pmids
#' #    <chr>   <chr>      <chr>       <chr>   <chr> <chr> <chr> <chr>
#' #  1 P04217  A1BG       protein     Schizo… 0.3   0.857 0.172 1
#' #  2 P04217  A1BG       protein     Hepato… 0.3   0.857 0.172 1
#' #  3 P01023  A2M        protein     alpha-… 0.31  0.564 0.724 0
#' #  4 P01023  A2M        protein     Fibros… 0.3   0.564 0.724 1
#' #  5 P01023  A2M        protein     Hepato… 0.3   0.564 0.724 1
#' # # . with 119,541 more rows, and 2 more variables: nof_snps <chr>,
#' # #   source <chr>
#'
#' # multiple resources: the result is a list
#' annotations <- import_omnipath_annotations(
#'     resources = c('DisGeNet', 'SignaLink_function', 'DGIdb', 'kinase.com')
#' )
#' annotations <- pivot_annotations(annotations)
#' names(annotations)
#' # [1] "DGIdb"              "DisGeNet"           "kinase.com"
#' # [4] "SignaLink_function"
#' annotations$kinase.com
#' # # A tibble: 825 x 6
#' #    uniprot genesymbol entity_type group family subfamily
#' #    <chr>   <chr>      <chr>       <chr> <chr>  <chr>
#' #  1 P31749  AKT1       protein     AGC   Akt    NA
#' #  2 P31751  AKT2       protein     AGC   Akt    NA
#' #  3 Q9Y243  AKT3       protein     AGC   Akt    NA
#' #  4 O14578  CIT        protein     AGC   DMPK   CRIK
#' #  5 Q09013  DMPK       protein     AGC   DMPK   GEK
#' # # . with 815 more rows
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select pull group_split
#' @importFrom purrr map
#' @importFrom rlang set_names
#'
#' @seealso \code{\link{import_omnipath_annotations}}
pivot_annotations <- function(annotations){

    # NSE vs. R CMD check workaround
    record_id <- NULL

    if(annotations %>% pull(source) %>% unique %>% length %>% `>`(1)){

        annotations %>%
        group_split(source) %>%
        set_names(annotations %>% pull(source) %>% unique %>% sort) %>%
        map(pivot_annotations)

    }else{

        (
            annotations %>%
            tidyr::pivot_wider(
                id_cols = c(
                    'record_id',
                    'uniprot',
                    'genesymbol',
                    'entity_type'
                ),
                names_from = 'label',
                values_from = 'value'
            ) %>%
            select(-record_id)
        )

    }

}

########## ########## ########## ##########
########## Intercell             ##########
########## ########## ########## ##########

#' Imports OmniPath intercell annotations
#'
#' Imports the OmniPath intercellular communication role annotation database
#' from \url{https://omnipathdb.org/intercell}. It provides information
#' on the roles in inter-cellular signaling. E.g. if a protein is
#' a ligand, a receptor, an extracellular matrix (ECM) component, etc.
#'
#' @return A dataframe cotaining information about roles in intercellular
#' signaling.
#'
#' @param categories vector containing the categories to be retrieved.
#'     All the genes belonging to those categories will be returned. For
#'     further information about the categories see
#'     code{\link{get_intercell_categories}}.
#' @param parent vector containing the parent classes to be retrieved.
#'     All the genes belonging to those classes will be returned. For
#'     furter information about the main classes see
#'     \code{\link{get_intercell_categories}}.
#' @param resources limit the query to certain resources; see the available
#'     resources by \code{\link{get_intercell_resources}}.
#' @param scope either `specific` or `generic`
#' @param aspect either `locational` or `functional`
#' @param source either `resource_specific` or `composite`
#' @param transmitter logical, include only transmitters i.e. proteins
#'     delivering signal from a cell to its environment.
#' @param receiver logical, include only receivers i.e. proteins delivering
#'     signal to the cell from its environment.
#' @param plasma_membrane_peripheral logical, include only plasma membrane
#'     peripheral membrane proteins.
#' @param plasma_membrane_transmembrane logical, include only plasma membrane
#'     transmembrane proteins.
#' @param secreted logical, include only secreted proteins
#' @param proteins limit the query to certain proteins
#' @param topology topology categories: one or more of `secreted` (sec),
#'     `plasma_membrane_peripheral` (pmp), `plasma_membrane_transmembrane`
#'     (pmtm) (both short or long notation can be used).
#' @param causality `transmitter` (trans), `receiver` (rec) or `both` (both
#'     short or long notation can be used).
#' @param consensus_percentile Numeric: a percentile cut off for the
#'     consensus score of generic categories. The consensus score is the
#'     number of resources supporting the classification of an entity into a
#'     category based on combined information of many resources. Here you can
#'     apply a cut-off, keeping only the annotations supported by a higher
#'     number of resources than a certain percentile of each category. If
#'     \code{NULL} no filtering will be performed. The value is either in the
#'     0-1 range, or will be divided by 100 if greater than 1. The
#'     percentiles will be calculated against the generic composite
#'     categories and then will be applied to their resource specific
#'     annotations and specific child categories.
#' @param loc_consensus_percentile Numeric: similar to
#'     code{consensus_percentile} for major localizations. For example, with
#'     a value of 50, the secreted, plasma membrane transmembrane or
#'     peripheral attributes will be true only where at least 50 percent
#'     of the resources support these.
#' @param ... Additional optional arguments, ignored.
#'
#' @examples
#' intercell <- import_omnipath_intercell(categories = 'ecm')
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_intercell_categories}}}
#'     \item{\code{\link{get_intercell_generic_categories}}}
#'     \item{\code{\link{import_intercell_network}}}
#'     \item{\code{\link{intercell_consensus_filter}}}
#' }
#'
#' @aliases import_Omnipath_intercell import_OmniPath_intercell
import_omnipath_intercell <- function(
    categories = NULL,
    resources = NULL,
    parent = NULL,
    scope = NULL,
    aspect = NULL,
    source = NULL,
    transmitter = NULL,
    receiver = NULL,
    secreted = NULL,
    plasma_membrane_peripheral = NULL,
    plasma_membrane_transmembrane = NULL,
    proteins = NULL,
    topology = NULL,
    causality = NULL,
    consensus_percentile = NULL,
    loc_consensus_percentile = NULL,
    ...
){

    args <- c(as.list(environment()), list(...))
    args$query_type <- 'intercell'
    args$logicals <- c(
        'transmitter',
        'receiver',
        'secreted',
        'plasma_membrane_peripheral',
        'plasma_membrane_transmembrane'
    )
    args$consensus_percentile <- NULL
    args$loc_consensus_percentile <- NULL

    result <-
        do.call(import_omnipath, args) %>%
        intercell_consensus_filter(
            consensus_percentile,
            loc_consensus_percentile
        )

    return(result)

}


#' Quality filter for intercell annotations
#'
#' @param data A data frame with intercell annotations, as provided by
#'     \code{\link{import_omnipath_intercell}}.
#' @param percentile Numeric: a percentile cut off for the consensus score
#'     of composite categories. The consensus score is the number of
#'     resources supporting the classification of an entity into a category
#'     based on combined information of many resources. Here you can apply
#'     a cut-off, keeping only the annotations supported by a higher number
#'     of resources than a certain percentile of each category. If
#'     \code{NULL} no filtering will be performed. The value is either in the
#'     0-1 range, or will be divided by 100 if greater than 1. The
#'     percentiles will be calculated against the generic composite
#'     categories and then will be applied to their resource specific
#'     annotations and specific child categories.
#' @param loc_percentile Numeric: similar to \code{percentile} for major
#'     localizations. For example, with a value of 50, the secreted, plasma
#'     membrane transmembrane or peripheral attributes will be \code{TRUE}
#'     only where at least 50 percent of the resources support these.
#'
#' @return The data frame in \code{data} filtered by the consensus scores.
#'
#' @examples
#' intercell <- import_omnipath_intercell(parent = c('ligand', 'receptor'))
#' nrow(intercell)
#' # [1] 50174
#' intercell_q50 <- intercell_consensus_filter(intercell, 50)
#' nrow(intercell_q50)
#' # [1] 42863
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr group_by filter ungroup bind_rows
#' @importFrom dplyr select distinct inner_join pull
#' @importFrom stats quantile
#' @importFrom rlang !! := sym
#' @importFrom purrr reduce
#' @export
intercell_consensus_filter <- function(
    data,
    percentile = NULL,
    loc_percentile = NULL
){

    # NSE vs. R CMD check workaround
    scope <- source <- parent <- consensus_score <- NULL

    percentile %<>%
        if_null(0L) %>%
        {`if`(. > 1, . / 100, .)}

    thresholds <-
        data %>%
        filter(scope == 'generic' & source == 'composite') %>%
        group_by(parent) %>%
        filter(
            consensus_score >= quantile(consensus_score, percentile)
        ) %>%
        ungroup %>%
        select(parent, uniprot) %>%
        distinct

    composite_parents <-
        data %>%
        filter(source == 'composite') %>%
        pull(parent) %>%
        unique

    data %<>%
        inner_join(thresholds, by = c('parent', 'uniprot')) %>%
        bind_rows(
            data %>%
            filter(!parent %in% composite_parents)
        )

    if(!is.null(loc_percentile)){

        major_locations <- c(
            'secreted',
            'plasma_membrane_transmembrane',
            'plasma_membrane_peripheral'
        )

        locations <- import_omnipath_intercell(
            aspect = 'locational',
            parent = major_locations,
            consensus_percentile = loc_percentile
        )

        data %<>%
        {reduce(
            major_locations,
            function(data, loc){

                in_location <-
                    locations %>%
                    filter(!!sym(loc)) %>%
                    pull(uniprot) %>%
                    unique

                data %>%
                mutate(!!sym(loc) := uniprot %in% in_location)

            },
            .init = .
        )}

    }

    return(data)

}


# Aliases (old names) to be deprecated
#' @rdname import_omnipath_intercell
#' @param ... Passed to \code{import_omnipath_intercell}.
#' @export
#'
#' @noRd
import_Omnipath_intercell <- function(...){
    .Deprecated("import_omnipath_intercell")
    import_omnipath_intercell(...)
}


#' @rdname import_omnipath_intercell
#' @param ... Passed to \code{import_omnipath_intercell}.
#' @export
#'
#' @noRd
import_OmniPath_intercell <- function(...){
    .Deprecated("import_omnipath_intercell")
    import_omnipath_intercell(...)
}


#' Retrieves a list of intercellular communication resources available in
#' OmniPath
#'
#' Retrieves a list of the databases from
#' \url{https://omnipath.org/intercell}.
#'
#' @param dataset ignored at this query type
#'
#' @return character vector with the names of the databases
#'
#' @examples
#' get_intercell_resources()
#'
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_omnipath_intercell}}}
#' }
get_intercell_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'intercell', datasets = dataset))

}


#' Intercellular communication network
#'
#' Imports an intercellular network by combining intercellular annotations
#' and protein interactions. First imports a network of protein-protein
#' interactions. Then, it retrieves annotations about the proteins
#' intercellular communication roles, once for the transmitter (delivering
#' information from the expressing cell) and second, the receiver (receiving
#' signal and relaying it towards the expressing cell) side. These 3 queries
#' can be customized by providing parameters in lists which will be passed to
#' the respective methods (\code{\link{import_omnipath_interactions}} for
#' the network and \code{\link{import_omnipath_intercell}} for the
#' annotations). Finally the 3 data frames combined in a way that the source
#' proteins in each interaction annotated by the transmitter, and the target
#' proteins by the receiver categories. If undirected interactions present
#' (these are disabled by default) they will be duplicated, i.e. both
#' partners can be both receiver and transmitter.
#'
#' @return A dataframe containing information about protein-protein
#' interactions and the inter-cellular roles of the protiens involved in those
#' interactions.
#'
#' @param interactions_param a list with arguments for an interactions query:
#'     \code{\link{import_omnipath_interactions}},
#'     \code{\link{import_pathwayextra_interactions}},
#'     \code{\link{import_kinaseextra_interactions}},
#'     \code{\link{import_ligrecextra_interactions}}
#' @param transmitter_param a list with arguments for
#'     \code{\link{import_omnipath_intercell}}, to define the transmitter side
#'     of intercellular connections
#' @param receiver_param a list with arguments for
#'     \code{\link{import_omnipath_intercell}}, to define the receiver side
#'     of intercellular connections
#' @param resources A character vector of resources to be applied to
#'     both the interactions and the annotations. For example, \code{resources
#'     = 'CellChatDB'} will download the transmitters and receivers defined by
#'     CellChatDB, connected by connections from CellChatDB.
#' @param entity_types Character, possible values are "protein", "complex" or
#'     both.
#' @param ligand_receptor Logical. If \code{TRUE}, only \emph{ligand} and
#'     \emph{receptor} annotations will be used instead of the more generic
#'     \emph{transmitter} and \emph{receiver} categories.
#' @param high_confidence Logical: shortcut to do some filtering in order to
#'     include only higher confidence interactions. The intercell database
#'     of OmniPath covers a very broad range of possible ways of cell to cell
#'     communication, and the pieces of information, such as localization,
#'     topology, function and interaction, are combined from many, often
#'     independent sources. This unavoidably result some weird and unexpected
#'     combinations which are false positives in the context of intercellular
#'     communication. This option sets some minimum criteria to remove most
#'     (but definitely not all!) of the wrong connections. These criteria
#'     are the followings: 1) the receiver must be plasma membrane
#'     transmembrane; 2) the curation effort for interactions must be larger
#'     than one; 3) the consensus score for annotations must be larger than
#'     the 50 percentile within the generic category (you can override this
#'     by \code{consensus_percentile}). 4) the transmitter must be secreted
#'     or exposed on the plasma membrane. 5) The major localizations have
#'     to be supported by at least 30 percent of the relevant resources (
#'     you can override this by \code{loc_consensus_percentile}). 6) The
#'     datasets with lower level of curation (\emph{kinaseextra} and \emph{
#'     pathwayextra}) will be disabled. These criteria are of medium
#'     stringency, you can always tune them to be more relaxed or stringent
#'     by filtering manually, using \code{\link{filter_intercell_network}}.
#' @param simplify Logical: keep only the most often used columns. This
#'     function combines a network data frame with two copies of the
#'     intercell annotation data frames, all of them already having quite
#'     some columns. With this option we keep only the names of the
#'     interacting pair, their intercellular communication roles, and the
#'     minimal information of the origin of both the interaction and
#'     the annotations.
#' @param unique_pairs Logical: instead of having separate rows for each
#'     pair of annotations, drop the annotations and reduce the data frame to
#'     unique interacting pairs. See \code{\link{unique_intercell_network}}
#'     for details.
#' @param consensus_percentile Numeric: a percentile cut off for the consensus
#'     score of generic categories in intercell annotations. The consensus
#'     score is the number of resources supporting the classification of an
#'     entity into a category based on combined information of many resources.
#'     Here you can apply a cut-off, keeping only the annotations supported
#'     by a higher number of resources than a certain percentile of each
#'     category. If \code{NULL} no filtering will be performed. The value is
#'     either in the 0-1 range, or will be divided by 100 if greater than 1.
#'     The percentiles will be calculated against the generic composite
#'     categories and then will be applied to their resource specific
#'     annotations and specific child categories.
#' @param loc_consensus_percentile Numeric: similar to
#'     \code{consensus_percentile} for major localizations. For example, with
#'     a value of 50, the secreted, plasma membrane transmembrane or
#'     peripheral attributes will be \code{TRUE} only where at least 50
#'     percent of the resources support these.
#' @param omnipath Logical: shortcut to include the \emph{omnipath} dataset
#'     in the interactions query.
#' @param ligrecextra Logical: shortcut to include the \emph{ligrecextra}
#'     dataset in the interactions query.
#' @param kinaseextra Logical: shortcut to include the \emph{kinaseextra}
#'     dataset in the interactions query.
#' @param pathwayextra Logical: shortcut to include the \emph{pathwayextra}
#'     dataset in the interactions query.
#' @param ... If \code{simplify} or \code{unique_pairs} is \code{TRUE},
#'     additional column  names can be passed here to \code{dplyr::select}
#'     on the final data frame. Otherwise ignored.
#'
#' @details
#' By default this function creates almost the largest possible network of
#' intercellular interactions. However, this might contain a large number
#' of false positives. Please refer to the documentation of the arguments,
#' especially \code{high_confidence}, and the \code{
#' \link{filter_intercell_network}} function. Note: if you restrict the query
#' to certain intercell annotation resources or small categories, it's not
#' recommended to use the \code{consensus_percentile} or
#' \code{high_confidence} options, instead filter the network with \code{
#' \link{filter_intercell_network}} for more consistent results.
#'
#' @examples
#' intercell_network <- import_intercell_network(
#'     interactions_param = list(datasets = 'ligrecextra'),
#'     receiver_param = list(categories = c('receptor', 'transporter')),
#'     transmitter_param = list(categories = c('ligand', 'secreted_enzyme'))
#' )
#'
#' @importFrom dplyr rename bind_rows filter inner_join distinct group_by
#' @importFrom dplyr summarize_all first
#' @importFrom rlang %||%
#' @importFrom magrittr %>% %<>%
#' @importFrom utils modifyList
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_intercell_categories}}}
#'     \item{\code{\link{get_intercell_generic_categories}}}
#'     \item{\code{\link{import_omnipath_intercell}}}
#'     \item{\code{\link{import_omnipath_interactions}}}
#'     \item{\code{\link{import_pathwayextra_interactions}}}
#'     \item{\code{\link{import_kinaseextra_interactions}}}
#'     \item{\code{\link{import_ligrecextra_interactions}}}
#'     \item{\code{\link{unique_intercell_network}}}
#'     \item{\code{\link{simplify_intercell_network}}}
#'     \item{\code{\link{filter_intercell_network}}}
#' }
import_intercell_network <- function(
    interactions_param = list(),
    transmitter_param = list(),
    receiver_param = list(),
    resources = NULL,
    entity_types = NULL,
    ligand_receptor = FALSE,
    high_confidence = FALSE,
    simplify = FALSE,
    unique_pairs = FALSE,
    consensus_percentile = NULL,
    loc_consensus_percentile = NULL,
    omnipath = TRUE,
    ligrecextra = TRUE,
    kinaseextra = !high_confidence,
    pathwayextra = !high_confidence,
    ...
){

    # NSE vs. R CMD check workaround
    parent <- secreted <- plasma_membrane_transmembrane <-
    plasma_membrane_peripheral <- NULL

    datasets <-
        environment() %>%
        select_interaction_datasets

    # retrieving interactions
    interactions_param <- list(
            query_type = 'interactions',
            datasets = datasets,
            fields = 'datasets'
        ) %>%
        insert_if_not_null(
            resources = resources,
            entity_types = entity_types
        ) %>%
        modifyList(interactions_param)

    interactions <- do.call(
        import_omnipath,
        interactions_param
    )
    interactions <- swap_undirected(interactions)

    # retrieving intercell annotations

    consensus_percentile %<>%
        {`if`(high_confidence, . %||% 50, .)}

    loc_consensus_percentile %<>%
        {`if`(high_confidence, . %||% 30, .)}

    transmitter_param <- list(
            causality = 'trans',
            scope = 'generic'
        ) %>%
        insert_if_not_null(
            resources = resources,
            entity_types = entity_types,
            consensus_percentile = consensus_percentile,
            loc_consensus_percentile = loc_consensus_percentile
        ) %>%
        {`if`(
            ligand_receptor,
            `[[<-`(., 'parent', 'ligand'),
            .
        )} %>%
        modifyList(transmitter_param)

    receiver_param <- list(
            causality = 'rec',
            scope = 'generic'
        ) %>%
        insert_if_not_null(
            resources = resources,
            entity_types = entity_types,
            consensus_percentile = consensus_percentile
        ) %>%
        {`if`(
            ligand_receptor,
            `[[<-`(., 'parent', 'receptor'),
            .
        )} %>%
        modifyList(receiver_param)

    intracell <- c('intracellular_intercellular_related', 'intracellular')
    transmitters <-
        do.call(import_omnipath_intercell, transmitter_param) %>%
        filter(!parent %in% intracell) %>%
        rename(category_source = source) %>%
        {`if`(
            high_confidence,
            filter(
                .,
                secreted |
                plasma_membrane_transmembrane |
                plasma_membrane_peripheral
            ),
            .
        )}
    receivers <-
        do.call(import_omnipath_intercell, receiver_param) %>%
        filter(!parent %in% intracell) %>%
        rename(category_source = source) %>%
        {`if`(
            high_confidence,
            filter(., plasma_membrane_transmembrane),
            .
        )}

    interactions %>%
    {`if`(
        high_confidence,
        filter(., curation_effort > 1),
        .
    )} %>%
    inner_join(
        transmitters,
        by = c('source' = 'uniprot')
    ) %>%
    group_by(
        category, parent, source, target
    ) %>%
    mutate(
        database = paste(database, collapse = ';')
    ) %>%
    summarize_all(first) %>%
    inner_join(
        receivers,
        by = c('target' = 'uniprot'),
        suffix = c('_intercell_source', '_intercell_target')
    ) %>%
    group_by(
        category_intercell_source,
        parent_intercell_source,
        source,
        target,
        category_intercell_target,
        parent_intercell_target
    ) %>%
    mutate(
        database_intercell_target = paste(
            database_intercell_target,
            collapse = ';'
        )
    ) %>%
    summarize_all(first) %>%
    ungroup() %>%
    {`if`(
        simplify,
        simplify_intercell_network(., ...),
        .
    )} %>%
    {`if`(
        unique_pairs,
        unique_intercell_network(., ...),
        .
    )}

}


#' Create a vector with dataset names from an environment with logical
#' variables.
#'
#' @param envir Environment from the calling function where dataset names
#'     present as logical variables.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr keep
#'
#' @noRd
select_interaction_datasets <- function(envir){

    envir %>%
    as.list %>%
    `[`(
        c(
            'omnipath',
            'pathwayextra',
            'kinaseextra',
            'ligrecextra'
        )
    ) %>%
    keep(identity) %>%
    names

}


#' Quality filter an intercell network
#'
#' The intercell database  of OmniPath covers a very broad range of possible
#' ways of cell to cell communication, and the pieces of information, such as
#' localization, topology, function and interaction, are combined from many,
#' often independent sources. This unavoidably result some weird and
#' unexpected combinations which are false positives in the context of
#' intercellular communication. \code{\link{import_intercell_network}}
#' provides a shortcut (\code{high_confidence}) to do basic quality filtering.
#' For custom filtering or experimentation with the parameters we offer this
#' function.
#'
#' @param network An intercell network data frame, as provided by
#'     \code{\link{import_intercell_network}}, without \code{simplify}.
#' @param transmitter_topology Character vector: topologies allowed for the
#'     entities in transmitter role. Abbreviations allowed: "sec", "pmtm"
#'     and "pmp".
#' @param receiver_topology Same as \code{transmitter_topology} for the
#'     entities in the receiver role.
#' @param min_curation_effort Numeric: a minimum value of curation effort
#'     (resource-reference pairs) for network interactions. Use zero to
#'     disable filtering.
#' @param min_resources Numeric: minimum number of resources for
#'     interactions. The value 1 means no filtering.
#' @param min_references Numeric: minimum number of references for
#'     interactions. Use zero to disable filtering.
#' @param min_provenances Numeric: minimum number of provenances (either
#'     resources or references) for interactions. Use zero or one to
#'     disable filtering.
#' @param consensus_percentile Numeric: percentile threshold for the consensus
#'     score of generic categories in intercell annotations. The consensus
#'     score is the number of resources supporting the classification of an
#'     entity into a category based on combined information of many resources.
#'     Here you can apply a cut-off, keeping only the annotations supported
#'     by a higher number of resources than a certain percentile of each
#'     category. If \code{NULL} no filtering will be performed. The value is
#'     either in the 0-1 range, or will be divided by 100 if greater than 1.
#'     The percentiles will be calculated against the generic composite
#'     categories and then will be applied to their resource specific
#'     annotations and specific child categories.
#' @param loc_consensus_percentile Numeric: similar to
#'     \code{consensus_percentile} for major localizations. For example, with
#'     a value of 50, the secreted, plasma membrane transmembrane or
#'     peripheral attributes will be \code{TRUE} only where at least 50
#'     percent of the resources support these.
#' @param ligand_receptor Logical. If \code{TRUE}, only \emph{ligand} and
#'     \emph{receptor} annotations will be used instead of the more generic
#'     \emph{transmitter} and \emph{receiver} categories.
#' @param simplify Logical: keep only the most often used columns. This
#'     function combines a network data frame with two copies of the
#'     intercell annotation data frames, all of them already having quite
#'     some columns. With this option we keep only the names of the
#'     interacting pair, their intercellular communication roles, and the
#'     minimal information of the origin of both the interaction and
#'     the annotations.
#' @param unique_pairs Logical: instead of having separate rows for each
#'     pair of annotations, drop the annotations and reduce the data frame to
#'     unique interacting pairs. See \code{\link{unique_intercell_network}}
#'     for details.
#' @param omnipath Logical: shortcut to include the \emph{omnipath} dataset
#'     in the interactions query.
#' @param ligrecextra Logical: shortcut to include the \emph{ligrecextra}
#'     dataset in the interactions query.
#' @param kinaseextra Logical: shortcut to include the \emph{kinaseextra}
#'     dataset in the interactions query.
#' @param pathwayextra Logical: shortcut to include the \emph{pathwayextra}
#'     dataset in the interactions query.
#' @param ... If \code{simplify} or \code{unique_pairs} is \code{TRUE},
#'     additional column  names can be passed here to \code{dplyr::select}
#'     on the final data frame. Otherwise ignored.
#'
#' @return An intercell network data frame filtered.
#'
#' @examples
#' icn <- import_intercell_network()
#' icn_f <- filter_intercell_network(
#'     icn,
#'     consensus_percentile = 75,
#'     min_provenances = 3,
#'     simplify = TRUE
#' )
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct filter inner_join left_join
#' @importFrom rlang !!! parse_expr exprs syms
#' @importFrom logger log_warn
#' @importFrom purrr walk
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_intercell_network}}}
#'     \item{\code{\link{unique_intercell_network}}}
#'     \item{\code{\link{simplify_intercell_network}}}
#' }
filter_intercell_network <- function(
    network,
    transmitter_topology = c(
        'secreted',
        'plasma_membrane_transmembrane',
        'plasma_membrane_peripheral'
    ),
    receiver_topology = 'plasma_membrane_transmembrane',
    min_curation_effort = 2,
    min_resources = 1,
    min_references = 0,
    min_provenances = 1,
    consensus_percentile = 50,
    loc_consensus_percentile = 30,
    ligand_receptor = FALSE,
    simplify = FALSE,
    unique_pairs = FALSE,
    omnipath = TRUE,
    ligrecextra = TRUE,
    kinaseextra = FALSE,
    pathwayextra = FALSE,
    ...
){

    # NSE vs. R CMD check workaround
    parent <- curation_effort <- n_resources <- n_references <- NULL

    major_locations <- c(
        'secreted',
        'plasma_membrane_transmembrane',
        'plasma_membrane_peripheral'
    )

    consensus <-
        import_omnipath_intercell(
            consensus_percentile = consensus_percentile,
            loc_consensus_percentile = loc_consensus_percentile
        )

    consensus_annot <-
        consensus %>%
        select(uniprot, parent) %>%
        distinct

    consensus_loc <-
        consensus %>%
        select(uniprot, !!!syms(major_locations)) %>%
        distinct

    topology_short <-
        major_locations %>%
        set_names(c('sec', 'pmtm', 'pmp'))
    topologies <- unlist(topology_short)
    check_topo <- function(x){
        if(!x %in% topologies){
            log_warn('Unknown topology: %s', x)
        }
    }

    datasets <-
        environment() %>%
        select_interaction_datasets

    missing_datasets <- datasets %>% setdiff(colnames(network))

    if(length(missing_datasets)){

        msg <- sprintf(
            paste0(
                'filter_intercell_network: cannot select %s %s, ',
                'apparently %s %s not included in the original ',
                'download.'
            ),
            missing_datasets %>%
            plural('dataset'),
            missing_datasets %>%
            pretty_list,
            missing_datasets %>%
            plural('this', 'these'),
            missing_datasets %>%
            plural('was', 'were')
        )
        log_warn(msg)
        warning(msg)

    }

    transmitter_topology %<>%
        recode(!!!topology_short) %>%
        walk(check_topo) %>%
        intersect(topologies) %>%
        sprintf('%s_intercell_source', .) %>%
        paste(collapse = ' | ')

    receiver_topology %<>%
        recode(!!!topology_short) %>%
        walk(check_topo) %>%
        intersect(topologies) %>%
        sprintf('%s_intercell_target', .) %>%
        paste(collapse = ' | ')

    datasets %<>%
        intersect(colnames(network)) %>%
        paste(collapse = ' | ')

    network %>%
    {`if`(
        is.null(loc_consensus_percentile),
        .,
        select(
            .,
            -(!!!exprs(sprintf('%s_intercell_source', major_locations))),
            -(!!!exprs(sprintf('%s_intercell_target', major_locations)))
        ) %>%
        left_join(consensus_loc, by = c('source' = 'uniprot')) %>%
        left_join(consensus_loc, by = c('target' = 'uniprot'),
            suffix = c('_intercell_source', '_intercell_target')
        )
    )} %>%
    filter(eval(parse_expr(receiver_topology))) %>%
    filter(eval(parse_expr(transmitter_topology))) %>%
    filter(eval(parse_expr(datasets))) %>%
    filter(
        curation_effort >= min_curation_effort &
        n_resources >= min_resources &
        n_references >= min_references &
        (
            n_resources >= min_provenances |
            n_references >= min_provenances
        )
    ) %>%
    inner_join(
        consensus_annot,
        by = c(
            'parent_intercell_source' = 'parent',
            'source' = 'uniprot'
        )
    ) %>%
    inner_join(
        consensus_annot,
        by = c(
            'parent_intercell_target' = 'parent',
            'target' = 'uniprot'
        )
    ) %>%
    {`if`(
        ligand_receptor,
        filter(
            .,
            parent_intercell_source == 'ligand' &
            parent_intercell_target == 'receptor'
        ),
        .
    )} %>%
    {`if`(
        simplify,
        simplify_intercell_network(., ...),
        .
    )} %>%
    {`if`(
        unique_pairs,
        unique_intercell_network(., ...),
        .
    )}

}


#' Simplify an intercell network
#'
#' The intercellular communication network data frames, created by
#' \code{\link{import_intercell_network}}, are combinations of a network data
#' frame with two copies of the intercell annotation data frames, all of them
#' already having quite some columns. Here we keep only the names of the
#' interacting pair, their intercellular communication roles, and the minimal
#' information of the origin of both the interaction and the annotations.
#' Optionally further columns can be selected.
#'
#' @param network An intercell network data frame, as provided by
#'     \code{\link{import_intercell_network}}.
#' @param ... Optional, further columns to select.
#'
#' @return An intercell network data frame with some columns removed.
#'
#' @examples
#' icn <- import_intercell_network()
#' icn_s <- simplify_intercell_network(icn)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom rlang ensyms !!!
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_intercell_network}}}
#'     \item{\code{\link{unique_intercell_network}}}
#'     \item{\code{\link{filter_intercell_network}}}
#' }
simplify_intercell_network <- function(network, ...){

    # NSE vs. R CMD check workaround
    source <- target <- source_genesymbol <- target_genesymbol <-
    category_intercell_source <- database_intercell_source <-
    category_intercell_target <- database_intercell_target <-
    is_directed <- is_stimulation <- is_inhibition <-
    sources <- references <- NULL

    simplify_cols <-
        alist(
            source,
            target,
            source_genesymbol,
            target_genesymbol,
            category_intercell_source,
            database_intercell_source,
            category_intercell_target,
            database_intercell_target,
            is_directed,
            is_stimulation,
            is_inhibition,
            sources,
            references
        ) %>%
        c(ensyms(...)) %>%
        unique

    network %>%
    select(!!!simplify_cols)

}


#' Unique intercellular interactions
#'
#' In the intercellular network data frames produced by \code{
#' \link{import_intercell_network}}, by default each pair of annotations for
#' an interaction is represented in a separate row. This function drops the
#' annotations and keeps only the distinct interacting pairs.
#'
#' @param network An intercellular network data frame as produced by
#'     \code{\link{import_intercell_network}}.
#' @param ... Additional columns to keep. Note: if these have multiple
#'     values for an interacting pair, only the first row will be
#'     preserved.
#'
#' @return A data frame with interacting pairs and interaction attributes.
#'
#' @examples
#' icn <- import_intercell_network()
#' icn_unique <- unique_intercell_network(icn)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct
#' @importFrom rlang !!!
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_intercell_network}}}
#'     \item{\code{\link{simplify_intercell_network}}}
#'     \item{\code{\link{filter_intercell_network}}}
#' }
unique_intercell_network <- function(network, ...){

    # NSE vs. R CMD check workaround
    source <- target <- source_genesymbol <- target_genesymbol <-
    is_directed <- is_stimulation <- is_inhibition <-
    sources <- references <- NULL

    cols <-
        alist(
            source,
            target,
            source_genesymbol,
            target_genesymbol,
            is_directed,
            is_stimulation,
            is_inhibition,
            sources,
            references
        ) %>%
        c(ensyms(...)) %>%
        unique

    network %>%
    select(!!!cols) %>%
    distinct(source, target)

}


#' Categories in the intercell database of OmniPath
#'
#' Retrieves a list of categories from \url{https://omnipath.org/intercell}.
#'
#' @return character vector with the different intercell categories
#' @export
#'
#' @examples
#' get_intercell_categories()
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_omnipath_intercell}}}
#'     \item{\code{\link{get_intercell_generic_categories}}}
#' }
get_intercell_categories <- function(){

    return(
        unique(
            import_omnipath('intercell_summary', license = NA)$category
        )
    )

}


#' Full list of intercell categories and resources
#'
#' @return A data frame of categories and resources.
#'
#' @examples
#' ic_cat <- intercell_categories()
#' ic_cat
#' # # A tibble: 1,125 x 3
#' #    category                parent                  database
#' #    <chr>                   <chr>                   <chr>
#' #  1 transmembrane           transmembrane           UniProt_location
#' #  2 transmembrane           transmembrane           UniProt_topology
#' #  3 transmembrane           transmembrane           UniProt_keyword
#' #  4 transmembrane           transmembrane_predicted Phobius
#' #  5 transmembrane_phobius   transmembrane_predicted Almen2009
#' # # . with 1,120 more rows
#'
#' @export
intercell_categories <- function(){

    import_omnipath('intercell_summary', license = NA)

}

#' Retrieves a list of the generic categories in the intercell database
#' of OmniPath
#'
#' Retrieves a list of the generic categories from
#' \url{https://omnipath.org/intercell}.
#'
#' @return character vector with the different intercell main classes
#' @export
#'
#' @examples
#' get_intercell_generic_categories()
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_omnipath_intercell}}}
#'     \item{\code{\link{get_intercell_categories}}}
#' }
#'
#' @aliases get_intercell_classes
get_intercell_generic_categories <- function(){

    return(
        unique(
            import_omnipath('intercell_summary', license = NA)$parent
        )
    )
}


# Aliases (old names) to be deprecated
#' @rdname get_intercell_generic_categories
#' @param ... Passed to \code{get_intercell_generic_categories}.
#' @export
#'
#' @noRd
get_intercell_classes <- function(...){
    .Deprecated("get_intercell_generic_categories")
    get_intercell_generic_categories(...)
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

        # unfortunately the column title is different across the various
        # query types, so we need to guess
        for(field in c('sources', 'database', 'source')){

            if(field %in% names(data)){

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
                break

            }

        }

        after <- nrow(data)

        logger::log_success(
            'Filtering by resources: removed %d records.',
            before - after
        )

    }

    return(data)
}


#' Alias for old function name.
#' @param ... Passed to \code{filter_by_resource}.
#'
#' @noRd
filter_sources <- function(...){
    .Deprecated("filter_by_resource")
    filter_by_resource(...)
}


## Filtering intercell records according to the categories and/or classes
## selected
#TODO: actually this we could export as it might be useful for users
## Filters an intercell data table according to various criteria
#' Filters intercell annotations
#'
#' Filters a data frame retrieved by \code{\link{import_omnipath_intercell}}.
#'
#' @importFrom dplyr recode rename_all
#' @importFrom magrittr %>%
#' @importFrom rlang .data set_names
#'
#'
#' @noRd
filter_intercell <- function(
    data,
    categories = NULL,
    resources = NULL,
    parent = NULL,
    scope = NULL,
    aspect = NULL,
    source = NULL,
    transmitter = NULL,
    receiver = NULL,
    secreted = NULL,
    plasma_membrane_peripheral = NULL,
    plasma_membrane_transmembrane = NULL,
    proteins = NULL,
    causality = NULL,
    topology = NULL,
    ...
){

    before <- nrow(data)

    topology <-
        data.frame(set_names(
            as.list(topology),
            rep(TRUE, length(topology))
        )) %>%
        rename_all(
            recode,
            secreted = 'sec',
            plasma_membrane_peripheral = 'pmp',
            plasma_membrane_transmembrane = 'pmtm'
        )

    if('both' %in% causality){
        causality <- c('transmitter', 'receiver')
    }
    causality <-
        data.frame(set_names(
            as.list(causality),
            rep(TRUE, length(causality))
        )) %>%
        rename_all(
            recode,
            transmitter = 'trans',
            receiver = 'rec'
        )

    data <-
        data %>%
        filter(
            (is.null(categories) | category %in% categories) &
            (is.null(parent) | .data$parent %in% parent) &
            (is.null(scope) | .data$scope %in% scope) &
            (is.null(aspect) | .data$aspect %in% aspect) &
            (is.null(source) | .data$source %in% source) &
            (is.null(transmitter) | .data$transmitter) &
            (is.null(receiver) | .data$receiver) &
            (is.null(secreted) | .data$secreted) &
            (
                is.null(plasma_membrane_peripheral) |
                .data$plasma_membrane_peripheral
            ) &
            (
                is.null(plasma_membrane_transmembrane) |
                .data$plasma_membrane_transmembrane
            ) &
            (
                is.null(proteins) |
                uniprot %in% proteins |
                genesymbol %in% proteins
            )
        ) %>%
        {`if`(
            ncol(topology) > 0,
            inner_join(., topology, by = names(topology)),
            .
        )} %>%
        {`if`(
            ncol(causality) > 0,
            inner_join(., causality, by = names(causality)),
            .
        )}

    after <- nrow(data)

    message(
        sprintf(
            'Removed %d records from intercell data.',
            before - after
        )
    )

    return(data)

}


#' Downloader dedicated to OmniPath web service URLs
#'
#' Just a thin wrapper around \code{download_base}.
#'
#' @param url Character: the URL to download.
#' @param fun The downloader function. Should be able to accept \code{url}
#'     as its first argument.
#' @param ... Passed to \code{fun}.
#'
#' @noRd
omnipath_download <- function(url, fun, ...) {

    from_cache <- omnipath_cache_load(url = url)

    if(!is.null(from_cache)){

        log_info('Loaded from cache: `%s`', url)
        return(from_cache)

    }

    result <- download_base(url, fun, ...)

    omnipath_cache_save(data = result, url = url)

    return(result)

}