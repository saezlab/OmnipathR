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
    from_cache_file = 'cache_file',
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
    'causality'
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
    directed = 'directed',
    entity_type = 'entity_types'
)

## Downloads data from the OmniPath web service
## Generic method for retrieval of a table and creating a data frame.
## All methods specific for certain query types or datasets use this function
## to manage the download.
## Not exported.
import_omnipath <- function(
    query_type,
    organism = 9606,
    resources = NULL,
    datasets = NULL,
    cache_file = NULL,
    genesymbols = 'yes',
    fields = NULL,
    default_fields = TRUE,
    silent = FALSE,
    logicals = NULL,
    download_args = list(),
    references_by_resource = TRUE,
    add_counts = TRUE,
    ...
){

    param <- c(as.list(environment()), list(...))
    param <- omnipath_check_param(param)

    if(!is.null(cache_file) && file.exists(cache_file)){
        loaded <- load(cache_file)
        if(length(loaded) > 0){
            result <- get(loaded[1])
        }else{
            stop(sprintf('Cache file `%s` yielded no data.', cache_file))
        }
        result <- filter_by_resource(result, resources = resources)
        msg <- 'Loaded %d %s from cache.'
    }else{

        url <- omnipath_url(param)
        download_args_defaults <- list(
            URL = url
        )
        dataframe_defaults <- list(
            FUN = read.table,
            header = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE,
            quote = ''
        )
        json_defaults <- list(
            FUN = jsonlite::fromJSON
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
        if(!is.null(cache_file)){
            save(result, file = cache_file)
        }
        msg <- 'Downloaded %d %s.'

    }

    result <- cast_logicals(result, logicals)
    result <- strip_resource_labels(result, references_by_resource)
    if(param$query_type %in% c('interactions', 'enzsub') && add_counts){
        result <- count_references(result)
        result <- count_resources(result)
    }


    if(!silent){
        message(sprintf(
            msg,
            `if`(
                is.data.frame(result),
                nrow(result),
                length(result)
            ),
            param$qt_message)
        )
    }

    return(result)

}


## Check the arguments of \link{import_omnipath}, corrects some easy to
## confuse or deprecated synonyms and selects the message printed by
## the download function.
## Not exported.
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
        warning('DoRothEA confidence levels available are A, B, C and D.')
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

    return(param)

}


## Constructs the URL by creating a base URL according to the query type and
## adding all user or package defined query string parameters.
## Not exported.
omnipath_url <- function(param){

    baseurl <- options('omnipath.url')
    baseurl <- sprintf('%s%s', baseurl, param$query_type)

    url <- Reduce(
        function(url, key){
            omnipath_url_add_param(url, key, param[[key]])
        },
        .omnipath_querystring_param,
        init = baseurl
    )

    return(url)

}


## Appends a query string parameter to the URL.
## Not exported, used internally for assembling the URLs.
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
        is.null(values),
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


## Checks whether the response is real data or an error message.
## In case of error stops the execution and prints the URL and the message
## from the server.
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

## Makes sure the boolean variables, listed in argument `logicals`, are of
## R logical type. Converts various string and numeric representations.
## Checks only for TRUE values, whatever does not match remains FALSE.
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


## Removes the resource labels from references (PubMed IDs) in the
## interactions and enzyme-substrate data frames.
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


## For a character vector splits each element and re-joins sorted unique
## values.
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


## For a character vector splits each element and applies a method for
## each sub vector.
split_apply <- function(
    x,
    method,
    sep = ';',
    ...
){
    return(
        sapply(
            strsplit(x, sep),
            function(x){method(x, ...)}
        )
    )
}


## For an interactions or enzyme-substrate data frame adds a column
## `n_resources` with the number of resources for each record.
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


## For an interactions or enzyme-substrate data frame adds a column
## `n_references` with the number of references for each record.
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

## For each undirected interaction adds a duplicate with the source and
## target nodes swapped
swap_undirected <- function(data){

    data <- data %>%
        dplyr::filter(is_directed == 0) %>%
        dplyr::rename(
            source = target,
            target = source,
            source_genesymbol = target_genesymbol,
            target_genesymbol = source_genesymbol
        ) %>%
        {`if`(
            'ncbi_tax_id_source' %in% names(.),
            dplyr::rename(
                .,
                ncbi_tax_id_source = ncbi_tax_id_target,
                ncbi_tax_id_target = ncbi_tax_id_source
            ),
            .
        )} %>%
        dplyr::bind_rows(data)

    return(data)

}

########## ########## ########## ##########
########## Enzyme-substrate      ##########
########## ########## ########## ##########

#' Imports enzyme-substrate relationships from OmniPath
#'
#' Imports the enzyme-substrate (more exactly, enzyme-PTM) relationship
#' database from \url{http://omnipathdb.org/enzsub}
#'
#' @return A data frame containing the information about ptms
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
#' @param resources PTMs not reported in these databases are
#' removed. See \code{\link{get_ptms_databases}} for more information
#' @param organism PTMs are available for human, mouse and rat.
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
#' ptms = import_omnipath_enzsub(
#'     resources = c('PhosphoSite', 'SIGNOR'),
#'     organism = 9606
#' )
#'
#' @seealso \code{\link{get_ptms_databases},
#'   \link{import_omnipath_interactions}}
#'
#' @aliases import_Omnipath_PTMS import_OmniPath_PTMS
import_omnipath_enzsub <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'enzsub',
        cache_file = cache_file,
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
import_Omnipath_PTMS <- function(...){
    .Deprecated("import_omnipath_enzsub")
    import_omnipath_enzsub(...)
} 
# Aliases (old names) to be Deprecated
#' @rdname import_omnipath_enzsub
#' @param ... Passed to \code{import_omnipath_enzsub}.
#' @export
import_OmniPath_PTMS <- function(...){
    .Deprecated("import_omnipath_enzsub")
    import_omnipath_enzsub(...)
}


#' Retrieve a list of enzyme-substrate resources available in OmniPath
#'
#' get the names of the enzyme-substrate relationship resources available
#' in \url{http://omnipath.org/enzsub}
#'
#' @param dataset ignored for this query type
#' @return character vector with the names of the enzyme-substrate resources
#' @export
#' @importFrom utils read.table
#'
#' @examples
#' get_enzsub_resources()
#'
#' @seealso  \code{\link{get_resources},
#' \link{import_omnipath_enzsub}}
#'
#' @aliases get_ptms_databases
get_enzsub_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'enzsub', datasets = dataset))

}

# Aliases (old names) to be deprecated
#' @rdname get_enzsub_resources
#' @param ... Passed to \code{get_enzsub_resources}.
#' @export
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
#' imports the database from \url{http://omnipathdb.org/interactions}, which
#' contains only interactions supported by literature references.
#' This part of the interaction database compiled a similar way as it has
#' been presented in the first paper describing OmniPath (Turei et al. 2016).
#'
#' @return A dataframe of protein-protein interactions
#' @export
#' @importFrom utils read.table
#'
#' @param cache_file path to an earlier data file
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param datasets Names of the interaction datasts to download: omnipath 
#' (by default). Other possiblites are: pathwayextra, kinaseextra, ligrecextra,
#' dorothea,tf_target, mirnatarget, tf_mirna, lncrna_mrna. The user can select 
#' multiple datasets as for example: c('omnipath', 'pathwayextra','kinaseextra')
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
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
#'
#' @aliases import_Omnipath_Interactions import_OmniPath_Interactions
import_omnipath_interactions <- function(
    cache_file = NULL,
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
        cache_file = cache_file,
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
import_Omnipath_Interactions <- function(...){
    .Deprecated("import_omnipath_interactions")
    import_omnipath_interactions(...)
}
    
# Aliases (old names) to be deprecated
#' @rdname import_omnipath_interactions
#' @param ... Passed to \code{import_omnipath_interactions}.
#' @export
import_OmniPath_Interactions <- function(...){
    .Deprecated("import_omnipath_interactions")
    import_omnipath_interactions(...)
}


#' Imports interactions from the `pathway extra` dataset of Omnipath
#'
#' Imports the dataset from:
#' \url{http://omnipathdb.org/interactions?datasets=pathwayextra},
#' which contains activity flow interactions without literature reference.
#' The activity flow interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @return A dataframe containing activity flow interactions between proteins
#' without literature reference
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
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
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
#'
#' @aliases import_PathwayExtra_Interactions
import_pathwayextra_interactions <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
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
#' @rdname import_pathwayextra_interactions
#' @param ... Passed to \code{import_pathwayextra_interactions}.
#' @export
import_PathwayExtra_Interactions <- function(...){
    .Deprecated("import_pathwayextra_interactions")
    import_pathwayextra_interactions(...)
}

#' Imports interactions from the `kinase extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{http://omnipathdb.org/interactions?datasets=kinaseextra},
#' which contains enzyme-substrate interactions without literature reference.
#' The enzyme-substrate interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @return A dataframe containing enzyme-substrate interactions without
#' literature reference
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
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
#'    import_kinaseextra_interactions(
#'        resources = c('PhosphoPoint', 'PhosphoSite'),
#'        organism = 9606
#'    )
#'
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
#'
#' @aliases import_KinaseExtra_Interactions
import_kinaseextra_interactions <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    fields = NULL, 
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
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
import_KinaseExtra_Interactions <- function(...){
    .Deprecated("import_kinaseextra_interactions")
    import_kinaseextra_interactions(...)
}

#' Imports interactions from the `ligrec extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{http://omnipathdb.org/interactions?datasets=ligrecextra},
#' which contains ligand-receptor interactions without literature reference.
#' The ligand-receptor interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @return A dataframe containing ligand-receptor interactions including
#' the ones without literature references
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
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
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
#'
#' @aliases import_LigrecExtra_Interactions
import_ligrecextra_interactions <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    fields = NULL, 
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
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
import_LigrecExtra_Interactions <- function(...){
    .Deprecated("import_ligrecextra_interactions")
    import_ligrecextra_interactions(...)
}

#' Imports all post-translational interactions from OmniPath
#'
#' Imports the dataset from all post-translational datasets of OmniPath
#'
#' @return A dataframe containing post-translational interactions
#' @export
#' @importFrom utils read.table
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
#' @examples
#' interactions <-
#'     import_post_translational_interactions(
#'         resources = c('BioGRID')
#'     )
#'
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
import_post_translational_interactions <- function(
    resources = NULL,
    organism = 9606,
    exclude = NULL,
    references_by_resource = TRUE,
    ...
){

    datasets <- c('omnipath', 'pathwayextra', 'kinaseextra', 'ligrecextra')
    datasets <- setdiff(datasets, exclude)


    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = datasets,
        references_by_resource = references_by_resource,
        ...
    )

    return(result)

}


#' Imports from Omnipath webservice the interactions from
#' Dorothea dataset
#'
#' Imports the dataset from:
#' \url{http://omnipathdb.org/interactions?datasets=dorothea}
#' which contains transcription factor (TF)-target interactions from DoRothEA
#' \url{https://github.com/saezlab/DoRothEA}
#'
#' @return A dataframe containing TF-target interactions from DoRothEA
#' @export
#' @importFrom utils read.table
#'
#' @param cache_file path to an earlier data file
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
#'     resources = c('DoRothEA_A', 'ARACNe-GTEx_DoRothEA'),
#'     organism = 9606
#' )
#'
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
#'
#' @aliases import_TFregulons_Interactions import_tfregulons_interactions
import_dorothea_interactions <- function(
    cache_file = NULL,
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
        cache_file = cache_file,
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
import_TFregulons_Interactions <- function(...){
    .Deprecated("import_dorothea_interactions")
    import_dorothea_interactions(...)
}
    
#' @rdname import_dorothea_interactions
#' @param ... Passed to \code{import_dorothea_interactions}.
#' @export
import_tfregulons_interactions <- function(...){
    .Deprecated("import_dorothea_interactions")
    import_dorothea_interactions(...)
}
    
#' Imports interactions from the TF-target dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{http://omnipathdb.org/interactions?datasets=tf_target},
#' which contains transcription factor-target protein coding gene
#' interactions. Note: this is not the only TF-target dataset in OmniPath,
#' `dorothea` is the other one and the `tf_mirna` dataset provides
#' TF-miRNA gene interactions.
#'
#' @return A dataframe containing TF-target interactions
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
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
#'     import_tf_target_interactions(
#'         resources = c('DoRothEA_A', 'SIGNOR')
#'     )
#'
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
#'
import_tf_target_interactions <- function(
    cache_file = NULL,
    resources = NULL, 
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
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
#' \url{http://omnipathdb.org/interactions?datasets=tf_target,dorothea},
#' which contains transcription factor-target protein coding gene
#' interactions.
#'
#' @return A dataframe containing TF-target interactions
#' @export
#' @importFrom utils read.table
#' @importFrom dplyr %>% mutate select
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
#'         resources = c('PAZAR', 'ORegAnno', 'DoRothEA_A')
#'     )
#'
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
import_transcriptional_interactions <- function(
    resources = NULL, 
    organism = 9606,
    dorothea_levels = c('A', 'B'),
    references_by_resource = TRUE,
    ...
){

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
            dplyr::mutate(dorothea_level = "") %>% 
            dplyr::select(source, target, source_genesymbol, target_genesymbol, 
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
#' \url{http://omnipathdb.org/interactions?datasets=mirnatarget},
#' which contains miRNA-mRNA interactions
#'
#' @return A dataframe containing miRNA-mRNA interactions
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
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
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
#'
#' @aliases import_miRNAtarget_Interactions
import_mirnatarget_interactions <- function(
    cache_file = NULL,
    resources = resources,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
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
import_miRNAtarget_Interactions <- function(...){
    .Deprecated("import_mirnatarget_interactions")
    import_mirnatarget_interactions(...)
}

#' Imports interactions from the TF-miRNA dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{http://omnipathdb.org/interactions?datasets=tf_mirna},
#' which contains transcription factor-miRNA gene interactions
#'
#' @return A dataframe containing TF-miRNA interactions
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
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
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
import_tf_mirna_interactions <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    fields = NULL, 
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
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
#' \url{http://omnipathdb.org/interactions?datasets=lncrna_mrna},
#' which contains lncRNA-mRNA interactions
#'
#' @return A dataframe containing lncRNA-mRNA interactions
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
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
#' @seealso \code{\link{get_interaction_resources},
#'   \link{import_all_interactions}}
import_lncrna_mrna_interactions <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
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
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param dorothea_levels The confidence levels of the dorothea 
#' interactions (TF-target) which range from A to D. Set to A and B by default. 
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
#' @seealso \code{\link{get_interaction_resources}}
#'
#' @aliases import_AllInteractions
import_all_interactions <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    dorothea_levels = c('A', 'B'),
    exclude = NULL,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    all_datasets <- jsonlite::fromJSON(
        txt = 'http://omnipathdb.org/queries/interactions?format=json'
    )$datasets

    all_datasets <- setdiff(all_datasets, exclude)

    # it does not make sense without the type field
    fields <- unique(c(fields), 'type')

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
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
#' gets the names of the resources from \url{http://omnipath.org/interactions}
#'
#' @param dataset a dataset within the interactions query type. Currently
#' available datasets are `omnipath`, `kinaseextra`, `pathwayextra`,
#' `ligrecextra`, `dorothea`, `tf_target`, `tf_mirna`, `mirnatarget` and
#' `lncrna_mrna`
#'
#' @return character vector with the names of the interaction databases
#' @export
#' @importFrom utils read.table
#'
#' @examples
#' get_interaction_resources()
#'
#' @seealso \code{\link{get_resources},
#' \link{import_all_interactions},
#' \link{import_omnipath_interactions}, \link{import_pathwayextra_interactions},
#' \link{import_kinaseextra_interactions},
#' \link{import_ligrecextra_interactions},
#' \link{import_mirnatarget_interactions},
#' \link{import_dorothea_interactions}}
#'
#' @aliases get_interaction_databases
get_interaction_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'interactions', datasets = dataset))

}

# Aliases (old names) to be deprecated
#' @rdname get_interaction_resources
#' @param ... Passed to \code{get_interaction_resources}.
#' @export
get_interaction_databases <- function(...){
    .Deprecated("get_interaction_resources")
    get_interaction_resources(...)
}

#' Retrieve the available resources for a given query type
#'
#' collects the names of the resources available in OmniPath for a certain
#' query type and optionally for a dataset within that.
#'
#' @param query_type one of the query types `interactions`, `enz_sub`,
#' `complexes`, `annotations` or `intercell`
#' @param datasets currently within the `interactions` query type only,
#' multiple datasets are available: `omnipath`, `kinaseextra`, `pathwayextra`,
#' `ligrecextra`, `dorothea`, `tf_target`, `tf_mirna`, `mirnatarget` and
#' `lncrna_mrna`
#' @param generic_categories for the `intercell` query type, restrict the
#' search for some generic categories e.g. `ligand` or `receptor`
#'
#' @return a character vector with resource names
#' @export
#'
#' @import jsonlite
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

    resources <- jsonlite::fromJSON(txt = 'http://omnipathdb.org/resources')

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

#' Imports protein complexes from Omnipath
#'
#' imports the complexes stored in Omnipath database from
#' \url{http://omnipathdb.org/complexes}
#'
#' @return A dataframe containing information about complexes
#' @export
#' @importFrom utils read.table
#'
#' @param cache_file path to an earlier data file
#' @param resources complexes not reported in these databases are
#' removed. See \code{\link{get_complexes_databases}} for more information.
#' @param ... optional additional arguments 
#'
#' @examples
#' complexes = import_omnipath_complexes(
#'     resources = c('CORUM', 'hu.MAP')
#' )
#'
#' @seealso \code{\link{get_complexes_databases}}
#'
#' @aliases import_Omnipath_complexes import_OmniPath_complexes
import_omnipath_complexes <- function(
    cache_file = NULL,
    resources = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'complexes',
        cache_file = cache_file,
        resources = resources,
        ...
    )

    return(result)

}

# Aliases (old names) to be deprecated
#' @rdname import_omnipath_complexes
#' @param ... Passed to \code{import_omnipath_complexes}.
#' @export
import_Omnipath_complexes <- function(...){
    .Deprecated("import_omnipath_complexes")
    import_omnipath_complexes(...)
}

#' @rdname import_omnipath_complexes
#' @param ... Passed to \code{import_omnipath_complexes}.
#' @export
import_OmniPath_complexes <- function(...){
    .Deprecated("import_omnipath_complexes")
    import_omnipath_complexes(...)    
}

#' Retrieve a list of complex resources available in Omnipath
#'
#' get the names of the resources from \url{http://omnipath.org/complexes}
#' @param dataset ignored for this query type
#' @return character vector with the names of the databases
#' @export
#' @importFrom utils read.csv
#'
#' @examples
#' get_complex_resources()
#'
#' @seealso \code{\link{get_resources},
#' \link{import_omnipath_complexes}}
#'
#' @aliases get_complexes_databases
get_complex_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'complexes', datasets = dataset))

}

# Aliases (old names) to be deprecated
#' @rdname get_complex_resources
#' @param ... Passed to \code{import_omnipath_enzsub}.
#' @export
get_complexes_databases <- function(...){
    .Deprecated("get_complex_resources")
    get_complex_resources(...)
}

########## ########## ########## ##########
########## Annotations           ##########
########## ########## ########## ##########

#' Imports annotations from OmniPath
#'
#' imports protein annotations about function, localization, expression,
#' structure and other properties of proteins from OmniPath
#' \url{http://omnipathdb.org/annotations}.
#' Note: there might be also a few miRNAs annotated; a vast majority of
#' protein complex annotations are inferred from the annotations of the
#' members: if all members carry the same annotation the complex inherits.
#'
#' @return A data.frame containing different gene/complex annotations
#' @export
#' @importFrom utils read.csv
#'
#' @param cache_file Path to an earlier data file
#' @param proteins Vector containing the genes or proteins for whom
#' annotations will be retrieved (UniProt IDs or HGNC Gene Symbols or
#' miRBase IDs). It is also possible to donwload annotations for protein
#' complexes. To do so, write 'COMPLEX:' right before the genesymbols of
#' the genes integrating the complex. Check the vignette for examples.
#' @param resources Load the annotations only from these databases.
#' See \code{\link{get_annotation_resources}} for possible values.
#' @param force_full_download Force the download of the entire annotations
#' dataset. This is disabled by default because the size of this data is
#' around 1GB. We recommend to retrieve the annotations for a set of proteins
#' or only from a few resources, depending on your interest.
#' @param ... Additional arguments.
#' @examples
#' annotations = import_omnipath_annotations(
#'     proteins = c('TP53', 'LMNA'),
#'     resources = c('HPA_subcellular')
#' )
#' @seealso \code{\link{get_annotation_databases}}
#'
#' @aliases import_Omnipath_annotations import_OmniPath_annotations
import_omnipath_annotations <- function(
    cache_file = NULL,
    proteins = NULL,
    resources = NULL,
    force_full_download = FALSE,
    ...
){

    # account for the old argument name
    proteins <- c(proteins, list(...)$select_genes)

    if(
        !force_full_download &&
        is.null(cache_file) &&
        is.null(proteins) &&
        is.null(resources)
    ){

        stop(
            paste(
                'Downloading the entire annotations database is not allowed',
                'by default because of its huge size (>1GB). If you really',
                'want to do this use the `force_full_download` parameter.',
                'However we recommend to query a set of proteins or a few',
                'resources, depending on your interest.'
            )
        )

    }

    if(
        (
            !is.null(cache_file) &&
            file.exists(cache_file)
        ) ||
        length(proteins) < 600
    ){

        result <- import_omnipath(
            query_type = 'annotations',
            proteins = proteins,
            resources = resources,
            cache_file = cache_file,
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
                        query_type = annotations,
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

        result <- do.call(rbind, result)

        if(!is.null(cache_file)){
            save(result, cache_file)
        }

        message(sprintf(
            'Downloaded %d annotation records.',
            nrow(result)
        ))

    }

    return(result)

}

# Aliases (old names) to be deprecated
#' @rdname import_omnipath_annotations
#' @param ... Passed to \code{import_omnipath_annotations}.
#' @export
import_Omnipath_annotations <- function(...){
    .Deprecated("import_omnipath_annotations")
    import_omnipath_annotations(...)
}
#' @rdname import_omnipath_annotations
#' @param ... Passed to \code{import_omnipath_annotations}.
#' @export
import_OmniPath_annotations <- function(...){
    .Deprecated("import_omnipath_annotations")
    import_omnipath_annotations(...)
}
    

#' Get the resources available in the annotations database of OmniPath
#'
#' get the names of the resources from \url{http://omnipath.org/annotations}
#'
#' @return character vector with the names of the annotation resources
#' @export
#' @param dataset ignored for this query type
#' @param ... optional additional arguments 
#'
#' @examples
#' get_annotation_resources()
#'
#' @seealso \code{\link{get_resources},
#' \link{import_omnipath_annotations}}
#'
#' @aliases get_annotation_databases
get_annotation_resources <- function(dataset = NULL, ...){

    return(get_resources(query_type = 'annotations', datasets = dataset))

}

# Aliases (old names) to be deprecated
#' @rdname get_annotation_resources
#' @param ... Passed to \code{get_annotation_resources}.
#' @export
get_annotation_databases <- function(...){
    .Deprecated("get_annotation_resources")
    get_annotation_resources(...)
}

########## ########## ########## ##########
########## Intercell             ##########
########## ########## ########## ##########

#' Imports OmniPath Intercell Data
#'
#' imports the OmniPath intercellular communication role annotation database
#' from \url{http://omnipathdb.org/intercell}. It provides information
#' on the roles in inter-cellular signaling. E.g. if a protein is
#' a ligand, a receptor, an extracellular matrix (ECM) component, etc.
#'
#' @return A dataframe cotaining information about roles in intercellular
#' signaling.
#' @export
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @param cache_file path to an earlier data file
#' @param categories vector containing the categories to be retrieved.
#' All the genes belonging to those categories will be returned. For further
#' information about the categories see \code{\link{get_intercell_categories}}
#' @param parent vector containing the parent classes to be retrieved.
#' All the genes belonging to those classes will be returned. For furter
#' information about the main classes see
#' \code{\link{get_intercell_categories}}
#' @param resources limit the query to certain resources; see the available
#' resources by \code{\link{get_intercell_resources}}
#' @param scope either `specific` or `generic`
#' @param aspect either `locational` or `functional`
#' @param source either `resource_specific` or `composite`
#' @param transmitter logical, include only transmitters i.e. proteins
#' delivering signal from a cell to
#' its environment
#' @param receiver logical, include only receivers i.e. proteins delivering
#' signal to the cell from its environment
#' @param plasma_membrane_peripheral logical, include only plasma membrane
#' peripheral membrane proteins
#' @param plasma_membrane_transmembrane logical, include only plasma membrane
#' transmembrane proteins
#' @param secreted logical, include only secreted proteins
#' @param proteins limit the query to certain proteins
#' @param topology topology categories: one or more of `secreted` (sec),
#' `plasma_membrane_peripheral` (pmp), `plasma_membrane_transmembrane` (pmtm)
#' (both short or long notation can be used)
#' @param causality `transmitter` (trans), `receiver` (rec) or `both` (both
#' short or long notation can be used)
#' @param ... Additional optional arguments 
#'
#' @examples
#' intercell = import_omnipath_intercell(categories = c('ecm'))
#'
#' @seealso \code{\link{get_intercell_categories},
#' \link{get_intercell_generic_categories}, \link{import_intercell_network}}
#'
#' @aliases import_Omnipath_intercell import_OmniPath_intercell
import_omnipath_intercell <- function(
    cache_file = NULL,
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
    ...
){

    from_cache <- !is.null(cache_file) && file.exists(cache_file)
    args <- c(as.list(environment()), list(...))
    args$query_type <- 'intercell'
    args$logicals <- c(
        'transmitter',
        'receiver',
        'secreted',
        'plasma_membrane_peripheral',
        'plasma_membrane_transmembrane'
    )

    result <- do.call(import_omnipath, args)

    if(from_cache){
        args$data <- result
        result <- do.call(filter_intercell, args)
    }

    return(result)

}

# Aliases (old names) to be deprecated
#' @rdname import_omnipath_intercell
#' @param ... Passed to \code{import_omnipath_intercell}.
#' @export
import_Omnipath_intercell <- function(...){
    .Deprecated("import_omnipath_intercell")
    import_omnipath_intercell(...)
}
#' @rdname import_omnipath_intercell
#' @param ... Passed to \code{import_omnipath_intercell}.
#' @export
import_OmniPath_intercell <- function(...){
    .Deprecated("import_omnipath_intercell")
    import_omnipath_intercell(...)
}

#' Retrieve a list of intercellular communication resources available in
#' Omnipath
#'
#' get the names of the databases from \url{http://omnipath.org/intercell}
#' @return character vector with the names of the databases
#' @export
#' @importFrom utils read.csv
#' @param dataset ignored at this query type
#'
#' @examples
#' get_intercell_resources()
#'
#' @seealso \code{\link{get_resources},
#' \link{import_omnipath_intercell}}
get_intercell_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'intercell', datasets = dataset))

}


#' Imports an intercellular network combining annotations and interactions
#'
#' Imports an intercellular network by mapping intercellular annotations
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
#' If a cache file provided, its content will be returned without any further
#' filtering.
#'
#' @return A dataframe containing information about protein-protein
#' interactions and the inter-cellular roles of the protiens involved in those
#' interactions.
#' @export
#' @importFrom utils read.csv modifyList
#' @importFrom dplyr %>% rename bind_rows filter inner_join distinct group_by
#' summarize_all first
#'
#' @param cache_file path to an earlier data file; if exists, will be loaded
#' as it is, the further arguments have no effect; if does not exists, the
#' result will be dumped into this file.
#' @param interactions_param a list with arguments for an interactions query: 
#' \code{\link{import_omnipath_interactions}, 
#' \link{import_pathwayextra_interactions},
#' \link{import_kinaseextra_interactions},
#' \link{import_ligrecextra_interactions}}
#' @param transmitter_param a list with arguments for
#' \code{\link{import_omnipath_intercell}}, to define the transmitter side
#' of intercellular connections
#' @param receiver_param a list with arguments for
#' \code{\link{import_omnipath_intercell}}, to define the receiver side
#' of intercellular connections
#'
#' @examples
#' intercellNetwork <- import_intercell_network(
#'    interactions_param = list(datasets = 'ligrecextra'),
#'    receiver_param = list(categories = c('receptor', 'transporter')),
#'    transmitter_param = list(categories = c('ligand', 'secreted_enzyme')))
#'
#' @seealso \code{\link{get_intercell_categories},
#' \link{get_intercell_generic_categories}, 
#' \link{import_omnipath_intercell},
#' \link{import_omnipath_interactions}, 
#' \link{import_pathwayextra_interactions},
#' \link{import_kinaseextra_interactions}, 
#' \link{import_ligrecextra_interactions}}
import_intercell_network <- function(
    cache_file = NULL,
    interactions_param = list(),
    transmitter_param = list(),
    receiver_param = list()
){

    result <- NULL

    if(!is.null(cache_file) && file.exists(cache_file)){
        loaded <- load(cache_file)
        if(length(loaded) > 0){
            result <- get(loaded[1])
        }
    }

    if(is.null(result)){

        interactions_param_default <- list(
            query_type = 'interactions',
            datasets = c(
                'omnipath',
                'pathwayextra',
                'kinaseextra',
                'ligrecextra'
            )
        )
        interactions_param <- modifyList(
            interactions_param_default,
            interactions_param
        )
        interactions <- do.call(
            import_omnipath,
            interactions_param
        )
        interactions <- swap_undirected(interactions)

        transmitter_param_defaults <- list(
            causality = 'trans',
            scope = 'generic'
        )
        transmitter_param <- modifyList(
            transmitter_param_defaults,
            transmitter_param
        )

        receiver_param_defaults <- list(
            causality = 'rec',
            scope = 'generic'
        )
        receiver_param <- modifyList(
            receiver_param_defaults,
            receiver_param
        )

        intracell <- c('intracellular_intercellular_related', 'intracellular')
        transmitters <-
            do.call(import_omnipath_intercell, transmitter_param) %>%
            dplyr::filter(!parent %in% intracell) %>%
            dplyr::rename(category_source = source)
        receivers <-
            do.call(import_omnipath_intercell, receiver_param) %>%
            dplyr::filter(!parent %in% intracell) %>%
            dplyr::rename(category_source = source)

        result <-
            interactions %>%
            dplyr::inner_join(
                transmitters,
                by = c('source' = 'uniprot')
            ) %>%
            dplyr::group_by(
                category, parent, source, target
            ) %>%
            dplyr::mutate(
                database = paste(database, collapse = ';')
            ) %>%
            dplyr::summarize_all(first) %>%
            dplyr::inner_join(
                receivers,
                by = c('target' = 'uniprot'),
                suffix = c('_intercell_source', '_intercell_target')
            ) %>%
            dplyr::group_by(
                category_intercell_source,
                parent_intercell_source,
                source,
                target,
                category_intercell_target,
                parent_intercell_target
            ) %>%
            dplyr::mutate(
                database_intercell_target = paste(
                    database_intercell_target,
                    collapse = ';'
                )
            ) %>%
            dplyr::summarize_all(first) %>%
            dplyr::ungroup()

        if(!is.null(cache_file)){
            save(result, cache_file)
        }

    }

    return(result)

}


#' Retrieve a list of categories from the intercell database of OmniPath
#'
#' get the names of the categories from \url{http://omnipath.org/intercell}
#' @return character vector with the different intercell categories
#' @export
#' @importFrom utils read.csv
#' @examples
#' get_intercell_categories()
#' @seealso \code{\link{import_omnipath_intercell},
#' \link{get_intercell_classes}}
get_intercell_categories <- function(){

    return(
        unique(
            import_omnipath('intercell_summary')$category
        )
    )

}


#' Retrieve a list of the generic categories in the intercell database
#' of OmniPath
#'
#' get the names of the generic categories from
#' \url{http://omnipath.org/intercell}
#' @return character vector with the different intercell main classes
#' @export
#' @importFrom utils read.csv
#' @examples
#' get_intercell_generic_categories()
#' @seealso \code{\link{import_omnipath_intercell},
#' \link{get_intercell_categories}}
get_intercell_generic_categories <- function(){

    return(
        unique(
            import_omnipath('intercell_summary')$parent
        )
    )
}


# Aliases (old names) to be deprecated
#' @rdname get_intercell_generic_categories
#' @param ... Passed to \code{get_intercell_generic_categories}.
#' @export
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

## Filtering Interactions, PTMs and complexes
#TODO: actually this we could export as it might be useful for users
filter_by_resource <- function(data, resources = NULL){

    if(!is.null(resources)){

        before <- nrow(data)

        # unfortunately the column title is different across the various
        # query types, so we need to guess
        for(field in c('sources', 'database', 'source')){

            if(field %in% names(data)){

                data <-data[
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

        message(
            sprintf(
                'Filtering by resources: removed %d records.',
                before - after
            )
        )

    }

    return(data)
}

# synonym (old name)
filter_sources <- filter_by_resource


## Filtering intercell records according to the categories and/or classes
## selected
#TODO: actually this we could export as it might be useful for users
## Filters an intercell data table according to various criteria
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
        data.frame(setNames(
            as.list(topology),
            rep(TRUE, length(topology))
        )) %>%
        dplyr::rename_all(
            dplyr::recode,
            secreted = 'sec',
            plasma_membrane_peripheral = 'pmp',
            plasma_membrane_transmembrane = 'pmtm'
        )

    if('both' %in% causality){
        causality <- c('transmitter', 'receiver')
    }
    causality <-
        data.frame(setNames(
            as.list(causality),
            rep(TRUE, length(causality))
        )) %>%
        dplyr::rename_all(
            dplyr::recode,
            transmitter = 'trans',
            receiver = 'rec'
        )

    data <-
        data %>%
        dplyr::filter(
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
            dplyr::inner_join(., topology, by = names(topology)),
            .
        )} %>%
        {`if`(
            ncol(causality) > 0,
            dplyr::inner_join(., causality, by = names(causality)),
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

########## ########## ########## ##########
########## Resource Queries      ##########
########## ########## ########## ##########
## This function is convenient for appropriate resource retrieval. Following:
## http://bioconductor.org/developers/how-to/web-query/
## It tries to retrieve the resource one or several times before failing.
omnipath_download <- function(URL, FUN, ..., N.TRIES = 1L) {
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

    while (N.TRIES > 0L) {
        result <- tryCatch(FUN(URL, ...), error = identity)
        if (!inherits(result, 'error'))
            break
            N.TRIES <- N.TRIES - 1L
        }

    if (N.TRIES == 0L) {
        stop(
            sprintf(
                'omnipath_download() failed:\n  URL: %s\n  error: %s',
                URL,
                conditionMessage(result)
            )
        )
    }

    return(result)
}
