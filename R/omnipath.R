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
    'qt_message',
    'exclude'
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
    ...
){

    datasets %<>% setdiff(exclude)
    resources %<>% setdiff(exclude)

    param <- c(as.list(environment()), list(...))
    param <- omnipath_check_param(param)

    url <- omnipath_url(param)
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

    result %<>% cast_logicals(logicals)
    result %<>% strip_resource_labels(references_by_resource)
    result %<>% apply_exclude(exclude)
    result %<>% deserialize_extra_attrs()

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


#' Converts the extra_attrs column from JSON encoded to list
#'
#' @param data A data frame from the OmniPath web service.
#'
#' @return The input data frame with the extra_attrs column converted
#'     to list.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom jsonlite fromJSON
#' @importFrom purrr map
#' @noRd
deserialize_extra_attrs <- function(data){

    data %>%
    {`if`(
        has_extra_attrs(.),
        mutate(., extra_attrs = map(extra_attrs, fromJSON)),
        .
    )}

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
                references,
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
#' download the full database from
#' \url{https://archive.omnipathdb.org/omnipath_webservice_annotations__recent.tsv}
#' using any standard R or \code{readr} method.
#'
#' @return A data frame containing different gene and complex annotations.
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
#' @importFrom logger log_fatal log_success
#' @export
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

        msg <- paste(
            'Downloading the entire annotations database is not allowed',
            'by default because of its huge size (>1GB). If you really',
            'want to do that, you find static files at',
            'https://archive.omnipathdb.org/.',
            'However we recommend to query a set of proteins or a few',
            'resources, depending on your interest.'
        )

        log_fatal(msg)
        stop(msg)

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

        log_success(
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
#' @importFrom magrittr %>% is_greater_than
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select pull group_split
#' @importFrom purrr map
#' @importFrom rlang set_names
#'
#' @seealso \code{\link{import_omnipath_annotations}}
pivot_annotations <- function(annotations){

    # NSE vs. R CMD check workaround
    record_id <- NULL

    more_than_one_resources <-
        annotations %>%
        pull(source) %>%
        unique %>%
        length %>%
        is_greater_than(1)

    if(more_than_one_resources){

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


#' Extra attribute names in an interaction data frame
#'
#' Interaction data frames might have an `extra_attrs` column if this field
#' has been requested in the query by passing the `fields = 'extra_attrs'
#' argument. This column contains resource specific attributes for the
#' interactions. The names of the attributes consist of the name of the
#' resource and the name of the attribute, separated by an underscore.
#' This function returns the names of the extra attributes available in
#' the provided data frame.
#'
#' @return Character: the names of the extra attributes in the data frame.
#'
#' @examples
#' i <- import_omnipath_interactions(fields = 'extra_attrs')
#' extra_attrs(i)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @export
extra_attrs <- function(data){

    data %>%
    {`if`(
        has_extra_attrs(data),
        pull(., extra_attrs) %>% map(names) %>% unlist() %>% unique(),
        character(0)
    )}

}


#' New columns from extra attributes
#'
#' @param data An interaction data frame.
#' @param ... The names of the extra attributes; NSE is supported.
#' @param flatten Logical: unnest the list column even if some records have
#'     multiple values for the attributes; these will yield multiple records
#'     in the resulted data frame.
#' @param keep_empty Logical: if `flatten` is `TRUE`, shall we keep the
#'     records which do not have the attribute?
#'
#' @return Data frame with the new column created; the new column is list
#'     type if one interaction might have multiple values of the attribute,
#'     or character type if
#'
#' @examples
#' i <- import_omnipath_interactions(fields = 'extra_attrs')
#' extra_attrs_to_cols(i, Cellinker_type, Macrophage_type)
#' extra_attrs_to_cols(
#'     i,
#'     Cellinker_type,
#'     Macrophage_type,
#'     flatten = TRUE,
#'     keep_empty = FALSE
#' )
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang enquos !!!
#' @importFrom purrr reduce map_chr
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{has_extra_attrs}}}
#'     \item{\code{\link{with_extra_attrs}}}
#' }
extra_attrs_to_cols <- function(
    data,
    ...,
    flatten = FALSE,
    keep_empty = TRUE
){

    if(!keep_empty){

        data %<>% with_extra_attrs(!!!enquos(...))
    }

    map_chr(enquos(...), .nse_ensure_str) %>%
    reduce(
        .extra_attr_to_col,
        .init = data,
        flatten = flatten
    )

}


#' New column from one extra attribute
#'
#' @param data An interaction data frame.
#' @param attr The name of an extra attribute; NSE is supported.
#' @param flatten Logical: unnest the list column even if some records have
#'     multiple values for the attributes; these will yield multiple records
#'     in the resulted data frame.
#' @param keep_empty Logical: if `flatten` is `TRUE`, shall we keep the
#'     records which do not have the attribute?
#'
#' @return Data frame with the new column created; the new column is list
#'     type if one interaction might have multiple values of the attribute,
#'     or character type if
#'
#' @importFrom magrittr %>% is_less_than
#' @importFrom rlang sym !! := enquo
#' @importFrom purrr map map_int pluck
#' @importFrom dplyr first mutate pull
#' @importFrom tidyr unnest
#' @noRd
.extra_attr_to_col <- function(
    data,
    attr,
    flatten = FALSE,
    keep_empty = TRUE
){

    # NSE vs. R CMD check workaround
    extra_attrs <- NULL

    attr_str <- .nse_ensure_str(!!enquo(attr))
    attr <- as.symbol(attr_str)

    data %>%
    {`if`(
        has_extra_attrs(data),
        mutate(
            .,
            !!attr := map(extra_attrs, pluck, attr_str)
        ) %>%
        {`if`(
            flatten,
            unnest(., !!attr, keep_empty = keep_empty),
            {`if`(
                pull(., !!attr) %>%
                    map_int(length) %>%
                    magrittr::is_less_than(2) %>%
                    all,
                mutate(
                    .,
                    !!attr := map(!!attr, first) %>%
                        null_to_na %>%
                        unlist
                ),
                .
            )}
        )},
        .
    )}

}


#' Tells if an interaction data frame has an extra_attrs column
#'
#' @param data An interaction data frame.
#'
#' @return Logical: TRUE if the data frame has the "extra_attrs" column.
#'
#' @examples
#' i <- import_omnipath_interactions(fields = 'extra_attrs')
#' has_extra_attrs(i)
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{extra_attrs_to_cols}}}
#'     \item{\code{\link{with_extra_attrs}}}
#' }
has_extra_attrs <- function(data){

    'extra_attrs' %in% colnames(data)

}


#' Interaction records having certain extra attributes
#'
#' @param data An interaction data frame.
#' @param ... The name(s) of the extra attributes; NSE is supported.
#'
#' @return The data frame filtered to the records having the extra attribute.
#'
#' @examples
#' i <- import_omnipath_interactions(fields = 'extra_attrs')
#' with_extra_attrs(i, Macrophage_type)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull
#' @importFrom purrr map_lgl map_chr
#' @importFrom rlang enquos
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{has_extra_attrs}}}
#'     \item{\code{\link{extra_attrs_to_cols}}}
#' }
with_extra_attrs <- function(data, ...){

    attrs_str <- map_chr(enquos(...), .nse_ensure_str)

    data %>%
    {`if`(
        has_extra_attrs(.),
        filter(
            .,
            pull(., 'extra_attrs') %>%
            map_lgl(
                function(x){
                    attrs_str %>%
                    intersect(names(x)) %>%
                    length %>%
                    as.logical
                }
            )
        ),
        .
    )}

}
