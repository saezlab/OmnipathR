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

ORGANISMS_SUPPORTED <- c(9606L, 10090L, 10116L)

.omnipath_qt_synonyms <- list(
    ptms = 'enzsub',
    enz_sub = 'enzsub',
    complex = 'complexes',
    interaction = 'interactions',
    network = 'interactions'
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
    'strict_evidences',
    'genesymbol_resource',
    'orthology_targets',
    'keep_evidences',
    'cache'
)


#' Download data from the OmniPath web service
#'
#' This is the most generic method for accessing data from the OmniPath web
#' service. All other functions retrieving data from OmniPath call this
#' function with various parameters. In general, every query can retrieve data
#' in tabular or JSON format, the tabular (data frame) being the default.
#'
#' @param query_type Character: "interactions", "enzsub", "complexes",
#'     "annotations", or "intercell".
#' @param organism Character or integer: name or NCBI Taxonomy ID of the
#'     organism. OmniPath is built of human data, and the web service provides
#'     orthology translated interactions and enzyme-substrate relationships for
#'     mouse and rat. For other organisms and query types, orthology
#'     translation will be called automatically on the downloaded human data
#'     before returning the result.
#' @param resources Character vector: name of one or more resources. Restrict
#'     the data to these resources. For a complete list of available resources,
#'     call the `get_<query_type>_resources` functions for the query type of
#'     interst.
#' @param datasets Character vector: name of one or more datasets. In the
#'     interactions query type a number of datasets are available. The default
#'     is caled "omnipath", and corresponds to the curated causal signaling
#'     network published in the 2016 OmniPath paper.
#' @param genesymbols Character or logical: TRUE or FALS or "yes" or "no".
#'     Include the `genesymbols` column in the results. OmniPath uses UniProt
#'     IDs as the primary identifiers, gene symbols are optional.
#' @param fields Character vector: additional fields to include in the result.
#'     For a list of available fields, call `query_info("interactions")`.
#' @param default_fields Logical: if TRUE, the default fields will be included.
#' @param silent Logical: if TRUE, no messages will be printed. By default a
#'     summary message is printed upon successful download.
#' @param logicals Character vector: fields to be cast to logical.
#' @param format Character: if "json", JSON will be retrieved and processed
#'     into a nested list; any other value will return data frame.
#' @param download_args List: parameters to pass to the download function,
#'     which is `readr::read_tsv` by default, and `jsonlite::safe_load`.
#' @param references_by_resource Logical: if TRUE,, in the `references`
#'     column the PubMed IDs will be prefixed with the names of the resources
#'     they are coming from. If FALSE, the `references` column will be a list
#'     of unique PubMed IDs.
#' @param add_counts Logical: if TRUE, the number of references and number of
#'     resources for each record will be added to the result.
#' @param license Character: license restrictions. By default, data from
#'     resources allowing "academic" use is returned by OmniPath. If you use
#'     the data for work in a company, you can provide "commercial" or
#'     "for-profit", which will restrict the data to those records which are
#'     supported by resources that allow for-profit use.
#' @param password Character: password for the OmniPath web service. You can
#'     provide a special password here which enables the use of `license =
#'     "ignore"` option, completely bypassing the license filter.
#' @param exclude Character vector: resource or dataset names to be excluded.
#'     The data will be filtered after download to remove records of the
#'     excluded datasets and resources.
#' @param json_param List: parameters to pass to the `jsonlite::fromJSON` when
#'     processing JSON columns embedded in the downloaded data. Such columns
#'     are "extra_attrs" and "evidences". These are optional columns which
#'     provide a lot of extra details about interactions.
#' @param strict_evidences Logical: reconstruct the "sources" and "references"
#'     columns of interaction data frames based on the "evidences" column,
#'     strictly filtering them to the queried datasets and resources. Without
#'     this, the "sources" and "references" fields for each record might
#'     contain information for datasets and resources other than the queried
#'     ones, because the downloaded records are a result of a simple filtering
#'     of an already integrated data frame.
#' @param genesymbol_resource Character: "uniprot" (default) or "ensembl". The
#'     OmniPath web service uses the primary gene symbols as provided by
#'     UniProt. By passing "ensembl" here, the UniProt gene symbols will be
#'     replaced by the ones used in Ensembl. This translation results in a loss
#'     of a few records, and multiplication of another few records due to
#'     ambiguous translation.
#' @param cache Logical: use caching, load data from and save to the. The cache
#'     directory by default belongs to the user, located in the user's default
#'     cache directory, and named "OmnipathR". Find out about it by
#'     \code{\link{omnipath_get_cachedir}}. Can be changed by
#'     \code{\link{omnipath_set_cachedir}}.
#' @param ... Additional parameters for the OmniPath web service. These
#'     parameters will be processed, validated and included in the query
#'     string. Many parameters are already explicitly set by the arguments
#'     above. A number of query type specific parameters are also available,
#'     learn more about these by the \code{\link{query_info}} function. For
#'     functions more specific than \code{\link{omnipath_query}}, arguments for
#'     all downstream functions are also passed here.
#'
#' @return Data frame (tibble) or list: the data returned by the OmniPath web
#'     service (or loaded from cache), after processing. Nested list if the
#'     "format" parameter is "json", otherwise a tibble.
#'
#' @examples
#' interaction_data <- omnipath_query("interaction", datasets = "omnipath")
#' interaction_data
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom tibble as_tibble
#' @importFrom readr read_tsv cols col_character
#' @importFrom utils modifyList
#' @importFrom rlang !!!
#' @export
omnipath_query <- function(
    query_type,
    organism = 9606L,
    resources = NULL,
    datasets = NULL,
    genesymbols = 'yes',
    fields = NULL,
    default_fields = TRUE,
    silent = FALSE,
    logicals = NULL,
    download_args = list(),
    format = 'data.frame',
    references_by_resource = TRUE,
    add_counts = TRUE,
    license = NULL,
    password = NULL,
    exclude = NULL,
    json_param = list(),
    strict_evidences = FALSE,
    genesymbol_resource = 'UniProt',
    cache = NULL,
    ...
){

    datasets %<>% setdiff(exclude)
    resources %<>% setdiff(exclude)
    cache %<>% use_cache

    param <-
        environment() %>%
        as.list %>%
        c(list(...)) %>%
        omnipath_check_param

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
        url = url,
        cache = cache
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
    download_args %<>%
        modifyList(
            `if`(
                !is_empty_2(param$format) && param$format == 'json',
                json_defaults,
                dataframe_defaults
            ),
            .
        ) %>%
        modifyList(download_args_defaults, .)

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
            resources = param$resources,
            .keep = param$keep_evidences
        )
    }

    if(param$query_type %in% c('interactions', 'enzsub') && add_counts){
        result %<>% count_references
        result %<>% count_resources
    }

    if(is.data.frame(result)){
        result %<>% as_tibble
    }

    if(!is_empty_2(param$orthology_targets)){
        result %<>% omnipath_orthology_translate(param)
    } else if(
        !is_empty_2(param$genesymbol_resource) &&
        str_to_lower(param$genesymbol_resource[1L]) != 'uniprot'
    ) {
        result %<>% update_genesymbols(param, param$organisms[1L])
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


#' Consolidate arguments for OmniPath calls
#'
#' Most importantly, control argument expansion by ensuring the priority of
#' the explicitely built in synonyms.
#'
#' @param dots List: dots from the parent.
#' @param ... Override arguments.
#'
#' @return A list of arguments for `import_omnipath`.
#'
#' @importFrom magrittr %>% inset2 extract
#' @importFrom stringr str_extract
#' @importFrom purrr map
#' @noRd
omnipath_args <- function(dots, ...) {

    calling_env <- parent.frame()
    override <- list(...) %>% qs_synonyms
    calling_fun <- sys.function(-1L)

    defaults <-
        calling_fun %>%
        formals %>%
        qs_synonyms

    calling_env %>%
    as.list %>%
    qs_synonyms %>%
    modifyList(defaults, .) %>%
    modifyList(dots %>% qs_synonyms) %>%
    modifyList(override) %>%
    inset2('...', NULL)

}


#' Checks the arguments of \link{import_omnipath}, corrects some easy to
#' confuse or deprecated synonyms and selects the message printed by
#' the download function.
#' Not exported.
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom logger log_warn
#' @importFrom purrr map_int
#' @importFrom dplyr first
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
    param %<>% qs_synonyms

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

    param$keep_evidences <- 'evidences' %in% param$fields

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

    organisms_supported <- `if`(
        param$query_type %in% c('interactions', 'enzsub'),
        ORGANISMS_SUPPORTED,
        9606L
    )
    # allow organism names
    param$organisms %<>% map_int(ncbi_taxid)
    param$orthology_targets <- param$organisms %>% setdiff(organisms_supported)
    param$organisms %<>% intersect(organisms_supported)
    log_trace('Organism(s): %s', enum_format(param$organisms))
    log_trace('Orthology targets: %s', enum_format(param$orthology_targets))

    if(length(param$orthology_targets > 0L)) {

        param$organisms %<>% if_null_len0(9606L)

        if(length(param$organisms) > 1L || param$organisms != 9606L) {

            log_warn(
                paste0(
                    'The result will be translated to `%s` by orthology; ',
                    'querying for human.'
                ),
                param$orthology_targets %>% first
            )
            param$organisms <- 9606L

        }

    }

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


#' Replace synonymous keys in query string parameters
#'
#' @importFrom magrittr %>% %<>% inset2
#' @noRd
qs_synonyms <- function(param) {

    # mapping the query string parameter synonyms
    for(name in names(param)){
        if(
            name %in% names(.omnipath_querystring_synonyms) &&
            !.omnipath_querystring_synonyms[[name]] %in% names(param)
        ){
            new_name <- .omnipath_querystring_synonyms[[name]]
            param %<>%
                inset2(new_name, param[[name]]) %>%
                inset2(name, NULL)
        }
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
#' @importFrom purrr reduce
#' @importFrom utils URLencode
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

    url <-
        .omnipath_querystring_param %>%
        reduce(
            function(url, key){
                omnipath_url_add_param(url, key, param[[key]])
            },
            .init = baseurl
        ) %>%
        URLencode

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


#' Translate an OmniPath data frame by orthology
#'
#' @param data A data frame from OmniPath, or any data frame with UniProt IDs
#'     in a column named "uniprot", or columns named "source" and "target", or
#'     "enzyme" and "substrate".
#' @param param List or character: OmniPath query parameters or the name of a
#'
#' @return Data frame: the input data frame with the UniProt IDs translated to
#'     the target organism by orthology.
#'
#' @importFrom magrittr %T>% %>% extract2
#' @importFrom dplyr first
#' @importFrom rlang !! sym
#' @importFrom purrr reduce
#' @importFrom logger log_warn log_info
#' @noRd
omnipath_orthology_translate <- function(data, param) {

    target <-
        param %>%
        {`if`(is.list(.), extract2(., 'orthology_targets'), .)} %T>%
        {`if`(
            !is.null(.) && length(.) > 1L,
            log_warn(
                paste0(
                       'Orthology translation works for only one organism in one ',
                       'query. Translating only to the first target organism: `%i`.'
                ),
                first(.)
            ),
            NA
        )} %>%
        first

    if(!is.na(target)) {

        log_info(
            'Translating to `%s` by orthology, %i records before translation',
            target,
            nrow(data)
        )

        genesymbol_resource <-
            param %>%
            {`if`(is.list(.), extract2(., 'genesymbol_resource'), 'uniprot')}

        data <-
            reduce(
                uniprot_columns(data),
                ~orthology_translate_column(
                    .x,
                    !!sym(.y),
                    source_organism = 9606L,
                    target_organism = target,
                    replace = TRUE
                ),
                .init = data
            ) %T>%
            {log_info('%i records after orthology translation.', nrow(.))} %>%
            update_genesymbols(genesymbol_resource, organism = target)

    }

    return(data)

}


#' Update gene symbols in the data frame
#'
#' @param data A data frame from OmniPath, or any data frame with UniProt IDs
#'     in a column named "uniprot", or columns named "source" and "target", or
#'     "enzyme" and "substrate".
#' @param param List or character: OmniPath query parameters or the name of a
#'     Gene Symbol resource, either "uniprot" or "ensembl".
#' @param organism Character or integer: name or NCBI Taxonomy ID of the
#'     organism the UniProt IDs in the data frame belong to.
#'
#' @return Data frame: the input data frame with the gene symbol columns
#'     updated.
#'
#' @importFrom magrittr %>% extract2 is_in
#' @importFrom logger log_warn log_info
#' @importFrom purrr reduce
#' @importFrom stringr str_to_lower
#' @importFrom rlang !! !!! sym :=
#' @noRd
update_genesymbols <- function(data, param, organism = 9606L) {

    up_cols <- data %>% uniprot_columns
    gs_cols <- data %>% genesymbol_columns

    if(is_empty_2(up_cols)) {

        warn <-
            paste0(
               'No columns with UniProt IDs found, ',
               'not updating gene symbols.'
            )
        log_warn(warn)
        warning(warn)
        return(data)

    }

    resource <-
        param %>%
        {`if`(is.list(.), extract2(., 'genesymbol_resource'), .)} %>%
        str_to_lower

    if(resource %>% is_in(c('ensembl', 'uniprot')) %>% not) {

        log_warn(
            paste0(
               'Unknown genesymbol resource: `%s`. ',
               'Using UniProt, the default one, instead.'
            ),
            resource
        )
        resource <- 'uniprot'

    }

    resource %T>%
    {log_info('Setting gene symbols from `%s`.', .)} %>%
    {reduce(
        up_cols,
        ~translate_ids(
            .x,
            !!sym(.y) := uniprot,
            !!sym(
                `if`(
                    .y == 'uniprot',
                    'genesymbol',
                    sprintf('%s_genesymbol', .y)
                )
            ) := genesymbol,
            organism = organism,
            ensembl = resource == 'ensembl',
            keep_untranslated = FALSE
        ),
        .init = data
    )} %T>%
    {log_info('%i records after setting gene symbols.', nrow(.))}

}

#' Names of UniProt columns in an OmniPath data frame
#'
#' @param data A data frame from OmniPath.
#'
#' @importFrom magrittr %>% is_in
#' @noRd
uniprot_columns <- function(data) {

    UP_COLS <- list(
        'uniprot',
        c('enzyme', 'substrate'),
        c('source', 'target')
    )

    for(cols in UP_COLS) {

        if(cols %>% is_in(data %>% colnames) %>% all) {
            return(cols)
        }

    }

}


#' Names of Gene Symbol columns in an OmniPath data frame
#'
#' @param data A data frame from OmniPath.
#'
#' @importFrom magrittr %>%
#' @noRd
genesymbol_columns <- function(data) {

    data %>%
    uniprot_columns %>%
    {`if`(
        is.null(.) || length(.) == 0L,
        .,
        `if`(
             .[1L] == 'uniprot',
             'genesymbol',
             map_chr(., ~sprintf('%s_genesymbol', .))
        )
    )} %>%
    intersect(data %>% colnames)

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
#' @param cache Logical: use the cache.
#' @param ... Passed to the internal function \code{download_base} and
#'     from there ultimately to \code{fun}.
#'
#' @importFrom logger log_trace log_info log_error
#' @noRd
omnipath_download <- function(url, fun, cache = NULL, ...) {

    cache %<>% use_cache

    if(cache) {

        for(the_url in url) {

            log_trace('Looking up in cache: `%s`', the_url)
            from_cache <- omnipath_cache_load(url = the_url)

            if(!is.null(from_cache)){

                log_info('Loaded from cache: `%s`', the_url)
                attr(from_cache, 'url') <- the_url
                return(from_cache)

            }

        }

    }

    for(the_url in url) {

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
            if(cache) {
                omnipath_cache_save(data = result, url = the_url)
            }
            attr(result, 'url') <- the_url
            return(result)

        }

    }

}
