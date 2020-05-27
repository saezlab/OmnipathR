########## ########## ########## ##########
########## Generic non exported functions #
########## ########## ########## ##########

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
    enzsub = c('sources', 'references'),
    interactions = c('sources', 'references')
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
    'plasma_membrane_peripheral'
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

#' Downloads data from the OmniPath web service
#' Generic method for retrieval of a table and creating a data frame.
#' All methods specific for certain query types or datasets use this function
#' to manage the download.
#' Not exported.
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
        message(url)
        result <- omnipath_download(
            url,
            read.table,
            sep = '\t',
            header = TRUE,
            stringsAsFactors = FALSE,
            quote = ''
        )
        omnipath_check_result(result, url)
        if(!is.null(cache_file)){
            save(result, file = cache_file)
        }
        msg <- 'Downloaded %d %s.'
    }

    if(!silent){
        message(sprintf(msg, nrow(result), param$qt_message))
    }

    return(result)

}


#' Check the arguments of \link{import_omnipath}, corrects some easy to
#' confuse or deprecated synonyms and selects the message printed by
#' the download function.
#' Not exported.
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
        if(name %in% names(.omnipath_querystring_synonyms)){
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
                .omnipath_default_fields[[param$query_type]]
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


#' Constructs the URL by creating a base URL according to the query type and
#' adding all user or package defined query string parameters.
#' Not exported.
omnipath_url <- function(param){

    baseurl <- sprintf('http://omnipathdb.org/%s', param$query_type)

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


#' Checks whether the response is real data or an error message.
#' In case of error stops the execution and prints the URL and the message
#' from the server.
omnipath_check_result <- function(result, url){

    if(ncol(result) == 1){
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
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' ptms = import_omnipath_enzsub(
#'     resources = c('PhosphoSite', 'SIGNOR'),
#'     organism = 9606
#' )
#'
#' @seealso \code{\link{get_ptms_databases},
#'   \link{import_Omnipath_Interactions}}
#'
#' @aliases import_Omnipath_PTMS import_OmniPath_PTMS
import_omnipath_enzsub <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    fields = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'enzsub',
        cache_file = cache_file,
        resources = resources,
        organism = organism,
        fields = fields,
        ...
    )

    return(result)

}


# synonyms (old name)
import_Omnipath_PTMS <- import_omnipath_enzsub
import_OmniPath_PTMS <- import_omnipath_enzsub

#' Retrieve a list of enzyme-substrate resources available in OmniPath
#'
#' get the names of the enzyme-substrate relationship resources available
#' in \url{http://omnipath.org/enzsub}
#'
#' @return character vector with the names of the enzyme-substrate resources
#' @export
#' @importFrom utils read.table
#'
#' @examples
#' get_enzsub_resources()
#'
#' @seealso  \code{\link{get_resources},
#' \link{import_omnipath_enzsub}
#'
#' @aliases get_ptms_databases
get_enzsub_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'enzsub', dataset = dataset))

}

# synonym (old name)
get_ptms_databases <- get_enzsub_resources

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
#' @param cache_file path to an earlier data file
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions = import_omnipath_interactions(
#'     resources = c('SignaLink3'),
#'     organism = 9606
#' )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
#'
#' @aliases import_Omnipath_Interactions import_OmniPath_Interactions
import_omnipath_interactions <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    datasets = 'omnipath',
    fields = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        cache_file = cache_file,
        resources = resources,
        organism = organism,
        datasets = datasets,
        fields = fields,
        ...
    )

    return(result)

}

# synonyms (old name)
import_Omnipath_Interactions <- import_omnipath_interactions
import_OmniPath_Interactions <- import_omnipath_interactions


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
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose one of those: 9606 human (default), 10116 rat or 10090 Mouse
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <-
#'     import_pathwayextra_interactions(
#'         resources = c('BioGRID', 'IntAct'),
#'         organism = 9606
#'     )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
#'
#' @aliases import_PathwayExtra_Interactions
import_pathwayextra_interactions <- function(
    cache_file = NULL,
    organism = 9606,
    ...
){

    result <- import_omnipath_interactions(
        datasets = 'pathwayextra',
        cache_file = cache_file,
        organism = organism,
        ...
    )

    return(result)

}

# synonym (old name)
import_PathwayExtra_Interactions <- import_pathwayextra_interactions


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
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <-
#'    import_kinaseextra_interactions(
#'        resources = c('PhosphoPoint', 'PhosphoSite'),
#'        organism = 9606
#'    )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
#'
#' @aliases import_KinaseExtra_Interactions
import_kinaseextra_interactions <- function(
    cache_file = NULL,
    organism = 9606,
    ...
){

    result <- import_omnipath_interactions(
        datasets = 'kinaseextra',
        cache_file = cache_file,
        organism = organism,
        ...
    )

    return(result)

}

# synonym (old name)
import_KinaseExtra_Interactions <- import_kinaseextra_interactions


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
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <- import_ligrecextra_interactions(
#'     resources = c('HPRD', 'Guide2Pharma'),
#'     organism = 9606
#' )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
#'
#' @aliases import_LigrecExtra_Interactions
import_ligrecextra_interactions <- function(
    cache_file = NULL,
    organism = 9606,
    ...
){

    result <- import_omnipath_interactions(
        datasets = 'ligrecextra',
        cache_file = cache_file,
        organism = organism,
        ...
    )

    return(result)

}

# synonym (old name)
import_LigrecExtra_Interactions <- import_ligrecextra_interactions


#' Imports all post-translational interactions from OmniPath
#'
#' Imports the dataset from all post-translational datasets of OmniPath
#'
#' @return A dataframe containing post-translational interactions
#' @export
#' @importFrom utils read.table
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param exclude datasets to exclude
#'
#' @examples
#' interactions <-
#'     import_transcriptional_interactions(
#'         resources = c('PAZAR', 'ORegAnnO', 'DoRothEA')
#'     )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
import_post_translational_interactions <- function(
    organism = 9606,
    exclude = NULL,
    ...
){

    datasets <- c('omnipath', 'pathwayextra', 'kinaseextra', 'ligrecextra')
    datasets <- setdiff(datasets, exclude)


    result <- import_omnipath_interactions(
        organism = organism,
        datasets = datasets,
        ...
    )

    return(result)

}


#' Imports from Omnipath webservice the interactions from
#' Dorothea dataset
#'
#' Imports the dataset from:
#' \url{http://omnipathdb.org/interactions?datasets=tfregulons}
#' which contains transcription factor (TF)-target interactions from DoRothEA
#' \url{https://github.com/saezlab/DoRothEA}
#'
#' @return A dataframe containing TF-target interactions from DoRothEA
#' @export
#' @importFrom utils read.table
#' @param cache_file path to an earlier data file
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param confidence_level Vector detailing the confidence levels of the
#' interactions to be downloaded. In dorothea, every TF-target interaction
#' has a confidence score ranging from A to E, being A the most reliable
#' interactions.
#' By default we take A and B level interactions (\code{c(A, B)}).
#' It is to note that E interactions are not available in OmnipathR.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <- import_dorothea_interactions(
#'     resources = c('DoRothEA_A', 'ARACNe-GTEx_DoRothEA'),
#'     organism = 9606
#' )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
#'
#' @aliases import_TFregulons_Interactions import_tfregulons_interactions
import_dorothea_interactions <- function(
    cache_file = NULL,
    organism = 9606,
    confidence_levels = c('A', 'B'),
    ...
){

    result <- import_omnipath_interactions(
        datasets = 'dorothea',
        cache_file = cache_file,
        organism = organism,
        confidence_levels = confidence_levels,
        ...
    )

    return(result)

}


# synonym (old name)
import_TFregulons_Interactions <- import_dorothea_interactions
import_tfregulons_interactions <- import_dorothea_interactions


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
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <-
#'     import_tf_target_interactions(
#'         resources = c('miRTarBase', 'miRecords')
#'     )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
#'
#' @aliases import_miRNAtarget_Interactions
import_tf_target_interactions <- function(
    cache_file = NULL,
    organism = 9606,
    ...
){

    result <- import_omnipath_interactions(
        datasets = 'tf_target',
        cache_file = cache_file,
        organism = organism,
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
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <-
#'     import_transcriptional_interactions(
#'         resources = c('PAZAR', 'ORegAnnO', 'DoRothEA')
#'     )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
import_transcriptional_interactions <- function(
    organism = 9606,
    confidence_levels = c('A', 'B'),
    ...
){

    result <- rbind(
        import_dorothea_interactions(
            organism = organism,
            confidence_levels = confidence_levels,
            ...
        ),
        import_tf_target_interactions(organism = organism, ...)
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
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <-
#'     import_mirnatarget_interactions(
#'         resources = c('miRTarBase', 'miRecords')
#'     )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
#'
#' @aliases import_miRNAtarget_Interactions
import_mirnatarget_interactions <- function(
    cache_file = NULL,
    organism = 9606,
    ...
){

    result <- import_omnipath_interactions(
        datasets = 'mirnatarget',
        cache_file = cache_file,
        organism = organism,
        ...
    )

    return(result)

}

# synonym (old name)
import_miRNAtarget_Interactions <- import_mirnatarget_interactions


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
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <-
#'     import_tf_mirna_interactions(
#'         resources = c('TransmiR')
#'     )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
import_tf_mirna_interactions <- function(
    cache_file = NULL,
    organism = 9606,
    ...
){

    result <- import_omnipath_interactions(
        datasets = 'tf_mirna',
        cache_file = cache_file,
        organism = organism,
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
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#'
#' @examples
#' interactions <-
#'     import_lncrna_mrna_interactions(
#'         resources = c('ncRDeathDB')
#'     )
#'
#' @seealso \code{\link{get_interaction_databases},
#'   \link{import_all_interactions}}
import_lncrna_mrna_interactions <- function(
    cache_file = NULL,
    organism = 9606,
    ...
){

    result <- import_omnipath_interactions(
        datasets = 'lncrna_mrna',
        cache_file = cache_file,
        organism = organism,
        ...
    )

    return(result)

}


#' Imports all interaction datasets available in Omnipath
#'
#' Imports the dataset from:
#' \url{http://omnipathdb.org/interactions?datasets=omnipath, pathwayextra,
#' kinaseextra, ligrecextra, tfregulons, mirnatarget&fields=sources,
#' references&genesymbols=1},
#' which contains all the different interactions available in the webserver:
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
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param exclude datasets to exclude
#'
#' @examples
#' interactions <- import_all_interactions(
#'     resources = c('HPRD', 'BioGRID'),
#'     organism = 9606
#' )
#'
#' @seealso \code{\link{get_interaction_databases}}
#'
#' @aliases import_AllInteractions
import_all_interactions <- function(
    cache_file = NULL,
    resources = NULL,
    organism = 9606,
    dorothea_confidence_levels = c('A', 'B'),
    exclude = NULL,
    ...
){

    all_datasets <- jsonlite::fromJSON(
        txt = 'http://omnipathdb.org/queries/interactions?format=json'
    )$datasets

    all_datasets <- setdiff(all_datasets, exclude)

    result <- import_omnipath_interactions(
        cache_file = cache_file,
        datasets = all_datasets,
        organism = organism,
        confidence_levels = dorothea_confidence_levels,
        ...
    )

    return(result)

}

# synonym (old name)
import_AllInteractions <- import_all_interactions


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
#' get_interaction_databases()
#'
#' @seealso \code{\link{get_resources},
#' \link{import_all_interactions},
#' \link{import_omnipath_interactions}, \link{improt_pathwayextra_interactions},
#' \link{import_kinaseextra_interactions},
#' \link{import_ligrecextra_interactions},
#' \link{import_mirnatarget_interactions},
#' \link{import_dorothea_interactions}}
#'
#' @aliases get_interaction_databases
get_interaction_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'interactions', dataset = dataset))

}

# synonym (old name)
get_interaction_databases <- get_interaction_resources


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

# synonyms (old name)
import_Omnipath_complexes <- import_omnipath_complexes
import_OmniPath_complexes <- import_omnipath_complexes


#' Retrieve a list of complex resources available in Omnipath
#'
#' get the names of the resources from \url{http://omnipath.org/complexes}
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

    return(get_resources(query_type = 'complexes', dataset = dataset))

}

# synonym (old name)
get_complexes_databases <- get_complex_resources


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
#'
#' @examples
#' annotations = import_omnipath_annotations(
#'     select_genes = c('TP53', 'LMNA'),
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

        # account for the old argument name
        proteins <- c(proteins, list(...)$select_genes)

        if(!is.null(proteins)){
            result <- result[
                which(
                    result$uniprot %in% proteins ||
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

# synonyms (old name)
import_Omnipath_annotations <- import_omnipath_annotations
import_OmniPath_annotations <- import_omnipath_annotations

#' Get the resources available in the annotations database of OmniPath
#'
#' get the names of the resources from \url{http://omnipath.org/annotations}
#'
#' @return character vector with the names of the annotation resources
#' @export
#' @param dataset ignored for this query type
#'
#' @examples
#' get_annotation_resources()
#'
#' @seealso \code{\link{get_resources},
#' \link{import_omnipath_annotations}}
#'
#' @aliases get_annotation_databases
get_annotation_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'annotations', dataset = dataset))

}

# synonym (old name)
get_annotation_databases <- get_annotation_resources

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
#' @return A dataframe cotaining information about roles in inter-cellular
#' signaling.
#' @export
#' @importFrom utils read.csv
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
#'
#' @examples
#' intercell = import_omnipath_intercell(categories = c('ecm'))
#'
#' @seealso \code{\link{get_intercell_categories}}
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
    ...
){

    from_cache <- !is.null(cache_file) && file.exists(cache_file)
    args <- c(as.list(environment()), list(...))
    args$query_type <- 'intercell'

    result <- do.call(import_omnipath, args)

    if(from_cache){
        args$data <- result
        result <- do.call(filter_intercell, args)
    }

    return(result)

}

# synonyms (old name)
import_Omnipath_intercell <- import_omnipath_intercell
import_OmniPath_intercell <- import_omnipath_intercell


#' Retrieve a list of intercellular communication resources available in
#' Omnipath
#'
#' get the names of the databases from \url{http://omnipath.org/intercell}
#' @return character vector with the names of the databases
#' @export
#' @importFrom utils read.csv
#'
#' @examples
#' get_intercell_resources()
#'
#' @seealso \code{\link{get_resources},
#' \link{import_omnipath_intercell}}
get_intercell_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'intercell', dataset = dataset))

}


#' Imports an intercellular network combining annotations and interactions
#'
#' Imports an intercellular network by mapping intercellular annotations
#' and protein interactions. It first imports the PPI interactions from the
#' different datasets here described. Then, it takes proteins with the
#' intercellular roles defined by the user. Some proteins should be defined
#' as transmiters (eg. ligand) and other as receivers (receptor). We find the
#' interactions which source is a transmiter and its target a receiver.
#'
#' @return A dataframe containing information about protein-protein
#' interactions and the inter-cellular roles of the protiens involved in those
#' interactions.
#' @export
#' @importFrom utils read.csv
#' @param cache_file path to an earlier data file
#' @param resources vector containing interactions databases.
#' Interactions not reported in these databases are removed.
#' See \code{\link{get_interaction_databases}} for more information.
#' @param classes_source A list containing two vectors. The first one with
#' the main classes to be considered as transmiters and the second with the
#' main classes to be considered as receivers. For furter information
#' about the main classes see \code{\link{get_intercell_classes}}
#'
#' @examples
#' intercellNetwork <- import_intercell_network(
#' classes_source = list(transmiters = c('ligand'), receivers = c('receptor')))
#'
#' @seealso \code{\link{get_intercell_categories}}
import_intercell_network <- function(
    cache_file = NULL,
    resources = get_interaction_databases(),
    classes_source = list(transmiters = c('ligand'),
    receivers = c('receptor'))
){

    mainclass <- genesymbol <- NULL
    AllClasses <- unlist(classes_source)

    if (!all(AllClasses %in% get_intercell_classes())){
        stop('Some all the classes are not correct.
            Check get_intercell_classes()')
    }

    url_allinteractions_common <-
        paste0('http://omnipathdb.org/interactions?datasets=omnipath',
            ', pathwayextra,kinaseextra,ligrecextra',
            '&fields=sources,references')

    url_allinteractions <- organism_url(url_allinteractions_common, 9606)

    if(is.null(cache_file)){
        interactions <- omnipath_download(url_allinteractions, read.table, sep = '\t',
            header = TRUE, stringsAsFactors = FALSE)
        message('Downloaded ', nrow(interactions), ' interactions')
    } else {
        load(cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions, resources)

    intercellAnnotations <-
        import_Omnipath_intercell(select_classes = AllClasses)

    genesTransmiters <- intercellAnnotations %>%
        dplyr::filter(mainclass %in% classes_source$transmiters) %>%
        dplyr::distinct(genesymbol, mainclass)
    genesReceivers <- intercellAnnotations %>%
        dplyr::filter(mainclass %in% classes_source$receivers) %>%
        dplyr::distinct(genesymbol, mainclass)

    intercelNetwork <-
        dplyr::inner_join(filteredInteractions, genesTransmiters,
            by = c('source_genesymbol' = 'genesymbol')) %>%
        dplyr::rename(class_source = mainclass) %>%
        dplyr::inner_join(genesReceivers,
            by = c('target_genesymbol' = 'genesymbol')) %>%
        dplyr::rename(class_target = mainclass)

    return(intercelNetwork)
}

#' Get the different intercell categories described in Omnipath
#'
#' get the names of the categories from \url{http://omnipath.org/intercell}
#' @return character vector with the different intercell categories
#' @export
#' @importFrom utils read.csv
#' @examples
#' get_intercell_categories()
#' @seealso \code{\link{import_Omnipath_intercell},
#' \link{get_intercell_classes}}
get_intercell_categories <- function(){

    url_intercell <- 'http://omnipathdb.org/intercell_summary'
    intercell <- omnipath_download(url_intercell, read.csv, sep = '\t', header = TRUE,
        stringsAsFactors = FALSE)

    return(unique(intercell$category))
}

#' Get the different intercell main classes described in Omnipath
#'
#' get the names of the main classes from \url{http://omnipath.org/intercell}
#' @return character vector with the different intercell main classes
#' @export
#' @importFrom utils read.csv
#' @examples
#' get_intercell_classes()
#' @seealso \code{\link{import_Omnipath_intercell},
#' \link{get_intercell_categories}}
get_intercell_classes <- function(){

    url_intercell <- 'http://omnipathdb.org/intercell_summary'
    intercell <- omnipath_download(url_intercell, read.csv, sep = '\t', header = TRUE,
        stringsAsFactors = FALSE)

    return(unique(intercell$mainclass))
}

########## ########## ########## ##########
########## RESOURCE FILTERING      ########
########## ########## ########## ##########
## Non exported functions (package internal functions) to filter PTMs,
## interactions, complexes and annotations according to the databases passed
## to the main functions

## Filtering Interactions, PTMs and complexes
filter_by_resource <- function(data, resources = NULL){

    if(!is.null(resources)){

        before <- nrow(data)

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

## Filtering Annotations
filter_sources_annotations <- function(annotations, databases){
## takes annotations and removes those which are
## not reported by the given databases.

    if(is.null(databases)){
        return(annotations)
    }

    nAnnot <- nrow(annotations)
    subsetAnnotations <- dplyr::filter(annotations, source %in% databases)
    nAnnotPost <- nrow(subsetAnnotations)

    message(
        sprintf(
            'Removed %d annotations during database filtering.',
            nAnnot - nAnnotPost
        )
    )

    if (nAnnotPost > 0){
        return(subsetAnnotations)
    } else {
        return(NULL)
    }
}

## Filtering intercell records according to the categories and/or classes
## selected
filter_intercell <- function(intercell, categories, classes){
## takes intercell removes and removes those not reported by the given
## databases
    nIntercell = nrow(intercell)
    subsetIntercell <- dplyr::filter(intercell, .data$category %in% categories)
    subsetIntercell <-
        dplyr::filter(subsetIntercell, .data$mainclass %in% classes)

    nIntercellPost = nrow(subsetIntercell)

    message(
        sprintf(
            'Removed %d intercell records during category/class filtering.',
            nIntercell - nIntercellPost,
        )
    )

    if (nIntercellPost > 0){
        return(subsetIntercell)
    } else {
        return(NULL)
    }
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

########## ########## ########## ##########
########## Queries Format        ##########
########## ########## ########## ##########
## This function format de url for the queries to the Omnipath webserver
## according to the selected organism


########## ########## ########## ########## ##########
########## Format and filter of interactions #########
########## ########## ########## ########## ##########
## This function calls to the filtering functions and gives format
## to the data frames containing the interactions (For instance, it generates
## a new field with the number of sources/references reporting a given
## interaction)

filter_format_inter <- function(interaction_df, databases){

    if(!is.null(databases)){
        interaction_df <-
            filter_sources(interaction_df, databases = databases)
    } else {
        interaction_df <- interaction_df
    }

    if (nrow(interaction_df) ==0){
        stop('Try another database filtering: No records found')
    }

    if ('residue_offset' %in% colnames(interaction_df)){
        interaction_df$residue_offset <-
            as.character(as.numeric(interaction_df$residue_offset))
    }

    interaction_df$sources <- as.character(interaction_df$sources)
    interaction_df$nsources <-
        unlist(lapply(strsplit(interaction_df$sources, split = ';'), length))

    if ('references' %in% colnames(interaction_df)){
        interaction_df$references <- as.character(interaction_df$references)
        ## we remove references mentioned multiple times:
        interaction_df$references <-
            unlist(lapply(strsplit(interaction_df$references, split = ';'),
            function(x)paste(unique(x), collapse = ';')))
        interaction_df$nrefs <-
            unlist(lapply(strsplit(interaction_df$references, split = ';'),
            length))
    }
    return(interaction_df)
}
