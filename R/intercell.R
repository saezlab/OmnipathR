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
