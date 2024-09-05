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


TOPOLOGIES <-
    c(
        'secreted',
        'plasma_membrane_transmembrane',
        'plasma_membrane_peripheral'
    )

TOPOLOGIES_SHORT <-
    rlang::set_names(TOPOLOGIES, c('sec', 'pmtm', 'pmp'))


#' Cell-cell communication roles from OmniPath
#'
#' Roles of proteins in intercellular communication from the
#' \url{https://omnipathdb.org/intercell} endpoint of the OmniPath web service.
#' It provides information on the roles in inter-cellular signaling. E.g. if
#' a protein is a ligand, a receptor, an extracellular matrix (ECM) component,
#' etc.
#'
#' @return A data frame of intercellular communication roles.
#'
#' @param categories vector containing the categories to be retrieved.
#'     All the genes belonging to those categories will be returned. For
#'     further information about the categories see
#'     \code{\link{get_intercell_categories}}.
#' @param parent vector containing the parent classes to be retrieved.
#'     All the genes belonging to those classes will be returned. For
#'     furter information about the main classes see
#'     \code{\link{get_intercell_categories}}.
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
#'     \code{consensus_percentile} for major localizations. For example, with
#'     a value of 50, the secreted, plasma membrane transmembrane or
#'     peripheral attributes will be true only where at least 50 percent
#'     of the resources support these.
#' @param ... Further arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -query_type -datasets -types -genesymbols
#'     -json_param -references_by_resource -add_counts
#'
#' @examples
#' ecm_proteins <- intercell(categories = "ecm")
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom purrr reduce
#' @importFrom rlang exec !!!
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{intercell_network}}}
#'     \item{\code{\link{intercell_consensus_filter}}}
#'     \item{\code{\link{filter_intercell}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell_resources}}}
#'     \item{\code{\link{intercell_summary}}}
#'     \item{\code{\link{intercell_network}}}
#' }
#' @aliases import_omnipath_intercell
intercell <- function(
    categories = NULL,
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

    topology %<>%
        topology_long %>%
        {reduce(
            TOPOLOGIES,
            function(topos, topo){
                arg_value <- get(topo)
                `if`(
                    is.null(arg_value),
                    topos,
                    `if`(
                        arg_value,
                        union(topos, topo),
                        setdiff(topos, topo)
                    )
                )
            },
            .init = .
        )}

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
        exec(omnipath_query, !!!args) %>%
        intercell_consensus_filter(
            percentile = consensus_percentile,
            loc_percentile = loc_consensus_percentile,
            topology = topology
        )

    return(result)

}


# Aliases (old names) to be Deprecated
#' @rdname complexes
#' @param ... Passed to \code{intercell}.
#' @export
#'
#' @noRd
import_omnipath_intercell <- function(...){
    .Deprecated('import_omnipath_intercell')
    intercell(...)
}


#' Quality filter for intercell annotations
#'
#' @param data A data frame with intercell annotations, as provided by
#'     \code{\link{intercell}}.
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
#' @param topology Character vector: list of allowed topologies, possible
#'     values are *"secreted"*, *"plasma_membrane_peripheral"* and
#'     *"plasma_membrane_transmembrane"*.
#'
#' @return The data frame in \code{data} filtered by the consensus scores.
#'
#' @examples
#' ligand_receptor <- intercell(parent = c("ligand", "receptor"))
#' nrow(ligand_receptor)
#' # [1] 50174
#' lr_q50 <- intercell_consensus_filter(ligand_receptor, 50)
#' nrow(lr_q50)
#' # [1] 42863
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr group_by filter ungroup bind_rows
#' @importFrom dplyr select distinct inner_join pull
#' @importFrom stats quantile
#' @importFrom rlang !! := sym eval_tidy parse_expr
#' @importFrom purrr reduce
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{resources}}}
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{filter_intercell}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell_resources}}}
#'     \item{\code{\link{intercell_summary}}}
#'     \item{\code{\link{intercell_network}}}
#' }
intercell_consensus_filter <- function(
    data,
    percentile = NULL,
    loc_percentile = NULL,
    topology = NULL
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

        locations <- intercell(
            aspect = 'locational',
            parent = TOPOLOGIES,
            consensus_percentile = loc_percentile
        )

        data %<>%
        {reduce(
            TOPOLOGIES,
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
        )} %>%
        {`if`(
            is.null(topology),
            .,
            filter(
                .,
                eval_tidy(parse_expr(paste(topology, collapse = ' | ')))
            )
        )}

    }

    return(data)

}


#' Retrieves a list of intercellular communication resources available in
#' OmniPath
#'
#' Retrieves a list of the databases from
#' \url{https://omnipathdb.org/intercell}.
#'
#' @param dataset ignored at this query type
#'
#' @return character vector with the names of the databases
#'
#' @examples
#' intercell_resources()
#'
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{resources}}}
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{filter_intercell}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell_summary}}}
#'     \item{\code{\link{intercell_network}}}
#' }
#' @aliases get_intercell_resources
intercell_resources <- function(dataset = NULL){

    return(resources(query_type = 'intercell', datasets = dataset))

}


# Aliases (old names) to be Deprecated
#' @rdname complexes
#' @param ... Passed to \code{intercell_resources}.
#' @export
#'
#' @noRd
get_intercell_resources <- function(...){
    .Deprecated('get_intercell_resources')
    intercell_resources(...)
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
#' the respective methods (\code{\link{omnipath_interactions}} for
#' the network and \code{\link{intercell}} for the
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
#' @param interactions_param a list with arguments for an interactions query;
#'     \code{\link{omnipath-interactions}}.
#' @param transmitter_param a list with arguments for
#'     \code{\link{intercell}}, to define the transmitter side
#'     of intercellular connections
#' @param receiver_param a list with arguments for
#'     \code{\link{intercell}}, to define the receiver side
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
#' intercell_network <- intercell_network(
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
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{intercell_summary}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{omnipath}}}
#'     \item{\code{\link{pathwayextra}}}
#'     \item{\code{\link{kinaseextra}}}
#'     \item{\code{\link{ligrecextra}}}
#'     \item{\code{\link{unique_intercell_network}}}
#'     \item{\code{\link{simplify_intercell_network}}}
#'     \item{\code{\link{filter_intercell_network}}}
#' }
#' @aliases import_intercell_network
intercell_network <- function(
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

    .slow_doctest()

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
        omnipath_query,
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
        do.call(intercell, transmitter_param) %>%
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
        do.call(intercell, receiver_param) %>%
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


# Aliases (old names) to be Deprecated
#' @rdname complexes
#' @param ... Passed to \code{intercell_network}.
#' @export
#'
#' @noRd
import_intercell_network <- function(...){
    .Deprecated('import_intercell_network')
    intercell_network(...)
}


#' Filter intercell annotations
#'
#' Filters a data frame retrieved by \code{\link{intercell}}.
#'
#' @param data An intercell annotation data frame as provided by
#'     \code{\link{intercell}}.
#' @param categories Character: allow only these values in the \code{category}
#'     column.
#' @param resources Character: allow records only from these resources.
#' @param parent Character: filter for records with these parent categories.
#' @param scope Character: filter for records with these annotation scopes.
#'     Possible values are \code{generic} and \code{specific}.
#' @param aspect Character: filter for records with these annotation aspects.
#'     Possible values are \code{functional} and \code{locational}.
#' @param source Character: filter for records with these annotation sources.
#'     Possible values are \code{composite} and \code{resource_specific}.
#' @param transmitter Logical: if \code{TRUE} only transmitters, if
#'     \code{FALSE} only non-transmitters will be selected, if \code{NULL}
#'     it has no effect.
#' @param receiver Logical: works the same way as \code{transmitters}.
#' @param secreted Logical: works the same way as \code{transmitters}.
#' @param plasma_membrane_peripheral Logical: works the same way as
#'     \code{transmitters}.
#' @param plasma_membrane_transmembrane Logical: works the same way as
#'     \code{transmitters}.
#' @param proteins Character: filter for annotations of these proteins.
#'     Gene symbols or UniProt IDs can be used.
#' @param causality Character: filter for records with these causal roles.
#'     Possible values are \code{transmitter} and \code{receiver}. The filter
#'     applied simultaneously to the \code{transmitter} and \code{receiver}
#'     arguments, it's just a different notation for the same thing.
#' @param topology Character: filter for records with these localization
#'     topologies. Possible values are \code{secreced},
#'     \code{plasma_membrane_peripheral} and
#'     \code{plasma_membrane_transmembrane}; the shorter notations \code{sec},
#'     \code{pmp} and \code{pmtm} can be used. Has the same effect as the
#'     logical type arguments, just uses a different notation.
#' @param ... Ignored.
#'
#' @return The intercell annotation data frame filtered according to the
#'     specified conditions.
#'
#' @examples
#' ic <- intercell()
#' ic <- filter_intercell(
#'     ic,
#'     transmitter = TRUE,
#'     secreted = TRUE,
#'     scope = "specific"
#' )
#'
#' @importFrom dplyr recode rename_all
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang .data set_names
#' @importFrom stringr str_detect
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell_summary}}}
#'     \item{\code{\link{intercell_network}}}
#' }
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

    args <- environment() %>% as.list()
    args$cats <- args$categories

    before <- nrow(data)

    topology %<>%
        {data.frame(set_names(
            as.list(.),
            rep(TRUE, length(.))
        ))} %>%
        rename_all(
            recode,
            secreted = 'sec',
            plasma_membrane_peripheral = 'pmp',
            plasma_membrane_transmembrane = 'pmtm'
        )

    causality %<>%
        {`if`('both' %in% ., c('transmitter', 'receiver'), .)} %>%
        {data.frame(set_names(
            as.list(.),
            rep(TRUE, length(.))
        ))} %>%
        rename_all(
            recode,
            transmitter = 'trans',
            receiver = 'rec'
        )

    data <-
        data %>%
        filter(
            (is.null(args$cats) | category %in% args$cats) &
            (is.null(args$parent) | .data$parent %in% args$parent) &
            (is.null(args$scope) | .data$scope %in% args$scope) &
            (is.null(args$aspect) | .data$aspect %in% aspect) &
            (is.null(args$source) | .data$source %in% source) &
            (is.null(args$transmitter) | .data$transmitter) &
            (is.null(args$receiver) | .data$receiver) &
            (is.null(args$secreted) | .data$secreted) &
            (
                is.null(args$plasma_membrane_peripheral) |
                .data$plasma_membrane_peripheral
            ) &
            (
                is.null(args$plasma_membrane_transmembrane) |
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
            'Removed %d and kept %d records of intercell data.',
            before - after,
            after
        )
    )

    return(data)

}


#' Quality filter an intercell network
#'
#' The intercell database  of OmniPath covers a very broad range of possible
#' ways of cell to cell communication, and the pieces of information, such as
#' localization, topology, function and interaction, are combined from many,
#' often independent sources. This unavoidably result some weird and
#' unexpected combinations which are false positives in the context of
#' intercellular communication. \code{\link{intercell_network}}
#' provides a shortcut (\code{high_confidence}) to do basic quality filtering.
#' For custom filtering or experimentation with the parameters we offer this
#' function.
#'
#' @param network An intercell network data frame, as provided by
#'     \code{\link{intercell_network}}, without \code{simplify}.
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
#' icn <- intercell_network()
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
#'     \item{\code{\link{intercell_network}}}
#'     \item{\code{\link{unique_intercell_network}}}
#'     \item{\code{\link{simplify_intercell_network}}}
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell_summary}}}
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

    .slow_doctest()

    # NSE vs. R CMD check workaround
    parent <- curation_effort <- n_resources <- n_references <- NULL

    consensus <-
        intercell(
            consensus_percentile = consensus_percentile,
            loc_consensus_percentile = loc_consensus_percentile
        )

    consensus_annot <-
        consensus %>%
        select(uniprot, parent) %>%
        distinct

    consensus_loc <-
        consensus %>%
        select(uniprot, !!!syms(TOPOLOGIES)) %>%
        distinct

    topologies <- unlist(TOPOLOGIES_SHORT)

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
        topology_long %>%
        walk(check_topo) %>%
        intersect(topologies) %>%
        sprintf('%s_intercell_source', .) %>%
        paste(collapse = ' | ')

    receiver_topology %<>%
        topology_long %>%
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
            -(!!!exprs(sprintf('%s_intercell_source', TOPOLOGIES))),
            -(!!!exprs(sprintf('%s_intercell_target', TOPOLOGIES)))
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
#' \code{\link{intercell_network}}, are combinations of a network data
#' frame with two copies of the intercell annotation data frames, all of them
#' already having quite some columns. Here we keep only the names of the
#' interacting pair, their intercellular communication roles, and the minimal
#' information of the origin of both the interaction and the annotations.
#' Optionally further columns can be selected.
#'
#' @param network An intercell network data frame, as provided by
#'     \code{\link{intercell_network}}.
#' @param ... Optional, further columns to select.
#'
#' @return An intercell network data frame with some columns removed.
#'
#' @examples
#' icn <- intercell_network()
#' icn_s <- simplify_intercell_network(icn)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom rlang ensyms !!!
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{intercell_network}}}
#'     \item{\code{\link{filter_intercell_network}}}
#'     \item{\code{\link{unique_intercell_network}}}
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell_summary}}}
#' }
simplify_intercell_network <- function(network, ...){

    .slow_doctest()

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
#' \link{intercell_network}}, by default each pair of annotations for
#' an interaction is represented in a separate row. This function drops the
#' annotations and keeps only the distinct interacting pairs.
#'
#' @param network An intercellular network data frame as produced by
#'     \code{\link{intercell_network}}.
#' @param ... Additional columns to keep. Note: if these have multiple
#'     values for an interacting pair, only the first row will be
#'     preserved.
#'
#' @return A data frame with interacting pairs and interaction attributes.
#'
#' @examples
#' icn <- intercell_network()
#' icn_unique <- unique_intercell_network(icn)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct
#' @importFrom rlang !!!
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{intercell_network}}}
#'     \item{\code{\link{simplify_intercell_network}}}
#'     \item{\code{\link{filter_intercell_network}}}
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell_summary}}}
#' }
unique_intercell_network <- function(network, ...){

    .slow_doctest()

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
#' Retrieves a list of categories from \url{https://omnipathdb.org/intercell}.
#'
#' @return character vector with the different intercell categories
#' @export
#'
#' @examples
#' intercell_categories()
#'
#' @seealso \itemize{
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{intercell_generic_categories}}}
#'     \item{\code{\link{intercell_summary}}}
#' }
#' @aliases get_intercell_categories
intercell_categories <- function(){

    return(
        unique(
            omnipath_query('intercell_summary', license = NA)$category
        )
    )

}


# Aliases (old names) to be Deprecated
#' @rdname complexes
#' @param ... Passed to \code{intercell_categories}.
#' @export
#'
#' @noRd
get_intercell_categories <- function(...){
    .Deprecated('get_intercell_categories')
    intercell_categories(...)
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
intercell_summary <- function(){

    omnipath_query('intercell_summary', license = NA)

}

#' Retrieves a list of the generic categories in the intercell database
#' of OmniPath
#'
#' Retrieves a list of the generic categories from
#' \url{https://omnipathdb.org/intercell}.
#'
#' @return character vector with the different intercell main classes
#' @export
#'
#' @examples
#' intercell_generic_categories()
#'
#' @seealso \itemize{
#'     \item{\code{\link{intercell}}}
#'     \item{\code{\link{intercell_categories}}}
#'     \item{\code{\link{intercell_summary}}}
#' }
#' @aliases get_intercell_generic_categories
intercell_generic_categories <- function(){

    return(
        unique(
            omnipath_query('intercell_summary', license = NA)$parent
        )
    )
}


# Aliases (old names) to be Deprecated
#' @rdname complexes
#' @param ... Passed to \code{intercell_generic_categories}.
#' @export
#'
#' @noRd
get_intercell_generic_categories <- function(...){
    .Deprecated('get_intercell_generic_categories')
    intercell_generic_categories(...)
}


#' Short to long topology names
#'
#' @param topologies Character vector of short topology names. Long names and
#'     any other strings will be left intact.
#'
#' @importFrom rlang !!!
#' @importFrom dplyr recode
#' @importFrom magrittr %>%
#' @noRd
topology_long <- function(topologies){

    topologies %>%
    {`if`(is.null(.), ., recode(., !!!TOPOLOGIES_SHORT))}

}
