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
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


LR_RESOURCES <- c(
    'CellPhoneDB',
    'Cellinker',
    'CellTalkDB',
    'CellChatDB',
    'CellCall',
    'connectomeDB2020',
    'Guide2Pharma',
    'Baccin2019',
    'Kirouac2010',
    'Ramilowski2015',
    'scConnect',
    'talklr',
    'ICELLNET',
    'EMBRACE',
    'LRdb',
    'iTALK',
    'SignaLink',
    'HPMR',
    'Wojtowicz2020'
)


#' Curated ligand-receptor interactions
#'
#' The OmniPath \emph{intercell} database annotates individual proteins and
#' complexes, and we combine these annotations with network interactions
#' on the client side, using \code{\link{import_intercell_network}}. The
#' architecture of this database is complex, aiming to cover a broad range
#' of knowledge on various levels of details and confidence. We can use the
#' \code{\link{intercell_consensus_filter}} and
#' \code{\link{filter_intercell_network}} functions for automated, data
#' driven quality filtering, in order to enrich the cell-cell communication
#' network in higher confidence interactions. However, for many users, a
#' simple combination of the most established, expert curated ligand-receptor
#' resources, provided by this function, fits better their purpose.
#'
#' @param curated_resources Character vector of the resource names which
#'     are considered to be expert curated. You can include any
#'     post-translational network resource here, but if you include non
#'     ligand-receptor or non curated resources, the result will not
#'     fulfill the original intention of this function.
#' @param cellphonedb Logical: include the curated interactions from
#'     \emph{CellPhoneDB} (not the whole \emph{CellPhoneDB} but a subset
#'     of it).
#' @param cellinker Logical: include the curated interactions from
#'     \emph{Cellinker} (not the whole \emph{Cellinker} but a subset of it).
#' @param talklr Logical: include the curated interactions from
#'     \emph{talklr} (not the whole \emph{talklr} but a subset of it).
#' @param signalink Logical: include the ligand-receptor interactions
#'     from \emph{SignaLink.} These are all expert curated.
#' @param ... Passed to \code{\link{import_post_translational_interactions}}:
#'     further parameters for the interaction data. Should not contain
#'     `resources` argument as that would interfere with the downstream calls.
#'
#' @details
#' Some resources are a mixture of curated and bulk imported interactions,
#' and sometimes it's not trivial to separate these, we take care of these
#' here. This function does not use the \emph{intercell} database of
#' OmniPath, but retrieves and filters a handful of network resources. The
#' returned data frame has the layout of \emph{interactions} (network) data
#' frames, and the \emph{source} and \emph{target} partners implicitly
#' correspond to \emph{ligand} and \emph{receptor.} The data frame shows
#' all resources and references for all interactions, but each interaction
#' is supported by at least one ligand-receptor resource which is supposed
#' to based on expert curation in a ligand-receptor context.
#'
#' @return A data frame similar to \emph{interactions} (network) data frames,
#'     the \emph{source} and \emph{target} partners being ligand and
#'     receptor, respectively.
#'
#' @examples
#' lr <- curated_ligand_receptor_interactions()
#' lr
#'
#' @importFrom magrittr %>% %<>% equals
#' @importFrom dplyr filter select bind_rows distinct across
#' @importFrom purrr reduce discard
#' @importFrom tidyselect everything
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_intercell_network}}}
#'     \item{\code{\link{filter_intercell_network}}}
#'     \item{\code{\link{annotated_network}}}
#'     \item{\code{\link{import_post_translational_interactions}}}
#'     \item{\code{\link{import_ligrecextra_interactions}}}
#'     \item{\code{\link{curated_ligrec_stats}}}
#' }
curated_ligand_receptor_interactions <- function(
    curated_resources = c(
        'Guide2Pharma', 'HPMR', 'ICELLNET', 'Kirouac2010',
        'CellTalkDB', 'CellChatDB', 'connectomeDB2020'
    ),
    cellphonedb = TRUE,
    cellinker = TRUE,
    talklr = TRUE,
    signalink = TRUE,
    ...
){

    cellchatdb <- 'CellChatDB' %in% curated_resources

    curated_resources %<>% discard(equals, 'CellChatDB')

    curated_resources %>%
    {`if`(
        length(.) > 0L,
        import_post_translational_interactions(resources = ., ...) %>%
        with_references(resources = curated_resources),
        NULL
    )} %>%
    {reduce(
        c('cellphonedb', 'cellinker', 'talklr'),
        function(d, resource){
            `if`(
                get(resource),
                bind_rows(d, get(sprintf('%s_curated', resource))(...)),
                d
            )
        },
        .init = .
    )} %>%
    {`if`(
        cellchatdb,
        bind_rows(., cellchatdb_ligrec(...)),
        .
    )} %>%
    {`if`(
        signalink,
        bind_rows(., signalink_ligrec(...)),
        .
    )} %>%
    distinct(across(everything()))

}


#' Curated interactions from CellPhoneDB
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom rlang exec !!!
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select
#'
#' @noRd
cellphonedb_curated <- function(...){

    # NSE vs. R CMD check workaround
    sources <- CellPhoneDB_type <- extra_attrs <- NULL

    args <- list(...)
    keep_eattrs <- 'extra_attrs' %in% names(args)
    args$fields %<>% c('extra_attrs') %>% unique

    exec(
        import_post_translational_interactions,
        resources = 'CellPhoneDB',
        !!!args
    ) %>%
    extra_attrs_to_cols(CellPhoneDB_type) %>%
    filter(
        # these are the interactions which are curated by
        # the CellPhoneDB team (labelled as "curated" in
        # cellphonedb/src/core/data/interaction_input.csv)
        !str_detect(sources, '_CellPhoneDB') &
        CellPhoneDB_type == 'ligand-receptor'
    ) %>%
    select(-CellPhoneDB_type) %>%
    {`if`(keep_eattrs, ., select(., -extra_attrs))}

}


#' Ligand-receptor interactions from CellChatDB
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom rlang exec !!!
#' @importFrom dplyr filter select
#' @noRd
cellchatdb_ligrec <- function(curated = TRUE, ...){

    # NSE vs. R CMD check workaround
    CellChatDB_category <- extra_attrs <- NULL

    args <- list(...)
    keep_eattrs <- 'extra_attrs' %in% names(args)
    args$fields %<>% c('extra_attrs') %>% unique

    exec(
        import_post_translational_interactions,
        # fully curated, ligand-receptor only resources
        resources = 'CellChatDB',
        !!!args
    ) %>%
    extra_attrs_to_cols(CellChatDB_category) %>%
    filter(
        CellChatDB_category %in% c(
            NA,
            'Cell-Cell Contact',
            'Secreted Signaling'
        )
    ) %>%
    select(-CellChatDB_category) %>%
    {`if`(keep_eattrs, ., select(., -extra_attrs))} %>%
    {`if`(curated, with_references(., resources = 'CellChatDB'), .)}

}


#' Curated interactions from Cellinker
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom rlang exec !!!
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select
#'
#' @noRd
cellinker_curated <- function(...){

    # NSE vs. R CMD check workaround
    sources <- Cellinker_type <- extra_attrs <- NULL

    lr_types <- c(
        'Cytokine-cytokine receptor interaction',
        'Secreted protein to receptor interaction'
    )

    args <- list(...)
    keep_eattrs <- 'extra_attrs' %in% names(args)
    args$fields %<>% c('extra_attrs') %>% unique

    exec(
        import_post_translational_interactions,
        resources = 'Cellinker',
        !!!args
    ) %>%
    extra_attrs_to_cols(Cellinker_type) %>%
    filter(
        # these are the interactions which are curated by
        # the CellPhoneDB team (labelled as "curated" in
        # cellphonedb/src/core/data/interaction_input.csv)
        !str_detect(sources, '_Cellinker') &
        Cellinker_type %in% lr_types
    ) %>%
    select(-Cellinker_type) %>%
    {`if`(keep_eattrs, ., select(., -extra_attrs))} %>%
    with_references(resources = 'Cellinker')

}


#' Curated interactions from talklr
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom rlang exec !!!
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select
#'
#' @noRd
talklr_curated <- function(...){

    # NSE vs. R CMD check workaround
    references <- extra_attrs <- NULL

    args <- list(...)
    keep_eattrs <- 'extra_attrs' %in% names(args)
    args$fields %<>% c('extra_attrs') %>% unique

    exec(
        import_post_translational_interactions,
        resources = 'talklr',
        !!!args
    ) %>%
    filter_extra_attrs(talklr_putative = 0) %>%
    filter(str_detect(references, 'talklr:')) %>%
    {`if`(keep_eattrs, ., select(., -extra_attrs))}

}



#' Ligand-receptor interactions from SignaLink
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select
#'
#' @noRd
signalink_ligrec <- function(...){

    # NSE vs. R CMD check workaround
    function_source <- function_target <- NULL

    annotated_network(
        network = 'SignaLink3',
        annot = 'SignaLink_function',
        network_args = list(...)
    ) %>%
    filter(
        # these interactions are directly expert curated in a
        # ligand-receptor context, and contain some well established
        # ligand inputs of the most important canonical pathways
        function_source == 'Ligand' &
        function_target == 'Receptor'
    ) %>%
    select(
        -function_source,
        -function_target
    )

}


#' Statistics about literature curated ligand-receptor interactions
#'
#' @param ... Passed to \code{\link{curated_ligand_receptor_interactions}},
#'     determines the set of all curated L-R interactions which will be
#'     compared against each of the individual resources.
#'
#' @return A data frame with estimated counts of curated ligand-receptor
#'     interactions for each L-R resource.
#'
#' @details
#' The data frame contains the total number of interactions, the number
#' of interactions which overlap with the set of curated interactions
#' \emph{(curated_overlap),} the number of interactions with literature
#' references from the given resource \emph{(literature)} and the number
#' of interactions which are curated by the given resource
#' \emph{(curated_self)}. This latter we defined according to our best
#' knowledge, in many cases it's not possible to distinguish curated
#' interactions). All these numbers are also presented as a percent of
#' the total. Importantly, here we consider interactions curated only if
#' they've been curated in a cell-cell communication context.
#'
#' @examples
#' clr <- curated_ligrec_stats()
#' clr
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom rlang set_names
#' @importFrom tidyr unnest_wider
#' @importFrom tibble tibble
#' @export
#' @seealso \code{\link{curated_ligand_receptor_interactions}}
curated_ligrec_stats <- function(...){

    # NSE vs. R CMD check workaround
    resource <- data <- NULL

    curated <- curated_ligand_receptor_interactions(...)

    LR_RESOURCES %>%
    set_names(
        map(., stats_one_resource, curated),
        .
    ) %>%
    tibble(
        resource = names(.),
        data = .
    ) %>%
    unnest_wider(col = data)

}


#' Statistics about curated interactions in one resource
#'
#' @importFrom magrittr %>% inset2
#' @importFrom purrr map
#' @importFrom rlang set_names
#' @noRd
stats_one_resource <- function(resource, curated){

    c('total', 'curated_overlap', 'literature', 'curated_self') %>%
    set_names(
        map(
            .,
            function(field){
                get(sprintf('%s_one_resource', field))(resource, curated)
            }
        ),
        .
    ) %>%
    map(nrow) %>%
    inset2(
        'curated_overlap_pct',
        .$curated_overlap / .$total * 100L
    ) %>%
    inset2(
        'literature_pct',
        .$literature / .$total * 100L
    ) %>%
    inset2(
        'curated_self_pct',
        .$curated_self / .$total * 100L
    )

}


#' Overlap of a resource with a curated
#'
#' @param resource Character: name of a ligand-receptor resource.
#' @param curated Data frame with all the curated L-R interactions, as
#'     provided by \code{\link{curated_ligand_receptor_interactions}}.
#'
#' @return Data frame with the interactions from the resource of interest
#'     in overlap with the curated set
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr inner_join select distinct
#' @noRd
curated_overlap_one_resource <- function(resource, curated){

    # NSE vs. R CMD check workaround
    source <- target <- NULL

    total_one_resource(resource) %>%
    inner_join(curated, by = c('source', 'target')) %>%
    select(source, target) %>%
    distinct()

}


#' Curated interactions from one resource
#'
#' @param resource Character: name of a single resource
#'
#' @noRd
curated_self_one_resource <- function(resource, ...){

    # some resources have their dedicated function:
    suffixes <- c('curated', 'ligrec')

    for(suf in suffixes){

        func_name <- sprintf('%s_%s', tolower(resource), suf)

        func <- tryCatch(
            get(func_name, envir = asNamespace('OmnipathR')),
            error = function(cond){NULL}
        )

        if(!is.null(func)){

            return(func())

        }

    }

    # at the end we fall back to the generic method:
    # interactions with literature references from the resource
    literature_one_resource(resource = resource)

}


#' All interactions from one resource
#'
#' @importFrom magrittr %<>%
#' @noRd
total_one_resource <- function(resource, ...){

    resource %<>% network_resource_name

    import_post_translational_interactions(resources = resource)

}


#' Interactions with literature references from one resource
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @noRd
literature_one_resource <- function(resource, ...){

    # NSE vs. R CMD check workaround
    references <- NULL

    resource %<>% network_resource_name

    import_post_translational_interactions(resources = resource) %>%
    with_references(resources = resource)

}


#' Network resource name
#'
#' @importFrom magrittr %>%
#' @noRd
network_resource_name <- function(resource){

    resource %>%
    {`if`(. == 'SignaLink', 'SignaLink3', .)}

}
