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


#' Enzyme-substrate graph
#'
#' Transforms the a data frame with enzyme-substrate relationships
#' (obtained by \code{\link{import_omnipath_enzsub}}) to an igraph
#' graph object.
#'
#' @param enzsub Data frame created by \code{\link{import_omnipath_enzsub}}
#'
#' @return An igraph directed graph object.
#'
#' @examples
#' enzsub <- import_omnipath_enzsub(resources = c('PhosphoSite', 'SIGNOR'))
#' enzsub_g <- enzsub_graph(enzsub = enzsub)
#'
#' @export
#' @importFrom magrittr %>%
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_omnipath_enzsub}}}
#'     \item{\code{\link{giant_component}}}
#'     \item{\code{\link{find_all_paths}}}
#' }
#' @aliases ptms_graph
enzsub_graph <- function(enzsub){
    # This is a gene_name based conversion to igraph, i.e. the vertices are
    # identified by genenames, and not by uniprot IDs.
    # This might cause issue when a gene name encodes multiple uniprot IDs.

    # We check that the input dataframe contain the required columns.
    if(!all(
        c(
            'enzyme', 'substrate', 'enzyme_genesymbol',
            'substrate_genesymbol','sources'
        ) %in% colnames(enzsub)
        )
    ){
        stop('The input data frame does not contain the required columns')
    }

    enzsub %>%
    # We give the proper format to the edges
    # by calling to the function below
    format_graph_edges(flag = 'enzsub_dataset')

}


# Aliases (old names) to be deprecated
#' @rdname enzsub_graph
#' @param ... Passed to \code{enzsub_graph}.
#' @export
#' @importFrom rlang %||%
#'
#' @noRd
ptms_graph <- function(...){
    .Deprecated('enzsub_graph')
    args <- list(...)
    enzsub <- args$enzsub %||% args$ptms
    enzsub_graph(enzsub)
}


#' Build Omnipath interaction graph
#'
#' Transforms the interactions data frame to an igraph graph object.
#'
#' @param interactions data.frame created by \itemize{
#'     \item{\code{\link{import_omnipath_enzsub}}}
#'     \item{\code{\link{import_omnipath_interactions}}}
#'     \item{\code{\link{import_pathwayextra_interactions}}}
#'     \item{\code{\link{import_kinaseextra_interactions}}}
#'     \item{\code{\link{import_ligrecextra_interactions}}}
#'     \item{\code{\link{import_post_translational_interactions}}}
#'     \item{\code{\link{import_dorothea_interactions}}}
#'     \item{\code{\link{import_tf_target_interactions}}}
#'     \item{\code{\link{import_transcriptional_interactions}}}
#'     \item{\code{\link{import_mirnatarget_interactions}}}
#'     \item{\code{\link{import_all_interactions}}}}
#'
#' @return An igraph graph object.
#'
#' @export
#' @importFrom dplyr select group_by summarise ungroup
#' @importFrom magrittr %>%
#'
#' @examples
#' interactions <- import_omnipath_interactions(resources = c('SignaLink3'))
#' g <- interaction_graph(interactions)
#'
#' @seealso \itemize{
#'     \item{\code{\link{graph_interaction}}}
#'     \item{\code{\link{import_omnipath_interactions}}}
#'     \item{\code{\link{import_pathwayextra_interactions}}}
#'     \item{\code{\link{import_kinaseextra_interactions}}}
#'     \item{\code{\link{import_ligrecextra_interactions}}}
#'     \item{\code{\link{import_dorothea_interactions}}}
#'     \item{\code{\link{import_mirnatarget_interactions}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{giant_component}}}
#'     \item{\code{\link{find_all_paths}}}
#' }
interaction_graph <- function(interactions = interactions){
    # This is a gene_name based conversion to igraph, i.e. the vertices are
    # identified by genenames, and not by uniprot IDs.
    # This might cause issue when a gene name encodes multiple uniprot IDs.

    # We check that the input dataframe contain the required columns.
    if (!all(
        c(
            'source', 'target', 'source_genesymbol','target_genesymbol',
            'sources'
        ) %in% colnames(interactions)
        )
    ){
        stop('The input data frame does not contain the required columns')
    }

    # We give the proper format to the edges by calling to the function below
    format_graph_edges(interactions, flag = 'interactions_dataset')

}


#' Non exported function to give the proper format to the edges in order to
#' transform the interactions o PTMs data frame into a graph. The input
#' parameters are the data frame containing interactions or PTMs and a flag
#' indicating if we are dealing with the interactions dataset or the PTMs
#' dataset.
#'
#' @importFrom igraph graph_from_data_frame E ecount
#' @importFrom dplyr select
#' @importFrom rlang .data
#'
#' @noRd
format_graph_edges <- function(df_interact, flag){

    if(flag == 'enzsub_dataset'){
        # keep only edge attributes
        edges <- df_interact %>%
            select(-c(.data$enzyme, .data$substrate))

        # build vertices: gene_names and gene_uniprotIDs
        nodesA <-
            select(
                df_interact,
                c(.data$enzyme_genesymbol, .data$enzyme)
            )
        nodesB <-
            select(
                df_interact,
                c(.data$substrate_genesymbol, .data$substrate)
            )
    }else{
        if(flag == 'interactions_dataset'){
            # keep only edge attributes
            edges <- df_interact %>%
                select(- c(.data$source, .data$target))
            # build vertices: gene_names and gene_uniprotIDs
            nodesA <-
                select(df_interact,
                c(.data$source_genesymbol,.data$source))
            nodesB <-
                select(df_interact,
                c(.data$target_genesymbol,.data$target))
        }else{
            stop('Incorrect Input Flag')
        }
    }

    colnames(nodesA) <- colnames(nodesB) <- c('genesymbol', 'up_id')
    nodes <- rbind(nodesA,nodesB)
    nodes <- unique(nodes)
    nodes <- nodes %>% group_by(.data$genesymbol) %>%
        summarise('up_ids' = paste0(.data$up_id,collapse=',')) %>%
        ungroup()

    op_dfs <- list(edges = edges, nodes = nodes)
    directed <- TRUE
    op_g <- igraph::graph_from_data_frame(
        d = op_dfs$edges,
        directed = directed,
        vertices = op_dfs$nodes
    )

    has_edges <- ecount(op_g) > 0L

    if(has_edges) {

        igraph::E(op_g)$sources <- strsplit(igraph::E(op_g)$sources, ';')

        if('references' %in% colnames(df_interact)) {

            igraph::E(op_g)$references <- strsplit(
                igraph::E(op_g)$references, ';'
            )
        }

    }

    return(op_g)
}


#' All paths between two groups of vertices
#'
#' Finds all paths up to length `maxlen` between specified groups of
#' vertices. This function is needed only becaues igraph`s
#' `all_shortest_paths` finds only the shortest, not any
#' path up to a defined length.
#'
#' @usage
#' find_all_paths(
#'     graph,
#'     start,
#'     end,
#'     attr = NULL,
#'     mode = 'OUT',
#'     maxlen = 2,
#'     progress = TRUE
#' )
#'
#' @param graph An igraph graph object.
#' @param start Integer or character vector with the indices or names
#'     of one or more start vertices.
#' @param end Integer or character vector with the indices or names
#'     of one or more end vertices.
#' @param attr Character: name of the vertex attribute to identify the
#'     vertices by. Necessary if `start` and `end` are not igraph vertex ids
#'     but for example vertex names or labels.
#' @param mode Character: IN, OUT or ALL. Default is OUT.
#' @param maxlen Integer: maximum length of paths in steps, i.e. if
#'     maxlen = 3, then the longest path may consist of 3 edges and 4 nodes.
#' @param progress Logical: show a progress bar.
#'
#' @return List of vertex paths, each path is a character or integer vector.
#'
#' @importFrom igraph ego vertex_attr vertex_attr_names vcount
#' @importFrom purrr map cross2 map2 transpose
#' @importFrom magrittr %>% %<>%
#' @importFrom progress progress_bar
#' @importFrom rlang set_names
#' @export
#'
#' @examples
#' interactions <- import_omnipath_interactions()
#' graph <- interaction_graph(interactions)
#' paths <- find_all_paths(
#'     graph = graph,
#'     start = c('EGFR', 'STAT3'),
#'     end = c('AKT1', 'ULK1'),
#'     attr = 'name'
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{enzsub_graph}}}
#'     \item{\code{\link{giant_component}}}
#' }
find_all_paths <- function(
    graph,
    start,
    end,
    attr = NULL,
    mode = 'OUT',
    maxlen = 2,
    progress = TRUE
){

    find_all_paths_aux <- function(start, end, path = NULL){

        path %<>% append(start)

        if(start == end) return(list(path))

        paths <- list()

        if(length(path) <= maxlen){

            paths <- adjlist[[start]] %>%
                setdiff(path) %>%
                map(find_all_paths_aux, end = end, path = path) %>%
                unlist(recursive = FALSE)

        }

        return(paths)

    }

    adjlist <- graph %>% ego(mode = mode) %>% map(as.numeric)

    if(!is.null(attr)){

        if(!attr %in% vertex_attr_names(graph)){

            stop(sprintf('No such vertex attribute: `%s`.', attr))

        }

        attr_to_id <- graph %>%
            vertex_attr(attr) %>%
            set_names(graph %>% vcount %>% seq, .)
        start <- attr_to_id[start]
        end <- attr_to_id[end]

    }

    if(progress){
        pbar <- progress_bar$new(
            format = 'Finding paths [:bar] :percent eta: :eta',
            total = length(start * length(end))
        )
        fun <- {pbar$tick(); find_all_paths_aux}
    }else{
        fun <- find_all_paths_aux
    }

    paths <- cross2(start, end) %>%
        transpose() %>%
        c(fun) %>%
        do.call(map2, .) %>%
        unlist(recursive = FALSE)

    if(!is.null(attr)){

        paths %<>% map(
            function(path){vertex_attr(graph, attr)[path]}
        )

    }

    return(paths)

}


#' Giant component of a graph
#'
#' For an igraph graph object returns its giant component.
#'
#' @param graph An igraph graph object.
#'
#' @return An igraph graph object containing only the giant component.
#'
#' @examples
#' interactions <- import_post_translational_interactions()
#' graph <- interaction_graph(interactions)
#' graph_gc <- giant_component(graph)
#'
#' @importFrom magrittr %>%
#' @importFrom igraph components induced_subgraph
#' @export
giant_component <- function(graph){

    graph %>%
    components() %>%
    {induced_subgraph(graph, which(.$members == which.max(.$csize)))}

}


#' Extract a custom subnetwork from a large network
#'
#' @param network Either an OmniPath interaction data frame, or an igraph
#'     graph object.
#' @param nodes Character or integer vector: names, identifiers or indices
#'     of the nodes to build the subnetwork around.
#' @param order Integer: order of neighbourhood around nodes; i.e., number
#'     of steps starting from the provided nodes.
#' @param mode Character: "all", "out" or "in". Follow directed edges from
#'     the provided nodes in any, outbound or inbound direction, respectively.
#' @param mindist Integer: The minimum distance to include the vertex in the
#'     result.
#' @param return_df Logical: return an interaction data frame instead of an
#'     igraph object.
#'
#' @return A network data frame or an igraph object, depending on the
#'     ``return_df`` parameter.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom igraph ego induced_subgraph
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{graph_interaction}}}
#'     \item{\code{\link{show_network}}}
#' }
subnetwork <- function(
    network,
    nodes = NULL,
    order = 1L,
    mode = 'all',
    mindist = 0L,
    return_df = TRUE
) {

    # NSE vs R CMD check workaround
    igraph_vs <- NULL

    network %<>% ensure_igraph
    nodes %<>% igraph_vs(network)

    network %>%
    ego(nodes = nodes, mode = mode, order = order, mindist = mindist) %>%
    unlist %>%
    unique %>%
    induced_subgraph(network, vids = .) %>%
    {`if`(return_df, graph_interaction(.), .)}

}


#' Converts a network to igraph object unless it is already one
#'
#' @param network Either an OmniPath interaction data frame, or an igraph
#'     graph object.
#'
#' @return An igraph graph object.
#'
#' @importFrom magrittr %>%
#' @importFrom igraph is_igraph
#' @export
ensure_igraph <- function(network) {

    network %>% {`if`(is_igraph(.), ., interaction_graph(.))}

}


#' Interaction data frame from igraph graph object
#'
#' Convert an igraph graph object to interaction data frame. This is the
#' reverse of the operation done by thje \code{\link{interaction_graph}}
#' function. Networks can be easily converted to igraph objects, then
#' you can make use of all igaph methods, and at the end, get back the
#' interactions in a data frame, along with all new edge and node attributes.
#'
#' @param graph An igraph graph object created formerly from an OmniPath
#'     interactions data frame.
#' @param implode Logical: restore the original state of the list type
#'     columns by imploding them to character vectors, subitems separated
#'     by semicolons.
#'
#' @return An interaction data frame.
#'
#' @importFrom igraph as_data_frame
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate across rename inner_join select
#' @importFrom purrr map_chr
#' @importFrom tidyselect where everything
#' @export
#' @seealso \code{\link{interaction_graph}}
graph_interaction <- function(graph, implode = FALSE) {

    # NSE vs. R CMD check workaround
    from <- to <- source <- target <- up_ids_source <- up_ids_target <- NULL

    graph %>%
    as_data_frame(what = 'both') %>%
    {inner_join(
        inner_join(.$edges, .$vertices, by = c('from' = 'name')),
        .$vertices,
        by = c('to' = 'name'),
        suffix = c('_source', '_target')
    )} %>%
    as_tibble %>%
    rename(
        source = up_ids_source,
        target = up_ids_target,
        source_genesymbol = from,
        target_genesymbol = to
    ) %>%
    select(source, target, everything()) %>%
    {`if`(
        implode,
        mutate(., across(where(is.list), map_chr, paste, collapse = ';')),
        .
    )}

}
