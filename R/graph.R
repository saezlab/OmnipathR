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
#' @seealso \code{\link{import_omnipath_enzsub}}
#' @aliases enzsub_graph
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
#'     \item{\code{\link{import_post_translationsl_interactions}}}
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
#'     \item{\code{\link{import_omnipath_interactions}}}
#'     \item{\code{\link{import_pathwayextra_interactions}}}
#'     \item{\code{\link{import_kinaseextra_interactions}}}
#'     \item{\code{\link{import_ligrecextra_interactions}}}
#'     \item{\code{\link{import_dorothea_interactions}}}
#'     \item{\code{\link{import_mirnatarget_interactions}}}
#'     \item{\code{\link{import_all_interactions}}}
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
#' @importFrom igraph graph_from_data_frame E
#' @importFrom dplyr select
#'
#' @noRd
format_graph_edges <- function(df_interact, flag){

    if(flag == 'enzsub_dataset'){
        # keep only edge attributes
        edges <- df_interact %>%
            dplyr::select(-c(.data$enzyme, .data$substrate))

        # build vertices: gene_names and gene_uniprotIDs
        nodesA <-
            dplyr::select(
                df_interact,
                c(.data$enzyme_genesymbol, .data$enzyme)
            )
        nodesB <-
            dplyr::select(
                df_interact,
                c(.data$substrate_genesymbol, .data$substrate)
            )
    }else{
        if(flag == 'interactions_dataset'){
            # keep only edge attributes
            edges <- df_interact %>%
                dplyr::select(- c(.data$source, .data$target))
            # build vertices: gene_names and gene_uniprotIDs
            nodesA <-
                dplyr::select(df_interact,
                c(.data$source_genesymbol,.data$source))
            nodesB <-
                dplyr::select(df_interact,
                c(.data$target_genesymbol,.data$target))
        }else{
            stop('Incorrect Input Flag')
        }
    }

    colnames(nodesA) <- colnames(nodesB) <- c('genesymbol', 'up_id')
    nodes <- rbind(nodesA,nodesB)
    nodes <- unique(nodes)
    nodes <- nodes %>% dplyr::group_by(.data$genesymbol) %>%
        dplyr::summarise('up_ids' = paste0(.data$up_id,collapse=',')) %>%
        dplyr::ungroup()

    op_dfs <- list(edges = edges, nodes = nodes)
    directed <- TRUE
    op_g <- igraph::graph_from_data_frame(
        d = op_dfs$edges,
        directed = directed,
        vertices = op_dfs$nodes
    )

    igraph::E(op_g)$sources <- strsplit(igraph::E(op_g)$sources, ';')

    if('references' %in% colnames(df_interact)){
        igraph::E(op_g)$references <- strsplit(
            igraph::E(op_g)$references, ';'
        )
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
#' @param progress Logical: show a progress bar. Default is FALSE.
#'
#' @return List of vertex paths, each path is a character or integer vector.
#'
#' @importFrom igraph ego vertex_attr vertex_attr_names vcount
#' @importFrom purrr map cross2 map2 transpose
#' @importFrom magrittr %>% %<>%
#' @importFrom progress progress_bar
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
            setNames(graph %>% vcount %>% seq, .)
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