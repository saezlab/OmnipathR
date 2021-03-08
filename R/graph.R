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


#' Post-translational modifications (PTMs) graph
#'
#' transforms the ptms interactions data.frame to igraph object
#'
#' @return An igraph object
#' @export
#' @import igraph
#' @importFrom dplyr %>% select group_by summarise ungroup
#' @importFrom rlang .data
#' @param ptms data.frame created by \code{\link{import_omnipath_enzsub}}
#' @examples
#' \donttest{
#' ptms = import_omnipath_enzsub(resources=c('PhosphoSite', 'SIGNOR'))
#' ptms_g = ptms_graph(ptms = ptms)
#' }
#' @seealso  \code{\link{import_omnipath_enzsub}}
ptms_graph <- function(ptms){
    # This is a gene_name based conversion to igraph, i.e. the vertices are
    # identified by genenames, and not by uniprot IDs.
    # This might cause issue when a gene name encodes multiple uniprot IDs.

    # We check that the input dataframe contain the required columns.
    if (!all(c('enzyme', 'substrate', 'enzyme_genesymbol',
        'substrate_genesymbol','sources') %in% colnames(ptms))) {
    stop('The input data frame does not contain the required columns')
    }

    # We give the proper format to the edges by calling to the function below
    output_graph <- format_graph_edges(ptms,flag = 'ptms_dataset')
    return(output_graph)
}


#' Build Omnipath interaction graph
#'
#' Transforms the interactions data frame to an igraph graph object.
#'
#' @param interactions data.frame created by
#'     \code{\link{import_omnipath_interactions}},
#'     \code{\link{import_pathwayextra_interactions}},
#'     \code{\link{import_kinaseextra_interactions}},
#'     \code{\link{import_ligrecextra_interactions}},
#'     \code{\link{import_dorothea_interactions}},
#'     \code{\link{import_mirnatarget_interactions}} or
#'     \code{\link{import_all_interactions}}
#'
#' @return An igraph object
#'
#' @export
#' @import igraph
#' @importFrom dplyr select group_by summarise ungroup
#' @importFrom magrittr %>%
#'
#' @examples
#' \donttest{
#' interactions = import_omnipath_interactions(resources=c('SignaLink3'))
#' OPI_g = interaction_graph(interactions)
#' }
#'
#' @seealso \code{\link{import_omnipath_interactions}},
#'     \code{\link{import_pathwayextra_interactions}},
#'     \code{\link{import_kinaseextra_interactions}},
#'     \code{\link{import_ligrecextra_interactions}},
#'     \code{\link{import_dorothea_interactions}},
#'     \code{\link{import_mirnatarget_interactions}} or
#'     \code{\link{import_all_interactions}}
interaction_graph <- function(interactions = interactions){
    # This is a gene_name based conversion to igraph, i.e. the vertices are
    # identified by genenames, and not by uniprot IDs.
    # This might cause issue when a gene name encodes multiple uniprot IDs.

    # We check that the input dataframe contain the required columns.
    if (!all(c('source', 'target', 'source_genesymbol','target_genesymbol',
        'sources') %in% colnames(interactions))) {
        stop('The input data frame does not contain the required columns')
    }

    # We give the proper format to the edges by calling to the function below
    output_graph <-
        format_graph_edges(interactions,flag = 'interactions_dataset')
    return(output_graph)
}


#' Non exported function to give the proper format to the edges in order to
#' transform the interactions o PTMs data frame into a graph. The input
#' parameters are the data frame containing interactions or PTMs and a flag
#' indicating if we are dealing with the interactions dataset or the PTMs
#' dataset.
#'
#' @importFrom dplyr select
#'
#' @noRd
format_graph_edges <- function(df_interact, flag){

    if(flag == 'ptms_dataset'){
        # keep only edge attributes
        edges <- df_interact %>%
            dplyr::select(- c(.data$enzyme, .data$substrate))

        # build vertices: gene_names and gene_uniprotIDs
        nodesA <-
            dplyr::select(df_interact, c(.data$enzyme_genesymbol, .data$enzyme))
        nodesB <-
            dplyr::select(df_interact, c(.data$substrate_genesymbol,
            .data$substrate))
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

    colnames(nodesA) = colnames(nodesB) = c('genesymbol', 'up_id')
    nodes <- rbind(nodesA,nodesB)
    nodes <- unique(nodes)
    nodes <- nodes %>% dplyr::group_by(.data$genesymbol) %>%
        dplyr::summarise('up_ids' = paste0(.data$up_id,collapse=',')) %>%
        dplyr::ungroup()

    op_dfs <- list(edges = edges, nodes = nodes)
    directed <- TRUE
    op_g <- igraph::graph_from_data_frame(d = op_dfs$edges,directed = directed,
        vertices = op_dfs$nodes)

    igraph::E(op_g)$sources    <- strsplit(igraph::E(op_g)$sources,    ';')
    if('references' %in% colnames(df_interact)){
        igraph::E(op_g)$references <- strsplit(igraph::E(op_g)$references, ';')
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
#' @param graph graph An igraph graph object.
#' @param start integer Numeric or character vector with the indices or names
#' of one or more start vertices.
#' @param end integer Numeric or character vector with the indices or names
#' of one or more end vertices.
#' @param attr character Name of the vertex attribute to identify the
#' vertices by. Necessary if `start` and `end` are not igraph vertex ids
#' but for example vertex names or labels.
#' @param mode character IN, OUT or ALL. Default is OUT.
#' @param maxlen integer Maximum length of paths in steps, i.e. if maxlen = 3,
#' then the longest path may consist of 3 edges and 4 nodes.
#' @param progress logical Show a progress bar. Default is FALSE.
#'
#' @return List of vertex paths
#'
#' @importFrom igraph ego vertex_attr vertex_attr_names
#' @importFrom purrr map cross2 map2 transpose
#' @importFrom magrittr %>% %<>%
#' @importFrom progress progress_bar
#' @export
#'
#' @examples
#' \donttest{
#' interactions <- import_omnipath_interactions()
#' graph <- interaction_graph(interaction)
#' paths <- find_all_paths(
#'     c('EGFR', 'STAT3'),
#'     c('AKT1', 'ULK1'),
#'     attr = 'name'
#' )
#' }
#'
#' @seealso \code{\link{interaction_graph}, \link{ptms_graph}}
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