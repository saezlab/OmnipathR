#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
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

#' Visualize node neighborhood with SigmaJS
#'
#' This function takes an OmniPath interaction data frame as input and
#' returns a sigmaJS object for the subgraph formed by the neighbors of a node
#' of interest.
#'
#' @param interactions An OmniPath interaction data frame.
#' @param node The node of interest.
#'
#' @return A sigmaJS object, check http://sigmajs.john-coene.com/index.html
#'     for further details and customization options.
#'
#' @export
#'
#' @examples
#' interactions <- tibble::tibble(
#'     source = c('P00533', 'P04626', 'P00533'),
#'     target = c('P04626', 'P62993', 'P62993'),
#'     source_genesymbol = c('EGFR', 'ERBB2', 'EGFR'),
#'     target_genesymbol = c('ERBB2', 'GRB2', 'GRB2'),
#'     sources = c('SIGNOR;SPIKE', 'SIGNOR', 'SPIKE'),
#'     consensus_stimulation = c(1L, 1L, 0L),
#'     curation_effort = c(4L, 2L, 1L)
#' )
#' if (requireNamespace('sigmajs', quietly = TRUE)) {
#'     network_plot <- show_network(interactions, node = 'EGFR')
#'     # open the plot in the viewer only in an interactive session
#'     if (interactive()) network_plot
#' }
#'
#' @importFrom igraph induced_subgraph neighborhood
#' @importFrom magrittr %>%
#' @export
show_network <- function(interactions, node = NULL) {

    # import the OmniPath network as an igraph object
    g <- interaction_graph(interactions)

    # if node is not in the network, return NULL and message
    if(!is.null(node) && !node %in% V(g)$name) {

        if(!node %in% V(g)$name) {

            warning(sprintf('Node `%s` not found in the network.', node))
            return(NULL)

        }

        # get neighborhood graph
        g <- induced_subgraph(g, neighborhood(g, nodes = node)[[1]])

    }

    V(g)$label <- V(g)$name
    E(g)$color <- ifelse(E(g)$consensus_stimulation == 1, 'red', 'blue')
    E(g)$size <- E(g)$curation_effort

    # create sigmaJS
    outp <-
        sigmajs::sigmajs() %>%
        sigmajs::sg_from_igraph(g) %>%
        sigmajs::sg_settings(drawLabels = TRUE, drawEdgeLabels = TRUE) %>%
        sigmajs::sg_drag_nodes()

    return(outp)

}
