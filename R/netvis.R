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
#' \dontrun{
#' # get interactions from omnipath
#' interactions <- omnipath()
#' # create and plot the network containing ATM neighbors
#' viz_sigmajs_neighborhood(interactions_df = interactions, int_node = "ATM")
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
