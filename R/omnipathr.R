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


#' The OmnipathR package
#'
#' @description
#' OmnipathR is an R package built to provide easy access to the data stored
#' in the OmniPath web service:
#'
#' \url{https://omnipathdb.org/}
#'
#' And a number of other resources, such as BioPlex, ConsensusPathDB, EVEX,
#' Guide to Pharmacology (IUPHAR/BPS), Harmonizome, HTRIdb, InWeb InBioMap,
#' KEGG Pathway, Pathway Commons, Ramilowski et al. 2015, RegNetwork, ReMap,
#' TF census, TRRUST and Vinayagam et al. 2011.
#'
#' The OmniPath web service implements a very simple REST style API. This
#' package make requests by the HTTP protocol to retreive the data. Hence,
#' fast Internet access is required for a propser use of OmnipathR.
#'
#' The package also provides some utility functions to filter, analyse and
#' visualize the data. Furthermore, OmnipathR features a close integration
#' with the NicheNet method for ligand activity prediction from
#' transcriptomics data, and its R implementation nichenetr (available in
#' CRAN).
#'
#' @examples
#' \dontrun{
#' # Download post-translational modifications:
#' enzsub <- enzyme_substrate(resources = c("PhosphoSite", "SIGNOR"))
#'
#' # Download protein-protein interactions
#' interactions <- omnipath(resources = "SignaLink3")
#'
#' # Convert to igraph objects:
#' enzsub_g <- enzsub_graph(enzsub = enzsub)
#' OPI_g <- interaction_graph(interactions = interactions)
#'
#' # Print some interactions:
#' print_interactions(head(enzsub))
#'
#' # interactions with references:
#' print_interactions(tail(enzsub), writeRefs = TRUE)
#'
#' # find interactions between kinase and substrate:
#' print_interactions(dplyr::filter(ptms,enzyme_genesymbol=="MAP2K1",
#'    substrate_genesymbol=="MAPK3"))
#'
#' # find shortest paths on the directed network between proteins
#' print_path_es(shortest_paths(OPI_g, from = "TYRO3", to = "STAT3",
#'    output = 'epath')$epath[[1]], OPI_g)
#'
#' # find all shortest paths between proteins
#' print_path_vs(
#'     all_shortest_paths(
#'         enzsub_g,
#'         from = "SRC",
#'         to = "STAT1"
#'     )$res,
#'     enzsub_g
#' )
#' }
#'
#' @author Alberto Valdeolivas <\email{alvaldeolivas@@gmail}> and
#' Denes Turei <\email{turei.denes@@gmail.com}> and Attila Gabor
#' <\email{gaborattila87@@gmail.com}>
#'
#' @name OmnipathR
#' @import logger
'_PACKAGE'

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
