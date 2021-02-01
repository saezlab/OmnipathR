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


#' The OmnipathR package
#'
#' @description
#' OmnipathR is an R package built to provide easy access to the data stored
#' in the OmniPath web service:
#'
#' \url{https://omnipathdb.org/}
#'
#' The web service implements a very simple REST style API. This package make
#' requests by the HTTP protocol to retreive the data. Hence, fast Internet
#' access is required for a propser use of OmnipathR.
#'
#' The package also provides some utility functions to filter, analyse and
#' visualize the data.
#'
#' @examples
#' # Download post-translational modifications:
#' ptms = import_omnipath_enzsub(resources=c("PhosphoSite", "SIGNOR"))
#'
#' # Download protein-protein interactions
#' interactions = import_omnipath_interactions(resources=c("SignaLink3"))
#'
#' # Convert to igraph objects:
#' ptms_g = ptms_graph(ptms = ptms )
#' OPI_g = interaction_graph(interactions = interactions )
#'
#' # Print some interactions:
#' print_interactions(head(ptms))
#'
#' # interactions with references:
#' print_interactions(tail(ptms),writeRefs=TRUE)
#'
#' # find interactions between kinase and substrate:
#' print_interactions(dplyr::filter(ptms,enzyme_genesymbol=="MAP2K1",
#'    substrate_genesymbol=="MAPK3"))
#'
#' # find shortest paths on the directed network between proteins
#' print_path_es(shortest_paths(OPI_g,from = "TYRO3",to = "STAT3",
#'    output = 'epath')$epath[[1]],OPI_g)
#'
#' # find all shortest paths between proteins
#' print_path_vs(
#'     all_shortest_paths(
#'         ptms_g,
#'         from = "SRC",
#'         to = "STAT1"
#'     )$res,
#'     ptms_g
#' )
#'
#' @author Alberto Valdeolivas <\email{alvaldeolivas@@gmail}> and
#' Denes Turei <\email{turei.denes@@gmail.com}> and Attila Gabor
#' <\email{gaborattila87@@gmail.com}>
#'
#' @docType package
#' @name OmnipathR
#' @import methods
#' @import logger
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
