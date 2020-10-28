#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2020
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


#' Builds all prior knowledge data required by NicheNet
#'
#' @param signaling_network A list of parameters for building the signaling
#'     network, passed to \code{\link{nichenet_signaling_network}}
#' @export
#'
#' @seealso \code{\link{nichenet_signaling_network}}
nichenet_prior_knowledge <- function(
    signaling_network = list()
){

    list(
        signaling_network = do.call(
            nichenet_signaling_network,
            signaling_network
        )
    )

}


#' Builds signaling network prior knowledge for NicheNet using multiple
#' resources
#'
#' @param omnipath List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_omnipath}}
#' @param pathwaycommons List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_pathwaycommons}}
#' @param harmonizome List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_harmonizome}}
#' @export
#'
#' @seealso \code{\link{nichenet_signaling_network_omnipath},
#'     \link{nichenet_signaling_network_pathwaycommons},
#'     \link{nichenet_signaling_network_harmonizome}}
nichenet_signaling_network <- function(
    omnipath = list(),
    pathwaycommons = NULL,
    harmonizome = NULL,
    vinayagam = NULL,
    cpdb = NULL,
    evex = NULL,
    inweb = NULL
){



}


#' Builds signaling network prior knowledge for NicheNet using OmniPath
#'
#' @param min_curation_effort Lower threshold for curation effort
#' @param ... Passed to \code{\link{import_post_translational_interactions}}
#' @importsFrom dplyr %>% mutate select
#' @export
#'
#' @seealso
nichenet_signaling_network_omnipath <- function(
    min_curation_effort = 0,
    ...
){

    import_post_translational_interactions(entity_types = 'protein', ...) %>%
    select(from = source_genesymbol, to = target_genesymbol, is_directed) %>%
    mutate(
        source = ifelse(
            is_directed,
            'omnipath_directed',
            'omnipath_undirected'
        ),
        database = 'omnipath'
    ) %>%
    select(-is_directed)

}

#' importsFrom dplyr %>% mutate
#' @importsFrom readr read_tsv
nichenet_signaling_network_pathwaycommons <- function(...){

    read_tsv()

}