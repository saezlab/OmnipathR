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


#' Protein complexes from OmniPath
#'
#' A comprehensive dataset of protein complexes from the
#' \url{https://omnipathdb.org/complexes} endpoint of the OmniPath web service.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -query_type -datasets -types -json_param
#'     -add_counts -references_by_resource
#'
#' @return A dataframe containing information about complexes
#'
#' @examples
#' complexes = import_omnipath_complexes(
#'     resources = c('CORUM', 'hu.MAP')
#' )
#'
#' @seealso \itemize{\item{\code{\link{get_complexes_databases}}}}
#'
#' @importFrom rlang exec !!!
#' @export
import_omnipath_complexes <- function(...){

    args <- omnipath_args(list(...), query_type = 'complexes')

    exec(omnipath_query, !!!args)

}


#' Retrieve a list of complex resources available in Omnipath
#'
#' Get the names of the resources from \url{https://omnipath.org/complexes}
#'
#' @param dataset ignored for this query type
#'
#' @examples
#' get_complex_resources()
#'
#' @return character vector with the names of the databases
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_omnipath_complexes}}}
#' }
#'
#' @export
get_complex_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'complexes', datasets = dataset))

}
