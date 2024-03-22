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


#' Imports protein complexes from OmniPath
#'
#' Imports the complexes stored in Omnipath database from
#' \url{https://omnipathdb.org/complexes}.
#'
#' @return A dataframe containing information about complexes
#' @export
#'
#' @param resources complexes not reported in these databases are
#' removed. See \code{\link{get_complexes_databases}} for more information.
#' @param ... optional additional arguments
#'
#' @examples
#' complexes = import_omnipath_complexes(
#'     resources = c('CORUM', 'hu.MAP')
#' )
#'
#' @seealso \itemize{\item{\code{\link{get_complexes_databases}}}}
#'
#' @aliases import_Omnipath_complexes import_OmniPath_complexes
import_omnipath_complexes <- function(
    resources = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'complexes',
        resources = resources,
        ...
    )

    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_omnipath_complexes
#' @param ... Passed to \code{import_omnipath_complexes}.
#' @export
#'
#' @noRd
import_Omnipath_complexes <- function(...){
    .Deprecated("import_omnipath_complexes")
    import_omnipath_complexes(...)
}


#' @rdname import_omnipath_complexes
#' @param ... Passed to \code{import_omnipath_complexes}.
#' @export
#'
#' @noRd
import_OmniPath_complexes <- function(...){
    .Deprecated("import_omnipath_complexes")
    import_omnipath_complexes(...)
}


#' Retrieve a list of complex resources available in Omnipath
#'
#' Get the names of the resources from \url{https://omnipath.org/complexes}
#'
#' @param dataset ignored for this query type
#' @return character vector with the names of the databases
#'
#' @examples
#' get_complex_resources()
#'
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_omnipath_complexes}}}
#' }
#'
#' @aliases get_complexes_databases
get_complex_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'complexes', datasets = dataset))

}


# Aliases (old names) to be deprecated
#' @rdname get_complex_resources
#' @param ... Passed to \code{import_omnipath_enzsub}.
#' @export
#'
#' @noRd
get_complexes_databases <- function(...){
    .Deprecated("get_complex_resources")
    get_complex_resources(...)
}
