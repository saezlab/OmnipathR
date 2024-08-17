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


#' Get all the molecular complexes for a given gene(s)
#'
#' This function returns all the molecular complexes where an input set
#' of genes participate. User can choose to retrieve every complex where
#' any of the input genes participate or just retrieve these complexes where
#' all the genes in input set participate together.
#'
#' @param complexes Data frame of protein complexes (obtained using
#'     \code{\link{complexes}}).
#' @param genes Character: search complexes where these genes present.
#' @param all_genes Logical: select only complexes where all of the genes
#'     present together. By default complexes where any of the genes can be
#'     found are returned.
#'
#' @export
#'
#' @return Data frame of complexes
#'
#' @examples
#' complexes <- complexes(
#'     filter_databases = c("CORUM", "hu.MAP")
#' )
#' query_genes <- c("LMNA", "BANF1")
#' complexes_with_query_genes <- complex_genes(complexes, query_genes)
#'
#' @seealso \code{\link{complexes}}
#' @aliases get_complex_genes
complex_genes <- function(
    complexes = complexes(),
    genes,
    all_genes = FALSE
){

    if(is.null(genes)){
        stop('A vector of genes should be provided')
    }

    if (!is.logical(all_genes)){
        stop('all_genes parameter should be logical')
    }

    if (all_genes){
        complexes_geneset <-
        complexes[
            which(unlist(lapply(
                strsplit(complexes$components_genesymbols, '_'),
                function(x){
                    sum(x %in% genes) == length(genes)
                }
            ))),
        ]
    } else {
        complexes_geneset <-
        complexes[
            which(unlist(lapply(
                strsplit(complexes$components_genesymbols, '_'),
                function(x){any(x %in% genes)}
            ))),
        ]
    }

    return(complexes_geneset)
}


# Aliases (old names) to be Deprecated
#' @rdname complex_genes
#' @param ... Passed to \code{complex_genes}.
#' @export
#'
#' @noRd
get_complex_genes <- function(...){
    .Deprecated('get_complex_genes')
    complex_genes(...)
}
