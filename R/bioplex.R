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


#' @param version Character: either 1.0, 2.0, 3.0 or HCT116_1.0
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv cols
#' @importFrom logger log_info
bioplex_download <- function(version){

    url <-
        'omnipath.bioplex_%s_url' %>%
        sprintf(version) %>%
        options() %>%
        `[[`(1)

    from_cache <- omnipath_cache_load(url = url)

    if(!is.null(from_cache)){

        logger::log_info('Loaded from cache: %s', url)
        return(from_cache)

    }

    result <- url %>% read_tsv(col_types = cols())

    omnipath_cache_save(data = result, url = url)

    names(result) <- c(
        'GeneA',
        'GeneB',
        'UniprotA',
        'UniprotB',
        'SymbolA',
        'SymbolB',
        'p_wrong',
        'p_no_interaction',
        'p_interaction'
    )

    return(result)

}


#' Downloads the BioPlex version 1.0 interaction dataset
#'
#' This dataset contains ~24,000 interactions detected in HEK293T cells
#' using 2,594 baits.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @export
#' @seealso \code{\link{bioplex2}, \link{bioplex3}, \link{bioplex_hct116_1}}
bioplex1 <- function(){

    bioplex_download(version = '1.0')

}


#' Downloads the BioPlex version 2.0 interaction dataset
#'
#' This dataset contains ~56,000 interactions detected in HEK293T cells
#' using 5,891 baits.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @export
#' @seealso \code{\link{bioplex1}, \link{bioplex3}, \link{bioplex_hct116_1}}
bioplex2 <- function(){

    bioplex_download(version = '2.0')

}


#' Downloads the BioPlex version 3.0 interaction dataset
#'
#' This dataset contains ~120,000 interactions detected in HEK293T cells
#' using 10,128 baits.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @export
#' @seealso \code{\link{bioplex1}, \link{bioplex2}, \link{bioplex_hct116_1}}
bioplex3 <- function(){

    bioplex_download(version = '3.0')

}


#' Downloads the BioPlex HCT116 version 1.0 interaction dataset
#'
#' This dataset contains ~71,000 interactions detected in HCT116 cells
#' using 5,522 baits.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @export
#' @seealso \code{\link{bioplex1}, \link{bioplex2}, \link{bioplex3}}
bioplex_hct116_1 <- function(){

    bioplex_download(version = 'HCT116_1.0')

}