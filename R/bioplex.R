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


#' Downloads a BioPlex interaction dataset
#'
#' @param version Character: either 1.0, 2.0, 3.0 or HCT116_1.0
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv cols
#'
#' @seealso \code{\link{bioplex1}, \link{bioplex2}, \link{bioplex3},
#'     \link{bioplex_hct116_1}}
bioplex_download <- function(version){

    result <-
        'omnipath.bioplex_%s_url' %>%
        generic_downloader(
            url_key_param = list(version),
            resource = sprintf('BioPlex %s', version)
        )

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

    result %T>% load_success

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


#' Downloads all BioPlex interaction datasets
#'
#' BioPlex provides four interaction datasets: version 1.0, 2.0, 3.0 and
#' HCT116 version 1.0. This function downloads all of them, merges them to
#' one data frame, removes the duplicates (based on unique pairs of UniProt
#' IDs) and separates the isoform numbers from the UniProt IDs.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @param unique Logical. Collapse the duplicate interactions into single
#'     rows or keep them as they are. In case of merging duplicate records
#'     the maximum p value will be choosen for each record.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate first summarize_all group_by bind_rows
#' @importFrom tidyr separate
#' @export
#' @seealso \code{\link{bioplex1}, \link{bioplex2}, \link{bioplex3},
#' \link{bioplex_hct116_1}}
bioplex_all <- function(unique = TRUE){

    bind_rows(
        bioplex1(),
        bioplex2(),
        bioplex3(),
        bioplex_hct116_1()
    ) %>%
    {`if`(
        unique,
        group_by(., UniprotA, UniprotB) %>%
        mutate(
            p_wrong = max(p_wrong),
            p_no_interaction = max(p_no_interaction),
            p_interaction = max(p_interaction)
        ) %>%
        summarize_all(first),
        .
    )} %>%
    separate(
        UniprotA,
        into = c('UniprotA', 'IsoformA'),
        sep = '-',
        fill = 'right',
        convert = TRUE
    ) %>%
    separate(
        UniprotB,
        into = c('UniprotB', 'IsoformB'),
        sep = '-',
        fill = 'right',
        convert = TRUE
    )

}