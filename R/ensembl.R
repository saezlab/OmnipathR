#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2022
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


#' Table of Ensembl organisms
#'
#' A table with various taxon IDs and metadata about related Ensembl database
#' contents, as shown at https://www.ensembl.org/info/about/species.html.
#' The "Taxon ID" column contains the NCBI Taxonomy identifiers.
#'
#' @return The table described above as a data frame.
#'
#' @examples
#' ens_org <- ensembl_organisms_raw()
#' ens_org
#'
#' @importFrom magrittr %>%
#' @importFrom rvest read_html html_element html_table
#' @export
ensembl_organisms_raw <- function(){

    'ensembl_organisms' %>%
    download_to_cache %>%
    read_html %>%
    html_element('table') %>%
    html_table()

}


#' Organism names and identifiers from Ensembl
#'
#' A table with various taxon names and identifiers: English common names,
#' latin (scientific) names, Ensembl organism IDs and NCBI taxonomy IDs.
#'
#' @return A data frame with the above mentioned columns.
#'
#' @examples
#' ens_org <- ensembl_organisms()
#' ens_org
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate recode
#' @importFrom purrr pmap_chr map_chr
#' @importFrom stringr str_extract_all str_extract str_sub
#' @importFrom rlang !!!
#' @export
ensembl_organisms <- function(){

    SPECIALS <- list(
        bixbtaurus = 'bihybrid'
    )

    ensembl_organisms_raw() %>%
    select(
        common_name = `Common name`,
        latin_name = `Scientific name`,
        ncbi_tax_id = `Taxon ID`
    ) %>%
    mutate(
        ensembl_id =
            pmap_chr(
                list(
                    str_extract_all(latin_name, '\\b[A-z]') %>%
                    map_chr(paste0, collapse = '') %>%
                    str_to_lower %>%
                    str_sub(end = -2L),
                    str_extract(latin_name, '\\w+$')
                ),
                paste0
            ),
        ensembl_id = recode(ensembl_id, !!!SPECIALS)
    )

}
