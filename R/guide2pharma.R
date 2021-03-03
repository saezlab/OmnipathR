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


#' Downloads interactions from the Guide to Pharmacology database
#'
#' Downloads ligand-receptor interactions from the Guide to Pharmacology
#' (IUPHAR/BPS) database.
#'
#' @return A tibble (data frame) of interactions as it is provided by the
#' database
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom readr cols col_character col_number
#' @importFrom tidyr separate_rows
guide2pharma_download <- function(){

    'omnipath.guide2pharma_url' %>%
    generic_downloader(
        reader_param = list(
            col_types = cols(
                ligand_context = col_character(),
                receptor_site = col_character(),
                target_ligand = col_character(),
                target_ligand_id = col_character(),
                target_ligand_pubchem_sid = col_number(),
                target_ligand_gene_symbol = col_character(),
                target_ligand_uniprot = col_character(),
                target_ligand_ensembl_gene_id = col_character()
            )
        ),
        resource = 'Guide to Pharmacology (IUPHAR/BPS)'
    ) %>%
    separate_rows(ligand_gene_symbol, sep = '\\|') %>%
    separate_rows(target_gene_symbol, sep = '\\|') %T>%
    load_success()

}