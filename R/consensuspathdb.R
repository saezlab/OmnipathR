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


#' Retrieves the ConsensusPathDB network
#'
#' Compiles a table of binary interactions from ConsensusPathDB and
#' translates the UniProtKB ACs to Gene Symbols.
#'
#' @importsFrom readr read_tsv
#' @importsFrom tidyr separate_rows
#' @importsFrom dplyr mutate select inner_join left_join pull filter rename
#' @importsFrom dplyr group_by ungroup
#' @importsFrom magrittr %>%
#' @export
consensuspathdb <- function(){

    cpdb_raw <-  consensuspathdb_raw_table()

    uniprot_genesymbol <- cpdb_raw %>%
        separate_rows(participants, sep = '[,.]') %>%
        pull(participants) %>%
        unique() %>%
        uniprot_id_mapping(from = 'ACC', to = 'GENENAME')

    cpdb_raw %>%
    mutate(
        record_id = n(),
        uniprot_b = participants
    ) %>%
    rename(uniprot_a = participants) %>%
    separate_rows(uniprot_a, sep = ',') %>%
    separate_rows(uniprot_b, sep = ',') %>%
    group_by(record_id) %>%
    mutate(in_complex = n() > 2) %>%
    ungroup() %>%
    separate_rows(uniprot_a, sep = '.') %>%
    separate_rows(uniprot_b, sep = '.') %>%
    filter(uniprot_a != uniprot_b) %>%
    left_join(uniprot_genesymbol, by = c('uniprot_a' = 'From')) %>%
    rename(genesymbol_a = To) %>%
    left_join(uniprot_genesymbol, by = c('uniprot_b' = 'From')) %>%
    rename(genesymbol_b = To)

}


#' Downloads interaction data from ConsensusPathDB
#'
#' @importsFrom readr read_tsv cols
#' @importsFrom magrittr %>%
#' @export
consensuspathdb_raw_table <- function(){

    'omnipath.cpdb_url' %>%
    options() %>%
    as.character() %>%
    read_tsv(
        col_names = c(
            'databases',
            'references',
            'participants',
            'confidence'
        ),
        col_types = cols(),
        skip = 2,
        progress = FALSE
    )

}