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


#' Retrieves the ConsensusPathDB network
#'
#' Compiles a table of binary interactions from ConsensusPathDB and
#' translates the UniProtKB ACs to Gene Symbols.
#'
#' @param complex_max_size Numeric: do not expand complexes with a higher
#'     number of elements than this. ConsensusPathDB does not contain
#'     conventional interactions but lists of participants, which might be
#'     members of complexes. Some records include dozens of participants and
#'     expanding them to binary interactions result thousands, sometimes
#'     hundreds of thousands of interactions from one single record. At the
#'     end, this process consumes >10GB of memory and results rather unusable
#'     data, hence it is recommended to limit the complex sizes at some low
#'     number.
#' @param min_score Numeric: each record in ConsensusPathDB comes with a
#'     confidence score, expressing the amount of evidences. The default
#'     value, a minimum score of 0.9 retains approx. the top 30% of the
#'     interactions.
#'
#' @importFrom readr read_tsv
#' @importFrom tidyr separate_rows
#' @importFrom dplyr mutate select inner_join left_join pull filter rename
#' @importFrom dplyr group_by ungroup
#' @importFrom magrittr %>% %T>%
#' @importFrom stringr str_count
#' @export
consensuspathdb_download <- function(
    complex_max_size = 4,
    min_score = .9
){

    cpdb_raw <- consensuspathdb_raw_table()

    uniprot_genesymbol <- cpdb_raw %>%
        separate_rows(participants, sep = '[,\\.]') %>%
        pull(participants) %>%
        unique() %>%
        uniprot_id_mapping_table(from = 'ACC', to = 'GENENAME')

    cpdb_raw %>%
    filter(
        str_count(participants, ',') <= complex_max_size,
        confidence >= min_score
    ) %>%
    mutate(
        record_id = seq(n()),
        uniprot_b = participants
    ) %>%
    rename(uniprot_a = participants) %>%
    separate_rows(uniprot_a, sep = ',') %>%
    separate_rows(uniprot_b, sep = ',') %>%
    group_by(record_id) %>%
    mutate(in_complex = n() > 2) %>%
    ungroup() %>%
    separate_rows(uniprot_a, sep = '\\.') %>%
    separate_rows(uniprot_b, sep = '\\.') %>%
    filter(!in_complex | uniprot_a != uniprot_b) %>%
    left_join(uniprot_genesymbol, by = c('uniprot_a' = 'From')) %>%
    rename(genesymbol_a = To) %>%
    left_join(uniprot_genesymbol, by = c('uniprot_b' = 'From')) %>%
    rename(genesymbol_b = To) %>%
    copy_source_attrs(cpdb_raw) %T>%
    load_success()

}


#' Downloads interaction data from ConsensusPathDB
#'
#' @importFrom readr read_tsv cols
#' @importFrom magrittr %>%
#' @export
consensuspathdb_raw_table <- function(){

    'omnipath.cpdb_url' %>%
    generic_downloader(
        reader_param = list(
            col_names = c(
                'databases',
                'references',
                'participants',
                'confidence'
            ),
            skip = 2,
            progress = FALSE
        ),
        resource = 'ConsensusPathDB'
    )

}