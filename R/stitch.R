#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Diego Mananes
#'                 Alberto Valdeolivas
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


#' Download STITCH link data frame from \url{http://stitch.embl.de/}
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Data frame of STITCH links.
#'
#' @importFrom magrittr %<>% %>% %T>%
#' @importFrom readr read_tsv
#'
#' @noRd
stitch_proteins <- function(organism = 'human') {

    .slow_doctest()

    organism %<>% organism_for('stitch')
    log_trace('Loading STITCH proteins.')

    'stitch_proteins' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = list(organism),
        reader_param = list(trim_ws = TRUE),
        resource = NULL,
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' Download STITCH actions data frame from \url{http://stitch.embl.de/}
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Data frame of STITCH actions.
#'
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom readr read_tsv
#'
#' @noRd
stitch_actions <- function(organism = 'human') {

    .slow_doctest()

    organism %<>% organism_for('stitch')
    log_trace('Loading STITCH actions.')

    'stitch_actions' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = list(organism),
        reader_param = list(trim_ws = TRUE),
        resource = NULL,
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' Chemical-protein interactions from STITCH
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#' @param threshold Confidence cutoff used for STITCH connections
#'     (700 by default).
#'
#' @return List containing PKN with COSMOS and OCEAN format.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter case_when mutate select
#' @importFrom tidyr unite
#'
#' @noRd
stitch_gem <- function(organism = 'human', min_score = 700L) {

    organism %<>% organism_for('stitch')

    actions <-
        stitch_actions(organism) %>%
        filter(
            combined_score >= threshold,
            experimental >= threshold | database >= threshold
        ) %>%
        select(
            item_id_a = chemical,
            item_id_b = protein
        ) %>%
        distinct %>%
        bind_rows(
            .,
            rename(
                item_id_a = item_id_b,
                item_id_b = item_id_a
            )
        )

    metabolites <-
        metaboliteIDMapping::metabolitesMapping %>%
        select(CID, HMDB) %>%
        mutate(across(CID, HMDB, ~sprintf('Metab__%s', .x)))


    stitch_actions(organism) %>%
    filter(
        mode == 'activation' |
        mode == 'inhibition',
        a_is_acting
    ) %>%
    inner_join(actions, by = c('item_id_a', 'item_id_b')) %>%
    select(-7L) %>%
    mutate(
        across(
            item_id_a,
            item_id_b,
            ~str_replace(.x, sprintf('%s\\.', organism), '')
        )
    ) %>%
    translate_ids(
        item_id_a = ensp,
        genesymbol_a = genesymbol,
        ensembl = TRUE,
        organism = organism
    ) %>%
    translate_ids(
        item_id_b = ensp,
        genesymbol_b = genesymbol,
        ensembl = TRUE,
        organism = organism
    ) %>%
    mutate(
        across(
            item_id_a,
            item_id_b,
            ~str_replace(.x, '^CID[a-z]0*', 'Metab__')
        ),
        sign = ifelse(action == 'inhibition', -1L, 1L)
    ) %>%
    select(
        source = item_id_a,
        target = item_id_b,
        sign
    ) %>%
    left_join(metabolites, by = c('source' = 'CID')) %>%
    left_join(metabolites, by = c('target' = 'CID')) %>%
    mutate(
        source = ifelse(is.na(HMDB.x), source, HMDB.x),
        target = ifelse(is.na(HMDB.y), target, HMDB.y)
    ) %>%
    select(-HMDB.x, -HMDB.y) %>%
    mutate(source = sprintf('%s_c', source))

}
