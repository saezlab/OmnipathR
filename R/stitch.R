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
stitch_links <- function(organism = 'human', prefixes = FALSE) {

    .slow_doctest()

    organism %<>% organism_for('stitch')
    log_trace('Loading STITCH protein-small molecule links.')

    'stitch_links' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = list(organism),
        reader_param = list(trim_ws = TRUE),
        resource = NULL,
        post = NULL,
        use_httr = FALSE
    ) %>%
    stitch_remove_prefixes(chemical, protein, remove = !prefixes) %T>%
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
stitch_actions <- function(organism = 'human', prefixes = FALSE) {

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
    ) %>%
    stitch_remove_prefixes(item_id_a, item_id_b, remove = !prefixes) %T>%
    load_success()

}


#' Remove the prefixes from STITCH identifiers
#'
#' STITCH adds the NCBI Taxonomy ID as a prefix to Ensembl protein identifiers,
#' e.g. "9606.ENSP00000170630", and "CID" followed by "s" or "m"
#' (stereospecific or merged, respectively) in front of PubChem Compound
#' Identifiers. It also pads the CID with zeros. This function removes these
#' prefixes, leaving only the identifiers.
#'
#' @param d Data frame, typically the output of \code{\link{stitch_links}} or
#'     \code{\link{stitch_actions}}.
#' @param ... Names of columns to remove prefixes from. NSE is supported.
#' @param remove Logical: remove the prefixes? If FALSE, this function does
#'     nothing.
#'
#' @examples
#' stitch_remove_prefixes(
#'     tibble(a = c('9606.ENSP00000170630', 'CIDs00012345')),
#'     a
#' )
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang enquos
#' @importFrom purrr map_chr
#' @importFrom dplyr mutate across
#' @export
stitch_remove_prefixes <- function(d, ..., remove = TRUE) {

    if(remove) {

        cols <- enquos(...) %>% map_chr(.nse_ensure_str)

        d %<>%
            mutate(
                across(
                    cols,
                    ~str_replace(.x, '^(\\d+\\.|CID[ms]0*)', '')
                )
            )

    }

    d

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

    .slow_doctest()

    # NSE vs. R CMD check workaround
    chemical <- protein <- item_id_a <- item_id_b <- CID <- HMDB <-
    a_is_acting <- ensp <- genesymbol <- pubchem <- hmdb <-
    hmdb_source <- hmdb_target <- NULL

    organism %<>% organism_for('stitch')

    links <-
        stitch_links(organism) %>%
        filter(
            combined_score >= min_score,
            experimental >= min_score | database >= min_score
        ) %>%
        select(
            item_id_a = chemical,
            item_id_b = protein
        ) %>%
        distinct %>%
        bind_rows(
            .,
            rename(
                .,
                item_id_a = item_id_b,
                item_id_b = item_id_a
            )
        )

    stitch_actions(organism) %>%
    filter(
        mode == 'activation' |
        mode == 'inhibition',
        a_is_acting
    ) %>%
    inner_join(links, by = c('item_id_a', 'item_id_b')) %>%
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
    translate_ids(item_id_a = pubchem, hmdb_a = hmdb, hmdb = TRUE) %>%
    translate_ids(item_id_b = pubchem, hmdb_b = hmdb, hmdb = TRUE) %>%
    mutate(
        item_id_a = coalesce(genesymbol_a, hmdb_a, item_id_a),
        item_id_b = coalesce(genesymbol_b, hmdb_b, item_id_b)
    ) %>%
    select(
        source = item_id_a,
        target = item_id_b,
        sign = action
    ) %>%
    mutate(
        across(
            c(source, target),
            ~str_replace(.x, '^(HMDB|\\d)', 'Metab__\\1')
        ),
        sign = ifelse(sign == 'inhibition', -1L, 1L)
    ) %>%
    mutate(source = sprintf('%s_c', source))

}
