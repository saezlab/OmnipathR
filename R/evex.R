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

#' Interactions from the EVEX database
#'
#' Downloads interactions from EVEX, a versatile text mining resource
#' (http://evexdb.org). Translates the Entrez Gene IDs to Gene Symbols and
#' combines the interactions and references into a single data frame.
#'
#' @param min_confidence Numeric: a threshold for confidence scores. EVEX
#'     confidence scores span roughly from -3 to 3. By providing a numeric
#'     value in this range the lower confidence interactions can be removed.
#' @param remove_negatives Logical: remove the records with the "negation"
#'     attribute set.
#'
#' @importsFrom magrittr %>% %T>%
#' @importsFrom readr read_tsv cols
#' @importsFrom dplyr left_join mutate group_by summarize_all first ungroup
#' @importsFrom dplyr rename filter
#' @export
#'
#' @examples
#' evex_interactions <- evex_download()
evex_download <- function(
    min_confidence = NULL,
    remove_negatives = TRUE
){

    relations <- archive_extractor(
        url_key = 'omnipath.evex_url',
        path = 'EVEX_relations_9606.tab',
        reader = read_tsv,
        reader_param = list(
            col_types = cols(
                source_entrezgene_id = col_character(),
                target_entrezgene_id = col_character()
            ),
            progress = FALSE
        ),
        resource = 'EVEX'
    )

    articles <- archive_extractor(
        url_key = 'omnipath.evex_url',
        path = 'EVEX_articles_9606.tab',
        reader = read_tsv,
        reader_param = list(
            col_types = cols(),
            progress = FALSE
        )
    )

    relations %>%
    {`if`(
        remove_negatives,
        filter(., negation == 0),
        .
    )} %>%
    {`if`(
        is.null(min_confidence),
        .,
        filter(., confidence >= min_confidence)
    )} %>%
    translate_ids(
        source_entrezgene_id,
        source_genesymbol,
        'entrez',
        'genesymbol',
        uploadlists = FALSE
    ) %>%
    translate_ids(
        target_entrezgene_id,
        target_genesymbol,
        'entrez',
        'genesymbol',
        uploadlists = FALSE
    ) %>%
    left_join(
        articles %>%
        rename(references = article_id),
        by = 'general_event_id'
    ) %>%
    mutate(references = sub('PMC?ID: ', '', references)) %>%
    group_by(general_event_id) %>%
    mutate(references = paste(references, sep = ',')) %>%
    summarize_all(first) %>%
    ungroup() %>%
    copy_source_attrs(relations) %T>%
    load_success()

}