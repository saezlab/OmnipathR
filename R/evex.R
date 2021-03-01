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
#' @importsFrom magrittr %>%
#' @importsFrom readr read_tsv cols
#' @importsFrom dplyr left_join mutate group_by summarize_all first ungroup
#' @importsFrom dplyr rename
#' @export
#'
#' @examples
#' evex_interactions <- evex_download()
evex_download <- function(...){

    tmp_tgz <- tempfile(fileext = '.tar.gz')
    tmpdir_ex <- tempdir()

    on.exit({
        unlink(tmpdir_ex)
        closeAllConnections()
    })

    'omnipath.evex_url' %>%
    options() %>%
    as.character() %>%
    download.file(destfile = tmp_tgz, quiet = TRUE)

    tmp_tgz %>%
    untar(exdir = tmpdir_ex)
    unlink(tmp_tgz)

    tmpdir_ex %>%
    file.path('EVEX_relations_9606.tab') %>%
    read_tsv(
        col_types = cols(
            source_entrezgene_id = col_character(),
            target_entrezgene_id = col_character()
        ),
        progress = FALSE
    ) %>%
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
        tmpdir_ex %>%
        file.path('EVEX_articles_9606.tab') %>%
        read_tsv(col_types = cols(), progress = FALSE) %>%
        rename(references = article_id),
        by = 'general_event_id'
    ) %>%
    mutate(references = sub('PMC?ID: ', '', references)) %>%
    group_by(general_event_id) %>%
    mutate(references = paste(references, sep = ',')) %>%
    summarize_all(first) %>%
    ungroup()

}