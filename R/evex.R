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
#' (\url{http://evexdb.org}). Translates the Entrez Gene IDs to Gene Symbols
#' and combines the interactions and references into a single data frame.
#'
#' @usage
#' evex_download(
#'     min_confidence = NULL,
#'     remove_negatives = TRUE,
#'     top_confidence = NULL
#' )
#'
#' @param min_confidence Numeric: a threshold for confidence scores. EVEX
#'     confidence scores span roughly from -3 to 3. By providing a numeric
#'     value in this range the lower confidence interactions can be removed.
#'     If NULL no filtering performed.
#' @param remove_negatives Logical: remove the records with the "negation"
#'     attribute set.
#' @param top_confidence Confidence cutoff as quantile (a number between
#'     0 and 1). If NULL no filtering performed.
#'
#' @return Data frame (tibble) with interactions.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv cols
#' @importFrom dplyr left_join mutate group_by summarize_all first ungroup
#' @importFrom dplyr rename filter
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' evex_interactions <- evex_download()
#' evex_interactions
#' # # A tibble: 368,297 x 13
#' #   general_event_id source_entrezge. target_entrezge. confidence negation
#' #               <dbl> <chr>            <chr>                 <dbl>    <dbl>
#' # 1               98 8651             6774                 -1.45         0
#' # 2              100 8431             6774                 -1.45         0
#' # 3              205 6261             6263                  0.370        0
#' # 4              435 1044             1045                 -1.09         0
#' # . with 368,287 more rows, and 8 more variables: speculation <dbl>,
#' #   coarse_type <chr>, coarse_polarity <chr>, refined_type <chr>,
#' #   refined_polarity <chr>, source_genesymbol <chr>,
#' #   target_genesymbol <chr>, references <chr>
evex_download <- function(
    min_confidence = NULL,
    remove_negatives = TRUE,
    top_confidence = NULL
){

    # NSE vs. R CMD check workaround
    negation <- confidence <- source_entrezgene_id <- target_entrezgene_id <-
        article_id <- general_event_id <- entrez <- NULL

    relations <- archive_extractor(
        url_key = 'evex',
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
        url_key = 'evex',
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
        is.null(top_confidence),
        .,
        filter(., confidence > quantile(confidence, top_confidence))
    )} %>%
    {`if`(
        is.null(min_confidence),
        .,
        filter(., confidence >= min_confidence)
    )} %>%
    translate_ids(
        source_entrezgene_id = entrez,
        source_genesymbol = genesymbol
    ) %>%
    translate_ids(
        target_entrezgene_id = entrez,
        target_genesymbol = genesymbol
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