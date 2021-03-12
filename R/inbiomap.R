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


#' Downloads network data from InWeb InBioMap
#'
#' Downloads the data from
#' \url{https://inbio-discover.com/map.html#downloads} in tar.gz format,
#' extracts the PSI MITAB table and returns it as a data frame.
#'
#' @param curl_verbose Logical. Perform CURL requests in verbose mode for
#'     debugging purposes.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom RCurl basicHeaderGatherer basicTextGatherer CFILE
#' @importFrom RCurl getCurlHandle curlSetOpt curlPerform close
#' @importFrom jsonlite parse_json
#' @importFrom readr read_tsv cols
#' @export
#'
#' @return A data frame (tibble) with the extracted interaction table.
#'
#' @seealso \code{\link{inbiomap_download}}
#'
#' @examples
#' inbiomap_psimitab <- inbiomap_raw()
inbiomap_raw <- function(curl_verbose = FALSE){

    url <- 'omnipath.inbiomap_url' %>% url_parser
    version <- omnipath_cache_latest_or_new(url = url, create = FALSE)
    token <- 'none'

    if(is.null(version) || version$status != CACHE_STATUS$READY){

        # acquiring an access token
        resp_headers <- basicHeaderGatherer()
        login_response <- basicTextGatherer()
        payload <- '{"ref":""}'
        login_curl <- getCurlHandle()
        opt_set <- curlSetOpt(
            curl = login_curl,
            writefunction = login_response$update,
            headerfunction = resp_headers$update,
            verbose = curl_verbose,
            post = TRUE,
            postfields = payload,
            followlocation = TRUE,
            httpheader = 'Content-Type: application/json;charset=utf-8',
            url = 'omnipath.inbiomap_login_url' %>% url_parser
        )
        success <- curlPerform(curl = login_curl)
        token <- (login_response$value() %>% parse_json())$token

    }

    # doing the actual download or reading from the cache
    archive_extractor(
        url_key = 'omnipath.inbiomap_url',
        path = 'InBio_Map_core_2016_09_12/core.psimitab',
        reader = read_tsv,
        reader_param = list(
            col_types = cols(),
            col_names = c(
                'id_a',
                'id_b',
                'id_alt_a',
                'id_alt_b',
                'synonyms_a',
                'synonyms_b',
                'detection_methods',
                'ref_first_authors',
                'references',
                'organism_a',
                'organism_b',
                'interaction_types',
                'databases',
                'database_ids',
                'score',
                'complex_expansion'
            ),
            progress = FALSE
        ),
        curl_verbose = curl_verbose,
        # this is passed to curlSetOpt
        httpheader = sprintf('Cookie: access_token=%s', token),
        resource = 'InWeb InBioMap'
    ) %T>%
    load_success()

}


#' Downloads and preprocesses network data from InWeb InBioMap
#'
#' Downloads the data by \code{\link{inbiomap_raw}}, extracts the
#' UniProt IDs, Gene Symbols and scores and removes the irrelevant columns.
#'
#' @param ... Passed to \code{\link{inbiomap_raw}}.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @importFrom tidyr separate
#' @export
#'
#' @return A data frame (tibble) of interactions.
#'
#' @seealso \code{\link{inbiomap_raw}}
#'
#' @examples
#' inbiomap_interactions <- inbiomap_download()
#' inbiomap_interactions
#' # # A tibble: 625,641 x 7
#' #    uniprot_a uniprot_b genesymbol_a genesymbol_b inferred score1 score2
#' #    <chr>     <chr>     <chr>        <chr>        <lgl>     <dbl>  <dbl>
#' #  1 A0A5B9    P01892    TRBC2        HLA-A        FALSE     0.417 0.458
#' #  2 A0AUZ9    Q96CV9    KANSL1L      OPTN         FALSE     0.155 0.0761
#' #  3 A0AV02    P24941    SLC12A8      CDK2         TRUE      0.156 0.0783
#' #  4 A0AV02    Q00526    SLC12A8      CDK3         TRUE      0.157 0.0821
#' #  5 A0AV96    P0CG48    RBM47        UBC          FALSE     0.144 0.0494
#' # # . with 625,631 more rows
inbiomap_download <- function(...){

    # NSE vs. R CMD check workaround
    score <- id_a <- id_b <- synonyms_a <- synonyms_b <- detection_methods <-
        score1 <- score2 <- uniprot_a <- uniprot_b <- genesymbol_a <-
        genesymbol_b <- inferred <- NULL

    .extract_id <- function(x){

        sub('[[:alnum:]]+:([[:alnum:]]+)', '\\1', x)

    }

    .extract_gs <- function(x){

        sub('.*:([[:alnum:]-]+)\\(gene name\\).*', '\\1', x)

    }

    inbiomap_raw(...) %>%
    separate(score, into = c('score1', 'score2'), sep = '\\|') %>%
    mutate(
        uniprot_a = .extract_id(id_a),
        uniprot_b = .extract_id(id_b),
        genesymbol_a = .extract_gs(synonyms_a),
        genesymbol_b = .extract_gs(synonyms_b),
        inferred = grepl('inference', detection_methods, fixed = TRUE),
        score1 = suppressWarnings(as.numeric(score1)),
        score2 = suppressWarnings(as.numeric(score2))
    ) %>%
    select(
        uniprot_a,
        uniprot_b,
        genesymbol_a,
        genesymbol_b,
        inferred,
        score1,
        score2
    )

}