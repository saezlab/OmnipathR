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


#' Downloads network data from InWeb InBioMap
#'
#' Downloads the data from
#' https://inbio-discover.intomics.com/map.html#downloads in tar.gz format,
#' extracts the PSI MITAB table and returns it as a data frame.
#'
#' @param .verbose Logical. Perform CURL requests in verbose mode for
#'     debugging purposes.
#'
#' @importsFrom magrittr %>%
#' @importsFrom RCurl basicHeaderGatherer basicTextGatherer CFILE
#' @importsFrom RCurl getCurlHandle curlSetOpt curlPerform close
#' @importsFrom jsonlite parse_json
#' @importsFrom readr read_tsv cols
#' @export
#'
#' @return A data frame (tibble) with the extracted interaction table.
#'
#' @seealso \code{\link{inbiomap}}
#'
#' @examples
#' inbiomap_psimitab <- inbiomap()
inbiomap_raw <- function(.verbose = FALSE){

    tmp_tgz <- tempfile(fileext = '.tar.gz')
    tmpdir_ex <- tempdir()

    on.exit({
        closeAllConnections()
        unlink(tmp_tgz)
        unlink(tmpdir_ex)
    })

    # acquiring an access token
    resp_headers <- basicHeaderGatherer()
    login_response <- basicTextGatherer()
    payload <- '{"ref":""}'
    login_curl <- getCurlHandle()
    opt_set <- curlSetOpt(
        curl = login_curl,
        writefunction = login_response$update,
        headerfunction = resp_headers$update,
        verbose = .verbose,
        post = TRUE,
        postfields = payload,
        httpheader = 'Content-Type: application/json;charset=utf-8',
        url = {
            'omnipath.inbiomap_login_url' %>%
            options() %>%
            as.character()
        }
    )
    success <- curlPerform(curl = login_curl)
    token <- (login_response$value() %>% parse_json())$token

    # downloading the data
    main_curl <- getCurlHandle()
    main_response <- CFILE(tmp_tgz, mode = 'wb')
    opt_set <- curlSetOpt(
        curl = main_curl,
        verbose = .verbose,
        httpheader = sprintf('Cookie: access_token=%s', token),
        url = {
            'omnipath.inbiomap_url' %>%
            options() %>%
            as.character()
        },
        writedata = main_response@ref
    )
    success <- curlPerform(curl = main_curl)
    RCurl::close(main_response)

    # extracting the archive
    tmp_tgz %>%
    untar(exdir = tmpdir_ex)

    tmpdir_ex %>%
    file.path('InBio_Map_core_2016_09_12', 'core.psimitab') %>%
    read_tsv(
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
    )

}


#' Downloads and preprocesses network data from InWeb InBioMap
#'
#' Downloads the data by \code{\link{inbiomap_raw}}, extracts the
#' UniProt IDs, Gene Symbols and scores and removes the irrelevant columns.
#'
#' @importsFrom magrittr %>%
#' @importsFrom dplyr mutate select
#' @importsFrom tidyr separate
#' @export
#'
#' @return A data frame (tibble) of interactions.
#'
#' @seealso \code{\link{inbiomap_raw}}
#'
#' @examples
#' inbiomap_interactions <- inbiomap()
inbiomap <- function(...){

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