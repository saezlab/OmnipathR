#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2026
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
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' List WikiPathways pathways
#'
#' Retrieves the full pathway list from WikiPathways and optionally filters
#' by organism.
#'
#' @param organism Character or integer: organism name, common name,
#'     latin name, or NCBI Taxonomy ID. If \code{NULL}, returns all
#'     organisms.
#'
#' @return A tibble with columns \code{pathway_id}, \code{pathway_name},
#'     \code{pathway_url}, \code{organism_species}.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr filter transmute
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom rlang %||% .data
#'
#' @export
#' @examples
#' \dontrun{
#' pathways <- wikipathways_pathways(organism = 9606)
#' head(pathways)
#' }
wikipathways_pathways <- function(organism = NULL) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    organism_species <- NULL

    data <- generic_downloader(
        url_key = 'wikipathways_list',
        reader = curl_read_json,
        reader_param = list(simplifyVector = FALSE),
        resource = 'WikiPathways'
    )

    organisms <- data$organisms

    if (is.null(organisms) || length(organisms) == 0L) {
        return(
            tibble(
                pathway_id = character(),
                pathway_name = character(),
                pathway_url = character(),
                organism_species = character()
            )
        )
    }

    result <- map_dfr(organisms, function(org) {

        pws <- org$pathways
        if (is.null(pws) || length(pws) == 0L) return(tibble())

        map_dfr(pws, function(pw) {
            tibble(
                pathway_id = pw$id %||% NA_character_,
                pathway_name = pw$name %||% NA_character_,
                pathway_url = pw$url %||% NA_character_,
                organism_species = pw$species %||% NA_character_
            )
        })

    })

    if (!is.null(organism)) {
        organism %<>% latin_name
        result %<>% filter(.data$organism_species == organism)
    }

    result %>%
        transmute(pathway_id, pathway_name, pathway_url, organism_species)

}


#' Normalize metabolite identifier
#'
#' @param db Character: database name.
#' @param id Character: identifier.
#'
#' @return Character: normalized \code{DB:ID} or \code{NA_character_} if
#'     missing.
#'
#' @importFrom stringr str_detect str_trim str_to_upper
#'
#' @noRd
norm_met_id <- function(db, id) {

    if (
        is.null(db) || is.null(id) ||
        is.na(db) || is.na(id) ||
        db == '' || id == ''
    ) {
        return(NA_character_)
    }

    db %<>% as.character %>% str_trim
    id %<>% as.character %>% str_trim

    if (str_detect(id, 'HMDB|CHEBI|ChEBI|PubChem-compound|Wikidata|CAS|:')) {
        str_to_upper(id)
    } else {
        paste0(str_to_upper(db), ':', id)
    }

}


#' Download GPML and return XML document
#'
#' @param wpid Character: WikiPathways ID.
#'
#' @return An \code{xml_document}.
#'
#' @importFrom xml2 read_xml
#'
#' @noRd
fetch_gpml <- function(wpid) {

    download_to_cache(
        url_key = 'wikipathways_gpml',
        url_param = list(wpid, wpid),
        ext = 'gpml'
    ) %>%
    read_xml

}


#' Build metabolite-pathway table from GPML files
#'
#' @param pathways_tbl Tibble: output of
#'     \code{\link{wikipathways_pathways}}.
#' @param out_path Character: optional CSV output path for the main table.
#' @param failures_path Character: optional CSV output path for failures.
#'
#' @return A tibble with columns \code{pathway_id}, \code{pathway_name},
#'     \code{pathway_url}, \code{metabolites}. If any pathway failed,
#'     the failures tibble is attached as \code{attr(out, "failures")}.
#'
#' @importFrom dplyr bind_rows
#' @importFrom logger log_warn log_trace
#' @importFrom magrittr %>% %<>%
#' @importFrom progress progress_bar
#' @importFrom purrr map map_chr discard compact
#' @importFrom readr write_csv
#' @importFrom stringr str_trunc str_pad
#' @importFrom tibble tibble
#' @importFrom xml2 xml_attr xml_find_all
#'
#' @noRd
wikipathways_metabolite_table <- function(
    pathways_tbl,
    out_path = NULL,
    failures_path = NULL
) {

    n <- nrow(pathways_tbl)

    if (n == 0L) {

        empty <- tibble(
            pathway_id = character(),
            pathway_name = character(),
            pathway_url = character(),
            metabolites = character()
        )

        if (!is.null(out_path)) write_csv(empty, out_path)

        return(empty)

    }

    pb <- progress_bar$new(
        total = n,
        format = paste0(
            '  WikiPathways GPML: :wpid ',
            '[:bar] :current/:total | failures: :nfail'
        )
    )

    rows <- vector('list', n)
    failures <- list()

    for (i in seq_len(n)) {

        wpid <- pathways_tbl$pathway_id[i]

        pb$tick(tokens = list(
            wpid = wpid %>% str_trunc(12L) %>% str_pad(12L, 'right'),
            nfail = length(failures)
        ))

        tryCatch({

            met_ids <-
                fetch_gpml(wpid) %>%
                xml_find_all(
                    "//*[local-name()='DataNode'][@Type='Metabolite']"
                ) %>%
                map(~xml_find_all(.x, "./*[local-name()='Xref']")) %>%
                map(function(xrefs) {
                    map_chr(xrefs, function(xr) {
                        norm_met_id(xml_attr(xr, 'Database'), xml_attr(xr, 'ID'))
                    })
                }) %>%
                unlist %>%
                discard(~is.na(.x) || !nzchar(.x)) %>%
                unique %>%
                sort

            rows[[i]] <- tibble(
                pathway_id = wpid,
                pathway_name = pathways_tbl$pathway_name[i],
                pathway_url = pathways_tbl$pathway_url[i],
                metabolites = met_ids %>% paste(collapse = '; ')
            )

        }, error = function(e) {

            failures[[length(failures) + 1L]] <<- tibble(
                pathway_id = wpid,
                pathway_name = pathways_tbl$pathway_name[i],
                error = conditionMessage(e)
            )
            rows[[i]] <<- NULL
            log_warn(
                'WikiPathways GPML parse failed for %s: %s',
                wpid,
                conditionMessage(e)
            )

        })

    }

    result <- bind_rows(rows)

    if (length(failures) > 0L) {
        failures_tbl <- bind_rows(failures)
        attr(result, 'failures') <- failures_tbl
        log_warn(
            'WikiPathways: %d pathway(s) failed.',
            nrow(failures_tbl)
        )
        if (!is.null(failures_path)) write_csv(failures_tbl, failures_path)
    }

    if (!is.null(out_path)) write_csv(result, out_path)

    result

}


#' WikiPathways metabolites per pathway (GPML)
#'
#' Retrieves WikiPathways pathways (optionally filtered by organism) and
#' parses GPML files to extract metabolite identifiers.
#'
#' @param organism Character or integer: organism name, common name, or
#'     NCBI Taxonomy ID. Defaults to human.
#' @param out_path Character: optional CSV output path for the main table.
#' @param failures_path Character: optional CSV output path for failures.
#'
#' @return A tibble with columns \code{pathway_id}, \code{pathway_name},
#'     \code{pathway_url}, \code{metabolites}. If any pathway failed,
#'     the failures tibble is attached as \code{attr(out, "failures")}.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom stringr str_trim
#'
#' @export
#' @examples
#' \dontrun{
#' df <- wikipathways_metabolites(organism = 9606)
#' head(df)
#' }
wikipathways_metabolites <- function(
    organism = 9606L,
    out_path = NULL,
    failures_path = NULL
) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    metabolites <- NULL

    organism %>%
        wikipathways_pathways(organism = .) %>%
        wikipathways_metabolite_table(
            pathways_tbl = .,
            out_path = out_path,
            failures_path = failures_path
        ) %>%
        filter(!is.na(metabolites), str_trim(metabolites) != '')

}


#' WikiPathways metabolites via SPARQL
#'
#' Queries the WikiPathways SPARQL endpoint to retrieve all metabolite
#' members for pathways of a given organism. Results are paginated and
#' combined into a single tibble.
#'
#' Metabolite identifiers are normalised to the following formats:
#'     \itemize{
#'         \item HMDB: \code{HMDB0001049}
#'         \item ChEBI: \code{CHEBI:16015}
#'         \item PubChem: \code{CID7098621}
#'         \item Wikidata: \code{Q715317}
#'         \item CAS: numeric or hyphenated form without prefix
#'     }
#'
#' The output contains one row per pathway, with metabolite identifiers
#' collapsed into a semicolon-separated string. The number of unique
#' metabolites per pathway is also reported.
#'
#' @param organism Character or integer: organism name, common name, or
#'     NCBI Taxonomy ID. Defaults to human. If \code{NULL}, returns all
#'     organisms.
#' @param page_size Integer: number of records retrieved per SPARQL page.
#'     Default is 50000.
#' @param max_pages Integer: maximum number of pages to retrieve.
#'     Default is 200.
#' @param max_retries Integer: number of retry attempts per page.
#'     Default is 4.
#' @param sleep_base Numeric: base sleep time in seconds for exponential
#'     backoff between retries. Default is 1.
#'
#' @return A tibble with columns:
#'     \describe{
#'         \item{pathway_id}{WikiPathways identifier (e.g. \code{"WP254"})}
#'         \item{pathway_name}{Pathway title}
#'         \item{pathway_url}{Direct URL to the pathway}
#'         \item{n_metabolites_in_pathway}{Number of unique metabolite IDs}
#'         \item{metabolites}{Semicolon-separated metabolite identifiers}
#'     }
#'
#' @details
#' This implementation avoids downloading individual GPML files and instead
#' performs a single paginated SPARQL query against the public WikiPathways
#' endpoint. Identifier IRIs from identifiers.org are normalised and
#' redundant prefixes are removed.
#'
#' @importFrom dplyr bind_rows case_when filter group_by left_join
#'     mutate n_distinct select summarise
#' @importFrom logger log_error log_trace
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr compact map_chr
#' @importFrom rlang %||%
#' @importFrom stringr str_detect str_replace str_replace_all
#'     str_to_upper str_trim
#' @importFrom tibble tibble
#'
#' @export
#' @examples
#' \dontrun{
#' df <- wikipathways_metabolites_sparql(organism = 9606)
#' head(df)
#' }
wikipathways_metabolites_sparql <- function(
    organism = 9606L,
    page_size = 50000L,
    max_pages = 200L,
    max_retries = 4L,
    sleep_base = 1
) {

    # NSE vs. R CMD check workaround
    pathway_id <- pathway_name <- met_db <- met_xref_raw <-
        met_xref_clean <- met_id <- pathway_url <- NULL

    organism_latin <- NULL

    if (!is.null(organism)) {
        organism_latin <- latin_name(organism)
    }

    endpoint <- url_parser(url_key = 'wikipathways_sparql')

    extract_binding_value <- function(bindings, field) {

        if (is.data.frame(bindings)) {

            if (
                field %in% names(bindings) &&
                is.data.frame(bindings[[field]]) &&
                'value' %in% names(bindings[[field]])
            ) {
                return(bindings[[field]]$value)
            }

            return(rep(NA_character_, nrow(bindings)))

        }

        map_chr(bindings, function(binding) {
            binding[[field]]$value %||% NA_character_
        })

    }

    normalize_pubchem <- function(xref) {

        xref %>%
            str_replace('^PUBCHEM:', '') %>%
            str_replace('^PUBCHEM\\.COMPOUND:', '') %>%
            {if_else(str_detect(., '^[0-9]+$'), paste0('CID', .), .)} %>%
            str_replace('^(CID[0-9]+).*$', '\\1')

    }

    species_clause <-
        if (is.null(organism_latin)) {
            ''
        } else {
            organism_latin %>%
                str_replace_all('"', '\\"') %>%
                sprintf('           wp:organismName "%s" ;\n', .)
        }

    all_pages <- list()

    for (page_idx in seq_len(max_pages)) {

        sparql_query <- sprintf(
            '
            PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>
            PREFIX dc:      <http://purl.org/dc/elements/1.1/>
            PREFIX dcterms: <http://purl.org/dc/terms/>

            SELECT DISTINCT
                ?pathway_id
                (STR(?pathway_name) AS ?pathway_name_str)
                ?met_db
                ?met_xref
            WHERE {
                ?pathway a wp:Pathway ;
                %s
                dcterms:identifier ?pathway_id ;
                dc:title ?pathway_name .

            ?metabolite a wp:Metabolite ;
                dcterms:isPartOf ?pathway .

            { ?metabolite wp:bdbHmdb ?met_xref . BIND("HMDB" AS ?met_db) }
            UNION
            { ?metabolite wp:bdbChEBI ?met_xref . BIND("CHEBI" AS ?met_db) }
            UNION
            { ?metabolite wp:bdbPubChem ?met_xref . BIND("PUBCHEM" AS ?met_db) }
            UNION
            { ?metabolite wp:bdbWikidata ?met_xref . BIND("WIKIDATA" AS ?met_db) }
            UNION
            { ?metabolite wp:bdbCas ?met_xref . BIND("CAS" AS ?met_db) }
            }
            LIMIT %d
            OFFSET %d
            ',
            species_clause,
            as.integer(page_size),
            as.integer((page_idx - 1L) * page_size)
        )

        response <- tryCatch(
            download_base(
                url = endpoint,
                fun = jsonlite::fromJSON,
                post = list(query = sparql_query),
                http_headers = list(
                    Accept = 'application/sparql-results+json'
                ),
                simplifyVector = FALSE
            ),
            error = function(e) {
                msg <- sprintf(
                    'WikiPathways SPARQL request failed (page %d): %s',
                    page_idx,
                    conditionMessage(e)
                )
                log_error(msg)
                stop(msg)
            }
        )

        bindings <- response$results$bindings

        empty <-
            is.null(bindings) ||
            (is.data.frame(bindings) && nrow(bindings) == 0L) ||
            (is.list(bindings) && length(bindings) == 0L)

        if (empty) break

        df_page <-
            tibble(
                pathway_id = extract_binding_value(bindings, 'pathway_id'),
                pathway_name = extract_binding_value(
                    bindings, 'pathway_name_str'
                ),
                met_db = bindings %>%
                    extract_binding_value('met_db') %>%
                    str_trim %>%
                    str_to_upper,
                met_xref_raw = bindings %>%
                    extract_binding_value('met_xref') %>%
                    str_trim %>%
                    str_to_upper
            ) %>%
            filter(
                !is.na(pathway_id),
                !is.na(met_db),
                !is.na(met_xref_raw),
                met_xref_raw != ''
            ) %>%
            mutate(
                met_xref_clean = met_xref_raw %>%
                    str_replace('^HTTPS?://IDENTIFIERS\\.ORG/', '') %>%
                    str_replace('^.*/', '') %>%
                    str_replace('^(CHEBI:)+', 'CHEBI:') %>%
                    str_replace('^(HMDB:)+', 'HMDB:') %>%
                    str_replace('^(WIKIDATA:)+', 'WIKIDATA:') %>%
                    str_replace('^(PUBCHEM:)+', 'PUBCHEM:')
            ) %>%
            mutate(
                met_id = case_when(
                    met_db == 'CHEBI' ~ met_xref_clean %>%
                        str_replace('^(CHEBI:[0-9]+).*$', '\\1'),
                    met_db == 'HMDB' ~ met_xref_clean %>%
                        str_replace('^HMDB:', '') %>%
                        str_replace('^(HMDB[0-9A-Z]+).*$', '\\1'),
                    met_db == 'WIKIDATA' ~ met_xref_clean %>%
                        str_replace('^WIKIDATA:', '') %>%
                        str_replace('^(Q[0-9]+).*$', '\\1'),
                    met_db == 'PUBCHEM' ~
                        normalize_pubchem(met_xref_clean),
                    met_db == 'CAS' ~ met_xref_clean %>%
                        str_replace('^CAS:', ''),
                    TRUE ~ met_xref_clean
                )
            ) %>%
            filter(!is.na(met_id), nzchar(met_id)) %>%
            select(pathway_id, pathway_name, met_id)

        all_pages[[length(all_pages) + 1L]] <- df_page

    }

    all_pages %<>%
        compact %>%
        bind_rows

    if (is.null(all_pages) || nrow(all_pages) == 0L) {
        return(
            tibble(
                pathway_id = character(),
                pathway_name = character(),
                pathway_url = character(),
                metabolites = character()
            )
        )
    }

    pws <-
        wikipathways_pathways(organism = organism) %>%
        select(pathway_id, pathway_url)

    all_pages %>%
        left_join(pws, by = 'pathway_id') %>%
        group_by(pathway_id, pathway_name, pathway_url) %>%
        summarise(
            n_metabolites_in_pathway = n_distinct(met_id),
            metabolites = met_id %>%
                unique %>%
                sort %>%
                paste(collapse = '; '),
            .groups = 'drop'
        )

}
