#!/usr/bin/env Rscript

#' List WikiPathways pathways
#'
#' Retrieves the full pathway list from WikiPathways and optionally filters by
#' species.
#'
#' @param species Character: filter by exact species name (e.g. "Homo sapiens").
#'   If `NULL`, returns all species.
#'
#' @return A tibble with columns `pathway_id`, `pathway_name`, `pathway_url`,
#'   `organism_species`.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter transmute
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom httr2 request req_perform resp_body_json
#' @importFrom rlang %||%
#'
#' @export
#' @examples
#' .slow_doctest()
#' \dontrun{
#' pathways <- get_wikipathways_pathways(species = "Homo sapiens")
#' head(pathways)
#' }
get_wikipathways_pathways <- function(species = NULL) {

    data <- generic_downloader(
        url_key = "wikipathways_list",
        reader = curl_read_json,
        reader_param = list(simplifyVector = FALSE),
        resource = "WikiPathways"
    )
    
    ## devel version with direct request instead of generic_downloader
    # data <- request("https://www.wikipathways.org/json/listPathways.json") %>%
    #     req_perform() %>%
    #     resp_body_json()
    
    
    organisms <- data$organisms
    if (is.null(organisms) || length(organisms) == 0) {
        return(
            tibble(
                pathway_id = character(),
                pathway_name = character(),
                pathway_url = character(),
                organism_species = character()
            )
        )
    }

    pathways_df <- map_dfr(organisms, function(org) {

        pws <- org$pathways
        if (is.null(pws) || length(pws) == 0) return(tibble())

        map_dfr(pws, function(pw) {
            tibble(
                pathway_id = pw$id %||% NA_character_,
                pathway_name = pw$name %||% NA_character_,
                pathway_url = pw$url %||% NA_character_,
                organism_species = pw$species %||% NA_character_
            )
        })

    })

    if (!is.null(species)) {
        pathways_df <- pathways_df %>% filter(organism_species == species)
    }

    pathways_df %>%
        transmute(pathway_id, pathway_name, pathway_url, organism_species)

}

#' Normalize metabolite identifiers
#'
#' @param db Character: database name.
#' @param id Character: identifier.
#'
#' @return Character: `db:id` or `NA_character_` if missing.
#'
#' @importFrom stringr fixed str_detect str_trim
#' @importFrom httr2 resp_body_json resp_body_string
#'
#' @noRd
norm_met_id <- function(db, id) {

    if (is.null(db) || is.null(id) || is.na(db) || is.na(id) || db == "" || id == "") {
        return(NA_character_)
    }

    db <- str_trim(as.character(db))
    id <- str_trim(as.character(id))

    if (str_detect(id, fixed("HMDB")) |
        str_detect(id, fixed("CHEBI")) |
        str_detect(id, fixed("ChEBI")) |
        str_detect(id, fixed("PubChem-compound")) |
        str_detect(id, fixed("Wikidata")) |
        str_detect(id, fixed("CAS")) |
        str_detect(id, fixed(":"))
    ){ 
        return(toupper(id))
    }

    paste0(toupper(db), ":", id)

}

#' Download GPML and return XML document
#'
#' @param wpid Character: WikiPathways ID.
#'
#' @return An `xml_document`.
#'
#' @importFrom magrittr %>%
#' @importFrom httr2 request req_perform resp_body_string
#' @importFrom xml2 read_xml
#'
#' @noRd
fetch_gpml <- function(wpid) {

    path <- download_to_cache(
        url_key = "wikipathways_gpml",
        url_param = list(wpid, wpid),
        ext = "gpml"
    )
    
    # ## devel version with direct request instead of generic_downloader
    # url = sprintf("https://www.wikipathways.org/wikipathways-assets/pathways/%s/%s.gpml", wpid, wpid)
    # resp <- request(url) %>%
    #     req_perform()

    read_xml(resp_body_string(resp))
    
    # read_xml(path)

}

#' Build metabolite-pathway table
#'
#' @param pathways_tbl Tibble: output of `get_wikipathways_pathways`.
#' @param out_path Character: optional CSV output path for the main table.
#' @param failures_path Character: optional CSV output path for failures.
#'
#' @return A tibble with columns `pathway_id`, `pathway_name`, `pathway_url`,
#'   `metabolites`. If any pathway failed, the failures tibble is attached as
#'   `attr(out, "failures")`.
#'
#' @importFrom dplyr bind_rows
#' @importFrom logger log_warn
#' @importFrom stringr str_pad
#' @importFrom tibble tibble
#' @importFrom utils flush.console write.csv
#' @importFrom xml2 xml_attr xml_find_all
#'
#' @noRd
get_metabolite_pathway_table <- function(
    pathways_tbl,
    out_path = NULL,
    failures_path = NULL
) {

    n <- nrow(pathways_tbl)
    if (n == 0) {
        empty <- tibble(
            pathway_id = character(),
            pathway_name = character(),
            pathway_url = character(),
            metabolites = character()
        )

        if (!is.null(out_path)) {
            write.csv(empty, out_path, row.names = FALSE, fileEncoding = "UTF-8")
        }

        return(empty)
    }

    rows <- vector("list", n)
    failures <- list()

    for (i in seq_len(n)) {

        wpid <- pathways_tbl$pathway_id[i]

        msg <- sprintf("[%d/%d] %s | failures: %d", i, n, wpid, length(failures))
        cat("\r", str_pad(msg, width = 60, side = "right"), sep = "")
        flush.console()

        tryCatch({

            doc <- fetch_gpml(wpid)

            dnodes <- xml_find_all(doc, "//*[local-name()='DataNode']")

            met_ids <- character(0)

            for (dn in dnodes) {

                dn_type <- xml_attr(dn, "Type")
                if (is.na(dn_type) || dn_type != "Metabolite") next

                xrefs <- xml_find_all(dn, "./*[local-name()='Xref']")
                if (length(xrefs) == 0) next

                for (xr in xrefs) {

                    db <- xml_attr(xr, "Database")
                    id <- xml_attr(xr, "ID")
                    met_id <- norm_met_id(db, id)

                    if (!is.na(met_id) && nzchar(met_id)) {
                        met_ids <- c(met_ids, met_id)
                    }

                }

            }

            met_ids <- sort(unique(met_ids))

            rows[[i]] <- tibble(
                pathway_id = wpid,
                pathway_name = pathways_tbl$pathway_name[i],
                pathway_url = pathways_tbl$pathway_url[i],
                metabolites = paste(met_ids, collapse = "; ")
            )

        }, error = function(e) {

            failures[[length(failures) + 1]] <<- tibble(
                pathway_id = wpid,
                pathway_name = pathways_tbl$pathway_name[i],
                error = conditionMessage(e)
            )
            rows[[i]] <<- NULL
            print(conditionMessage(e))
            print(e)
        })

        # Sys.sleep(0.1)
        
    }

    cat("\n")

    out <- bind_rows(rows)

    failures_tbl <- NULL
    if (length(failures) > 0) {
        failures_tbl <- bind_rows(failures)
        attr(out, "failures") <- failures_tbl
        msg <- sprintf("WikiPathways: %d pathway(s) failed.", nrow(failures_tbl))
        log_warn(msg)
    }

    if (!is.null(out_path)) {
        write.csv(out, out_path, row.names = FALSE, fileEncoding = "UTF-8")
    }

    if (!is.null(failures_path) && !is.null(failures_tbl)) {
        write.csv(
            failures_tbl,
            failures_path,
            row.names = FALSE,
            fileEncoding = "UTF-8"
        )
    }

    out

}

#' Fetch WikiPathways metabolites per pathway
#'
#' Retrieves WikiPathways pathways (optionally filtered by species) and parses
#' the GPML to extract metabolite identifiers.
#'
#' @param species Character: filter by exact species name (e.g. "Homo sapiens").
#' @param out_path Character: optional CSV output path for the main table.
#' @param failures_path Character: optional CSV output path for failures.
#'
#' @return A tibble with columns `pathway_id`, `pathway_name`, `pathway_url`,
#'   `metabolites`. If any pathway failed, the failures tibble is attached as
#'   `attr(out, "failures")`.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' 
#' @export
#' @examples
#' .slow_doctest()
#' \dontrun{
#' df <- get_wikipathways(species = "Homo sapiens")
#' head(df)
#' }
get_wikipathways <- function(
    species = "Homo sapiens",
    out_path = NULL,
    failures_path = NULL
) {

    pathways <- get_wikipathways_pathways(species = species)
    metabolite_tbl <- get_metabolite_pathway_table(
        pathways_tbl = pathways,
        out_path = out_path,
        failures_path = failures_path
    )
    
    # clean metabolite_table
    metabolite_tbl_clean <- metabolite_tbl %>%
        filter(!is.na(metabolites) & trimws(metabolites) != "")
        
    metabolite_tbl_clean

}

#' Retrieve WikiPathways metabolite members via SPARQL
#'
#' Queries the WikiPathways SPARQL endpoint to retrieve all metabolite members
#' for pathways of a given species. Results are paginated and combined into a
#' single tibble.
#'
#' Metabolite identifiers are normalised to the following formats:
#' \itemize{
#'   \item HMDB: \code{HMDB0001049}
#'   \item ChEBI: \code{CHEBI:16015}
#'   \item PubChem: \code{CID7098621}
#'   \item Wikidata: \code{Q715317}
#'   \item CAS: numeric or hyphenated form without prefix
#' }
#'
#' The output contains one row per pathway, with metabolite identifiers
#' collapsed into a semicolon separated string. The number of unique
#' metabolites per pathway is also reported.
#'
#' @param species Character scalar. Exact organism name as used by
#'   WikiPathways, for example \code{"Homo sapiens"}. If \code{NULL},
#'   do not filter by species.
#' @param page_size Integer. Number of records retrieved per SPARQL page.
#'   Larger values reduce the number of HTTP requests but may increase
#'   endpoint load. Default is \code{50000}.
#' @param max_pages Integer. Maximum number of pages to retrieve before
#'   stopping. Acts as a safety limit. Default is \code{200}.
#' @param max_retries Integer. Number of retry attempts per page in case of
#'   transient HTTP errors. Default is \code{4}.
#' @param sleep_base Numeric. Base sleep time in seconds for exponential
#'   backoff between retries. Default is \code{1}.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{pathway_id}{WikiPathways identifier, for example \code{"WP254"}}
#'   \item{pathway_name}{Pathway title}
#'   \item{pathway_url}{Direct URL to the pathway}
#'   \item{n_metabolites_in_pathway}{Number of unique metabolite identifiers}
#'   \item{metabolites}{Semicolon separated metabolite identifiers}
#' }
#'
#' @details
#' This implementation avoids downloading individual GPML files and instead
#' performs a single paginated SPARQL query against the public
#' WikiPathways endpoint. Identifier IRIs from identifiers.org are normalised
#' and redundant prefixes are removed.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr filter select bind_rows left_join group_by summarise n_distinct
#' 
#' @examples
#' \dontrun{
#' df <- get_wikipathways_metabolites_sparql(
#'   species = "Homo sapiens"
#' )
#' head(df)
#' }
#'
#' @export
get_wikipathways_metabolites_sparql <- function(
    species = "Homo sapiens",
    page_size = 50000,
    max_pages = 200,
    max_retries = 4,
    sleep_base = 1
) {
    if (!is.null(species)) {
        stopifnot(
            is.character(species), 
            length(species) == 1, 
            nzchar(species)
        )
    }

    endpoint <- url_parser(url_key = "wikipathways_sparql")

    extract_binding_value <- function(bindings, field) {
        ## parses sql field data
        
        if (is.data.frame(bindings)) {
            if (
                field %in% names(bindings) &&
                is.data.frame(bindings[[field]]) &&
                "value" %in% names(bindings[[field]])
            ) {
                return(bindings[[field]]$value)
            }
            
            return(
                rep(
                    NA_character_, 
                    nrow(bindings)
                )
            )
        }

        vapply(
            bindings,
            function(binding) {
                if (!is.null(binding[[field]]$value)) {
                    as.character(binding[[field]]$value)
                } else {
                    NA_character_
                }
            },
            character(1)
        )
    }

    normalize_pubchem <- function(xref) {
        # normalized pubchem meabolite id strings
        
        xref %>%
            str_replace(
                "^PUBCHEM:", ""
            ) %>%
            str_replace(
                "^PUBCHEM\\.COMPOUND:", ""
            ) %>%
            {
                ifelse(
                    grepl(
                        "^[0-9]+$", 
                        .
                    ), 
                    paste0(
                        "CID", 
                        .
                    ), 
                    .
                )
            } %>%
            str_replace(
                "^(CID[0-9]+).*$", "\\1"
            )
    }

    # optionally add species to sql query if provided by user
    if (is.null(species)) {
        species_clause <- ""
    } else {
        species_clause <- 
            sprintf(
                '           wp:organismName "%s" ;\n',
                gsub(
                    '"', '\\"', species, fixed = TRUE
                )
            )
    }
    
    
    all_pages <- list()
    
    # main loop over all pages
    # this code does not retrieve everything at once but rather iterates over
    # multiple pages, calling the API multiple times, doing:
    # - construct current SQL query
    # - retrieve response data
    # - parse response data
    # -> iteratively creating response object from multiple SQL queries
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

        # retrieve json response with sql query from endpoint
        response <- tryCatch(
            download_base(
                url = endpoint,
                fun = jsonlite::fromJSON,
                post = list(query = sparql_query),
                http_headers = list(
                    Accept = "application/sparql-results+json"
                ),
                simplifyVector = FALSE
            ),
            error = function(e) {
                msg <- sprintf(
                    "WikiPathways SPARQL request failed (page %d): %s",
                    page_idx,
                    conditionMessage(e)
                )
                logger::log_error(msg)
                stop(msg)
            }
        )

        bindings <- response$results$bindings

        empty <- is.null(bindings) ||
            (is.data.frame(bindings) && nrow(bindings) == 0) ||
            (is.list(bindings) && length(bindings) == 0)
        if (empty) break
        
        # extract all columns
        pathway_id <- extract_binding_value(bindings, "pathway_id")
        pathway_name <- extract_binding_value(bindings, "pathway_name_str")
        met_db <- extract_binding_value(bindings, "met_db")
        met_xref <- extract_binding_value(bindings, "met_xref")
        
        # construct dataframe from current SQL query response
        df_page <- tibble(
            pathway_id = pathway_id,
            pathway_name = pathway_name,
            met_db = toupper(trimws(met_db)),
            met_xref_raw = toupper(trimws(met_xref))
        ) %>%
        filter(
            !is.na(pathway_id),
            !is.na(met_db),
            !is.na(met_xref_raw),
            met_xref_raw != ""
        ) %>%
        mutate(
            met_xref_clean = met_xref_raw %>%
                str_replace(
                    "^HTTPS?://IDENTIFIERS\\.ORG/",
                    ""
                ) %>%
                str_replace("^.*/", "") %>%
                str_replace("^(CHEBI:)+", "CHEBI:") %>%
                str_replace("^(HMDB:)+", "HMDB:") %>%
                str_replace("^(WIKIDATA:)+", "WIKIDATA:") %>%
                str_replace("^(PUBCHEM:)+", "PUBCHEM:")
        ) %>%
        mutate(
            met_id = case_when(
                met_db == "CHEBI" ~ str_replace(
                    met_xref_clean,
                    "^(CHEBI:[0-9]+).*$",
                    "\\1"
                ),
                met_db == "HMDB" ~ met_xref_clean %>%
                    str_replace("^HMDB:", "") %>%
                    str_replace("^(HMDB[0-9A-Z]+).*$", "\\1"),
                met_db == "WIKIDATA" ~ met_xref_clean %>%
                    str_replace("^WIKIDATA:", "") %>%
                    str_replace("^(Q[0-9]+).*$", "\\1"),
                met_db == "PUBCHEM" ~ normalize_pubchem(met_xref_clean),
                met_db == "CAS" ~ str_replace(
                    met_xref_clean,
                    "^CAS:",
                    ""
                ),
                TRUE ~ met_xref_clean
            )
        ) %>%
        filter(!is.na(met_id), nzchar(met_id)) %>%
        select(pathway_id, pathway_name, met_id)

        
        # append parsed results of current query to all_pages
        all_pages[[length(all_pages) + 1L]] <- df_page
    }

    
    all_pages <- all_pages %>%
        Filter(f = Negate(is.null), .) %>%
        bind_rows()

    if (is.null(all_pages) || nrow(all_pages) == 0) {
        return(
            tibble(
                pathway_id = character(),
                pathway_name = character(),
                pathway_url = character(),
                metabolites = character()
            )
        )
    }
    
    pws <- get_wikipathways_pathways(species = species) %>%
        select(pathway_id, pathway_url)

    # return full result table
    all_pages %>%
        left_join(pws, by = "pathway_id") %>%
        group_by(pathway_id, pathway_name, pathway_url) %>%
        summarise(
            n_metabolites_in_pathway = n_distinct(met_id),
            metabolites = paste(sort(unique(met_id)), collapse = "; "),
            .groups = "drop"
        ) %>%
        as.data.frame()
    
}


