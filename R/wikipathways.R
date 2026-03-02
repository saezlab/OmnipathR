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
#' 
#' @importFrom httr2 resp_body_json resp_body_string
#'
#' @export
#' @examples
#' .slow_doctest()
#' \dontrun{
#' pathways <- get_wikipathways_pathways(species = "Homo sapiens")
#' head(pathways)
#' }
get_wikipathways_pathways <- function(species = NULL) {

    # data <- generic_downloader(
    #     url_key = "wikipathways_list",
    #     reader = curl_read_json,
    #     reader_param = list(simplifyVector = FALSE),
    #     resource = "WikiPathways"
    # )
    
    ## devel version with direct request instead of generic_downloader
    data <- request("https://www.wikipathways.org/json/listPathways.json") %>%
        req_perform() %>%
        resp_body_json()
    
    
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
                pathway_id = if_null(pw$id, NA_character_),
                pathway_name = if_null(pw$name, NA_character_),
                pathway_url = if_null(pw$url, NA_character_),
                organism_species = if_null(pw$species, NA_character_)
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
#' 
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
#' @importFrom xml2 read_xml
#'
#' @noRd
fetch_gpml <- function(wpid) {

    # path <- download_to_cache(
    #     url_key = "wikipathways_gpml",
    #     url_param = list(wpid, wpid),
    #     ext = "gpml"
    # )
    
    ## devel version with direct request instead of generic_downloader
    url = sprintf("https://www.wikipathways.org/wikipathways-assets/pathways/%s/%s.gpml", wpid, wpid)
    resp <- request(url) %>%
        req_perform()

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
