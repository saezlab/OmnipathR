#!/usr/bin/env Rscript

#' Reactome species table
#'
#' @return A tibble with `code` and `name`.
#'
#' @importFrom tibble tibble
#'
#' @noRd
reactome_species_table <- function() {

    tibble(
        code = c(
            "BTA", "CEL", "CFA", "DRE",
            "DDI", "DME", "GGA", "HSA",
            "MMU", "MTU", "PFA", "RNO",
            "SCE", "SPO", "SSC", "XTR"
        ),
        name = c(
            "Bos taurus",
            "Caenorhabditis elegans",
            "Canis familiaris",
            "Danio rerio",
            "Dictyostelium discoideum",
            "Drosophila melanogaster",
            "Gallus gallus",
            "Homo sapiens",
            "Mus musculus",
            "Mycobacterium tuberculosis",
            "Plasmodium falciparum",
            "Rattus norvegicus",
            "Saccharomyces cerevisiae",
            "Schizosaccharomyces pombe",
            "Sus scrofa",
            "Xenopus tropicalis"
        )
    )

}

#' Normalize Reactome species input
#'
#' @param species Character: species code or full species name.
#'
#' @return List with `code`, `name` and `valid`.
#'
#' @importFrom logger log_warn
#' @importFrom stringr str_trim
#'
#' @noRd
normalize_reactome_species <- function(species) {

    if (is.null(species)) {
        return(list(code = NULL, name = NULL, valid = TRUE))
    }

    if (length(species) == 0) {
        log_warn(sprintf("Reactome species is empty."))
        return(list(code = NA_character_, name = NA_character_, valid = FALSE))
    }

    if (length(species) > 1) {
        log_warn(sprintf(
            "Multiple Reactome species values provided; using first: %s",
            species[[1]]
        ))
        species <- species[[1]]
    }

    species <- str_trim(as.character(species))
    if (is.na(species) || species == "") {
        log_warn(sprintf("Reactome species is empty."))
        return(list(code = NA_character_, name = NA_character_, valid = FALSE))
    }

    table <- reactome_species_table()
    species_upper <- toupper(species)

    idx <- match(species_upper, table$code)
    if (!is.na(idx)) {
        return(list(code = table$code[idx], name = table$name[idx], valid = TRUE))
    }

    idx <- match(tolower(species), tolower(table$name))
    if (!is.na(idx)) {
        return(list(code = table$code[idx], name = table$name[idx], valid = TRUE))
    }

    log_warn(sprintf("Unknown Reactome species: %s", species))
    list(code = NA_character_, name = NA_character_, valid = FALSE)

}

#' Empty Reactome relations table
#'
#' @return A tibble with columns `Parent`, `Child`.
#'
#' @importFrom tibble tibble
#'
#' @noRd
empty_reactome_relations <- function() {

    tibble(Parent = character(), Child = character())

}

#' Empty Reactome pathways table
#'
#' @return A tibble with columns `pathway_id`, `pathway_name`, `species`.
#'
#' @importFrom tibble tibble
#'
#' @noRd
empty_reactome_pathways <- function() {

    tibble(
        pathway_id = character(),
        pathway_name = character(),
        species = character()
    )

}

#' Empty Reactome ChEBI mapping table
#'
#' @return A tibble with columns `pathway_id`, `pathway_name`, `pathway_url`,
#'   `species`, `chebi_ids`.
#'
#' @importFrom tibble tibble
#'
#' @noRd
empty_reactome_chebi_map <- function() {
    
    tibble(
        pathway_id = character(),
        pathway_name = character(),
        pathway_url = character(),
        species = character(),
        chebi_ids = character()
    )
    
}

#' Reactome pathway relations
#'
#' Retrieves parent-child relationships between Reactome pathways and
#' optionally filters by species.
#'
#' @param species Character: species code (e.g. "HSA") or full name
#'   (e.g. "Homo sapiens"). If `NULL`, returns all species.
#'
#' @return A tibble with columns `Parent`, `Child`.
#'
#' @importFrom dplyr filter transmute
#' @importFrom stringr fixed str_detect
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' .slow_doctest()
#' \dontrun{
#' relations <- get_reactome_pathway_relations("HSA")
#' head(relations)
#' }
get_reactome_pathway_relations <- function(species = NULL) {

    relations <- generic_downloader(
        url_key = "reactome_pathway_relations",
        reader_param = list(col_names = c("Parent", "Child")),
        resource = "Reactome"
    )

    if (is.null(species)) {
        return(relations %>% transmute(Parent, Child))
    }

    sp <- normalize_reactome_species(species)
    if (!sp$valid) {
        return(empty_reactome_relations())
    }

    pattern <- sprintf("R-%s-", sp$code)
    relations %>%
        filter(
            str_detect(Parent, fixed(pattern)),
            str_detect(Child, fixed(pattern))
        ) %>%
        transmute(Parent, Child)

}

#' Reactome pathways
#'
#' Retrieves Reactome pathway IDs, names and species. Optionally filters by
#' species.
#'
#' @param species Character: species code (e.g. "HSA") or full name
#'   (e.g. "Homo sapiens"). If `NULL`, returns all species.
#'
#' @return A tibble with columns `pathway_id`, `pathway_name`, `species`.
#'
#' @importFrom dplyr filter transmute
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' .slow_doctest()
#' \dontrun{
#' pathways <- get_reactome_pathways("Homo sapiens")
#' head(pathways)
#' }
get_reactome_pathways <- function(species = NULL) {

    pathways <- generic_downloader(
        url_key = "reactome_pathways",
        reader_param = list(
            col_names = c("pathway_id", "pathway_name", "species")
        ),
        resource = "Reactome"
    )

    if (is.null(species)) {
        return(pathways %>% transmute(pathway_id, pathway_name, species))
    }

    sp <- normalize_reactome_species(species)
    if (!sp$valid) {
        return(empty_reactome_pathways())
    }

    pathways %>%
        filter(species == sp$name) %>%
        transmute(pathway_id, pathway_name, species)

}

#' Reactome pathway participants (ChEBI mapping)
#'
#' Retrieves ChEBI identifiers mapped to Reactome pathways and optionally
#' filters by species and/or pathway IDs. Returns a per-pathway list of ChEBI
#' IDs.
#'
#' @param pathway_ids Character: Reactome pathway IDs to keep. If `NULL`, keeps
#'   all pathways (subject to species filter).
#' @param species Character: species code (e.g. "HSA") or full name
#'   (e.g. "Homo sapiens"). If `NULL`, returns all species.
#' @param out_path Character: optional CSV output path for the main table.
#'
#' @return A tibble with columns `pathway_id`, `pathway_name`, `pathway_url`,
#'   `species`, `chebi_ids`. When writing to CSV, the list
#'   column is serialized by `write.csv`.
#'
#' @importFrom dplyr filter group_by summarise
#' @importFrom logger log_warn
#' @importFrom magrittr %>%
#' @importFrom readr cols col_character
#' @importFrom stringr str_detect
#' @importFrom utils write.csv
#'
#' @export
#' @examples
#' .slow_doctest()
#' \dontrun{
#' tbl <- get_pathway_participants(species = "HSA")
#' head(tbl)
#' }
get_pathway_participants <- function(
    pathway_ids = NULL,
    species = NULL,
    out_path = NULL
) {

    mapping <- generic_downloader(
        url_key = "reactome_chebi2pathways",
        reader_param = list(
            col_names = c(
                "chebi_id",
                "pathway_id",
                "pathway_url",
                "pathway_name",
                "evidence",
                "species"
            ),
            col_types = cols(.default = col_character())
        ),
        resource = "Reactome"
    )

    if (!is.null(species)) {
        sp <- normalize_reactome_species(species)
        if (!sp$valid) {
            return(empty_reactome_chebi_map())
        }
        mapping <- mapping %>% filter(species == sp$name)
    }

    if (!is.null(pathway_ids)) {
        pathway_ids <- as.character(pathway_ids)
        invalid_ids <- pathway_ids[
            is.na(pathway_ids) |
            !nzchar(pathway_ids) |
            !str_detect(pathway_ids, "^R-[A-Z]{3}-")
        ]
        if (length(invalid_ids) > 0) {
            unique(invalid_ids) %>%
                sort() %>%
                lapply(function(id) {
                    log_warn(sprintf("Invalid Reactome pathway ID: %s", id))
                })
        }
        pathway_ids <- setdiff(pathway_ids, invalid_ids)
        mapping <- mapping %>% filter(pathway_id %in% pathway_ids)
    }

    mapping <- mapping %>%
        filter(!is.na(chebi_id) & chebi_id != "") %>%
        mutate(
            chebi_id = str_trim(chebi_id),
            chebi_id = ifelse(
                str_detect(chebi_id, "^CHEBI:"),
                chebi_id,
                paste0("CHEBI:", chebi_id)
            )
        )

    out <- mapping %>%
        group_by(pathway_id, pathway_name, pathway_url, species) %>%
        summarise(
            chebi_ids = paste0(sort(unique(chebi_id)), collapse = ", "),
            .groups = "drop"
        )
    
    if (nrow(out) == 0) {
        out <- empty_reactome_chebi_map()
    }

    if (!is.null(out_path)) {
        write.csv(out, out_path, row.names = FALSE, fileEncoding = "UTF-8")
    }

    out

}

#' Reactome ChEBI mapping per pathway
#'
#' Retrieves ChEBI identifiers mapped to Reactome pathways and returns a
#' per-pathway list of ChEBI IDs.
#'
#' @param species Character: species code (e.g. "HSA") or full name
#'   (e.g. "Homo sapiens"). If `NULL`, returns all species.
#' @param pathway_ids Character: optional list of pathway IDs to query.
#' @param out_path Character: optional CSV output path for the main table.
#'
#' @return A data.frame with columns `pathway_id`, `pathway_name`, `pathway_url`,
#'   `species`, `chebi_ids`.
#'
#' @export
#' @examples
#' .slow_doctest()
#' \dontrun{
#' df <- get_reactome(species = "Homo sapiens")
#' head(df)
#' }
get_reactome <- function(
    species = "Homo sapiens",
    pathway_ids = NULL,
    out_path = NULL
) {

    get_pathway_participants(
        pathway_ids = pathway_ids,
        species = species,
        out_path = out_path
    ) %>%
    as.data.frame()
    
}
