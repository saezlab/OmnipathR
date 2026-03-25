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


#' Reactome species codes
#'
#' Named character vector: latin names to three-letter Reactome codes.
#'
#' @noRd
REACTOME_SPECIES <- c(
    `Bos taurus` = 'BTA',
    `Caenorhabditis elegans` = 'CEL',
    `Canis familiaris` = 'CFA',
    `Danio rerio` = 'DRE',
    `Dictyostelium discoideum` = 'DDI',
    `Drosophila melanogaster` = 'DME',
    `Gallus gallus` = 'GGA',
    `Homo sapiens` = 'HSA',
    `Mus musculus` = 'MMU',
    `Mycobacterium tuberculosis` = 'MTU',
    `Plasmodium falciparum` = 'PFA',
    `Rattus norvegicus` = 'RNO',
    `Saccharomyces cerevisiae` = 'SCE',
    `Schizosaccharomyces pombe` = 'SPO',
    `Sus scrofa` = 'SSC',
    `Xenopus tropicalis` = 'XTR'
)


#' Reactome species code from organism identifier
#'
#' @param organism Character or integer: any organism name or identifier
#'     accepted by \code{\link{latin_name}}.
#'
#' @return Character: three-letter Reactome code, or \code{NA} if the
#'     organism is not in Reactome.
#'
#' @importFrom magrittr %>%
#' @importFrom logger log_warn
#'
#' @noRd
reactome_species_code <- function(organism) {

    organism %<>% latin_name

    code <- REACTOME_SPECIES[organism]

    if (is.na(code)) {
        log_warn('Unknown Reactome species: `%s`.', organism)
    }

    code

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
#' @return A tibble with columns `pathway_id`, `pathway_name`,
#'     `pathway_url`, `species`, `chebi_ids`.
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
#' optionally filters by organism.
#'
#' @param organism Character or integer: organism name, common name or
#'     NCBI Taxonomy ID. If \code{NULL}, returns all organisms.
#'
#' @return A tibble with columns \code{Parent}, \code{Child}.
#'
#' @importFrom dplyr filter transmute
#' @importFrom stringr fixed str_detect
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' \dontrun{
#' relations <- reactome_pathway_relations(9606)
#' head(relations)
#' }
reactome_pathway_relations <- function(organism = NULL) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    Parent <- Child <- NULL

    relations <-
        generic_downloader(
            url_key = 'reactome_pathway_relations',
            reader_param = list(col_names = c('Parent', 'Child')),
            resource = 'Reactome'
        )

    if (is.null(organism)) {

        relations %>% transmute(Parent, Child)

    } else {

        code <- reactome_species_code(organism)

        if (is.na(code)) {
            empty_reactome_relations()
        } else {
            pattern <- sprintf('R-%s-', code)
            relations %>%
                filter(
                    str_detect(Parent, fixed(pattern)),
                    str_detect(Child, fixed(pattern))
                ) %>%
                transmute(Parent, Child)
        }

    }

}


#' Reactome pathways
#'
#' Retrieves Reactome pathway IDs, names and species. Optionally filters
#' by organism.
#'
#' @param organism Character or integer: organism name, common name or
#'     NCBI Taxonomy ID. If \code{NULL}, returns all organisms.
#'
#' @return A tibble with columns \code{pathway_id}, \code{pathway_name},
#'     \code{species}.
#'
#' @importFrom dplyr filter transmute
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' \dontrun{
#' pathways <- reactome_pathways(9606)
#' head(pathways)
#' }
reactome_pathways <- function(organism = NULL) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    pathway_id <- pathway_name <- species <- NULL

    pathways <-
        generic_downloader(
            url_key = 'reactome_pathways',
            reader_param = list(
                col_names = c('pathway_id', 'pathway_name', 'species')
            ),
            resource = 'Reactome'
        )

    if (is.null(organism)) {

        pathways %>% transmute(pathway_id, pathway_name, species)

    } else {

        organism %<>% latin_name

        if (is.na(organism)) {
            empty_reactome_pathways()
        } else {
            pathways %>%
                filter(.data$species == organism) %>%
                transmute(pathway_id, pathway_name, species)
        }

    }

}


#' Reactome ChEBI pathway mapping
#'
#' Retrieves ChEBI identifiers mapped to Reactome pathways and optionally
#' filters by organism and/or pathway IDs. Returns a per-pathway list of
#' ChEBI IDs.
#'
#' @param pathway_ids Character: Reactome pathway IDs to keep. If
#'     \code{NULL}, keeps all pathways (subject to organism filter).
#' @param organism Character or integer: organism name, common name or
#'     NCBI Taxonomy ID. If \code{NULL}, returns all organisms.
#' @param out_path Character: optional CSV output path.
#'
#' @return A tibble with columns \code{pathway_id}, \code{pathway_name},
#'     \code{pathway_url}, \code{species}, \code{chebi_ids}.
#'
#' @importFrom dplyr filter group_by summarise mutate if_else
#' @importFrom logger log_warn
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr walk
#' @importFrom readr cols col_character write_csv
#' @importFrom stringr str_detect str_trim
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' \dontrun{
#' tbl <- reactome_chebi_pathways(organism = 9606)
#' head(tbl)
#' }
reactome_chebi_pathways <- function(
    pathway_ids = NULL,
    organism = NULL,
    out_path = NULL
) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    chebi_id <- pathway_id <- pathway_name <- pathway_url <-
        species <- chebi_ids <- NULL

    mapping <-
        generic_downloader(
            url_key = 'reactome_chebi2pathways',
            reader_param = list(
                col_names = c(
                    'chebi_id', 'pathway_id', 'pathway_url',
                    'pathway_name', 'evidence', 'species'
                ),
                col_types = cols(.default = col_character())
            ),
            resource = 'Reactome'
        )

    if (!is.null(organism)) {

        organism %<>% latin_name

        if (is.na(organism)) {
            return(empty_reactome_chebi_map())
        }

        mapping %<>% filter(.data$species == organism)

    }

    if (!is.null(pathway_ids)) {

        pathway_ids %<>% as.character

        invalid_ids <-
            pathway_ids[
                is.na(pathway_ids) |
                !nzchar(pathway_ids) |
                !str_detect(pathway_ids, '^R-[A-Z]{3}-')
            ]

        if (length(invalid_ids) > 0L) {
            invalid_ids %>%
                unique %>%
                sort %>%
                walk(~log_warn('Invalid Reactome pathway ID: %s', .x))
        }

        pathway_ids %<>% setdiff(invalid_ids)
        mapping %<>% filter(pathway_id %in% pathway_ids)

    }

    mapping %<>%
        filter(!is.na(chebi_id), chebi_id != '') %>%
        mutate(
            chebi_id = chebi_id %>%
                str_trim %>%
                {if_else(str_detect(., '^CHEBI:'), ., paste0('CHEBI:', .))}
        ) %>%
        group_by(pathway_id, pathway_name, pathway_url, species) %>%
        summarise(
            chebi_ids = chebi_id %>%
                unique %>%
                sort %>%
                paste(collapse = ', '),
            .groups = 'drop'
        )

    if (nrow(mapping) == 0L) {
        mapping <- empty_reactome_chebi_map()
    }

    if (!is.null(out_path)) {
        write_csv(mapping, out_path)
    }

    mapping

}


#' Reactome ChEBI mappings
#'
#' Retrieves ChEBI identifiers mapped to Reactome pathways and returns a
#' per-pathway list of ChEBI IDs as a data frame.
#'
#' @param organism Character or integer: organism name, common name or
#'     NCBI Taxonomy ID. Defaults to human.
#' @param pathway_ids Character: optional list of pathway IDs to query.
#' @param out_path Character: optional CSV output path.
#'
#' @return A data frame with columns \code{pathway_id},
#'     \code{pathway_name}, \code{pathway_url}, \code{species},
#'     \code{chebi_ids}.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' \dontrun{
#' df <- reactome_chebi(organism = 9606)
#' head(df)
#' }
reactome_chebi <- function(
    organism = 9606L,
    pathway_ids = NULL,
    out_path = NULL
) {

    .slow_doctest()

    reactome_chebi_pathways(
        pathway_ids = pathway_ids,
        organism = organism,
        out_path = out_path
    ) %>%
    as.data.frame

}
