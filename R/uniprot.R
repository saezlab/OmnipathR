#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
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


#' Try how non standard evaluation works
#'
#' @noRd
nse_test <- function(a, ..., b = FALSE){

    print(a)
    print(b)

    print(enquos(...) %>% map(.nse_ensure_str))

}


#' All UniProt ACs for one organism
#'
#' @param organism Character or integer: name or identifier of the organism.
#' @param reviewed Retrieve only reviewed (`TRUE`), only unreviewed (`FALSE`)
#'     or both (`NULL`).
#'
#' @return Character vector of UniProt accession numbers.
#'
#' @examples
#' human_swissprot_acs <- all_uniprot_acs()
#' human_swissprot_acs[1:5]
#' # [1] "P51451" "A6H8Y1" "O60885" "Q9Y3X0" "P22223"
#' length(human_swissprot_acs)
#' # [1] 20376
#' mouse_swissprot_acs <- all_uniprot_acs("mouse")
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom dplyr pull
#' @export
all_uniprot_acs <- function(organism = 9606, reviewed = TRUE){

    .slow_doctest()

    organism %<>% ncbi_taxid
    all_uniprots(organism = organism, reviewed = reviewed) %>% pull(1L)

}


#' A table with all UniProt records
#'
#' Retrieves a table from UniProt with all proteins for a certain organism.
#'
#' @param fields Character vector of fields as defined by UniProt. For
#'     possible values please refer to
#'     \url{https://www.uniprot.org/help/return_fields}
#' @param reviewed Retrieve only reviewed (`TRUE`), only unreviewed (`FALSE`)
#'     or both (`NULL`).
#' @param organism Character or integer: name or identifier of the organism.
#'
#' @return Data frame (tibble) with the requested UniProt entries and fields.
#'
#' @importFrom magrittr %<>% %T>% %>%
#' @importFrom rlang exec !!!
#' @importFrom logger log_trace
#' @export
#'
#' @examples
#' human_swissprot_entries <- all_uniprots(fields = 'id')
#' human_swissprot_entries
#' # # A tibble: 20,396 x 1
#' #    `Entry name`
#' #    <chr>
#' #  1 OR4K3_HUMAN
#' #  2 O52A1_HUMAN
#' #  3 O2AG1_HUMAN
#' #  4 O10S1_HUMAN
#' #  5 O11G2_HUMAN
#' # # . with 20,386 more rows
all_uniprots <- function(
    fields = 'accession',
    reviewed = TRUE,
    organism = 9606L
){

    organism %<>% ncbi_taxid

    fields <- fields %>% paste(collapse = ',')

    log_trace(
        paste0(
            'Loading all UniProt records for organism %d ',
            '(only reviewed: %s); fields: %s'
        ),
        organism, reviewed, fields
    )

    reviewed <- `if`(
        is.null(reviewed),
        '',
        sprintf(' AND reviewed:%s', `if`(reviewed, 'true', 'false'))
    )

    generic_downloader(
        url_key = 'new_uniprot',
        url_param = list(fields, organism, reviewed),
        reader_param = list(progress = FALSE),
        resource = 'UniProt'
    ) %T>%
    load_success()

}

#' UniProt taxonomy data
#'
#' @return A tibble with columns: \code{code}, \code{kingdom}, \code{ncbi_tax_id}, \code{latin_name}, \code{common_name}, \code{synonym}.
#'
#' @examples
#' uniprot_organisms()
#'
#' @importFrom stringr str_detect str_match str_remove str_trim
#' @importFrom purrr map2
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr unnest_wider
#' @export
uniprot_organisms <- function() {

    read_uniprot_taxonomy <- function(con, ...) {

        lines <- readLines(con)
        sep_idx <- which(str_detect(lines, "^_____"))[1]
        lines %<>% .[sep_idx + 1L:length(.)]
        entry_starts <- which(!str_detect(lines, "^\\s"))
        entry_ends <- c(entry_starts[-1L] - 1L, length(lines))

        # Parse each entry
        entries <- map2(
            entry_starts,
            entry_ends,
            function(start, end) {
                entry_lines <- lines[start:end]
                m <- str_match(entry_lines[1],
                "^([A-Z0-9]+)\\s+([A-Z])\\s+(\\d+):\\s+N=(.+)$")
                code <- m[2]
                kingdom <- m[3]
                ncbi_tax_id <- m[4]
                latin_name <- m[5]
                common_name <- NA_character_
                synonym <- NA_character_

                if (length(entry_lines) > 1L) {
                    for (l in entry_lines[-1L]) {
                        if (is.na(l)) {
                            break
                        } else if (str_detect(l, "C=")) {
                            common_name <- str_remove(l, "^\\s*C=") %>% str_trim()
                        } else if (str_detect(l, "S=")) {
                            synonym <- str_remove(l, "^\\s*S=") %>% str_trim()
                        }
                    }
                }

                list(
                    code = code,
                    kingdom = kingdom,
                    ncbi_tax_id = as.integer(ncbi_tax_id),
                    latin_name = latin_name,
                    common_name = common_name,
                    synonym = synonym
                )
            }
        )

        entries %>%
        tibble(entries = .) %>%
        unnest_wider(entries)

    }

    generic_downloader(
        url_key = 'uniprot_organisms',
        reader = read_uniprot_taxonomy,
        resource = 'UniProt'
    )

}

