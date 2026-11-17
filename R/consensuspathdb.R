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


#' Retrieves the ConsensusPathDB network
#'
#' Compiles a table of binary interactions from ConsensusPathDB
#' (\url{http://cpdb.molgen.mpg.de/}). The UniProtKB ACs and Gene Symbols
#' of the participants come from the ConsensusPathDB file.
#'
#' @param complex_max_size Numeric: do not expand complexes with a higher
#'     number of elements than this. ConsensusPathDB does not contain
#'     conventional interactions but lists of participants, which might be
#'     members of complexes. Some records include dozens of participants and
#'     expanding them to binary interactions result thousands, sometimes
#'     hundreds of thousands of interactions from one single record. At the
#'     end, this process consumes >10GB of memory and results rather unusable
#'     data, hence it is recommended to limit the complex sizes at some low
#'     number.
#' @param min_score Numeric: each record in ConsensusPathDB comes with a
#'     confidence score, expressing the amount of evidences. The default
#'     value, a minimum score of 0.9 retains approx. the top 30 percent
#'     of the interactions.
#'
#' @return Data frame (tibble) with interactions.
#'
#' @examples
#' cpdb_data <- consensuspathdb_download(
#'     complex_max_size = 1,
#'     min_score = .99
#' )
#' nrow(cpdb_data)
#' # [1] 63522
#' colnames(cpdb_data)
#' # [1] "databases"  "references" "uniprot_a"    "confidence"   "record_id"
#' # [6] "uniprot_b"  "in_complex" "genesymbol_a" "genesymbol_b"
#'
#' @importFrom readr read_tsv
#' @importFrom tidyr separate_rows
#' @importFrom dplyr mutate select inner_join left_join pull filter rename
#' @importFrom dplyr n
#' @importFrom magrittr %>% %T>%
#' @importFrom stringr str_count
#' @export
consensuspathdb_download <- function(
    complex_max_size = 4,
    min_score = .9
){

    .slow_doctest()

    # NSE vs. R CMD check workaround
    uniprots <- confidence <- uniprot_a <- uniprot_b <- record_id <-
        in_complex <- genesymbol <- databases <- references <- single <- NULL

    cpdb_raw <- consensuspathdb_raw_table()

    uniprot_genesymbol <- cpdb_raw %>% consensuspathdb_genesymbols

    cpdb_raw %>%
    filter(
        str_count(uniprots, ',') <= complex_max_size,
        confidence >= min_score
    ) %>%
    mutate(
        record_id = seq(n()),
        uniprot_b = uniprots,
        # more than two participants
        in_complex = str_count(uniprots, ',') > 1L,
        single = !grepl(',', uniprots, fixed = TRUE)
    ) %>%
    rename(uniprot_a = uniprots) %>%
    select(
        databases, references, uniprot_a, confidence, record_id, uniprot_b,
        in_complex, single
    ) %>%
    separate_rows(uniprot_a, sep = ',') %>%
    separate_rows(uniprot_b, sep = ',') %>%
    separate_rows(uniprot_a, sep = '\\.') %>%
    separate_rows(uniprot_b, sep = '\\.') %>%
    # self pairs only from the records of one participant
    filter(single | uniprot_a != uniprot_b) %>%
    select(-single) %>%
    left_join(uniprot_genesymbol, by = c('uniprot_a' = 'uniprot')) %>%
    rename(genesymbol_a = genesymbol) %>%
    left_join(uniprot_genesymbol, by = c('uniprot_b' = 'uniprot')) %>%
    rename(genesymbol_b = genesymbol) %>%
    copy_source_attrs(cpdb_raw) %T>%
    load_success()

}


#' Downloads interaction data from ConsensusPathDB
#'
#' @return Data frame (tibble) with interactions.
#'
#' @examples
#' cpdb_raw <- consensuspathdb_raw_table()
#'
#' @importFrom readr read_tsv cols
#' @importFrom magrittr %>%
#' @importFrom logger log_warn
#' @export
consensuspathdb_raw_table <- function(){

    .slow_doctest()

    col_names <- c(
        'databases',
        'references',
        'uniprot_entries',
        'uniprots',
        'genesymbols',
        'hgncs',
        'entrezs',
        'ensembls',
        'confidence'
    )

    download <- function() {

        'cpdb' %>%
        generic_downloader(
            reader_param = list(
                col_names = col_names,
                col_types = 'ccccccccd',
                skip = 2,
                progress = FALSE
            ),
            resource = 'ConsensusPathDB'
        )

    }

    result <- download()

    if (!identical(colnames(result), col_names)) {

        # cached by an older version, read with an outdated column layout
        log_warn(
            'ConsensusPathDB table in the cache has outdated columns, ',
            'downloading it again.'
        )
        omnipath_cache_remove(url = url_parser(url_key = 'cpdb')[1L])
        result <- download()

    }

    result

}


#' UniProt ACs and Gene Symbols from the ConsensusPathDB table
#'
#' The participant columns of the ConsensusPathDB table are parallel lists,
#' separated by commas. A participant can be a group of alternative
#' proteins, separated by dots, but the alternatives do not always come in
#' the same order in the UniProt and the Gene Symbol columns. Hence we take
#' the pairs only from the single protein participants.
#'
#' @param cpdb_raw The ConsensusPathDB table, as returned by
#'     \code{consensuspathdb_raw_table}.
#'
#' @return Data frame (tibble) with columns `uniprot` and `genesymbol`, one
#'     row for each UniProt AC.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter count arrange desc distinct select
#' @importFrom tibble tibble
#' @noRd
consensuspathdb_genesymbols <- function(cpdb_raw) {

    # NSE vs. R CMD check workaround
    uniprot <- genesymbol <- n <- NULL

    uniprots <- strsplit(cpdb_raw$uniprots, ',', fixed = TRUE)
    genesymbols <- strsplit(cpdb_raw$genesymbols, ',', fixed = TRUE)
    aligned <- lengths(uniprots) == lengths(genesymbols)

    tibble(
        uniprot = unlist(uniprots[aligned]),
        genesymbol = unlist(genesymbols[aligned])
    ) %>%
    filter(
        !is.na(uniprot) & !is.na(genesymbol) &
        !grepl('.', uniprot, fixed = TRUE) &
        !grepl('.', genesymbol, fixed = TRUE)
    ) %>%
    count(uniprot, genesymbol) %>%
    arrange(desc(n)) %>%
    distinct(uniprot, .keep_all = TRUE) %>%
    select(uniprot, genesymbol)

}
