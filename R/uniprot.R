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


#' Retrieves an identifier translation table from the UniProt uploadlists
#' service
#'
#' @param identifiers Character vector of identifiers
#' @param from Type of the identifiers provided. See Details for possible
#'     values.
#' @param to Identifier type to be retrieved from UniProt. See Details for
#'     possible values.
#'
#' @details
#' This function uses the uploadlists service of UniProt to obtain identifier
#' translation tables. The possible values for `from` and `to` are the
#' identifier type abbreviations used in the UniProt API, please refer to
#' the table here: https://www.uniprot.org/help/api_idmapping
#'
#' @importsFrom readr read_tsv cols
#' @importsFrom httr POST
#' @importsFrom magrittr %>%
#' @export
#'
#' @examples
#' uniprot_genesymbol <- uniprot_id_mapping_table(
#'     c('P00533', 'P23771'), 'ID', 'GENENAME'
#' )
uniprot_id_mapping_table <- function(identifiers, from, to){

    POST(
        url = 'https://www.uniprot.org/uploadlists/',
        body = list(
            from = from,
            to = to,
            format = 'tab',
            query = paste(identifiers, collapse = ' ')
        )
    ) %>%
    content(encoding = 'ASCII') %>%
    read_tsv(col_types = cols())

}


#' Translates a column of identifiers in a data frame by creating another
#' column with the target identifiers
#'
#' @param d A data frame
#' @param from_col Name of an existing column in the data frame containing
#'     the identifiers to be translated
#' @param to_col Name for a new column which will contain the target
#'     identifiers
#' @param from Identifier type for `from_col`. See Details for possible
#'     values.
#' @param to Identifier type for `to_col`. See Details for possible values.
#' @param keep_untranslated Keep the records where the source identifier
#'     could not be translated. At these records the target identifier will
#'     be NA.
#'
#' @details
#' This function uses the uploadlists service of UniProt to obtain identifier
#' translation tables. The possible values for `from` and `to` are the
#' identifier type abbreviations used in the UniProt API, please refer to
#' the table here: https://www.uniprot.org/help/api_idmapping
#' The mapping between identifiers can be ambiguous. In this case one row
#' in the original data frame yields multiple rows in the returned data
#' frame.
#'
#' @importsFrom rlang !! enquo := quo_text
#' @importsFrom magrittr %>%
#' @importsFrom dplyr pull left_join inner_join rename
#' @export
#'
#' @examples
#' d <- translate_ids(d, uniprot_id, genesymbol, 'ID', 'GENENAME')
#'
#' @seealso \code{\link{uniprot_id_mapping_table}}
translate_ids <- function(
    d, from_col, to_col, from, to,
    keep_untranslated = TRUE
){

    from_col <- enquo(from_col)
    to_col <- enquo(to_col)

    join_method <- `if`(keep_untranslated, left_join, inner_join)

    translation_table <- d %>%
    pull(!!from_col) %>%
    unique() %>%
    uniprot_id_mapping_table(from = from, to = to)

    d %>%
    join_method(
        translation_table,
        by = 'From' %>% setNames(from_col %>% quo_text)
    ) %>%
    mutate(!!to_col := To) %>%
    select(-To)

}


#' Retrieves a table from UniProt with all proteins for a certain organism
#'
#' @param fields Character vector of fields as defined by UniProt. For
#'     possible values please refer to
#'     https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames
#' @param reviewed Retrieve only reviewed (TRUE), only unreviewed (FALSE) or
#'     both (NULL).
#' @param organism Integer, NCBI Taxonomy ID of the organism (by default
#'     9606 for human).
#'
#' @importsFrom readr read_tsv cols
#' @importsFrom magrittr %>%
#' @export
#'
#' @examples
#' human_swissprot_ac <- all_uniprots(fields = 'entry name')
all_uniprots <- function(fields = 'id', reviewed = TRUE, organism = 9606){

    on.exit(closeAllConnections())

    fields <- fields %>% paste(collapse = ',')
    reviewed <- `if`(
        is.null(reviewed),
        '',
        sprintf(' AND reviewed:%s', `if`(reviewed, 'yes', 'no'))
    )

    'omnipath.all_uniprots_url' %>%
    options() %>%
    as.character() %>%
    sprintf(fields, organism, reviewed) %>%
    URLencode() %>%
    read_tsv(col_types = cols(), progress = FALSE)

}


#' Creates an ID translation table from UniProt data
#'
#' @param to Target ID type. See Details for possible values.
#' @param from Source ID type. See Details for possible values.
#' @param reviewed Retrieve only reviewed (TRUE), only unreviewed (FALSE) or
#'     both (NULL).
#' @param organism Integer, NCBI Taxonomy ID of the organism (by default
#'     9606 for human).
#'
#' @details
#' For both source and target ID type, this function accepts column codes
#' used by UniProt and some simple shortcuts defined here. For the UniProt
#' codes please refer to
#' https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames
#' The shortcuts are entrez, genesymbol, genesymbol_syn (synonym gene
#' symbols), hgnc, embl, refseqp (RefSeq protein), enst (Ensembl transcript),
#' uniprot_entry (UniProt AC, e.g. EGFR_HUMAN), protein_name (full name of
#' the protein).
#'
#' @importsFrom magrittr %>%
#' @importsFrom dplyr mutate rename
#' @importsFrom tidyr separate_rows
#' @importsFrom rlang !! enquo
#' @export
#'
#' @examples
#' uniprot_entrez <- all_uniprots_id_mapping_table(to = 'entrez')
all_uniprots_id_mapping_table <- function(
    to,
    from = 'id',
    reviewed = TRUE,
    organism = 9606
){

    id_types = list(
        entrez = c('database', 'GeneID'),
        genesymbol = c('genes', 'PREFERRED'),
        genesymbol_syn = c('genes', 'ALTERNATIVE'),
        hgnc = c('database', 'HGNC'),
        embl = c('database', 'embl'),
        entrez = c('database', 'geneid'),
        refseqp = c('database', 'refseq'),
        enst = c('database', 'ensembl'),
        uniprot_entry = c('entry name', NULL),
        protein_name = c('protein names', NULL)
    )

    get_field_name <- function(label){

        label %>%
        {`if`(
            . %in% names(id_types),
            sprintf(
                '%s%s',
                id_types[[.]][1],
                `if`(
                    is.null(id_types[[.]][2]),
                    '',
                    sprintf('(%s)', id_types[[.]][2])
                )
            ),
            .
        )}

    }

    strip_semicol <- function(v){sub(';$', '', v)}

    to <-
        .nse_ensure_str(!!enquo(to)) %>%
        get_field_name()
    from <-
        .nse_ensure_str(!!enquo(from)) %>%
        get_field_name()

    c(from, to) %>%
    all_uniprots(reviewed = reviewed, organism = organism) %>%
    rename(From = 1, To = 2) %>%
    mutate(
        From = strip_semicol(From),
        To = strip_semicol(To)
    ) %>%
    separate_rows(From, To, sep = ';')

}