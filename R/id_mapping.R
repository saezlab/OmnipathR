#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2023
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


#' ID translation data from UniProt Uploadlists
#'
#' Retrieves an identifier translation table from the UniProt Uploadlists
#' service.
#'
#' @param identifiers Character vector of identifiers
#' @param from Character or symbol: type of the identifiers provided.
#'     See Details for possible values.
#' @param to Character or symbol: identifier type to be retrieved from
#'     UniProt. See Details for possible values.
#' @param chunk_size Integer: query the identifiers in chunks of this size.
#'     If you are experiencing download failures, try lower values.
#'
#' @return A data frame (tibble) with columns `From` and `To`, the
#'     identifiers provided and the corresponding target IDs, respectively.
#'
#' @details
#' This function uses the uploadlists service of UniProt to obtain identifier
#' translation tables. The possible values for `from` and `to` are the
#' identifier type abbreviations used in the UniProt API, please refer to
#' the table here: \url{https://www.uniprot.org/help/api_idmapping} or
#' the table of synonyms supported by the current package:
#' \code{\link{translate_ids}}.
#' Note: if the number of identifiers is larger than the chunk size the log
#' message about the cache origin is not guaranteed to be correct (most
#' of the times it is still correct).
#'
#' @importFrom magrittr %<>% %>% %T>%
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom rlang !!
#' @export
#'
#' @examples
#' uniprot_genesymbol <- uniprot_id_mapping_table(
#'     c('P00533', 'P23771'), uniprot, genesymbol
#' )
#' uniprot_genesymbol
#' # # A tibble: 2 x 2
#' #   From   To
#' #   <chr>  <chr>
#' # 1 P00533 EGFR
#' # 2 P23771 GATA3
#'
#' @seealso \code{\link{translate_ids}}
uniprot_id_mapping_table <- function(
    identifiers,
    from,
    to,
    chunk_size = NULL
){

    .slow_doctest()

    from <- .nse_ensure_str(!!enquo(from))
    to <- .nse_ensure_str(!!enquo(to))

    chunk_size %<>% if_null(getOption('omnipath.uploadlists_chunk_size'))

    identifiers %>%
    unique %T>%
    {log_trace('Querying UniProt Uploadlists with %i IDs.', length(.))} %>%
    sort %>%
    chunks(chunk_size) %>%
    map(.uniprot_id_mapping_table, from, to) %>%
    bind_rows() %>%
    trim_and_distinct %T>%
    load_success()

}


#' R CMD check workaround, see details at \code{uniprot_id_mapping_table}
#'
#' @importFrom magrittr %<>%
#' @importFrom logger log_trace
#'
#' @noRd
.uniprot_id_mapping_table <- uniprot_domains %@% function(
    identifiers,
    from,
    to,
    .subdomain = 'www'
){

    from %<>% uploadlists_id_type
    to   %<>% uploadlists_id_type

    post <- list(
        from = from,
        to = to,
        format = 'tab',
        query = paste(identifiers, collapse = ' ')
    )

    log_trace(
        'UniProt uploadlists: querying `%s` to `%s`, %d identifiers.',
        from, to, length(identifiers)
    )

    generic_downloader(
        url_key = 'uniprot_uploadlists',
        url_param = list(.subdomain),
        post = post,
        content_param = list(encoding = 'ASCII'),
        resource = 'UniProt'
    )

}


#' Translate gene and protein identifiers
#'
#' Translates a vector of identifiers, resulting a new vector, or a column
#' of identifiers in a data frame by creating another column with the target
#' identifiers.
#'
#' @param d Character vector or data frame.
#' @param ... At least two arguments, with or without names. The first
#'     of these arguments describes the source identifier, the rest
#'     of them describe the target identifier(s). The values of all these
#'     arguments must be valid identifier types as shown in Details. The
#'     names of the arguments are column names. In case of the first
#'     (source) ID the column must exist. For the rest of the IDs new
#'     columns will be created with the desired names. For ID types provided
#'     as arguments without names, the name of the ID type will be used for
#'     column name.
#' @param uploadlists Force using the `uploadlists` service from UniProt.
#'     By default the plain query interface is used (implemented in
#'     \code{\link{uniprot_full_id_mapping_table}} in this package).
#'     If any of the provided ID types is only available in the uploadlists
#'     service, it will be automatically selected. The plain query interface
#'     is preferred because in the long term, with caching, it requires
#'     less download and data storage.
#' @param ensembl Logical: use data from Ensembl BioMart instead of UniProt.
#' @param keep_untranslated In case the output is a data frame, keep the
#'     records where the source identifier could not be translated. At
#'     these records the target identifier will be NA.
#' @param return_df Return a data frame even if the input is a vector.
#' @param reviewed Translate only reviewed (\code{TRUE}), only unreviewed
#'     (\code{FALSE}) or both (\code{NULL}) UniProt records. Matters only
#'     if \code{uploadlists} is \code{FALSE}.
#' @param organism Integer, NCBI Taxonomy ID of the organism (by default
#'     9606 for human). Matters only if \code{uploadlists} is \code{FALSE}.
#'
#' @return
#' \itemize{
#'     \item{Data frame: if the input is a data frame or the input is a
#'         vector and `return_df` is \code{TRUE}.}
#'     \item{Vector: if the input is a vector, there is only one target
#'         ID type and `return_df` is \code{FALSE}.}
#'     \item{List of vectors: if the input is a vector, there are more than
#'         one target ID types and `return_df` is \code{FALSE}. The names
#'         of the list will be ID types (as they were column names, see
#'         the description of the `...` argument), and the list will also
#'         include the source IDs.}
#' }
#'
#' @details
#' This function, depending on the `uploadlists` parameter, uses either
#' the uploadlists service of UniProt or plain UniProt queries to obtain
#' identifier translation tables. The possible values for `from` and `to`
#' are the identifier type abbreviations used in the UniProt API, please
#' refer to the table here: \url{https://www.uniprot.org/help/api_idmapping}.
#' In addition, simple synonyms are available which realize a uniform API
#' for the uploadlists and UniProt query based backends. These are the
#' followings:
#'
#' | **OmnipathR**  | **Uploadlists**      | **UniProt query**  | **Ensembl BioMart**        |
#' |----------------|----------------------|--------------------|----------------------------|
#' | uniprot        | ACC                  | id                 | uniprotswissprot           |
#' | uniprot_entry  | ID                   | entry name         |                            |
#' | trembl         | *reviewed = FALSE*   | *reviewed = FALSE* | uniprotsptrembl            |
#' | genesymbol     | GENENAME             | genes(PREFERRED)   | external_gene_name         |
#' | genesymbol_syn |                      | genes(ALTERNATIVE) | external_synonym           |
#' | hgnc           | HGNC_ID              | database(HGNC)     | hgnc_symbol                |
#' | entrez         | P_ENTREZGENEID       | database(GeneID)   |                            |
#' | ensembl        | ENSEMBL_ID           |                    | ensembl_gene_id            |
#' | ensg           | ENSEMBL_ID           |                    | ensembl_gene_id            |
#' | enst           | ENSEMBL_TRS_ID       | database(Ensembl)  | ensembl_transcript_id      |
#' | ensp           | ENSEMBL_PRO_ID       |                    | ensembl_peptide_id         |
#' | ensgg          | ENSEMBLGENOME_ID     |                    |                            |
#' | ensgt          | ENSEMBLGENOME_TRS_ID |                    |                            |
#' | ensgp          | ENSEMBLGENOME_PRO_ID |                    |                            |
#' | protein_name   |                      | protein names      |                            |
#' | pir            | PIR                  | database(PIR)      |                            |
#' | ccds           |                      | database(CCDS)     |                            |
#' | refseqp        | P_REFSEQ_AC          | database(refseq)   |                            |
#' | ipro           |                      |                    | interpro                   |
#' | ipro_desc      |                      |                    | interpro_description       |
#' | ipro_sdesc     |                      |                    | interpro_short_description |
#' | wikigene       |                      |                    | wikigene_name              |
#' | rnacentral     |                      |                    | rnacentral                 |
#' | gene_desc      |                      |                    | description                |
#' | wormbase       |                      | database(WormBase) |                            |
#' | flybase        |                      | database(FlyBase)  |                            |
#' | xenbase        |                      | database(Xenbase)  |                            |
#' | zfin           |                      | database(ZFIN)     |                            |
#' | pbd            | PBD_ID               | database(PDB)      | pbd                        |
#'
#' The mapping between identifiers can be ambiguous. In this case one row
#' in the original data frame yields multiple rows or elements in the
#' returned data frame or vector(s).
#'
#' @examples
#' d <- data.frame(uniprot_id = c('P00533', 'Q9ULV1', 'P43897', 'Q9Y2P5'))
#' d <- translate_ids(d, uniprot_id = uniprot, genesymbol)
#' d
#' #   uniprot_id genesymbol
#' # 1     P00533       EGFR
#' # 2     Q9ULV1       FZD4
#' # 3     P43897       TSFM
#' # 4     Q9Y2P5    SLC27A5
#'
#' @importFrom rlang !! !!! enquo := enquos quo_text set_names sym
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr pull left_join inner_join rename
#' @importFrom purrr map reduce2
#' @importFrom logger log_fatal
#' @importFrom utils tail
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#'     \item{\code{\link{uniprot_full_id_mapping_table}}}
#'     \item{\code{\link{ensembl_id_mapping_table}}}
#' }
#' @md
translate_ids <- function(
    d,
    ...,
    uploadlists = FALSE,
    ensembl = FALSE,
    keep_untranslated = TRUE,
    return_df = FALSE,
    organism = 9606,
    reviewed = TRUE
){

    .slow_doctest()

    # NSE vs. R CMD check workaround
    To <- NULL

    ids <-
        enquos(...) %>%
        map(.nse_ensure_str) %>%
        set_names(names(.) %||% unlist(.)) %>%
        set_names(ifelse(nchar(names(.)), names(.), unlist(.)))

    organism <- .nse_ensure_str(!!enquo(organism)) %>% ncbi_taxid
    id_cols <- names(ids)
    id_types <- unlist(ids)
    from_col <- id_cols[1]
    from_type <- id_types[1]
    to_cols <- id_cols %>% tail(-1)
    to_types <- id_types %>% tail(-1)

    use_vector <- !inherits(d, 'data.frame')

    if(use_vector){

        d %<>% list %>% set_names(from_col) %>% tibble(!!!.)

    }

    if(!from_col %in% names(d)){

        msg <- sprintf('translate_ids: no column named `%s`.', from_col)
        log_fatal(msg)
        stop(msg)

    }

    join_method <- `if`(keep_untranslated, left_join, inner_join)

    d %<>%
    reduce2(
        to_cols,
        to_types,
        function(d, to_col, to_type){

            translation_table <- id_translation_table(
                !!sym(from_type),
                !!sym(to_type),
                ensembl = ensembl,
                uploadlists = uploadlists,
                identifiers = d %>% pull(!!sym(from_col)),
                organism = organism,
                reviewed = reviewed
            )

            d %>%
            join_method(
                translation_table,
                by = 'From' %>% set_names(from_col)
            ) %>%
            rename(!!sym(to_col) := To)

        },
        .init = .
    )

    if(use_vector && !return_df){

        # convert output to a list
        d %<>%
            as.list %>%
            map(discard, is.na) %>%
            {`if`(length(.) == 2, .[[to_cols]], .)}

    }

    return(d)

}


#' Creates an ID translation table from UniProt data
#'
#' @param to Character or symbol: target ID type. See Details for possible
#'     values.
#' @param from Character or symbol: source ID type. See Details for possible
#'     values.
#' @param reviewed Translate only reviewed (\code{TRUE}), only unreviewed
#'     (\code{FALSE}) or both (\code{NULL}) UniProt records.
#' @param organism Integer, NCBI Taxonomy ID of the organism (by default
#'     9606 for human).
#'
#' @return A data frame (tibble) with columns `From` and `To`, UniProt IDs
#'     and the corresponding foreign IDs, respectively.
#'
#' @details
#' For both source and target ID type, this function accepts column codes
#' used by UniProt and some simple shortcuts defined here. For the UniProt
#' codes please refer to
#' https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames
#' The shortcuts are entrez, genesymbol, genesymbol_syn (synonym gene
#' symbols), hgnc, embl, refseqp (RefSeq protein), enst (Ensembl transcript),
#' uniprot_entry (UniProtKB AC, e.g. EGFR_HUMAN), protein_name (full name of
#' the protein), uniprot (UniProtKB ID, e.g. P00533). For a complete table
#' please refer to \code{\link{translate_ids}}.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate rename filter
#' @importFrom tidyr separate_rows
#' @importFrom rlang !! enquo
#' @importFrom logger log_trace
#' @importFrom stringr str_extract
#' @export
#'
#' @examples
#' uniprot_entrez <- uniprot_full_id_mapping_table(to = 'entrez')
#' uniprot_entrez
#' # # A tibble: 20,723 x 2
#' #    From   To
#' #    <chr>  <chr>
#' #  1 Q96R72 NA
#' #  2 Q9UKL2 23538
#' #  3 Q9H205 144125
#' #  4 Q8NGN2 219873
#' #  5 Q8NGC1 390439
#' # # . with 20,713 more rows
#'
#' @seealso \itemize{
#'     \item{\code{\link{translate_ids}}}
#'     \item{\code{\link{ensembl_id_mapping_table}}}
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#' }
uniprot_full_id_mapping_table <- function(
    to,
    from = 'id',
    reviewed = TRUE,
    organism = 9606
){

    .slow_doctest()

    # NSE vs. R CMD check workaround
    From <- To <- NULL

    strip_semicol <- function(v){sub(';$', '', v)}

    ids <- c(
        .nse_ensure_str(!!enquo(to)),
        .nse_ensure_str(!!enquo(from))
    )

    reviewed <- `if`(
        'trembl' %in% ids,
        FALSE,
        `if`('swissprot' %in% ids, TRUE, reviewed)
    )

    to   <- .nse_ensure_str(!!enquo(to))   %>% uniprot_id_type
    from <- .nse_ensure_str(!!enquo(from)) %>% uniprot_id_type

    reens <- 'ENS[A-Z]+\\d+'
    to_ens <- to == 'database(Ensembl)'
    from_ens <- from == 'database(Ensembl)'

    log_trace(
        paste0(
            'Creating ID mapping table from `%s` to `%s`, ',
            'for organism %d (only reviewed: %s)'
        ),
        from, to, organism, reviewed
    )

    c(from, to) %>%
    all_uniprots(reviewed = reviewed, organism = organism) %>%
    rename(From = 1, To = 2) %>%
    mutate(
        From = strip_semicol(From),
        To = strip_semicol(To)
    ) %>%
    separate_rows(From, sep = ';') %>%
    separate_rows(To, sep = ';') %>%
    filter(!is.na(From) & !is.na(To)) %>%
    {`if`(from_ens, mutate(., From = str_extract(From, reens)), .)} %>%
    {`if`(to_ens, mutate(., To = str_extract(To, reens)), .)} %>%
    trim_and_distinct

}


#' Trim padding whitespace and makes the records unique
#'
#' @param d A data frame with character columns.
#'
#' @return The same data frame with trailing and leading whitespaces removed
#'     and unique records.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate across
#' @importFrom tidyselect everything
#' @importFrom stringr str_trim
#' @noRd
trim_and_distinct <- function(d){

    d %>%
    mutate(across(everything(), str_trim)) %>%
    distinct

}


#' Choose an ID translation table
#'
#' @param from Character or symbol: the source identifier type.
#' @param to Character or symbol: the target identifier type.
#' @param uploadlists Logical: force to use the uploadlists service.
#' @param identifiers Character vector: query these identifiers from
#'     the uploadlists service.
#' @param organism Integer: NCBI Taxonomy ID of the organism.
#' @param reviewed Logical: use only the reviewed (SwissProt) records
#'     or only the unreviewed, or both (NULL).
#'
#' @return A data frame which is a translation table between the identifiers
#'     `from` and `to`.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang %||% !! sym enquo
#' @importFrom logger log_trace
#' @importFrom dplyr pull rename relocate
#' @noRd
#'
#' @seealso \itemize{
#'     \item{\code{\link{translate_ids}}}
#'     \item{\code{\link{uniprot_full_id_mapping_table}}}
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#'     \item{\code{\link{ensembl_id_mapping_table}}}
#' }
id_translation_table <- function(
    from,
    to,
    uploadlists = FALSE,
    ensembl = FALSE,
    identifiers = NULL,
    organism = 9606L,
    reviewed = TRUE
){

    # NSE vs. R CMD check workaround
    To <- From <- NULL

    from <- .nse_ensure_str(!!enquo(from))
    to <- .nse_ensure_str(!!enquo(to))

    if(ensembl){

        log_trace(
            'ID translation table: from `%s` to `%s`, using Ensembl BioMart.',
            from, to
        )

        result <-
            ensembl_id_mapping_table(
                to = !!sym(to),
                from = !!sym(from),
                organism = organism
            )

    }else if(
        uploadlists || (
            (
                !id_type_in(from, 'uniprot') ||
                !id_type_in(to, 'uniprot')
            ) &&
            id_type_in(from, 'uploadlists') &&
            id_type_in(to, 'uploadlists')
        )
    ){

        log_trace(
            'ID translation table: from `%s` to `%s`, using `uploadlists`.',
            from, to
        )

        swap <-
            (
                length(identifiers) > 10000L ||
                is.null(identifiers)
            ) &&
            to == 'uniprot'

        if(swap){
            to <- from
            from <- 'uniprot'
            identifiers <- NULL
            log_trace(
                'Querying Uploadlists by all UniProts, ',
                'swapping the table later.'
            )
        }

        result <-
            identifiers %>%
            {
                . %||% (
                    all_uniprots(organism = organism, reviewed = reviewed) %>%
                    pull(1L)
                )
            } %>%
            unique %>%
            uniprot_id_mapping_table(!!sym(from), !!sym(to)) %>%
            {`if`(
                swap,
                rename(., From = To, To = From) %>%
                relocate(From, .before = To),
                .
            )}

    }else{

        log_trace(
            'ID translation table: from `%s` to `%s`, using `uniprot`.',
            from, to
        )

        result <-
            uniprot_full_id_mapping_table(
                to = !!sym(to),
                from = !!sym(from),
                reviewed = reviewed,
                organism = organism
            )

    }

    result %>%
    trim_and_distinct

}


#' Is the ID type known to be available by a service
#'
#' @param id_type Character: name of the ID type.
#' @param service Character: name of the service; either "uniprot" or
#'     "uploadlists" (which is UniProt too, but another API).
#'
#' @noRd
id_type_in <- function(id_type, service){

    id_type %in% names(omnipath.env$id_types[[service]]) ||
    id_type %in% omnipath.env$id_types[[service]]

}


#' Read ID type information
#'
#' @importFrom magrittr %>%
#'
#' @noRd
.load_id_types <- function(pkgname){

    omnipath.env$id_types <-
        system.file(
            'internal',
            'id_types.json',
            package = pkgname,
            mustWork = TRUE
        ) %>%
        safe_json()

}


#' Identifier translation table from Ensembl
#'
#' @param to Character or symbol: target ID type. See Details for possible
#'     values.
#' @param from Character or symbol: source ID type. See Details for possible
#'     values.
#' @param organism Character or integer: NCBI Taxonomy ID or name of the
#'     organism (by default 9606 for human).
#'
#' @return A data frame (tibble) with columns `From` and `To`.
#'
#' @details The arguments \code{to} and \code{from} can be provided either
#' as character or as symbol (NSE). Their possible values are either Ensembl
#' attribute names or synonyms listed at \code{\link{translate_ids}}.
#'
#' @examples
#' ensp_up <- ensembl_id_mapping_table("ensp")
#' ensp_up
#' # # A tibble: 119,129 × 2
#' #    From   To
#' #    <chr>  <chr>
#' #  1 P03886 ENSP00000354687
#' #  2 P03891 ENSP00000355046
#' #  3 P00395 ENSP00000354499
#' #  4 P00403 ENSP00000354876
#' #  5 P03928 ENSP00000355265
#' # # . with 119,124 more rows
#'
#' @importFrom dplyr recode distinct
#' @importFrom magrittr %>%
#' @importFrom rlang enquo !! !!! set_names
#' @importFrom logger log_warn log_trace
#' @importFrom tibble tibble
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{translate_ids}}}
#'     \item{\code{\link{uniprot_full_id_mapping_table}}}
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#' }
ensembl_id_mapping_table <- function(
    to,
    from = 'uniprot',
    organism = 9606
){

    get_field_name <- function(label){

        label %>% recode(!!!omnipath.env$id_types$ensembl)

    }

    to <-
        .nse_ensure_str(!!enquo(to)) %>%
        ensembl_id_type()
    from <-
        .nse_ensure_str(!!enquo(from)) %>%
        ensembl_id_type()

    organism %<>% ensembl_name_warn

    if(is.na(organism)){

        log_warn('Returning empty table.')

        return(tibble(From = character(0), To = character(0)))

    }

    log_trace(
        paste0(
            'Creating ID mapping table from `%s` to `%s`, ',
            'for organism %s'
        ),
        from, to, organism
    )

    dataset <- ensembl_dataset(organism)

    biomart_query(
        attrs = c(from, to),
        dataset = dataset
    ) %>%
    set_names(c('From', 'To')) %>%
    trim_and_distinct

}


#' Ensembl identifier type label
#'
#' @param label Character: an ID type label, as shown in the table at
#'     \code{\link{translate_ids}}
#'
#' @return Character: the Ensembl specific ID type label, or the input
#'     unchanged if it could not be translated (still might be a valid
#'     identifier name). These labels should be valid Ensembl attribute
#'     names, directly usable in Ensembl queries.
#'
#' @examples
#' ensembl_id_type("uniprot")
#' # [1] "uniprotswissprot"
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{uniprot_id_type}}}
#'     \item{\code{\link{uploadlists_id_type}}}
#' }
ensembl_id_type <- function(label){

    resource_id_type(label, 'ensembl')

}



#' UniProt identifier type label
#'
#' @param label Character: an ID type label, as shown in the table at
#'     \code{\link{translate_ids}}
#'
#' @return Character: the UniProt specific ID type label, or the input
#'     unchanged if it could not be translated (still might be a valid
#'     identifier name). This is the label that one can use in UniProt
#'     REST queries.
#'
#' @examples
#' ensembl_id_type("entrez")
#' # [1] "database(GeneID)"
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{ensembl_id_type}}}
#'     \item{\code{\link{uploadlists_id_type}}}
#' }
uniprot_id_type <- function(label){

    resource_id_type(label, 'uniprot')

}


#' UniProt Uploadlists identifier type label
#'
#' @param label Character: an ID type label, as shown in the table at
#'     \code{\link{translate_ids}}
#'
#' @return Character: the UniProt Uploadlists specific ID type label, or
#'     the input unchanged if it could not be translated (still might be
#'     a valid identifier name). This is the label that one can use in
#'     UniProt Uploadlists (ID Mapping) queries.
#'
#' @examples
#' ensembl_id_type("entrez")
#' # [1] "P_ENTREZGENEID"
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{ensembl_id_type}}}
#'     \item{\code{\link{uniprot_id_type}}}
#' }
uploadlists_id_type <- function(label){

    resource_id_type(label, 'uploadlists')

}


#' Resource specific identifier type label
#'
#' @param label Character: an ID type label, as shown in the table at
#'     \code{\link{translate_ids}}
#'
#' @return Character: the resource specific ID type label, or the input
#'     unchanged if it could not be translated (still might be a valid)
#'     identifier name).
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !!!
#' @importFrom dplyr recode
#' @noRd
resource_id_type <- function(label, resource){

    label %>% recode(!!!omnipath.env$id_types[[resource]])

}


#' Is it a valid identifier type?
#'
#' @param label Character: identifier type label.
#'
#' @return Logical: true if label is a registered identifier type.
#'
#' @importFrom magrittr %>% is_in
#' @importFrom purrr map
#' @noRd
is_id_type <- function(label){

    omnipath.env$id_types %>%
    map(names) %>%
    unlist %>%
    is_in(label, .)

}

#' Looks like a UniProt ID?
#'
#' This function checks only the format of the IDs, no guarantee that these
#' IDs exist in UniProt.
#'
#' @param identifiers Character: one or more identifiers (typically a single
#'     string, a vector or a data frame column).
#'
#' @return Logical: true if all elements in the input (except NAs) looks like
#'     valid UniProt IDs. If the input is not a character vector, `FALSE`
#'     is returned.
#'
#' @examples
#' is_uniprot(all_uniprot_acs())
#' # [1] TRUE
#' is_uniprot("P00533")
#' # [1] TRUE
#' is_uniprot("pizza")
#' # [1] FALSE
#'
#' @importFrom magrittr %>%
#' @importFrom purrr discard
#' @importFrom stringr str_detect
#' @export
is_uniprot <- function(identifiers){

    reuni <- paste0(
        '^[OPQ][0-9][A-Z0-9]{3}[0-9]|',
        '[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$'
    )

    identifiers %>%
    {`if`(
        is.character(.),
        discard(., is.na) %>%
        str_detect(reuni) %>%
        all,
        FALSE
    )}

}


#' Attept to find valid SwissProt IDs
#'
#' In preparation, not yet functional
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom rlang enquos !! sym :=
#' @importFrom purrr reduce map
#' @importFrom dplyr left_join rename pull first mutate n
#' @importFrom utils tail
#' @importFrom tidyselect ends_with
#' @noRd
uniprot_cleanup <- function(d, ..., organism = 9606){

    # NSE vs. R CMD check workaround
    input <- output <- NULL

    SUFFIX <- '__uniprot_cleanup'

    organism %<>% ncbi_taxid

    columns <-
        enquos(...) %>%
        map(.nse_ensure_str)

    result_col <- columns %>% first
    columns %<>% tail(-1L)

    swissprots <- get_db('swissprots', organism = organism)
    trembls <- get_db('trembls', organism = organism)

    columns %>%
    reduce(
        function(d, col){
            d %>%
            left_join(
                d %>%
                pull(col) %>%
                uniprot_genesymbol_cleanup(organism = organism) %>%
                rename(!!sym(sprintf('%s%s', col, SUFFIX)) := output),
                by = c(!!sym(col) := input)
            )
        },
        .init = d %>% mutate(record_id = 1L:n())
    ) %>%
    mutate(
        !!sym(result_col) := c(ends_with(SUFFIX))
    )

}


#' TrEMBL to SwissProt by gene names
#'
#' @param uniprots Character vector possibly containing TrEMBL IDs.
#' @param organism Character or integer: organism name or identifier.
#' @param only_trembls Attempt to convert only known TrEMBL IDs of the
#'     organism. This is the recommended practice.
#'
#' @return Data frame with two columns: "input" and "output". The first one
#'     contains all identifiers from the input vector `uniprots`. The second
#'     one has the corresponding identifiers which are either SwissProt IDs
#'     with gene names identical to the TrEMBL IDs in the input, or if no
#'     such records are available, the output has the input items unchanged.
#'
#' @details
#' Sometimes one gene or protein is represented by multiple identifiers in
#' UniProt. These are typically slightly different isoforms, some of them
#' having TrEMBL IDs, some of the SwissProt. For the purposes of most systems
#' biology application, the most important is to identify the protein or gene
#' in a way that we can recognize it in other datasets. Unfortunately UniProt
#' or Ensembl do not seem to offer solution for this issue. Hence, if we find
#' that a TrEMBL ID has a gene name which is also associated with a SwissProt
#' ID, we replace this TrEMBL ID by that SwissProt. There might be a minor
#' difference in their sequence, but most of the omics analyses do not even
#' consider isoforms. And it is quite possible that later UniProt will convert
#' the TrEMBL record to an isoform within the SwissProt record. Typically
#' this translation is not so important (but still beneficial) for human,
#' but for other organisms it is critical especially when translating from
#' foreign identifiers.
#'
#' This function accepts a mixed input of UniProt IDs and provides a distinct
#' translation table that you can use to translate your data.
#'
#' @examples
#' \dontrun{
#' uniprot_genesymbol_cleanup('Q6PB82', organism = 10090)
#' # # A tibble: 1 × 2
#' #   input  output
#' #   <chr>  <chr>
#' # 1 Q6PB82 O70405
#' }
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom tibble tibble
#' @importFrom dplyr mutate select distinct full_join
#' @export
uniprot_genesymbol_cleanup <- function(
    uniprots,
    organism = 9606,
    only_trembls = TRUE
){

    # NSE vs. R CMD check workaround
    trembl <- swissprot <- input <- output <- NULL

    organism %<>% ncbi_taxid

    uniprots %<>% unique

    uniprots %>%
    {`if`(only_trembls, trembls_only(., organism = organism), .)} %>%
    tibble(input = .) %>%
    translate_ids(input = trembl, genesymbol, organism = organism) %>%
    translate_ids(genesymbol, swissprot, organism = organism) %>%
    full_join(
        uniprots %>%
        tibble(input = .),
        by = c('input')
    ) %>%
    mutate(output = ifelse(is.na(swissprot), input, swissprot)) %>%
    select(input, output) %>%
    distinct

}


#' Retain only TrEMBL IDs
#'
#' @param uniprots Character vector of UniProt IDs.
#' @param organism Character or integer: name or identifier of the organism.
#'
#' @return Character vector with only TrEMBL IDs.
#'
#' @examples
#' trembls_only(c("Q05BL1", "A0A654IBU3", "P00533"))
#' # [1] "Q05BL1" "A0A654IBU3"
#'
#' @importFrom magrittr %<>% %>% is_in
#' @importFrom purrr keep
#' @export
trembls_only <- function(uniprots, organism = 9606){

    uniprots %>%
    only_id_type('trembls', organism = organism)

}


#' Retain only SwissProt IDs
#'
#' @param uniprots Character vector of UniProt IDs.
#' @param organism Character or integer: name or identifier of the organism.
#'
#' @return Character vector with only SwissProt IDs.
#'
#' @examples
#' swissprots_only(c("Q05BL1", "A0A654IBU3", "P00533"))
#' # [1] "P00533"
#'
#' @importFrom magrittr %<>% %>% is_in
#' @importFrom purrr keep
#' @export
swissprots_only <- function(uniprots, organism = 9606){

    uniprots %>%
    only_id_type('swissprots', organism = organism)

}


#' @importFrom magrittr %>%
#' @importFrom purrr keep
#' @noRd
only_id_type <- function(identifiers, id_type, organism = 9606){

    identifiers %>%
    keep(is_id_of_type, id_type = id_type, organism = organism)
}


#' Check for SwissProt IDs
#'
#' @param uniprots Character vector of UniProt IDs.
#' @param organism Character or integer: name or identifier of the organism.
#'
#' @return Logical vector TRUE for SwissProt IDs and FALSE for any other
#'     element.
#'
#' @examples
#' is_swissprot(c("Q05BL1", "A0A654IBU3", "P00533"))
#' # [1] FALSE FALSE TRUE
#'
#' @importFrom magrittr %<>% %>% is_in
#' @importFrom purrr keep
#' @export
is_swissprot <- function(uniprots, organism = 9606){

    uniprots %>%
    is_id_of_type('swissprots', organism = organism)

}


#' Check for TrEMBL IDs
#'
#' @param uniprots Character vector of UniProt IDs.
#' @param organism Character or integer: name or identifier of the organism.
#'
#' @return Logical vector TRUE for TrEMBL IDs and FALSE for any other
#'     element.
#'
#' @examples
#' is_trembl(c("Q05BL1", "A0A654IBU3", "P00533"))
#' # [1] TRUE TRUE FALSE
#'
#' @importFrom magrittr %<>% %>% is_in
#' @importFrom purrr keep
#' @export
is_trembl <- function(uniprots, organism = 9606){

    uniprots %>%
    is_id_of_type('trembls', organism = organism)

}


#' @importFrom magrittr %<>% %>% is_in
#' @noRd
is_id_of_type <- function(identifiers, id_type, organism = 9606){

    organism %<>% ncbi_taxid

    all_ids <- get_db(id_type, organism = organism)

    identifiers %>%
    is_in(all_ids)

}


#' Ensembl-UniProt identifier table
#'
#' In preparation, not yet functional
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang !! enquo
#' @importFrom dplyr bind_rows select distinct mutate group_by filter ungroup
#' @noRd
ensembl_uniprot <- function(ens_id_type = 'ensg', organism = 9606L){

    # NSE vs. R CMD check workaround
    ensembl <- NULL

    dataset <- organism %>% ensembl_dataset

    ens_id_type <- .nse_ensure_str(!!enquo(ens_id_type))

    biomart_query(
        attrs = c('uniprotswissprot', 'uniprotsptrembl'),
        gene = ens_id_type %in% c('gene', 'ensg'),
        transcript = ens_id_type %in% c('transcript', 'enst'),
        peptide = ens_id_type %in% c('protein', 'peptide', 'ensp'),
        dataset = dataset
    ) %>%
    {bind_rows(
        select(., ensembl = 1L, uniprot = 2L),
        select(., ensembl = 1L, uniprot = 3L)
    )} %>%
    distinct %>%
    mutate(is_trembl = is_trembl(uniprot)) %>%
    group_by(ensembl) %>%
    filter(all(is_trembl) | !is_trembl) %>%
    ungroup %>%
    filter(!is.na(uniprot))


}
