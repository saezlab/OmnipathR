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


#' ID translation data from UniProt ID Mapping
#'
#' Retrieves an identifier translation table from the UniProt ID Mapping
#' service (https://www.uniprot.org/help/id_mapping).
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
#' the table here: \code{\link{uniprot_idmapping_id_types}} or
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

    chunk_size %<>% if_null(getOption('omnipathr.uploadlists_chunk_size'))

    identifiers %>%
    unique %T>%
    {log_trace('Querying UniProt ID Mapping with %i IDs.', length(.))} %>%
    sort %>%
    chunks(chunk_size) %>%
    map(.uniprot_id_mapping_table, from, to) %>%
    bind_rows() %>%
    trim_and_distinct %T>%
    load_success()

}


#' R CMD check workaround, see details at \code{uniprot_id_mapping_table}
#'
#' @importFrom magrittr %<>% extract2 %>%
#' @importFrom httr2 req_body_form req_headers req_perform resp_body_json
#' @importFrom httr2 request
#' @importFrom logger log_trace
#' @importFrom stringr str_replace
#' @importFrom readr read_tsv cols
#' @importFrom rlang !!!
#'
#' @noRd
.uniprot_id_mapping_table <- function(
    identifiers,
    from,
    to
){

    from %<>% uploadlists_id_type
    to   %<>% uploadlists_id_type

    post <- list(
        from = from,
        to = to,
        ids = paste(identifiers, collapse = ' ')
    )

    log_trace(
        'UniProt id-mapping: querying `%s` to `%s`, %d identifiers.',
        from, to, length(identifiers)
    )

    run_url <- url_parser(url_key = 'uniprot_idmapping_run')

    version <- omnipath_cache_latest_or_new(
        url = run_url,
        post = post,
        ext = 'tsv'
    )

    from_cache <- version$status == CACHE_STATUS$READY

    if(!from_cache){

        run_result <-
            request(run_url) %>%
            req_body_form(!!!post) %>%
            req_headers('Accept' = 'application/json') %>%
            req_perform() %>%
            resp_body_json()

        if(!is.null(run_result$messages)) {
            msg <- run_result$messages %>%
                paste(collapse = ' ') %>%
                sprintf('Error at querying UniProt ID mapping: %s', .)
            log_error_with_info(msg)
            stop(msg)
        }

        jobid <- run_result$jobId

        poll_interval <- getOption('omnipathr.uniprot_idmapping_poll_interval')
        timeout <- getOption('omnipathr.uniprot_idmapping_timeout')
        max_polls <- ceiling(timeout / poll_interval)

        for(i in 1L:max_polls) {

            status_url <-
                url_parser(
                    url_key = 'uniprot_idmapping_poll',
                    url_param = list(jobid)
                )

            log_trace('Polling `%s`', status_url)

            status <-
                request(status_url) %>%
                req_headers('Accept' = 'application/json') %>%
                req_perform() %>%
                resp_body_json()

            if(!is.null(status$results) || !is.null(status$failedIds)){
                break
            }else if(!is.null(status$messages)){
                msg <- sprintf(
                    'Error at querying UniProt ID mapping: %s',
                    status$messages
                )
                log_error_with_info(msg)
                stop(msg)

            }

            Sys.sleep(poll_interval)

        }

        log_trace(
            'UniProt ID mapping job is ready, getting results URL: `%s`',
            jobid
        )

        result_url <-
            url_parser(
                url_key = 'uniprot_idmapping_details',
                url_param = list(jobid)
            ) %>%
            request() %>%
            req_headers('Accept' = 'application/json') %>%
            req_perform() %>%
            resp_body_json() %>%
            extract2('redirectURL') %>%
            str_replace('/idmapping/results/', '/idmapping/stream/') %>%
            str_replace('/results/', '/results/stream/') %>%
            sprintf('%s?format=tsv', .)

        log_trace(
            'Retrieving UniProt ID mapping results from: `%s`',
            result_url
        )

        path <- download_base(url = result_url, path = version$path)

        omnipath_cache_download_ready(version)

    }

    version$path %>%
    read_tsv(col_types = cols(), progress = FALSE) %>%
    origin_cache(from_cache) %>%
    source_attrs('UniProt', run_url)


}


#' ID types available in the UniProt ID Mapping service
#'
#' @return A data frame listing the ID types.
#'
#' @examples
#' uniprot_idmapping_id_types()
#'
#' @importFrom magrittr %>% %T>% extract2
#' @importFrom tidyr unnest_wider unnest_longer
#' @importFrom tibble tibble
#' @importFrom jsonlite fromJSON
#' @export
uniprot_idmapping_id_types <- function() {

    # NSE vs. R CMD check workaround
    groups <- items <- NULL

    url_key <- 'uniprot_idmapping_id_types'

    generic_downloader(
        url_key = url_key,
        reader = fromJSON,
        reader_param = list(simplifyDataFrame = FALSE)
    ) %>%
    extract2('groups') %>%
    tibble(groups = .) %>%
    unnest_wider(col = groups) %>%
    unnest_longer(col = items) %>%
    unnest_wider(col = items) %>%
    source_attrs('UniProt', get_url(url_key)) %T>%
    load_success()

}


#' Translate gene, protein and small molecule identifiers
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
#' @param hmdb Logical: use HMDB ID translation data.
#' @param ramp Logical: use RaMP ID translation data.
#' @param chalmers Logical: use ID translation data from Chalmers Sysbio GEM.
#' @param entity_type Character: "gene" and "smol" are short symbols for
#'     proteins, genes and small molecules respectively. Several other synonyms
#'     are also accepted.
#' @param keep_untranslated In case the output is a data frame, keep the
#'     records where the source identifier could not be translated. At
#'     these records the target identifier will be NA.
#' @param return_df Return a data frame even if the input is a vector.
#' @param reviewed Translate only reviewed (\code{TRUE}), only unreviewed
#'     (\code{FALSE}) or both (\code{NULL}) UniProt records. Matters only
#'     if \code{uploadlists} is \code{FALSE}.
#' @param organism Character or integer, name or NCBI Taxonomy ID of the
#'     organism (by default 9606 for human). Matters only if
#'     \code{uploadlists} is \code{FALSE}.
#' @param complexes Logical: translate complexes by their members. Only
#'     complexes where all members can be translated will be included in the
#'     result. If \code{NULL}, the option `omnipathr.complex_translation` will
#'     be used.
#' @param complexes_one_to_many Logical: allow combinatorial expansion or
#'     use only the first target identifier for each member of each complex.
#'     If \code{NULL}, the option `omnipathr.complex_translation_one_to_many`
#'     will be used.
#' @param track Logical: Track the records (rows) in the input data frame by
#'     adding a column `record_id` with the original row numbers.
#' @param quantify_ambiguity Logical or character: inspect the mappings for each
#'     ID for ambiguity. If TRUE, for each translated column, two new columns
#'     will be created with numeric values, representing the ambiguity of the
#'     mapping on the "from" and "to" side of the translation, respectively.
#'     If a character value provided, it will be used as a column name suffix
#'     for the new columns.
#' @param qualify_ambiguity Logical or character: inspect the mappings for each
#'     ID for ambiguity. If TRUE, for each translated column, a new column
#'     will be inculded with values `one-to-one`, `one-to-many`, `many-to-one`
#'     or `many-to-many`. If a character value provided, it will be used as a
#'     column name suffix for the new column.
#' @param ambiguity_groups Character vector: additional column names to group by
#'     during inspecting ambiguity. By default, the identifier columns (from
#'     and to) will be used to determine the ambiguity of mappings.
#' @param expand Logical: if \code{TRUE}, ambiguous (to-many) mappings will be
#'     expanded to multiple rows, resulting character type columns; if
#'     \code{FALSE}, the original rows will be kept intact, and the target
#' @param ambiguity_global Logical or character: if `ambiguity_groups` are
#'      provided, analyse ambiguity also globally, across the whole data frame.
#'      Character value provides a custom suffix for the columns quantifying
#'      and qualifying global ambiguity.
#' @param ambiguity_summary Logical: generate a summary about the ambiguity of the
#'      translation and make it available as an attribute.
#'      columns will be lists of character vectors.
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
#' For a complete list of ID types and their synonyms, including metabolite and
#' chemical ID types which are not shown here, see \code{\link{id_types}}.
#'
#' The mapping between identifiers can be ambiguous. In this case one row
#' in the original data frame yields multiple rows or elements in the
#' returned data frame or vector(s).
#'
#' The columns in the translation must be character type. Some ID types are
#' numeric, such as the ones from NCBI, these are sometimes present in data
#' frames as double or integer type. This function will convert those columns
#' to character.
#'
#' @examples
#' d <- data.frame(
#'     uniprot_id = c(
#'         'P00533', 'Q9ULV1', 'P43897', 'Q9Y2P5',
#'         'P01258', 'P06881', 'P42771', 'Q8N726'
#'     )
#' )
#' d <- translate_ids(d, uniprot_id = uniprot, genesymbol)
#' d
#' #   uniprot_id genesymbol
#' # 1     P00533       EGFR
#' # 2     Q9ULV1       FZD4
#' # 3     P43897       TSFM
#' # 4     Q9Y2P5    SLC27A5
#'
#' @importFrom rlang !! !!! enquo := enquos quo_text set_names sym syms
#' @importFrom magrittr %>% %<>% or
#' @importFrom tibble tibble
#' @importFrom dplyr pull left_join inner_join mutate row_number
#' @importFrom dplyr relocate rename bind_rows group_by summarize
#' @importFrom purrr map reduce2 discard
#' @importFrom logger log_fatal
#' @importFrom utils tail
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{translate_ids_multi}}}
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#'     \item{\code{\link{uniprot_full_id_mapping_table}}}
#'     \item{\code{\link{ensembl_id_mapping_table}}}
#'     \item{\code{\link{hmdb_id_mapping_table}}}
#'     \item{\code{\link{id_types}}}
#'     \item{\code{\link{ensembl_id_type}}}
#'     \item{\code{\link{uniprot_id_type}}}
#'     \item{\code{\link{uploadlists_id_type}}}
#'     \item{\code{\link{hmdb_id_type}}}
#'     \item{\code{\link{chalmers_gem_id_type}}}
#' }
#' @md
translate_ids <- function(
    d,
    ...,
    uploadlists = FALSE,
    ensembl = FALSE,
    hmdb = FALSE,
    ramp = FALSE,
    chalmers = FALSE,
    entity_type = NULL,
    keep_untranslated = TRUE,
    return_df = FALSE,
    organism = 9606,
    reviewed = TRUE,
    complexes = NULL,
    complexes_one_to_many = NULL,
    track = FALSE,
    quantify_ambiguity = FALSE,
    qualify_ambiguity = FALSE,
    ambiguity_groups = NULL,
    ambiguity_global = FALSE,
    ambiguity_summary = FALSE,
    expand = TRUE
){

    # NSE vs. R CMD check workaround
    To <- From <- NULL

    complexes %<>% if_null(getOption('omnipathr.complex_translation'))
    organism %<>% ncbi_taxid
    entity_type %<>% ensure_entity_type
    ramp %<>% or(entity_type == 'small_molecule' && !hmdb && !chalmers)

    ids <- ellipsis_to_char(...)
    id_cols <- names(ids)
    id_types <- unlist(ids)
    from_col <- id_cols[1L]
    from_type <- id_types[1L]
    to_cols <- id_cols %>% tail(-1L)
    to_types <- id_types %>% tail(-1L)
    track %<>% toggle_label('record_id')

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
    ensure_character(!!sym(from_col)) %>%
    {`if`(track$toggle, mutate(., !!sym(track$label) := row_number()), .)} %>%
    reduce2(
        to_cols,
        to_types,
        function(d, to_col, to_type){

            translation_table <-
                id_translation_table(
                    !!sym(from_type),
                    !!sym(to_type),
                    ensembl = ensembl,
                    uploadlists = uploadlists,
                    hmdb = hmdb,
                    ramp = ramp,
                    chalmers = chalmers,
                    entity_type = entity_type,
                    identifiers = d %>% pull(!!sym(from_col)),
                    organism = organism,
                    reviewed = reviewed
                ) %>%
                {`if`(
                    complexes,
                    bind_rows(
                        .,
                        translate_complexes(
                            d,
                            !!sym(from_col),
                            mapping = .,
                            one_to_many = complexes_one_to_many
                        )
                    ),
                    .
                )} %>%
                ensure_character(From, To) %>%
                {`if`(expand, ., group_by(., From) %>% summarize(To = list(To)))}

            log_trace(
                '%i rows before translation, %i %s IDs in column `%s`.',
                nrow(d),
                d %>% pull(!!sym(from_col)) %>% n_distinct,
                from_type,
                from_col
            )

            d %<>%
            join_method(
                translation_table,
                by = 'From' %>% set_names(from_col)
            ) %>%
            {`if`(keep_untranslated, ., filter(., !is.na(To)))} %>%
            relocate(To, .after = last_col()) %>%
            mutate(!!sym(to_col) := To) %>%
            select(-To) %>%
            {`if`(
                bool(quantify_ambiguity) || bool(qualify_ambiguity),
                ambiguity(
                    .,
                    !!sym(from_col),
                    !!sym(to_col),
                    groups = ambiguity_groups,
                    quantify = quantify_ambiguity,
                    qualify = qualify_ambiguity,
                    summary = ambiguity_summary,
                    global = ambiguity_global
                ),
                .
            )} %>%
            {`if`(
                expand,
                .,
                mutate(., !!sym(to_col) := na_to_len0(!!sym(to_col)))
            )}

            log_trace(
                paste0(
                    '%i rows after translation; translated %i `%s` ',
                    'IDs in column `%s` to %i `%s` IDs in column `%s`.'
                ),
                nrow(d),
                d %>% pull(!!sym(from_col)) %>% unlist %>% n_distinct,
                from_type,
                from_col,
                d %>% pull(!!sym(to_col)) %>% unlist %>% n_distinct,
                to_type,
                to_col
            )

            return(d)

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


#' Inspect the ambiguity of a mapping
#'
#' @param d Data frame: a data frame with two columns to be inspected. It might
#'     contain arbitrary other columns. Existing grouping will be removed.
#' @param from_col Character: column name of the "from" side of the mapping.
#' @param to_col Character: column name of the "to" side of the mapping.
#' @param groups Character vector of column names. Inspect ambiguity within
#'     these groups; by default, ambiguity is determined across all rows.
#' @param quantify Logical or character: inspect the mappings for each
#'     ID for ambiguity. If TRUE, for each translated column, two new columns
#'     will be created with numeric values, representing the ambiguity of the
#'     mapping on the "from" and "to" side of the translation, respectively.
#'     If a character value provided, it will be used as a column name suffix
#'     for the new columns.
#' @param qualify Logical or character: inspect the mappings for each
#'     ID for ambiguity. If TRUE, for each translated column, a new column
#'     will be inculded with values `one-to-one`, `one-to-many`, `many-to-one`
#'     or `many-to-many`. If a character value provided, it will be used as a
#'     column name suffix for the new column.
#' @param expand Logical: override the expansion of target columns, including
#'     `to_col`: by default, this function expands data into multiple rows if
#'     the `to_col` has already been expanded. Using this argument, the
#'     `to_col` and other target columns will be lists of vectors for `expand =
#'     FALSE`, and simple vectors for `expand = TRUE`.
#' @param global Logical or character: if `groups` are provided, analyse
#'      ambiguity also globally, across the whole data frame. Character value
#'      provides a custom suffix for the columns quantifying and qualifying
#'      global ambiguity.
#' @param summary Logical: generate a summary about the ambiguity of the
#'      translation and make it available as an attribute.
#'
#' @return A data frame (tibble) with ambiguity information added in new
#'     columns, as described at the "quantify" and "qualify" arguments.
#'
#' @importFrom rlang !! !!! enquo := sym syms set_names
#' @importFrom magrittr %>% %<>% or extract not
#' @importFrom dplyr pull mutate select row_number cur_group_id n_distinct
#' @importFrom dplyr group_by summarize reframe first across relocate
#' @importFrom dplyr case_when
#' @importFrom tidyselect where all_of any_of last_col
#' @importFrom tidyr unnest_longer
#' @export
ambiguity <- function(
        d,
        from_col,
        to_col,
        groups = NULL,
        quantify = TRUE,
        qualify = TRUE,
        expand = NULL,
        global = FALSE,
        summary = FALSE
    ) {

    # NSE vs. R CMD check workaround
    To <- NULL

    if (!bool(quantify) && !bool(qualify)) {
        return(d)
    }

    from_col <- .nse_ensure_str(!!enquo(from_col))
    to_col <- .nse_ensure_str(!!enquo(to_col))
    bygrp <- function(x){doif(x, !is.null(groups), ~sprintf('%s_bygroup', x))}
    rm_bygrp <- function(x){doif(x, !is.null(groups), ~sub('_.*', '', x))}
    quantify %<>% toggle_label('ambiguity')
    qualify %<>% toggle_label('ambiguity')
    quantify$label %<>% bygrp
    qualify$label %<>% bygrp
    qn_from <- sprintf('%s_%s_from_%s', from_col, to_col, quantify$label)
    qn_to <- sprintf('%s_%s_to_%s', from_col, to_col, quantify$label)

    record_id <- '__omnipathr_tmp_record_id'
    global %<>% toggle_label('global')
    global$toggle %<>% and(is.null(groups) %>% not)
    summary %<>% toggle_label(sprintf('ambiguity_%s_%s', from_col, to_col))
    collapsed <- d %>% pull(to_col) %>% is.list
    expand %<>% if_null(collapsed %>% not)
    collapse <- collapsed %>% or(expand) %>% not

    qn_cols <- c(qn_from, qn_to)
    ql_col <- qualify$label %>% sprintf('%s_%s_%s', from_col, to_col, .)
    q_cols <- c(qn_cols, ql_col)
    lst_cols <- c(to_col, qn_from, ql_col)

    distinct_ids <- function(l) {
        unlist(l, recursive = FALSE) %>%
        n_distinct(na.rm = TRUE)
    }

    log_trace(
        'Translation ambiguity between columns `%s` and `%s`',
        from_col,
        to_col
    )

    d %>%
    # make sure we have NAs and not empty vectors for untranslated items
    doif(
        collapsed,
        ~mutate(.x, !!sym(to_col) := len0_to_na(!!sym(to_col)))
    ) %>%
    # quantify ambiguity on the "to" side (ID translates to N targets)
    group_by(., !!sym(from_col), !!!syms(groups)) %>%
    mutate(!!sym(qn_to) := distinct_ids(!!sym(to_col))) %>%
    # create the "record_id", so later we know how to restore the original rows
    doif(!collapse, ungroup) %>%
    mutate(!!sym(record_id) := `if`(collapse, cur_group_id(), row_number())) %>%
    doif(collapse, ungroup) %>%
    # make sure the "to" identifiers are expanded
    unnest_longer(!!sym(to_col)) %>%
    # quantify ambiguity on the "from" side (N IDs translate to one target)
    group_by(!!sym(to_col), !!!syms(groups)) %>%
    mutate(
        !!sym(qn_from) := first(ifelse(
            is.na(!!sym(to_col)),
            1L,  # NA is not a target ID, these should end up one-to-none
            distinct_ids(!!sym(from_col))
        ))
    ) %>%
    ungroup %>%
    # compile qualitative column based on the quantitative ones
    doif(
        qualify$toggle,
        ~mutate(
            .x,
            !!sym(ql_col) := sprintf(
                '%s-to-%s',
                one_many(!!sym(qn_from)),
                one_many(!!sym(qn_to))
            )
        )
    ) %>%
    # add summary if required
    doif(
        summary$toggle,
        ~exec(
            ambiguity_summarize,
            d = .x,
            from_col = from_col,
            to_col = to_col,
            label = summary$label,
            groups = groups,
            !!!q_cols
        )
    ) %>%
    # quantify columns: remove or move them to the end and name their values
    {`if`(
        quantify$toggle,
        rowwise(.) %>%
        mutate(across(all_of(qn_cols), ~set_names(.x, !!sym(to_col)))) %>%
        mutate(!!sym(qn_to) := map_int(!!sym(qn_to), ~extract(.x, 1L))) %>%
        ungroup %>%
        relocate(all_of(qn_cols), .after = last_col()),
        select(., -all_of(qn_cols))
    )} %>%
    # qualify column: remove or move it to the end and name its values
    {`if`(
        qualify$toggle,
        rowwise(.) %>%
        mutate(!!sym(ql_col) := set_names(!!sym(ql_col), !!sym(to_col))) %>%
        ungroup %>%
        relocate(!!sym(ql_col), .after = last_col()),
        select(., -(!!!syms(ql_col)))
    )} %>%
    # expand new list cols or restore original rows based on "record_id"
    {`if`(
        expand,
        mutate(., across(any_of(q_cols), ~unname(unlist(.x)))),
        group_by(., !!sym(record_id)) %>%
        {list(
            cols = colnames(.),
            df = reframe(
                .,
                across(any_of(lst_cols), ~list(.x)),
                across(where(is.list), ~extract(.x, 1L)),
                across(where(~!is.list(.x)), first)
            )
        )} %>%
        {relocate(.$df, all_of(.$cols))}
    )} %>%
    select(-!!sym(record_id)) %>%
    # do also the global ambiguity, if required
    doif(
        global$toggle,
        ~ambiguity(
            .x,
            to_col = !!sym(to_col),
            from_col = !!sym(from_col),
            quantify = quantify %>% label_toggle %>% rm_bygrp,
            qualify = qualify %>% label_toggle %>% rm_bygrp,
            summary = summary %>% label_toggle,
            expand = expand
        )
    )

}


#' Summarise the results of an ID translation ambiguity analysis
#'
#' @param d A data frame with all three output columns of the ambiguity
#'     analysis (quantify_to, quantify_from and qualify) available under the
#'     names provided in the rest of the arguments.
#' @param from_col Character: name of the column with the original IDs.
#' @param to_col Character: name of the column with the translated IDs.
#' @param qn_from Character: name of the column with the numeric ambiguity on
#'     the "from" side.
#' @param qn_to Character: name of the column with the numeric ambiguity on
#'     the "to" side.
#' @param qualify Character: name of the column with the qualitative ambiguity
#'     on the "to" side.
#' @param label Character: name of the attribute to store the summary in.
#'
#' @return The input data frame with the summary appended to a data frame under
#'     the attribute `label`.
#'
#' @importFrom tidyr unnest_longer nest unnest_wider
#' @importFrom dplyr group_split group_by mutate summarise
#' @importFrom dplyr case_when pull bind_cols first
#' @importFrom tibble as_tibble
#' @importFrom purrr map reduce
#' @importFrom rlang set_names sym syms !! !!! := exec is_integer
#' @importFrom stringr str_replace_all str_c
#' @importFrom magrittr %<>% %>% set_attr extract2
#' @noRd
ambiguity_summarize <- function(
        d,
        from_col,
        to_col,
        qn_from,
        qn_to,
        qualify,
        label = NULL,
        groups = NULL
    ) {

    # NSE vs. R CMD check workaround
    data <- NULL

    count <- function(
            d,
            grp_cols,
            item_col,
            transform = identity,
            distrib = FALSE
        ) {

        prefix <-
            grp_cols %>%
            {case_when(
                identical(., to_col) ~ 'from_',
                identical(., from_col) ~ 'to_',
                TRUE ~ ''
            )}

        .transform <- function(x) {

            x %>% transform %>% str_replace_all('-', '_') %>% str_c(prefix, .)
        }

        total_col <- str_c(prefix, 'total')
        distrib_col <- str_c(prefix, 'distrib')

        d %>%
        doif(
            !distrib,
            ~mutate(.x, !!sym(item_col) := .transform(!!sym(item_col)))
        ) %>%
        group_by(!!!syms(grp_cols)) %>%
        summarise(
            !!sym(item_col) := first(!!sym(item_col)),
            .groups = 'drop'
        ) %>%
        pull(!!sym(item_col)) %>%
        table %>%
        as.list %>%
        {`if`(
              distrib,
              tibble(!!sym(distrib_col) := .),
              as_tibble(.) %>%
              mutate(!!sym(total_col) := sum(.))
        )}

    }

    counters <- list(
        list(grp_cols = c(from_col, to_col), item_col = qualify),
        list(grp_cols = from_col, item_col = qn_to, transform = one_many),
        list(grp_cols = to_col, item_col = qn_from, transform = one_many)
# This is an attempt to keep the ambiuity count distribution in the table
# but so far it fails on the recursive nature of the downstream unnest
#        list(grp_cols = from_col, item_col = qn_to, distrib = TRUE),
#        list(grp_cols = to_col, item_col = qn_from, distrib = TRUE)
    )

    count_all <- function(d) {
        counters %>%
        map(~exec(count, d, !!!.x)) %>%
        bind_cols %>%
        mutate(scope = scope)
    }

    log_trace(
        'Summarizing translation ambiguity for columns `%s`, `%s` and `%s`',
        qn_from, qn_to, qualify
    )

    scope <- `if`(is.null(groups), 'global', 'group')

    d %>%
    unnest_longer(c(from_col, to_col, qualify, qn_from, qn_to)) %>%
    nest(.by = groups) %>%
    mutate(summary = map(data, count_all)) %>%
    select(-data) %>%
    unnest_wider(summary) %>%
    mutate(across(where(is_integer), ~replace_na(.x, 0L))) %>%
    bind_rows(attr(d, label), .) %>%
    relocate(scope, .after = last_col()) %>%
    arrange(scope) %>%
    set_attr(d, label, .) %>%
    tbl_attrs

}

#' Translate gene, protein and small molecule identifiers from multiple columns
#'
#' Especially when translating network interactions, where two ID columns exist
#' (source and target), it is convenient to call the same ID translation on
#' multiple columns. The \code{\link{translate_ids}} function is already able
#' to translate to multiple ID types in one call, but is able to work only from
#' one source column. Here too, multiple target IDs are supported. The source
#' columns can be listed explicitely, or they might share a common stem, in
#' this case the first element of \code{...} will be used as stem, and the
#' column names will be created by adding the \code{suffixes}. The
#' \code{suffixes} are also used to name the target columns. If no
#' \code{suffixes} are provided, the name of the source columns will be added
#' to the name of the target columns. ID types can be defined the same way as
#' for \code{\link{translate_ids}}. The only limitation is that, if the source
#' columns are provided as stem+suffixes, they must be the same ID type.
#'
#' @param d A data frame.
#' @param ... At least two arguments, with or without names. These arguments
#'     describe identifier columns, either the ones we translate from (source),
#'     or the ones we translate to (target). Columns existing in the data frame
#'     will be used as source columns. All the rest will be considered target
#'     columns. Alternatively, the source columns can be defined as a stem and
#'     a vector of suffixes, plus a separator between the stem and suffix. In
#'     this case, the source columns will be the ones that exist in the data
#'     frame with the suffixes added. The values of all these
#'     arguments must be valid identifier types as shown at
#'     \code{\link{translate_ids}}. If ID type is provided only for the first
#'     source column, the rest of the source columns will be assumed to have
#'     the same ID type. For the target identifiers new columns will be
#'     created with the desired names, with the suffixes added. If no suffixes
#'     provided, the names of the source columns will be used instead.
#' @param suffixes Column name suffixes in case the names should be composed of
#'     stem and suffix.
#' @param suffix_sep Character: separator between the stem and suffixes.
#' @param uploadlists Force using the `uploadlists` service from UniProt.
#'     By default the plain query interface is used (implemented in
#'     \code{\link{uniprot_full_id_mapping_table}} in this package).
#'     If any of the provided ID types is only available in the uploadlists
#'     service, it will be automatically selected. The plain query interface
#'     is preferred because in the long term, with caching, it requires
#'     less download and data storage.
#' @param ensembl Logical: use data from Ensembl BioMart instead of UniProt.
#' @param hmdb Logical: use HMDB ID translation data.
#' @param ramp Logical: use RaMP ID translation data.
#' @param chalmers Logical: use ID translation data from Chalmers Sysbio GEM.
#' @param entity_type Character: "gene" and "smol" are short symbols for
#'     proteins, genes and small molecules respectively. Several other synonyms
#'     are also accepted.
#' @param keep_untranslated In case the output is a data frame, keep the
#'     records where the source identifier could not be translated. At
#'     these records the target identifier will be NA.
#' @param organism Character or integer, name or NCBI Taxonomy ID of the
#'     organism (by default 9606 for human). Matters only if
#'     \code{uploadlists} is \code{FALSE}.
#' @param reviewed Translate only reviewed (\code{TRUE}), only unreviewed
#'     (\code{FALSE}) or both (\code{NULL}) UniProt records. Matters only
#'     if \code{uploadlists} is \code{FALSE}.
#' @param expand Logical: if \code{TRUE}, ambiguous (to-many) mappings will be
#'     expanded to multiple rows, resulting character type columns; if
#'     \code{FALSE}, the original rows will be kept intact, and the target
#'
#' @return
#' A data frame with all source columns translated to all target identifiers.
#' The number of new columns is the product of source and target columns. The
#' target columns are distinguished by the suffexes added to their names.
#'
#' @examples
#' ia <- omnipath()
#' translate_ids_multi(ia, source = uniprot, target, ensp, ensembl = TRUE)
#'
#' @importFrom magrittr %>% %<>% extract
#' @importFrom purrr map map_chr pluck reduce discard
#' @importFrom rlang set_names !!! !! := enquos
#' @seealso \code{\link{translate_ids}}
#' @export
translate_ids_multi <- function(
    d,
    ...,
    suffixes = NULL,
    suffix_sep = '_',
    uploadlists = FALSE,
    ensembl = FALSE,
    hmdb = FALSE,
    ramp = FALSE,
    chalmers = FALSE,
    entity_type = NULL,
    keep_untranslated = TRUE,
    organism = 9606,
    reviewed = TRUE,
    expand = TRUE
) {

    ids <- ellipsis_to_char(...)
    raw_ids <- enquos(...) %>% map(.nse_ensure_str) %>% unlist
    source_cols <- intersect(names(ids), colnames(d))

    if (length(source_cols) == 0L && length(suffixes > 0L)) {

        source_cols <-
            names(ids)[1L] %>%
            paste(suffixes, sep = suffix_sep) %>%
            intersect(colnames(d))

        target_cols <- names(ids) %>% tail(-1L)

        from_types <-
            source_cols %>%
            length %>%
            rep(ids[1L], .) %>%
            unlist %>% unname

    } else {

        target_cols <- names(ids) %>% setdiff(colnames(d))

        default_id <-
            raw_ids %>%
            extract(source_cols) %>%
            discard(is.na) %>%
            first

        from_types <-
            source_cols %>%
            map_chr(~pluck(raw_ids, .x, .default = default_id))

    }

    if(length(target_cols) == 0L) {

        return(d)

    }

    if(length(source_cols) == 0L) {

        msg <- 'translate_ids_multi: no source column provided.'
        log_error_with_info(msg)
        stop(msg)

    }

    to_types <- ids %>% extract(target_cols) %>% unlist %>% unname
    suffixes %<>% if_null(`if`(length(source_cols) == 1L, '', source_cols))
    suffix_sep %<>% {`if`(length(suffixes) == 1L && suffixes == '', '', .)}

    if(length(suffixes) != length(source_cols)) {

        msg <-
            paste0(
                'translate_ids_multi: number of suffixes (%i) does not ',
                'match number of source columns (%i).'
            ) %>%
            sprintf(length(suffixes), length(source_cols))

        log_error_with_info(msg)
        stop(msg)

    }

    source_cols %>%
    seq_along %>%
    reduce(
        ~translate_ids(
            .x,
            !!sym(source_cols[.y]) := !!sym(from_types[.y]),
            !!!syms(
                to_types %>%
                set_names(paste0(target_cols, suffix_sep, suffixes[.y]))
            ),
            uploadlists = uploadlists,
            ensembl = ensembl,
            hmdb = hmdb,
            ramp = ramp,
            chalmers = chalmers,
            entity_type = entity_type,
            keep_untranslated = keep_untranslated,
            organism = organism,
            reviewed = reviewed,
            expand = expand
        ),
        .init = d
    )

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
    from = 'accession',
    reviewed = TRUE,
    organism = 9606
){

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
    to_ens <- to == 'xref_ensembl'
    from_ens <- from == 'xref_ensembl'

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
    separate_rows(From, sep = '[; ]') %>%
    separate_rows(To, sep = '[; ]') %>%
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
#' @param ensembl Logical: use data from Ensembl BioMart instead of UniProt.
#' @param hmdb Logical: use HMDB ID translation data.
#' @param ramp Logical: use RaMP ID translation data.
#' @param chalmers Logical: use ID translation data from Chalmers Sysbio GEM.
#' @param entity_type Character: "gene" and "smol" are short symbols for
#'     proteins, genes and small molecules respectively. Several other synonyms
#'     are also accepted.
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
    hmdb = FALSE,
    ramp = FALSE,
    chalmers = FALSE,
    entity_type = NULL,
    identifiers = NULL,
    organism = 9606L,
    reviewed = TRUE
){

    # NSE vs. R CMD check workaround
    To <- From <- NULL

    from <- .nse_ensure_str(!!enquo(from))
    to <- .nse_ensure_str(!!enquo(to))

    if(ensembl) {

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

    } else if(chalmers) {

        log_trace(
            'ID translation table: from `%s` to `%s`, using Chalmers GEM.',
            from, to
        )

        result <-
            chalmers_gem_id_mapping_table(
                to = !!sym(to),
                from = !!sym(from),
                organism = organism
            )

    } else if(hmdb) {

        log_trace(
            'ID translation table: from `%s` to `%s`, using HMDB.',
            from, to
        )

        result <-
            hmdb_id_mapping_table(
                to = !!sym(to),
                from = !!sym(from),
                entity_type = entity_type
            )

    } else if(ramp) {

        log_trace(
            'ID translation table: from `%s` to `%s`, using RaMP',
            from, to
        )

        result <-
            ramp_id_mapping_table(
                to = !!sym(to),
                from = !!sym(from)
            )

    } else if(
        uploadlists || (
            (
                !id_type_in(from, 'uniprot') ||
                !id_type_in(to, 'uniprot')
            ) &&
            id_type_in(from, 'uploadlists') &&
            id_type_in(to, 'uploadlists')
        )
    ) {

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

    } else {

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
    hmdb_old_ids(from, to) %>%
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

    id_type %in% names(omnipathr.env$id_types[[service]]) ||
    id_type %in% omnipathr.env$id_types[[service]]

}


#' Read ID type information
#'
#' @importFrom magrittr %>%
#'
#' @noRd
.load_id_types <- function(pkgname){

    omnipathr.env$id_types <-
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
#'     \item{\code{\link{hmdb_id_mapping_table}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#' }
ensembl_id_mapping_table <- function(
    to,
    from = 'uniprot',
    organism = 9606
){

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


#' Identifier translation table from HMDB
#'
#' @param to Character or symbol: target ID type. See Details for possible
#'     values.
#' @param from Character or symbol: source ID type. See Details for possible
#'     values.
#' @param entity_type Character: "gene" and "smol" are short symbols for
#'     proteins, genes and small molecules respectively. Several other synonyms
#'     are also accepted.
#'
#' @return A data frame (tibble) with columns `From` and `To`.
#'
#' @details The arguments \code{to} and \code{from} can be provided either
#' as character or as symbol (NSE). Their possible values are either HMDB XML
#' tag names or synonyms listed at \code{\link{id_types}}.
#'
#' @examples
#' hmdb_kegg <- hmdb_id_mapping_table("kegg", "hmdb")
#' hmdb_kegg
#'
#' @importFrom dplyr recode
#' @importFrom magrittr %>%
#' @importFrom rlang enquo !! !!!
#' @importFrom tidyselect everything
#' @importFrom tidyr unnest_longer
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{translate_ids}}}
#'     \item{\code{\link{id_types}}}
#'     \item{\code{\link{hmdb_table}}}
#'     \item{\code{\link{uniprot_full_id_mapping_table}}}
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#'     \item{\code{\link{ensembl_id_mapping_table}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#' }
hmdb_id_mapping_table <- function(to, from, entity_type = 'metabolite') {

    .slow_doctest()

    DATASETS <- c(
        small_molecule = 'metabolites',
        protein = 'proteins'
    )

    to <-
        .nse_ensure_str(!!enquo(to)) %>%
        hmdb_id_type()
    from <-
        .nse_ensure_str(!!enquo(from)) %>%
        hmdb_id_type()

    hmdb_id_mapping_table_impl <- function(from, to) {
        entity_type %>%
        ensure_entity_type %>%
        recode(!!!DATASETS) %>%
        hmdb_table(fields = c(from, to)) %>%
        set_names(c('From', 'To')) %>%
        unnest_longer(everything()) %>%
        trim_and_distinct
    }

    with_cache(
        name = 'hmdb_id_mapping_table',
        args = list(from, to),
        callback = hmdb_id_mapping_table_impl
    )

}


#' Take care about 5 and 7 digit HMDB IDs
#'
#' @importFrom dplyr mutate bind_rows
#' @importFrom magrittr %>% extract
#' @importFrom stringr str_to_lower str_detect str_match
#' @noRd
hmdb_old_ids <- function(d, from, to){

    is_hmdb <- function(i) {i %>% str_to_lower %>% str_detect('^hmdb')}

    fmt_hmdb <- function(v, fmt = 7L) {
        fmt <- 'HMDB%%0%ii' %>% sprintf(fmt)

        v %>%
        str_match('HMDB(\\d+)') %>%
        extract(,2L) %>%
        as.integer %>%
        sprintf(fmt, .)
    }

    d %>%
    doif(
        is_hmdb(from),
        ~bind_rows(
            mutate(.x, From = fmt_hmdb(From)),
            mutate(.x, From = fmt_hmdb(From, 5L))
        )
    ) %>%
    doif(
        is_hmdb(to),
        ~mutate(.x, To = fmt_hmdb(To))
    )

}


#' Pairwise ID translation table from RaMP database

#' @param from Character or Symbol. Name of an identifier type.
#' @param to Character or Symbol. Name of an identifier type.
#' @param version Character. The version of RaMP to download.
#'
#' @return Dataframe of pairs of identifiers.
#'
#' @examples
#' ramp_id_mapping_table('hmdb', 'kegg')
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !! enquo
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{ramp_sqlite}}}
#'     \item{\code{\link{ramp_tables}}}
#'     \item{\code{\link{ramp_table}}}
#'     \item{\code{\link{translate_ids}}}
#'     \item{\code{\link{id_types}}}
#'     \item{\code{\link{hmdb_table}}}
#'     \item{\code{\link{uniprot_full_id_mapping_table}}}
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#'     \item{\code{\link{ensembl_id_mapping_table}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#' }
ramp_id_mapping_table <- function(from, to, version = '2.5.4') {

    .slow_doctest()

    from <- .nse_ensure_str(!!enquo(from)) %>% ramp_id_type
    to <- .nse_ensure_str(!!enquo(to)) %>% ramp_id_type

    with_cache(
        name = 'ramp_id_mapping_table',
        args = list(from, to, version),
        callback = ramp_id_mapping_table_impl
    )

}


#' Metabolite ID translation tables from Chalmers Sysbio
#'
#' @param to Character: type of ID to translate to, either label used
#'   internally in this package, or a column name from "metabolites.tsv"
#'   distributed by Chalmers Sysbio. NSE is supported.
#' @param from Character: type of ID to translate from, same format as "to".
#' @param organism Character or integer: name or identifier of the organism.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Tibble with two columns, "From" and "To", with the corresponding ID
#'   types.
#'
#' @examples
#' chalmers_gem_id_mapping_table('metabolicatlas', 'hmdb')
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !! enquo sym
#' @importFrom dplyr select distinct filter mutate if_any across
#' @importFrom tidyselect everything
#' @importFrom tidyr separate_longer_delim
#' @export
chalmers_gem_id_mapping_table <- function(
        to,
        from = 'metabolicatlas',
        organism = 'Human'
    ) {

    .slow_doctest()

    to <- .nse_ensure_str(!!enquo(to)) %>% chalmers_gem_id_type
    from <- .nse_ensure_str(!!enquo(from)) %>% chalmers_gem_id_type

    organism %>%
    chalmers_gem_metabolites %>%
    select(From = !!sym(from), To = !!sym(to)) %>%
    filter(!if_any(everything(), is.na)) %>%
    mutate(across(everything(), as.character)) %>%
    separate_longer_delim(everything(), delim = ';') %>%
    distinct

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
#'     \item{\code{\link{chalmers_gem_id_type}}}
#'     \item{\code{\link{hmdb_id_type}}}
#' }
ensembl_id_type <- function(label){

    resource_id_type(label, 'ensembl')

}


#' HMDB identifier type label
#'
#' @param label Character: an ID type label, as shown in the table at
#'     \code{\link{translate_ids}}
#'
#' @return Character: the HMDB specific ID type label, or the input
#'     unchanged if it could not be translated (still might be a valid
#'     identifier name). These labels should be valid HMDB field
#'     names, as used in HMDB XML files.
#'
#' @examples
#' hmdb_id_type("hmdb")
#' # [1] "accession"
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_id_type}}}
#'     \item{\code{\link{uniprot_id_type}}}
#'     \item{\code{\link{ensembl_id_type}}}
#'     \item{\code{\link{uploadlists_id_type}}}
#' }
hmdb_id_type <- function(label){

    resource_id_type(label, 'hmdb')

}


#' RaMP identifier type label
#'
#' @param label Character: an ID type label, as shown in the table returned
#'     by \code{\link{id_types}}
#'
#' @return Character: the RaMP specific ID type label, or the input
#'     unchanged if it could not be translated (still might be a valid
#'     identifier name). These labels should be valid value
#'     names, as used in RaMP SQL database.
#'
#' @examples
#' ramp_id_type("rhea")
#' # [1] "rhea-comp"
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_id_type}}}
#'     \item{\code{\link{uniprot_id_type}}}
#'     \item{\code{\link{ensembl_id_type}}}
#'     \item{\code{\link{uploadlists_id_type}}}
#' }
ramp_id_type <- function(label){

    resource_id_type(label, 'ramp')

}


#' Metabolite identifier type label used in Chalmers Sysbio GEM
#'
#' @param label Character: an ID type label, as shown in the table at
#'     \code{\link{translate_ids}}
#'
#' @return Character: the Chalmers GEM specific ID type label, or the input
#'     unchanged if it could not be translated (still might be a valid
#'     identifier name). These labels should be column names from the
#'     "metabolites.tsv" distributed with the GEMs.
#'
#' @examples
#' chalmers_gem_id_type("metabolicatlas")
#' # [1] "metsNoComp"
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{hmdb_id_type}}}
#'     \item{\code{\link{uniprot_id_type}}}
#'     \item{\code{\link{ensembl_id_type}}}
#'     \item{\code{\link{uploadlists_id_type}}}
#' }
chalmers_gem_id_type <- function(label){

    resource_id_type(label, 'chalmers_gem')

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
#' @param side Character: either "from" or "to": direction of the mapping.
#'
#' @return Character: the UniProt Uploadlists specific ID type label, or
#'     the input unchanged if it could not be translated (still might be
#'     a valid identifier name). This is the label that one can use in
#'     UniProt Uploadlists (ID Mapping) queries.
#'
#' @examples
#' ensembl_id_type("entrez")
#' # [1] "GeneID"
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr recode
#' @importFrom rlang !!!
#' @seealso \itemize{
#'     \item{\code{\link{ensembl_id_type}}}
#'     \item{\code{\link{uniprot_id_type}}}
#'     \item{\code{\link{hmdb_id_type}}}
#'     \item{\code{\link{chalmers_gem_id_type}}}
#' }
uploadlists_id_type <- function(label, side = 'from'){

    label %>%
    recode(!!!omnipathr.env$id_types[[sprintf('uniprot_%s', side)]]) %>%
    resource_id_type('uploadlists')

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

    label %>% recode(!!!omnipathr.env$id_types[[resource]])

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

    omnipathr.env$id_types %>%
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


#' ID types and synonyms in identifier translation
#'
#' @return Data frame with 4 columns: the ID type labels in the resource, their
#'     synonyms in OmniPath (this package), the name of the ID translation
#'     resource, and the entity type.
#'
#' @examples
#' id_types()
#'
#' @importFrom magrittr %>% is_in
#' @importFrom dplyr mutate case_when
#' @importFrom tidyr unnest_longer separate_longer_delim
#' @importFrom tibble tibble
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{translate_ids}}}
#'     \item{\code{\link{translate_ids_multi}}}
#'     \item{\code{\link{ensembl_id_mapping_table}}}
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#'     \item{\code{\link{hmdb_id_mapping_table}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#'     \item{\code{\link{uniprot_full_id_mapping_table}}}
#'     \item{\code{\link{ensembl_id_type}}}
#'     \item{\code{\link{uniprot_id_type}}}
#'     \item{\code{\link{hmdb_id_type}}}
#'     \item{\code{\link{chalmers_gem_id_type}}}
#' }
id_types <- function() {

    # NSE vs. R CMD check workaround
    idtype <- entity_type <- NULL

    hmdb_met <- hmdb_metabolite_fields()
    hmdb_pro <- hmdb_protein_fields()
    hmdb_com <- hmdb_met %>% intersect(hmdb_pro)
    ramp_met <- RAMP_METABOLITE_ID_TYPES
    ramp_pro <- RAMP_GENE_ID_TYPES
    ramp_com <- ramp_met %>% intersect(ramp_pro)

    omnipathr.env$id_types %>%
    {tibble(idtype = list(.))} %>%
    unnest_longer(idtype, indices_to = 'resource') %>%
    unnest_longer(
        idtype,
        indices_to = 'in_omnipath',
        values_to = 'in_resource'
    ) %>%
    mutate(
        entity_type = case_when(
            resource == 'hmdb' & in_resource %in% hmdb_com ~ 'metabolite;protein',
            resource == 'hmdb' & in_resource %in% hmdb_met ~ 'metabolite',
            resource == 'ramp' & in_resource %in% ramp_com ~ 'metabolite;protein',
            resource == 'ramp' & in_resource %in% ramp_met ~ 'metabolite',
            resource == 'chalmers_gem' ~ 'metabolite',
            TRUE ~ 'protein'
        )
    ) %>%
    separate_longer_delim(entity_type, delim = ';')

}


#' Translate identifiers within complexes
#'
#' Here we attempt to provide ID translation for complexes by translating the
#' identifiers of their members one by one. Importantly, iff any of the members
#' is not possible to translate (i.e. missing from the translation table), the
#' translation of the whole complex will be omitted, in order to avoid creating
#' complexes with missing members. One to many mappings of identifiers result a
#' combinatorial expansion, which can dramatically inflate the number of
#' records. By default this expansion is disabled, and only the first (randomly
#' choosen) target identifier will be used for each member. This function is
#' not aware of any parameter of the translation, it is applicable for any
#' organism, entity type, and even for translation by orthologous gene pairs.
#'
#' @param d A data frame.
#' @param ... Names of columns containing complex identifiers.
#' @param mapping Data frame: the ID translation table with two columns, the
#'     source and target identifiers.
#' @param one_to_many Logical or numeric: allow combinatorial expansion or use
#'     only the first target identifier for each member of each complex. If
#'     numeric, it will be the maximum number of complex. The value `TRUE`
#'     corresponds to 12, i.e. one complex in `d` can be translated to a
#'     maximum of 12 complexes, which is the default.  If \code{NULL}, the
#'     option `omnipathr.complex_translation_one_to_many` will be used.
#'
#' @param Data frame: the ID translation table with translation for complexes
#'     added to it.
#'
#' @importFrom magrittr %<>% %>% not
#' @importFrom dplyr select mutate left_join group_by filter row_number
#' @importFrom dplyr summarize first across pull
#' @importFrom stringr str_replace
#' @importFrom tibble tibble
#' @importFrom tidyr separate_longer_delim unnest_longer unite expand_grid
#' @importFrom purrr map
#' @importFrom tidyselect everything
#' @noRd
translate_complexes <- function(d, ..., mapping, one_to_many = NULL) {

    # NSE vs. R CMD check workaround
    from <- to <- From <- To <- ids <- original <-
    complex_id <- member_id <- comp <- NULL

    one_to_many %<>%
        if_null(
            getOption('omnipathr.complex_translation_one_to_many')
        ) %>%
        {`if`(identical(., TRUE), 12L, .)} %>%
        {`if`(identical(., FALSE), ., 1L)}

    mapping %<>% set_names(c('from', 'to'))

    d %>%
    select(...) %>%
    unlist %>%
    unname %>%
    unique %>%
    keep(str_starts(., 'COMPLEX:')) %>%
    str_replace('^COMPLEX:', '') %>%
    tibble(ids = .) %T>%
    {log_trace('Translating complexes: %i complexes in data.', nrow(.))} %>%
    mutate(
        original = ids,
        complex_id = row_number()
    ) %>%
    separate_longer_delim(ids, '_') %>%
    mutate(member_id = row_number()) %>%
    left_join(mapping, by = c('ids' = 'from')) %>%
    group_by(complex_id) %>%
    filter(to %>% is.na %>% any %>% not) %>%
    filter(table(member_id) %>% prod <= one_to_many) %T>%
    {log_trace(
        paste0(
           '%i complexes after removing the ones mapping to ',
           'more than %i items in target identifier space.'
        ),
        n_distinct(.$original),
        one_to_many
    )} %>%
    summarize(
        to =
            split(to, member_id) %>%
            do.call(expand_grid, .) %>%
            list,
        from =
            split(ids, member_id) %>%
            do.call(expand_grid, .) %>%
            list,
        original = first(original)
    ) %>%
    mutate(
        to = map(to, ~unite(.x, comp, sep = '_') %>% pull(comp)),
        from = map(from, ~unite(.x, comp, sep = '_') %>% pull(comp))
    ) %>%
    unnest_longer(to, values_to = 'To') %>%
    unnest_longer(from, values_to = 'From') %>%
    select(From, To) %>%
    mutate(across(everything(), ~sprintf('COMPLEX:%s', .x))) %T>%
    {log_trace(
        'Translated %i complexes to %i.',
        n_distinct(.$From),
        n_distinct(.$To)
    )}

}


#' List available ID translation resources
#'
#' @return A character vector with the names of the available ID translation
#'     resources.
#'
#' @examples
#' id_translation_resources()
#'
#' @export
id_translation_resources <- function() {

    names(omnipathr.env$id_types)

}
