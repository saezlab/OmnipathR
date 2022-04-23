#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2022
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


#' Decorator for trying UniProt subdomains
#'
#' This has any relevance only in rare cases with OS networking issues.
#'
#' @noRd
uniprot_domains <- decorator %@% function(FUN){

    function(...){

        for(subd in c('www', 'pir3')){

            result <- tryCatch(
                FUN(..., .subdomain = subd),
                error = identity
            )

            if(!inherits(result, 'error')){

                break

            }

        }

        if(inherits(result, 'error')){

            stop(conditionMessage(result))

        }

        return(result)

    }

}


#' Retrieves an identifier translation table from the UniProt uploadlists
#' service
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
#' @importFrom magrittr %>%
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
    chunk_size = 5000
){

    from <- .nse_ensure_str(!!enquo(from))
    to <- .nse_ensure_str(!!enquo(to))

    identifiers %>%
    sort %>%
    chunks(5000) %>%
    map(.uniprot_id_mapping_table, from, to) %>%
    bind_rows() %T>%
    load_success()

}


#' R CMD check workaround, see details at \code{uniprot_id_mapping_table}
#'
#' @importFrom magrittr %T>% %<>%
#' @importFrom rlang !!!
#' @importFrom dplyr recode
#' @importFrom logger log_trace
#'
#' @noRd
.uniprot_id_mapping_table <- uniprot_domains %@% function(
    identifiers,
    from,
    to,
    .subdomain = 'www'
){

    id_types <- omnipath.env$id_types$uploadlists

    from %<>% recode(!!!id_types)
    to %<>% recode(!!!id_types)

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
#'     | OmnipathR      | Uploadlists          | UniProt query           |
#'     | -------------- | -------------------- | ----------------------- |
#'     | uniprot        | ACC                  | id                      |
#'     | uniprot_entry  | ID                   | entry name              |
#'     | genesymbol     | GENENAME             | genes(PREFERRED)        |
#'     | genesymbol_syn |                      | genes(ALTERNATIVE)      |
#'     | hgnc           | HGNC_ID              | database(HGNC)          |
#'     | entrez         | P_ENTREZGENEID       | database(geneid)        |
#'     | ensg           | ENSEMBLGENOME_ID     |                         |
#'     | enst           | ENSEMBL_TRS_ID       | database(ensembl)       |
#'     | ensp           | ENSEMBL_PRO_ID       |                         |
#'     | ensgt          | ENSEMBLGENOME_TRS_ID |                         |
#'     | ensgp          | ENSEMBLGENOME_PRO_ID |                         |
#'     | ensembl        | ENSEMBL_ID           |                         |
#'     | protein_name   |                      | protein names           |
#'     | refseqp        | P_REFSEQ_AC          | database(refseq)        |
#'     | refseqn        | REFSEQ_NT_ID         |                         |
#'     | embl           | EMBL                 | database(embl)          |
#'     | embl_id        | EMBL_ID              |                         |
#'     | gi             | P_GI                 |                         |
#'     | pir            | PIR                  |                         |
#'     | pdb            | PDB_ID               |                         |
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

    # NSE vs. R CMD check workaround
    To <- NULL

    ids <-
        enquos(...) %>%
        map(.nse_ensure_str) %>%
        set_names(names(.) %||% unlist(.)) %>%
        set_names(ifelse(nchar(names(.)), names(.), unlist(.)))

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


#' Try how non standard evaluation works
#'
#' @noRd
nse_test <- function(a, ..., b = FALSE){

    print(a)
    print(b)

    print(enquos(...) %>% map(.nse_ensure_str))

}


#' A table with all UniProt IDs
#'
#' Retrieves a table from UniProt with all proteins for a certain organism.
#'
#' @param fields Character vector of fields as defined by UniProt. For
#'     possible values please refer to
#'     \url{https://www.uniprot.org/help/uniprotkb\%5Fcolumn\%5Fnames}
#' @param reviewed Retrieve only reviewed (`TRUE`), only unreviewed (`FALSE`)
#'     or both (`NULL`).
#' @param organism Integer, NCBI Taxonomy ID of the organism (by default
#'     9606 for human).
#'
#' @return Data frame (tibble) with the requested UniProt entries and fields.
#'
#' @importFrom rlang exec !!!
#' @export
#'
#' @examples
#' human_swissprot_ac <- all_uniprots(fields = 'entry name')
#' human_swissprot_ac
#' # # A tibble: 20,396 x 1
#' #    `Entry name`
#' #    <chr>
#' #  1 OR4K3_HUMAN
#' #  2 O52A1_HUMAN
#' #  3 O2AG1_HUMAN
#' #  4 O10S1_HUMAN
#' #  5 O11G2_HUMAN
#' # # . with 20,386 more rows
all_uniprots <- function(fields = 'id', reviewed = TRUE, organism = 9606){

    exec(.all_uniprots, !!!as.list(environment()))

}


#' R CMD check workaround, see details at \code{all_uniprots}
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom logger log_trace
#'
#' @noRd
.all_uniprots <- uniprot_domains %@% function(
    fields = 'id',
    reviewed = TRUE,
    organism = 9606,
    .subdomain = 'www'
){

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
        sprintf(' AND reviewed:%s', `if`(reviewed, 'yes', 'no'))
    )

    generic_downloader(
        url_key = 'all_uniprots',
        url_param = list(.subdomain, fields, organism, reviewed),
        reader_param = list(progress = FALSE),
        resource = 'UniProt'
    ) %T>%
    load_success()

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
#' @importFrom dplyr mutate rename filter recode
#' @importFrom tidyr separate_rows
#' @importFrom rlang !! !!! enquo
#' @importFrom logger log_trace
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
#' @seealso \code{\link{translate_ids}}
uniprot_full_id_mapping_table <- function(
    to,
    from = 'id',
    reviewed = TRUE,
    organism = 9606
){

    # NSE vs. R CMD check workaround
    From <- To <- NULL

    get_field_name <- function(label){

        label %>% recode(!!!omnipath.env$id_types$uniprot)

    }

    strip_semicol <- function(v){sub(';$', '', v)}

    to <-
        .nse_ensure_str(!!enquo(to)) %>%
        get_field_name()
    from <-
        .nse_ensure_str(!!enquo(from)) %>%
        get_field_name()

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
    filter(!is.na(From) & !is.na(To))

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
#' @importFrom dplyr pull
#'
#' @noRd
id_translation_table <- function(
    from,
    to,
    uploadlists = FALSE,
    ensembl = FALSE,
    identifiers = NULL,
    organism = 9606L,
    reviewed = TRUE
){

    from <- .nse_ensure_str(!!enquo(from))
    to <- .nse_ensure_str(!!enquo(to))

    if(ensembl){

        log_trace(
            'ID translation table: from `%s` to `%s`, using Ensembl BioMart.',
            from, to
        )

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

        result <-
            identifiers %>%
            {
                . %||%
                all_uniprots(organism = organism, reviewed = reviewed) %>%
                pull(1)
            } %>%
            uniprot_id_mapping_table(!!sym(from), !!sym(to))

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

    return(result)

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
#' @importFrom jsonlite fromJSON
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
        fromJSON()

}
