#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
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


KEGG_DATABASES <- c(
    'pathway',
    'brite',
    'module',
    'ko',
    'genome',
    'vg',
    'vp',
    'ag',
    'compound',
    'glycan',
    'reaction',
    'rclass',
    'enzyme',
    'disease',
    'drug',
    'dgroup',
    'genes',
    'ligand',
    'kegg',
    'network',
    'variant',
    'organism'
)

KEGG_DB_GENES <- c(
    'genes'
)

KEGG_OUTSIDE_DB_GENES <- c(
    'ncbi-geneid',
    'ncbi-proteinid',
    'uniprot'
)

KEGG_DB_COMPOUNDS <- c(
    'compound',
    'glycan',
    'drug'
)

KEGG_OUTSIDE_DB_COMPOUNDS <- c(
    'pubchem',
    'chebi'
)

KEGG_OUTSIDE_DB <- c(
    'pubmed',
    'atc',
    'jtc',
    'ndc',
    'yk'
)

KEGG_DB_DRUGS <- c(
    'drug'
)

KEGG_OUTSIDE_DB_DRUGS <- c(
    'ndc',
    'yj'
)

KEGG_OUTSIDE_DB_DRUG_TC <- c(
    'atc',
    'jtc'
)

KEGG_OUTSIDE_DB_ALL <- c(
    KEGG_OUTSIDE_DB,
    KEGG_OUTSIDE_DB_GENES,
    KEGG_OUTSIDE_DB_COMPOUNDS,
    KEGG_OUTSIDE_DB_DRUGS
)

KEGG_DATABASES_JA <- c(
    'disease_ja',
    'compound_ja',
    'dgroup_ja',
    'drug_ja'
)

KEGG_GET_FORMATS <- c(
    'aaseq',
    'ntseq',
    'mol',
    'kcf',
    'image',
    'conf',
    'kgml',
    'json'
)

KEGG_LINK_FORMATS <- c(
    'n-triple',
    'turtle'
)

KEGG_ORGANISM <- '<org>'

KEGG_DB_ORG <- KEGG_DATABASES %>% c(KEGG_ORGANISM)

KEGG_DB_GLKO <- c(
    'genes',
    'ligand',
    'kegg',
    'organism'
)

KEGG_DB_WO_GLKO <- KEGG_DB_ORG %>% setdiff(KEGG_DB_GLKO)

KEGG_COLUMNS <- list(
    list = c('id', 'name'),
    find = c('id', 'value'),
    conv = c('id_a', 'id_b'),
    link = c('id_a', 'id_b'),
    ddi = c('drug_a', 'drug_b', 'interaction', 'mechanism')
)

KEGG_READERS <- list()

KEGG_KID_PREFIXES <- list(
    pathway = c('map', 'ko', 'ec', 'rn', '<org>'),
    brite = c('br', 'jp', 'ko', '<org>'),
    module = c('M', '<org>_M'),
    orthology = 'K',
    genome = 'T',
    compound = 'C',
    glycan = 'G',
    reaction = 'R',
    rclass = 'RC',
    network = 'N',
    disease = 'H',
    drug = 'D',
    dgroup = 'DG'
)

KEGG_API <- list(
    # the order of templates sometimes matters, see `kegg_query` for details;
    # the more specific ones should come first, to avoid a less specific one
    # to match the query when it could match a more specific template.
    info = list(
        list(
            database = KEGG_DB_ORG
        )
    ),
    list = list(
        list(
            database = 'pathway',
            organism = KEGG_ORGANISM
        ),
        list(
            database = 'brite',
            option = c('br', 'jp', 'ko', KEGG_ORGANISM)
        ),
        list(
            database = KEGG_DB_ORG %>% setdiff(c('genes', 'ligand', 'kegg'))
        ),
        list(
            dbentries = KEGG_DB_WO_GLKO
        )
    ),
    find = list(
        list(
            database = c('compound', 'drug'),
            query = NULL,
            option = c('formula', 'exact_mass', 'mol_weight', 'nop')
        ),
        list(
            database = KEGG_DB_ORG %>% setdiff(c('kegg', 'organism')),
            query = NULL
        )
    ),
    get = list(
        list(
            dbentries = KEGG_DB_WO_GLKO %>% c(KEGG_DATABASES_JA),
            option = KEGG_GET_FORMATS
        )
    ),
    conv = list(
        list(
            target_db = KEGG_OUTSIDE_DB_GENES,
            source_db = KEGG_ORGANISM
        ),
        list(
            target_db = KEGG_OUTSIDE_DB_COMPOUNDS,
            source_db = KEGG_DB_COMPOUNDS
        ),
        list(
            target_db = c(KEGG_ORGANISM, KEGG_OUTSIDE_DB_GENES, KEGG_DB_GENES),
            dbentries = c(KEGG_ORGANISM, KEGG_OUTSIDE_DB_GENES, KEGG_DB_GENES)
        ),
        list(
            target_db = c(KEGG_OUTSIDE_DB_COMPOUNDS, KEGG_DB_COMPOUNDS),
            dbentries = c(KEGG_OUTSIDE_DB_COMPOUNDS, KEGG_DB_COMPOUNDS)
        )
    ),
    link = list(
        list(
            target_db = c(KEGG_DB_DRUGS, KEGG_OUTSIDE_DB_DRUG_TC),
            source_db = c(KEGG_DB_DRUGS, KEGG_OUTSIDE_DB_DRUG_TC),
            option = KEGG_LINK_FORMATS
        ),
        list(
            target_db = c(KEGG_DB_DRUGS, KEGG_OUTSIDE_DB_DRUG_TC),
            dbentries = c(KEGG_DB_DRUGS, KEGG_OUTSIDE_DB_DRUG_TC),
            option = KEGG_LINK_FORMATS
        ),
        list(
             target_db = KEGG_DB_WO_GLKO %>% c(KEGG_OUTSIDE_DB),
             source_db = KEGG_DB_WO_GLKO %>% c(KEGG_OUTSIDE_DB)
        ),
        list(
             target_db = KEGG_DB_WO_GLKO %>% c(KEGG_OUTSIDE_DB),
             dbentries = KEGG_DB_WO_GLKO %>% c(KEGG_OUTSIDE_DB)
        )
    ),
    ddi = list(
        list(
            dbentries = KEGG_DB_DRUGS %>% c(KEGG_OUTSIDE_DB_DRUGS)
        )
    )
)


#' List of databases (endpoints) in the KEGG REST API
#'
#' @return A character vector of KEGG databases.
#'
#' @examples
#' kegg_databases()
#'
#' @export
kegg_databases <- function() {

    KEGG_DATABASES

}


#' List of operations in the KEGG REST API
#'
#' @return A character vector of KEGG operations.
#'
#' @examples
#' kegg_operations()
#'
#' @export
kegg_operations <- function() {

    KEGG_API %>% names

}


#' Compile a query for the KEGG REST API
#'
#' @importFrom magrittr %>% %<>% extract2 extract
#' @importFrom rlang list2
#' @importFrom logger log_warn log_error
#' @noRd
kegg_query <- function(operation, ...) {

    op_templates <- KEGG_API %>% extract2(operation)

    if (is.null(op_templates)) {

        msg <- 'Invalid operation: %s' %>% sprintf(operation)
        log_error(msg)
        stop(msg)

    }

    bad_queries <- list()

    for (op_template in op_templates) {

        query <- kegg_query_match_template(op_template, operation, ...)

        if (query$complete && is.null(query$error)) {

            return(query)

        }

        bad_queries %<>% c(list(query))

    }

    errors <- c(
        sprintf(
            paste0(
                'kegg_query: failed to match arguments against ',
                'any query template for operation `%s`.'
            ),
            operation
        ),
        sprintf('The arguments provided: %s.', compact_repr(list2(...)))
    )

    for (query in bad_queries) {

        errors %<>% c(
            sprintf(
                'Errors for template `%s`:',
                paste0(query$names, collapse = '/')
            )
        )

        for (error in query$error) {

            errors %<>% c(error)

        }

    }

    msg <- errors %>% paste0(collapse = ' ')
    log_error(msg)
    stop(msg)

}


#' Attempts to compile a KEGG query based on a template
#'
#' @param template A KEGG API query template, as defined in the KEGG_API
#'     constant in this package.
#' @param operation Character: the KEGG API operation.
#' @param ... Arguments to the operation. Names are optional, unnamed elements
#'     must follow the order of the arguments in the template.
#'
#' @return A list with the following elements: \itemize{
#'     \item operation - The KEGG API operation.
#'     \item names - The names of the arguments.
#'     \item query - The values of the arguments.
#'     \item error - Error messages.
#'     \item complete - Whether the query has all mandatory arguments.
#' }
#'
#' @importFrom magrittr %>% %<>% inset extract extract2
#' @importFrom rlang list2 set_names
#' @importFrom stringr str_replace str_detect str_extract
#' @noRd
kegg_query_match_template <- function(template, operation, ...) {

    args <- list2(...)
    argnames <- args %>% (vctrs%:::%minimal_names)
    iarg_noname <- argnames %>% {nchar(.)  == 0L} %>% which
    template_names <- template %>% names
    optional <- template_names %>% str_detect('^\\.')
    template_names %<>% str_replace('^\\.', '')

    result <- list(
        operation = operation,
        names = NULL,
        query = NULL,
        error = NULL,
        complete = FALSE
    )

    args_noname <-
        template_names %>%
        setdiff(argnames) %>%
        extract(seq_along(iarg_noname))

    iargs <-
        argnames %>%
        inset(iarg_noname, args_noname) %>%
        set_names(seq_along(.), .)

    args_required <-
        template_names %>%
        extract(optional %>% not %>% which)

    args_missing <-
        args_required %>%
        setdiff(iargs %>% names)

    if (!is_empty(args_missing)) {

        msg <-
            'Missing mandatory argument(s) for KEGG operation `%s`: %s.' %>%
            sprintf(operation, args_missing %>% paste0(collapse = ', '))

        result$error %<>% c(msg)
        return(result)

    }

    for (i in seq_along(template)) {

        argname <- template_names[i]
        iarg <- iargs %>% extract(argname)
        arg_template <- template %>% extract2(i)

        if (is.na(iarg)) {

            if (optional[i]) {

                next

            } else {

                msg <-
                    'Missing argument for KEGG `%s` operation: %s.' %>%
                    sprintf(operation, argname)

                result$error %<>% c(msg)
                break

            }

        }

        arg <- args %>% extract2(iarg)

        if (argname == 'dbentries') {

            idtype <- argnames[iarg]

            if (!is_empty_2(idtype) && !idtype == 'dbentries') {

                arg %<>% str_replace('^(?:[-\\w]+:)?', sprintf('%s:', idtype))

            }

            result$error %<>% c(kegg_check_prefixes(arg, arg_template))
            arg %<>% paste0(collapse = '+')

        } else if (arg %>% kegg_valid_arg(arg_template) %>% not) {

            msg <-
                kegg_invalid_arg_msg(
                    operation,
                    arg,
                    argname,
                    arg_template
                )

            result$error %<>% c(msg)

        }

        result$query %<>% c(arg)
        result$names %<>% c(argname)

    }

    result$complete <- template_names %>% setdiff(result$names) %>% is_empty

    result

}


#' Checks whether all prefixes are valid in <dbentries>
#'
#' @param arg Character: the value of the argument.
#' @param arg_template Character: the template of the argument.
#'
#' @return Character or NULL: NULL if all prefixes are valid; an error message
#'     otherwise.
#'
#' @importFrom stringr str_detect str_extract str_replace
#' @importFrom magrittr %>% is_in
#' @importFrom purrr map
#' @noRd
kegg_check_prefixes <- function(arg, arg_template) {

    msg <- NULL

    valid_prefixes <-
        map(arg_template, ~KEGG_KID_PREFIXES[[.x]]) %>%
        unlist %>%
        {doif(
            str_detect(., KEGG_ORGANISM),
            map(
                .,
                ~str_replace(.x, KEGG_ORGANISM, kegg_organism_codes())
            )
        )} %>%
        unlist %>%
        c(arg_template %>% intersect(KEGG_OUTSIDE_DB_ALL))

    if (length(valid_prefixes) > 0L) {

        which_valid <-
            arg %>%
            str_extract('^[-\\A-z]+:?') %>%
            is_in(valid_prefixes) %>%
            which

        if (length(which_valid) > 0L) {

            noprefix <- arg %>% extract(which_valid %>% head)
            msg <- sprintf(
                paste0(
                    'In KEGG, <dbentries> must be prefixed with an ',
                    'ID type (database or KEGG organism code), which ',
                    'also can be provided here as an argument name.
                    In this case, valid prefixes are %s. A sample ',
                    'of the items missing this prefix: %s.'
                ),
                valid_prefixes %>% compact_repr,
                noprefix %>% compact_repr
            )

        }

    }

    msg

}


#' Match an argument against its template
#'
#' @param arg Character: the argument to match.
#' @param arg_template Character or NULL: the template to match against.
#'
#' @return Logical: whether the argument matches the template.
#' @noRd
kegg_valid_arg <- function(arg, arg_template) {

    is.null(arg_template) ||
    arg %in% arg_template ||
    (
        KEGG_ORGANISM %in% arg_template &&
        arg %in% kegg_organism_codes()
    )

}


#' Generate an invalid argument message
#'
#' @param operation Character: one of the KEGG REST API operations.
#' @param argname Character: the name of the argument.
#' @param arg Character: the value of the argument.
#' @param arg_template Character: the template of the argument.
#'
#' @return Character: an error message.
#' @noRd
kegg_invalid_arg_msg <- function(
        operation,
        argname,
        arg,
        arg_template
    ) {

    sprintf(
        paste0(
            '`%s` is an invalid value for argument <%s> in ',
            'KEGG operation `%s`. Valid values are %s%s.'
        ),
        arg,
        argname,
        operation,
        arg_template %>% paste0(collapse = ', '),
        `if`(
            KEGG_ORGANISM %in% arg_template,
            ' and KEGG organism codes',
            ''
        )
    )

}


#' Compile a KEGG REST API path
#'
#' @param operation Character: one of the KEGG REST API operations.
#' @param ... Further parameters to \code{\link{kegg_query}}.
#'
#' @return Character: a KEGG REST API path, without the base URL.
#' @noRd
kegg_api_path <- function(operation, ...) {

    query <- kegg_query(operation, ...)

    sprintf(
        '%s/%s',
        operation,
        query$query %>% paste0(collapse = '/')
    )

}


#' Perform a KEGG REST API request
#'
#' @param operation Character: one of the KEGG REST API operations.
#' @param ... Further parameters to \code{\link{kegg_request}}.
#'
#' @export
kegg_request <- function(operation, ...) {

    kegg_api_path(operation, ...) %>%
    {generic_downloader(url_key = 'kegg_rest', url_param = list(.))}

}



#' List of organisms in KEGG
#'
#' @return A data frame (tibble) with organism data.
#'
#' @examples
#' kegg_organisms()
#'
#' @export
kegg_organisms <- function() {

   get_db('kegg_organisms')

}


#' Load the KEGG organisms table
#'
#' @noRd
kegg_organisms_load <- function() {

    kegg_request('list', 'organism')

}


#' All 3 letter organism code from KEGG
#'
#' @return A character vector with all 3 letter codes.
#'
#' @examples
#' kegg_organism_codes()
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @export
kegg_organism_codes <- function() {

    kegg_organisms() %>% pull(2L)

}
