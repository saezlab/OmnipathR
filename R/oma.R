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


#' Organism identifiers from the Orthologous Matrix
#'
#' @return A data frame with organism identifiers.
#'
#' @examples
#' oma_organisms()
#'
#' @importFrom readr cols
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export
#' @seealso \code{\link{ensembl_organisms}}
oma_organisms <- function() {

    # NSE vs. R CMD check workaround
    ncbi_tax_id <- NULL

    'oma_species' %>%
    generic_downloader(
        reader_param = list(
            col_names = c(
                'oma_code',
                'ncbi_tax_id',
                'latin_name',
                'genome_source',
                'oma_version'
            ),
            col_types = cols(),
            skip = 3L
        )
    ) %>%
    mutate(
        gtdb_tax_id = ifelse(ncbi_tax_id < 0L, -ncbi_tax_id, NA),
        ncbi_tax_id = ifelse(ncbi_tax_id < 0L, NA, ncbi_tax_id)
    )

}


#' OMA codes of all supported organisms
#'
#' We need this wrapper only to help `organism_for`, i.e. to provide a callable
#' without arguments that returns a simple character vector.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @noRd
oma_supported_organisms <- function() {

    # NSE vs. R CMD check workaround
    oma_code <- NULL

    oma_organisms() %>%
    pull(oma_code) %>%
    unique

}


#' Orthologous gene pairs between two organisms
#'
#' From the web API of Orthologous Matrix (OMA). Items which could not be
#' translated to `id_type` (but present in the data with their internal OMA
#' IDs) are removed.
#'
#' @param organism_a Name or identifier of an organism.
#' @param organism_b Name or identifier of another organism.
#' @param id_type The gene or protein identifier to use in the table. For a
#'     list of supported ID types see `omnipathr.env$id_types$oma`. In addition,
#'     "genesymbol" is supported, in this case
#'     \code{\link{oma_pairwise_genesymbols}} will be called automatically.
#' @param mappings Character vector: control ambiguous mappings: \itemize{
#'         \item{1:1 - unambiguous}
#'         \item{1:m - one-to-many}
#'         \item{n:1 - many-to-one}
#'         \item{n:m - many-to-many}
#'     }
#' @param only_ids Logical: include only the two identifier columns, not the
#'     mapping type and the orthology group columns.
#'
#' @return A data frame with orthologous gene pairs.
#'
#' @examples
#' oma_pairwise("human", "mouse", "uniprot")
#' # # A tibble: 21,753 × 4
#' #    id_organism_a id_organism_b mapping oma_group
#' #    <chr>         <chr>         <chr>       <dbl>
#' #  1 Q15326        Q8R5C8        1:1       1129380
#' #  2 Q9Y2E4        B2RQ71        1:1        681224
#' #  3 Q92615        Q6A0A2        1:1       1135087
#' #  4 Q9BZE4        Q99ME9        1:1       1176239
#' #  5 Q9BXS1        Q8BFZ6        1:m            NA
#' # # … with 21,743 more rows
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom readr cols
#' @importFrom utils tail
#' @importFrom dplyr across mutate filter
#' @importFrom stringr str_extract str_detect
#' @importFrom tidyselect starts_with
#' @importFrom rlang exec !!!
#' @export
oma_pairwise <- function(
    organism_a = 'human',
    organism_b = 'mouse',
    id_type = 'uniprot',
    mappings = c('1:1', '1:m', 'n:1', 'n:m'),
    only_ids = TRUE
) {

    .slow_doctest()

    if (!id_type %in% names(omnipathr.env$id_types$oma)) {
        return(
            environment() %>%
            as.list %>%
            exec(oma_pairwise_translated, !!!.)
        )
    }

    # NSE vs. R CMD check workaround:
    id_organism_a <- id_organism_b <- mapping <- NULL

    organism_a %<>% organism_for('oma')
    organism_b %<>% organism_for('oma')
    id_type %<>% oma_id_type

    args <- match.call() %>% as.list %>% tail(-1L)

    for(arg in names(args)) {
        if(is.null(get(arg)) || any(is.na(get(arg)))) {
            msg <- sprintf(
                'Could not recognize %s: `%s`.',
                `if`(arg == 'id_type', 'ID type', 'organism'),
                args[[arg]]
            )
            log_error_with_info(msg)
            stop(msg)
        }
    }

    'oma_pairwise_new' %>%
    generic_downloader(
        reader_param = list(
            col_names = c(
                'oma_code_a',
                'id_organism_a',
                'oma_code_b',
                'id_organism_b',
                'mapping',
                'oma_group'
            ),
            col_types = cols(),
            comment = '#',
            skip = 1L
        ),
        url_param = list(organism_a, organism_b, id_type)
    ) %T>%
    {`if`(
        ncol(.) != 6L,
        {
            msg <- sprintf(
                paste0(
                    'Unexpected number of columns in OMA data: ',
                    'most likely the `id_type` (`%s`) is not correct.'
                ),
                id_type
            )
            log_error_with_info(msg)
            stop(msg)
        },
        .
    )} %>%
    filter(
        !str_detect(id_organism_a, sprintf('%s\\d+', organism_a)) &
        !str_detect(id_organism_b, sprintf('%s\\d+', organism_b)) &
        mapping %in% mappings
    ) %>%
    {`if`(only_ids, select(., id_organism_a, id_organism_b), .)}

}


#' Orthologous pairs between two organisms for ID types not supported by OMA
#'
#' The Orthologous Matrix (OMA), a resource of orthologous relationships
#' between genes, doesn't provide gene symbols, the identifier preferred in
#' many bioinformatics pipelines. Hence this function wraps
#' \code{\link{oma_pairwise}} by translating the identifiers used in OMA to
#' gene symbols. Items that can not be translated to `id_type` (but present
#' in the data with their internal OMA IDs) will be removed. Then,
#' in this function we translate the identifiers to the desired ID type.
#'
#' @param organism_a Name or identifier of an organism.
#' @param organism_b Name or identifier of another organism.
#' @param id_type The gene or protein identifier to use in the table. For a
#'     list of supported ID types see `omnipathr.env$id_types$oma`. These are
#'     the identifiers that will be translated to gene symbols.
#' @param oma_id_type Character: the gene or protein identifier to be queried
#'     from OMA. These IDs will be translated to `id_type`.
#' @param mappings Character vector: control ambiguous mappings: \itemize{
#'         \item{1:1 - unambiguous}
#'         \item{1:m - one-to-many}
#'         \item{n:1 - many-to-one}
#'         \item{n:m - many-to-many}
#'     }
#' @param only_ids Logical: include only the two identifier columns, not the
#'     mapping type and the orthology group columns.
#'
#' @return A data frame with orthologous gene pairs.
#'
#' @examples
#' oma_pairwise_translated("human", "mouse")
#'
#' @importFrom magrittr %>%
#' @importFrom rlang exec !!! !! := sym
#' @importFrom dplyr filter
#' @export
oma_pairwise_translated <- function(
    organism_a = 'human',
    organism_b = 'mouse',
    id_type = 'uniprot',
    oma_id_type = 'uniprot_entry',
    mappings = c('1:1', '1:m', 'n:1', 'n:m'),
    only_ids = TRUE
) {

    .slow_doctest()

    # NSE vs. R CMD check workaround:
    id_organism_a <- id_organism_b <- NULL

    environment() %>%
    as.list %>%
    extract(!startsWith(names(.), 'id_organism_')) %>%
    inset2('id_type', oma_id_type) %>%
    extract(names(.) != 'oma_id_type') %>%
    exec(oma_pairwise, !!!.) %>%
    translate_ids(
        id_organism_a := !!sym(oma_id_type),
        id_organism_a := !!sym(id_type),
        organism = organism_a
    ) %>%
    translate_ids(
        id_organism_b := !!sym(oma_id_type),
        id_organism_b := !!sym(id_type),
        organism = organism_b
    ) %>%
    filter(!is.na(id_organism_a) & !is.na(id_organism_b))

}


#' Orthologous pairs of gene symbols between two organisms
#'
#' The Orthologous Matrix (OMA), a resource of orthologous relationships
#' between genes, doesn't provide gene symbols, the identifier preferred in
#' many bioinformatics pipelines. Hence this function wraps
#' \code{\link{oma_pairwise}} by translating the identifiers used in OMA to
#' gene symbols. Items that can not be translated to `id_type` (but present
#' in the data with their internal OMA IDs) will be removed. Then,
#' in this function we translate the identifiers to gene symbols.
#'
#' @param organism_a Name or identifier of an organism.
#' @param organism_b Name or identifier of another organism.
#' @param oma_id_type Character: the gene or protein identifier to be queried
#'     from OMA. These IDs will be translated to `id_type`.
#' @param mappings Character vector: control ambiguous mappings: \itemize{
#'         \item{1:1 - unambiguous}
#'         \item{1:m - one-to-many}
#'         \item{n:1 - many-to-one}
#'         \item{n:m - many-to-many}
#'     }
#' @param only_ids Logical: include only the two identifier columns, not the
#'     mapping type and the orthology group columns.
#'
#' @return A data frame with orthologous gene pairs.
#'
#' @examples
#' oma_pairwise_genesymbols("human", "mouse")
#'
#' @importFrom magrittr %>%
#' @importFrom rlang exec !!!
#' @export
oma_pairwise_genesymbols <- function(
    organism_a = 'human',
    organism_b = 'mouse',
    oma_id_type = 'uniprot_entry',
    mappings = c('1:1', '1:m', 'n:1', 'n:m'),
    only_ids = TRUE
) {

    .slow_doctest()

    id_type <- 'genesymbol'

    environment() %>%
    as.list %>%
    exec(oma_pairwise_translated, !!!.)

}


#' OMA identifier type from synonyms and lower case version
#'
#' @param id_type Character: a synonym or a lower case version of an OMA
#'     identifier type.
#'
#' @return Character: an ID type symbol that can be used in OMA queries;
#'     NULL if id_type is not recognized.
#'
#' @importFrom magrittr %<>% %>% extract2
#' @importFrom rlang set_names
#' @importFrom stringr str_to_lower
#' @noRd
oma_id_type <- function(id_type) {

    id_type %<>% str_to_lower

    omnipathr.env$id_types$oma %>%
    extract2(id_type) %>%
    if_null(
        omnipathr.env$id_types$oma %>%
        set_names(., str_to_lower(.)) %>%
        extract2(id_type)
    ) %>%
    URLencode(reserved = TRUE)

}
