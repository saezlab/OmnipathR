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
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Preprocess the taxon name table from Ensembl
#'
#' @return A data frame with taxon names, each present also as all lower case,
#'     and the NCBI Taxonomy IDs converted to character.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate rename_with bind_cols across
#' @importFrom tidyselect vars_select_helpers everything
#' @importFrom stringr str_to_lower
#' @noRd
taxon_names_table <- function(){

    # NSE vs. R CMD check workaround
    ncbi_tax_id <- NULL

    ensembl_organisms() %>%
    {bind_cols(
        mutate(
            .,
            across(vars_select_helpers$where(is.character), str_to_lower),
            ncbi_tax_id = as.character(ncbi_tax_id)
        ) %>%
        rename_with(.fn = paste, .cols = everything(), 'l', sep = '_'),
        .
    )}

}

#' Translate between organism names
#'
#' Conversion between common English, latin (scientific), Ensembl and NCBI
#' taxon names and identifiers.
#'
#' @param name Character or integer: the name to be converted. Can be any
#'     of the types, if it's already the target type, it will be returned
#'     with no change. Case unsensitive.
#' @param name_type Character: type of the returned name or identifier.
#'     Many synonyms are accepted, the shortest ones: "latin", "ncbi",
#'     "common", "ensembl". Case unsensitive.
#'
#' @return Character or integer of length 1 or 0: the matched name of the
#'     requested name type.
#'
#' @examples
#' latin_name("human", "latin")
#' # [1] "Homo sapiens"
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom rlang enquo !! !!!
#' @importFrom dplyr filter first if_any recode pull
#' @importFrom tidyselect ends_with
#' @importFrom stringr str_to_lower
#' @noRd
taxon_name <- function(name, name_type){

    SYNONYMS <- list(
        latin = 'latin_name',
        scientific = 'latin_name',
        ncbi = 'ncbi_tax_id',
        taxid = 'ncbi_tax_id',
        common = 'common_name',
        english = 'common_name',
        ensembl = 'ensembl_id'
    )

    name_type <-
        .nse_ensure_str(!!enquo(name_type)) %>%
        str_to_lower %>%
        recode(!!!SYNONYMS)

    name %<>% as.character %>% str_to_lower

    get_db('organisms') %>%
    filter(
        if_any(ends_with('_l'), ~ .x == name)
    ) %>%
    pull(name_type) %>%
    first %>%
    if_null_len0(NA)

}


#' Latin (scientific) names of organisms
#'
#' @param name Vector with any kind of organism name or identifier, can be
#'     also mixed type.
#'
#' @return Character vector with latin (scientific) names, NA if a name
#'     in the input could not be found.
#'
#' @examples
#' latin_name(c(9606, "cat", "dog"))
#' # [1] "Homo sapiens" "Felis catus" "Canis lupus familiaris"
#' latin_name(c(9606, "cat", "doggy"))
#' # [1] "Homo sapiens" "Felis catus"  NA
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{ncbi_taxid}}}
#'     \item{\code{\link{common_name}}}
#'     \item{\code{\link{ensembl_name}}}
#' }
latin_name <- function(name){

    name %>%
    map_chr(taxon_name, 'latin')

}


#' NCBI Taxonomy IDs of organisms
#'
#' @param name Vector with any kind of organism name or identifier, can be
#'     also mixed type.
#'
#' @return Integer vector with NCBI Taxonomy IDs, NA if a name
#'     in the input could not be found.
#'
#' @examples
#' ncbi_taxid(c("Homo sapiens", "cat", "dog"))
#' # [1] 9606 9685 9615
#' ncbi_taxid(c(9606, "cat", "doggy"))
#' # [1] 9606 9685   NA
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_int
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{latin_name}}}
#'     \item{\code{\link{common_name}}}
#'     \item{\code{\link{ensembl_name}}}
#' }
ncbi_taxid <- function(name){

    name %>%
    map_int(taxon_name, 'ncbi') %>%
    as.integer

}


#' Ensembl identifiers of organisms
#'
#' @param name Vector with any kind of organism name or identifier, can be
#'     also mixed type.
#'
#' @return Character vector with Ensembl taxon names, NA if a name
#'     in the input could not be found.
#'
#' @examples
#' ensembl_name(c(9606, "cat", "dog"))
#' # [1] "hsapiens" "fcatus" "clfamiliaris"
#' ensembl_name(c("human", "kitten", "cow"))
#' # [1] "hsapiens" NA  "btaurus"
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{ncbi_taxid}}}
#'     \item{\code{\link{common_name}}}
#'     \item{\code{\link{latin_name}}}
#' }
ensembl_name <- function(name){

    name %>%
    map_chr(taxon_name, 'ensembl')

}


#' Common (English) names of organisms
#'
#' @param name Vector with any kind of organism name or identifier, can be
#'     also mixed type.
#'
#' @return Character vector with common (English) taxon names, NA if a name
#'     in the input could not be found.
#'
#' @examples
#' common_name(c(10090, "cjacchus", "Vicugna pacos"))
#' # [1] "Mouse" "White-tufted-ear marmoset" "Alpaca"
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{ncbi_taxid}}}
#'     \item{\code{\link{latin_name}}}
#'     \item{\code{\link{ensembl_name}}}
#' }
common_name <- function(name){

    name %>%
    map_chr(taxon_name, 'common')

}
