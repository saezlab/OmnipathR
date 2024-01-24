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


#' Preprocess the taxon name table from Ensembl
#'
#' @return A data frame with taxon names, each present also as all lower case,
#'     and the NCBI Taxonomy IDs converted to character.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate rename_with bind_cols across full_join select
#' @importFrom tidyselect vars_select_helpers everything
#' @importFrom stringr str_to_lower
#' @noRd
taxon_names_table <- function(){

    # NSE vs. R CMD check workaround
    ncbi_tax_id <- oma_version <- genome_source <-
    latin_name.x <- latin_name.y <- NULL

    ensembl_organisms() %>%
    full_join(
        oma_organisms() %>%
        select(-genome_source, -oma_version),
        by = 'ncbi_tax_id'
    ) %>%
    mutate(
        latin_name = ifelse(
            is.na(latin_name.y),
            latin_name.x,
            latin_name.y
        )
    ) %>%
    select(-latin_name.x, -latin_name.y) %>%
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
#' Conversion between common English, latin (scientific), Ensembl,
#' Orthologous Matrix (OMA) and NCBI taxon names and identifiers.
#'
#' @param name Character or integer: the name to be converted. Can be any
#'     of the types, if it's already the target type, it will be returned
#'     with no change. Case unsensitive.
#' @param name_type Character: type of the returned name or identifier.
#'     Many synonyms are accepted, the shortest ones: "latin", "ncbi",
#'     "common", "ensembl" and "oma". Case unsensitive.
#'
#' @return Character or integer of length 1 or 0: the matched name of the
#'     requested name type.
#'
#' @examples
#' taxon_name("human", "latin")
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
        ensembl = 'ensembl_id',
        oma = 'oma_code',
        oma_name = 'oma_code'
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


#' Orthologous Matrix (OMA) codes of organisms
#'
#' Note: OMA species codes are whenever possible identical to UniProt codes.
#'
#' @param name Vector with any kind of organism name or identifier, can be
#'     also mixed type.
#'
#' @return A character vector with the Orthologous Matrix (OMA) codes of the
#'     organisms.
#'
#' @examples
#' oma_code(c(10090, "cjacchus", "Vicugna pacos"))
#' # [1] "MOUSE" "CALJA" "VICPA"
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{ncbi_taxid}}}
#'     \item{\code{\link{latin_name}}}
#'     \item{\code{\link{ensembl_name}}}
#'     \item{\code{\link{common_name}}}
#' }
oma_code <- function(name) {

    name %>%
    map_chr(taxon_name, 'oma')

}


#' Reads a built-in table about organism support of resources
#'
#' @importFrom magrittr %>%
#' @importFrom yaml read_yaml
#'
#' @noRd
.load_organisms <- function(pkgname){

    omnipath.env$organisms <-
        system.file(
            'internal',
            'organisms.yaml',
            package = pkgname,
            mustWork = TRUE
        ) %>%
        read_yaml()

}

#' Make sure the resource supports the organism and it has the ID
#'
#' @param organism Character or integer: name or NCBI Taxonomy ID of the
#'     organism.
#' @param resource Charater: name of the resource.
#' @param error Logical: raise an error if the organism is not supported in the
#'     resource. Otherwise it only emits a warning.
#'
#' @return Character: the ID of the organism as it is used by the resource. NA
#'     if the organism can not be translated to the required identifier type.
#'
#' @examples
#' organism_for(10116, 'chalmers-gem')
#' # [1] "Rat"
#' organism_for(6239, 'chalmers-gem')
#' # [1] "Worm"
#' # organism_for('foobar', 'chalmers-gem')
#' # Error in organism_for("foobar", "chalmers-gem") :
#' # Organism `foobar` (common_name: `NA`; common_name: `NA`)
#' # is not supported by resource `chalmers-gem`. Supported organisms:
#' # Human, Mouse, Rat, Zebrafish, Drosophila melanogaster (Fruit fly),
#' # Caenorhabditis elegans (PRJNA13758).
#'
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom logger log_warn log_error
#' @export
organism_for <- function(organism, resource, error = TRUE) {

    resource_info <-
        omnipath.env$organisms %>%
        extract2(resource %>% str_to_lower)

    resource_info$supported %<>% maybe_call
    resource_info$check_by_id %<>% if_null(resource_info$id_type)
    organism_id_lookup <- get(resource_info$check_by_id)(organism)
    organism_id_resource <- get(resource_info$id_type)(organism)

    if (!organism_id_lookup %in% resource_info$supported) {

        supported <-
            resource_info$supported %>%
            {`if`(
                length(.) <= 9L,
                enum_format(resource_info$supported),
                sprintf(
                    paste0(
                        '%i taxons in total, see the built-in ',
                        '`organisms.yaml` for details'
                    ),
                    length(.)
                )
            )} %>%
            sprintf('Supported organisms: %s.', .)

        msg <- sprintf(
            paste0(
                'Organism `%s` (%s: `%s`; %s: `%s`) is not supported by ',
                'resource `%s`. %s'
            ),
            organism,
            resource_info$check_by_id,
            organism_id_lookup,
            resource_info$id_type,
            organism_id_resource,
            resource,
            supported
        )

        if (error) {

            log_error(msg)
            stop(msg)

        } else {

           log_warn(msg)
           warning(msg)

        }

    }

    organism_id_resource %>%
    if_null(extract2(resource_info$custom, .), .) %>%
    maybe_call(resource_info$finalize, .)

}
