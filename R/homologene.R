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


#' Orthology data from NCBI HomoloGene
#'
#' Retrieves NCBI HomoloGene data without any processing. Processed tables
#' are more useful for most purposes, see below other functions that provide
#' those. Genes of various organisms are grouped into homology groups
#' ("hgroup" column). Organisms are identified by NCBI Taxonomy IDs, genes
#' are identified by four different identifier types.
#'
#' @return A data frame as provided by NCBI HomoloGene.
#'
#' @examples
#' hg <- homologene_raw()
#' hg
#' # # A tibble: 275,237 × 6
#' #    hgroup ncbi_taxid entrez  genesymbol  gi        refseqp
#' #     <int>      <int> <chr>   <chr>       <chr>     <chr>
#' #  1      3       9606 34      ACADM       4557231   NP_000007.1
#' #  2      3       9598 469356  ACADM       160961497 NP_001104286.1
#' #  3      3       9544 705168  ACADM       109008502 XP_001101274.1
#' #  4      3       9615 490207  ACADM       545503811 XP_005622188.1
#' #  5      3       9913 505968  ACADM       115497690 NP_001068703.1
#' # # . with 275,232 more rows
#'
#' # which organisms are available?
#' common_name(unique(hg$ncbi_taxid))
#' #  [1] "Human" "Chimpanzee" "Macaque" "Dog" "Cow" "Mouse" "Rat" "Zebrafish"
#' #  [9] "D. melanogaster" "Caenorhabditis elegans (PRJNA13758)"
#' # [11] "Tropical clawed frog" "Chicken"
#' # ...and 9 more organisms with missing English names.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr cols col_integer col_character
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{homologene_download}}}
#' }
homologene_raw <- function(){

    .slow_doctest(value = list(ncbi_taxid = 'OmnipathR: no data'))

    hdr <- c('hgroup', 'ncbi_taxid', 'entrez', 'genesymbol', 'gi', 'refseqp')

    'homologene' %>%
    generic_downloader(
        reader_param = list(
            col_names = hdr,
            col_types = cols(
                hgroup = col_integer(),
                ncbi_taxid = col_integer(),
                entrez = col_character(),
                gi = col_character()
            )
        ),
        resource = 'NCBI HomoloGene'
    ) %T>%
    load_success()

}


#' Orthology table for a pair of organisms
#'
#' Orthologous pairs of genes for a pair of organisms from NCBI HomoloGene,
#' using one identifier type.
#'
#' @param target Character or integer: name or ID of the target organism.
#' @param source Character or integer: name or ID of the source organism.
#' @param id_type Symbol or character: identifier type, possible values are
#'     "genesymbol", "entrez", "refseqp" or "gi".
#' @param hgroup_size Logical: include a column with the size of the homology
#'     groups. This column distinguishes one-to-one and one-to-many or
#'     many-to-many mappings.
#'
#' @details
#' The operation of this function is symmetric, *source* and *target* are
#' interchangeable but determine the column layout of the output. The column
#' "hgroup" is a numberic identifier of the homology groups. Most of the
#' groups consist of one pair of orthologous genes (one-to-one mapping), and
#' a few of them multiple ones (one-to-many or many-to-many mappings).
#'
#' @return A data frame with orthologous identifiers between the two organisms.
#'
#' @examples
#' chimp_human <- homologene_download(chimpanzee, human, refseqp)
#' chimp_human
#' # # A tibble: 17,737 × 3
#' #    hgroup refseqp_source refseqp_target
#' #     <int> <chr>          <chr>
#' #  1      3 NP_000007.1    NP_001104286.1
#' #  2      5 NP_000009.1    XP_003315394.1
#' #  3      6 NP_000010.1    XP_508738.2
#' #  4      7 NP_001096.1    XP_001145316.1
#' #  5      9 NP_000014.1    XP_523792.2
#' # # . with 17,732 more rows
#'
#' @importFrom rlang !! enquo sym
#' @importFrom magrittr %>%
#' @importFrom dplyr inner_join filter select group_by mutate n ungroup
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{homologene_raw}}}
#'     \item{\code{\link{homologene_uniprot_orthology}}}
#' }
homologene_download <- function(
    target = 10090L,
    source = 9606L,
    id_type = 'genesymbol',
    hgroup_size = FALSE
){

    .slow_doctest()

    # NSE vs. R CMD check workaround
    hgroup <- NULL

    source <- .nse_ensure_str(!!enquo(source)) %>% ncbi_taxid
    target <- .nse_ensure_str(!!enquo(target)) %>% ncbi_taxid
    id_type <- .nse_ensure_str(!!enquo(id_type))

    homologene_raw() %>%
    {inner_join(
        filter(., ncbi_taxid == source) %>% select(hgroup, !!sym(id_type)),
        filter(., ncbi_taxid == target) %>% select(hgroup, !!sym(id_type)),
        by = 'hgroup',
        suffix = c('_source', '_target')
    )} %>%
    {`if`(
        hgroup_size,
        group_by(., hgroup) %>% mutate(hgroup_size = n()) %>% ungroup,
        .
    )}

}


#' Orthology table with UniProt IDs
#'
#' Orthologous pairs of UniProt IDs for a pair of organisms, based on NCBI
#' HomoloGene data.
#'
#' @param target Character or integer: name or ID of the target organism.
#' @param source Character or integer: name or ID of the source organism.
#' @param by Symbol or character: the identifier type in NCBI HomoloGene
#'     to use. Possible values are "refseqp", "entrez", "genesymbol", "gi".
#' @param ... Further arguments passed to \code{\link{translate_ids}}.
#'
#' @return A data frame with orthologous pairs of UniProt IDs.
#'
#' @examples
#' homologene_uniprot_orthology(by = genesymbol)
#' # # A tibble: 14,235 × 2
#' #    source target
#' #    <chr>  <chr>
#' #  1 P11310 P45952
#' #  2 P49748 P50544
#' #  3 P24752 Q8QZT1
#' #  4 Q04771 P37172
#' #  5 Q16586 P82350
#' # # . with 14,230 more rows
#'
#' @importFrom rlang !! enquo sym :=
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter distinct
#' @export
homologene_uniprot_orthology <- function(
    target = 10090L,
    source = 9606L,
    by = entrez,
    ...
){

    .slow_doctest()

    # NSE vs. R CMD check workaround
    entrez <- uniprot <- NULL

    by <- .nse_ensure_str(!!enquo(by))
    source <- .nse_ensure_str(!!enquo(source)) %>% ncbi_taxid
    target <- .nse_ensure_str(!!enquo(target)) %>% ncbi_taxid

    homologene_download(
        target = !!target,
        source = !!source,
        id_type = !!sym(by)
    ) %>%
    translate_ids(
        !!sym(sprintf('%s_source', by)) := !!sym(by),
        source = uniprot,
        organism = !!source,
        ...
    ) %>%
    translate_ids(
        !!sym(sprintf('%s_target', by)) := !!sym(by),
        target = uniprot,
        organism = !!target,
        ...
    ) %>%
    select(source, target) %>%
    filter(!is.na(source) & !is.na(target)) %>%
    distinct


}
