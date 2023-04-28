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

ORTHO_SOURCE_COL <- 'omnipathr_orthology_source'
ORTHO_TARGET_COL <- 'omnipathr_orthology_target'
ORTHO_GROUP_COL <- 'omnipatr_orthology_group'
ORTHO_COMP_OTM_COL <- 'omnipathr_complex_onetomany'
CPLEX_PREFIX <- 'COMPLEX:'


#' Translate a column of identifiers by orthologous gene pairs
#'
#' @param data A data frame with the column to be translated.
#' @param column Name of a character column with identifiers of the source
#'     organism of type `id_type`.
#' @param target_organism Name or NCBI Taxonomy ID of the target organism.
#' @param source_organism Name or NCBI Taxonomy ID of the source organism.
#' @param replace Logical or character: replace the column with the translated
#'     identifiers, or create a new column. If it is character, the
#' @param one_to_many Integer: maximum number of orthologous pairs for one
#'     gene of the source organism.
#' @param keep_untranslated Logical: keep records without orthologous pairs.
#'     If `replace` is TRUE, this option is ignoredx, and untranslated records
#'     will be dropped.
#' @param translate_complexes Logical: translate the complexes by translating
#'     their components.
#' @param uniprot_by_id_type Character: translate NCBI HomoloGene to UniProt
#'     by this ID type. One of "genesymbol", "entrez", "refseq" or "gi".
#'
#' @return The data frame with identifiers translated to other organism.
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom dplyr inner_join left_join relocate mutate pull
#' @importFrom dplyr filter rename group_by select ungroup
#' @importFrom rlang sym !! := enquo
#' @importFrom tidyselect any_of
#' @export
orthology_translate_column <- function(
        data,
        column,
        id_type,
        target_organism,
        source_organism = 9606L,
        replace = FALSE,
        one_to_many = NULL,
        keep_untranslated = FALSE,
        translate_complexes = FALSE,
        uniprot_by_id_type = 'entrez'
) {

    uniprot <- id_type == 'uniprot'
    column <- .nse_ensure_str(!!enquo(column))
    id_type <- .nse_ensure_str(!!enquo(id_type))
    unirprot_by_id_type <- .nse_ensure_str(!!enquo(uniprot_by_id_type))
    target_organism %<>% ncbi_taxid
    source_organism %<>% ncbi_taxid

    target_column <-
        column %>%
        {`if`(
            is.character(replace),
            replace,
            `if`(
                replace,
                column,
                sprintf('%s_%i', column, target_organism)
            )
        )}

    homologene_param <-
        list(
            source = source_organism,
            target = target_organism
        ) %>%
        c(
            `if`(
                uniprot,
                list(by = uniprot_by_id_type),
                list(id_type = id_type)
            )
        )

    db_name <- `if`(uniprot, 'homologene_uniprot', 'homologene')

    hg <-
        get_db(db_name, param = homologene_param) %>%
        select(-any_of('hgroup')) %>%
        set_names(c(ORTHO_SOURCE_COL, ORTHO_TARGET_COL)) %>%
        {`if`(
            translate_complexes,
            bind_rows(
                .,
                complex_orthology(
                    orthology = .,
                    identifiers = data %>% pull(!!sym(column)),
                    one_to_many = translate_complexes
                )
            ),
            .
        )}

    join <- `if`(
        !keep_untranslated || column == target_column,
        inner_join,
        left_join
    )

    data %>%
    mutate(!!sym(ORTHO_GROUP_COL) := 1L:n()) %>%
    join(hg, by = join_by(!!sym(column) == !!sym(ORTHO_SOURCE_COL))) %>%
    relocate(ORTHO_TARGET_COL, .after = column) %>%
    {`if`(
        column == target_column,
        select(., -!!sym(column)),
        .
    )} %>%
    rename(!!sym(target_column) := !!sym(ORTHO_TARGET_COL)) %>%
    {`if`(
        is.null(one_to_many),
        .,
        group_by(., !!sym(ORTHO_GROUP_COL)) %>%
        filter(n() <= one_to_many) %>%
        ungroup
    )} %>%
    select(-!!sym(ORTHO_GROUP_COL))

}


#' Translate complexes by their members
#'
#' @param orthology A data frame with the orthologous gene pairs.
#' @param identifiers A vector of identifiers, all or some of them are
#'     complexes.
#' @param one_to_many Integer: maximum number of orthologous complexes
#'     to derive from one complex in the source organism. To avoid
#'     unpleasant surprises, by default one-to-many mapping is disabled.
#'     To impose no limits on one-to-many mapping, set this argument to
#'     `FALSE`.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom tidyr separate_rows unnest_longer
#' @importFrom dplyr filter mutate left_join n pull group_by distinct
#' @importFrom dplyr ungroup mutate_all distinct c_across summarise_at
#' @importFrom tidyselect everything
#' @importFrom rlang set_names !! sym
#' @importFrom stringr str_detect str_replace
#' @noRd
complex_orthology <- function(
    orthology,
    identifiers,
    one_to_many = 1L
) {

    has_prefix <- identifiers %>% str_detect(CPLEX_PREFIX) %>% any

    comp_prod <- function(comp) {

        comp %>%
        expand.grid %>%
        mutate_all(as.character) %>%
        rowwise %>%
        mutate(
            complexes = paste(sort(c_across(everything())), collapse = '_')
        ) %>%
        pull(complexes)

    }

    cplex_table <<-
    identifiers %>%
    tibble(source = .) %>%
    filter(str_detect(source, '_')) %>%
    mutate(components = str_replace(source, CPLEX_PREFIX, '')) %>%
    separate_rows(components, sep = '_') %>%
    left_join(orthology, by = c('components' = ORTHO_SOURCE_COL)) %>%
    group_by(source) %>%
    filter(!any(is.na(!!sym(ORTHO_TARGET_COL)))) %>%
    ungroup %>%
    {`if`(
        one_to_many,
        group_by(., source, components) %>%
        mutate(!!sym(ORTHO_COMP_OTM_COL) := n()) %>%
        ungroup %>%
        group_by(source) %>%
        filter(prod(!!sym(ORTHO_COMP_OTM_COL)) <= one_to_many) %>%
        ungroup %>%
        select(-!!sym(ORTHO_COMP_OTM_COL)),
        .
    )} %>%
    group_by(source, components) %>%
    mutate(!!sym(ORTHO_TARGET_COL) := list(!!sym(ORTHO_TARGET_COL))) %>%
    summarise_at(.vars = vars(!!sym(ORTHO_TARGET_COL)), .funs = extract, 1L) %>%
    ungroup %>%
    group_by(source) %>%
    mutate(!!sym(ORTHO_TARGET_COL) := list(comp_prod(!!sym(ORTHO_TARGET_COL)))) %>%
    unnest_longer(!!sym(ORTHO_TARGET_COL)) %>%
    {`if`(
        has_prefix,
        mutate(
            .,
            !!sym(ORTHO_TARGET_COL) :=
            sprintf('%s%s', CPLEX_PREFIX, !!sym(ORTHO_TARGET_COL))
        ),
        .
    )} %>%
    select(-components) %>%
    set_names(c(ORTHO_SOURCE_COL, ORTHO_TARGET_COL)) %>%
    distinct

}
