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



#' Translate identifiers between organisms by orthologous gene pairs
#'
#' Translates identifiers between organisms using orthology data from NCBI
#' HomoloGene.
#'
#' @param data Data frame or character vector.
#' @param ... Column specification: from zero to up to three arguments, with
#'     or without names. NSE is supported. Arguments beyond the third one
#'     will be ignored.
#'     \itemize{
#'         \item{
#'             The name of the arguments should be column names, the
#'             values identifier types, either as character or as symbols.
#'         }
#'         \item{
#'             Arguments without names assumed to be both column names and
#'             identifier types, e.g. a column called "uniprot" containing
#'             UniProt IDs.
#'         }
#'         \item{
#'             The first column spefication describes the source column,
#'             with identifiers of the source organism. This column must
#'             exist in the data and this will be the input of the homology
#'             translation. This column will be removed from the returned
#'             data frame.
#'         }
#'         \item{
#'             In case of "uniprot", the source column name can be anything,
#'             if it contains only UniProt IDs it will be handled accordingly.
#'         }
#'         \item{
#'             In case of "genesymbol", is enough if the source column name
#'             contains the word "genesymbol", e.g. "ligand_genesymbol".
#'         }
#'         \item{
#'             The second column spefication describes the target column,
#'             with its name and identifier type. If not provided, both
#'             the column name and type will be the same as the source
#'         }
#'         \item{
#'             Optionally a third column can be specified with another
#'             identifier type. This is convenient if you want, for example
#'             also Gene Symbols along with UniProt IDs.
#'         }
#'         \item{
#'             If no specification provided, the input assumed to have
#'             a column named either "uniprot" or "genesymbol", or be
#'             a character vector of UniProt IDs or Gene Symbols.
#'         }
#'     }
#' @param target Character or integer: name or identifier of the target
#'     organism (the one we translate to). The default target organism is
#'     mouse.
#' @param source Character or integer: name of identifier of the source
#'     organism (the one the IDs in the input data belong to). The default
#'     source organism is human.
#'
#' @noRd
orthology_translate <- function(
    data,
    ...,
    target = 10090,
    source = 9606
) {

    # NSE vs R CMD check workaround
    d <- NULL

    UNIPROT_DEFAULTS <- c('uniprot', 'uniprot', 'genesymbol')
    GENESYMBOL_DEFAULTS <- c('genesymbol')
    HOMOLOGENE_ID_TYPES <- c('uniprot', 'genesymbol', 'entrez', 'refseq', 'gi')

    ids <-
        enquos(...) %>%
        map(.nse_ensure_str) %>%
        if_null_len0(
            `if`(
                'uniprot' %in% colnames(data) || is_uniprot(data),
                UNIPROT_DEFAULTS,
                GENESYMBOL_DEFAULTS
            )
        ) %>%
        set_names(names(.) %||% unlist(.)) %>%
        set_names(ifelse(nchar(names(.)), names(.), unlist(.)))

    id_cols <- names(ids)
    id_types <- unlist(ids)

    if(length(id_cols) == 1L || id_cols[2] == '.replace'){

        id_cols[2] <- id_cols[1]
        id_types[2] <- id_types[1]

    }

    target_col <- id_cols[2]
    target_id_type <- id_types[2]

    target2_col <- id_cols[3]
    target2_id_type <- sym(id_types[3])
    target2_col_tmp <- sym(sprintf('%s__tmp', target2_col))

    source_col <- id_cols[1]
    source_col_pos <- data %>% colnames %>% {which(. == source_col_str)}

    # to handle single vectors as inputs, we convert them to data frames
    vector_input <- !is.data.frame(data)

    if(vector_input){

        data <- tibble(!!sym(source_col) := data)

    }

    # have to keep these to move the new column
    # to the position of the original column
    before_col <-
        data %>% colnames %>%
        extract(source_col_pos - 1L) %>%
        {`if`(length(.) == 0L, NULL, .)}
    after_col <-
        data %>% colnames %>%
        extract(source_col_pos + 1L) %>%
        if_null_len0(NA) %>% # this happens if the input is vector
        {`if`(is.na(.) || !is.null(before_col), NULL, .)}

    # trying to be smart: if the name is not a registered ID type, but
    # "genesymbol" is in the name then we assume it's genesymbol, otherwise
    # if the column contains UniProt IDs, we know it's uniprot, and finally
    # we fall back to the provided value, it still might be a non registered
    # ID type
    source_id_type <-
        id_types[1] %>%
        {`if`(
            is_id_type(.),
            .,
            `if`(
                str_detect(., 'genesymbol'),
                'genesymbol',
                `if`(
                    is_uniprot(data$.),
                    'uniprot',
                    .
                )
            )
        )}

    source_ncbi <- ncbi_taxid(source)
    target_ncbi <- ncbi_taxid(target)

    # if we start from anything else than UniProt,
    # first we have to translate it to UniProt:
    if(source_id_type != 'uniprot'){

        source_uniprot <- sprintf('%s_uniprot', source_col_str)

        d %<>%
            translate_ids(
                !!source_col := source_id_type,
                !!sym(source_uniprot) := 'uniprot',
                organism = source_ncbi
            ) %>%
            select(-!!source_col)

        source_col_str <- source_uniprot
        source_col <- sym(source_col_str)

    }

}


#' Translate a column of identifiers by orthologous gene pairs
#'
#' @param data A data frame with the column to be translated.
#' @param column Name of a character column with identifiers of the source
#'     organism of type `id_type`.
#' @param id_type Type of identifiers in `column`. Available ID types include
#'     "uniprot", "entrez", "ensg", "refseq" and "swissprot" for OMA, and
#'     "uniprot", "entrez", "genesymbol", "refseq" and "gi" for NCBI
#'     HomoloGene. If you want to translate an ID type not directly available
#'     in your preferred resource, use first \code{\link{translate_ids}}
#'     to translate to an ID type directly available in the orthology resource.
#'     If not provided, it is assumed the column name is the ID type.
#' @param target_organism Name or NCBI Taxonomy ID of the target organism.
#' @param source_organism Name or NCBI Taxonomy ID of the source organism.
#' @param resource Character: source of the orthology mapping. Currently
#'     Orthologous Matrix (OMA) and NCBI HomoloGene are available, refer to
#'     them by "oma" and "homologene", respectively.
#' @param replace Logical or character: replace the column with the translated
#'     identifiers, or create a new column. If it is character, it will be
#'     used as the name of the new column.
#' @param one_to_many Integer: maximum number of orthologous pairs for one
#'     gene of the source organism. Genes mapping to higher number of
#'     orthologues will be dropped.
#' @param keep_untranslated Logical: keep records without orthologous pairs.
#'     If `replace` is TRUE, this option is ignored, and untranslated records
#'     will be dropped. Genes with more than `one_to_many` orthologues will
#'     always be dropped.
#' @param translate_complexes Logical: translate the complexes by translating
#'     their components.
#' @param uniprot_by_id_type Character: translate NCBI HomoloGene to UniProt
#'     by this ID type. One of "genesymbol", "entrez", "refseq" or "gi".
#'
#' @return The data frame with identifiers translated to other organism.
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom dplyr inner_join left_join relocate mutate pull
#' @importFrom dplyr filter rename group_by select ungroup n join_by
#' @importFrom rlang sym !! := enquo
#' @importFrom tidyselect any_of
#' @export
orthology_translate_column <- function(
        data,
        column,
        id_type = NULL,
        target_organism = 'mouse',
        source_organism = 'human',
        resource = 'oma',
        replace = FALSE,
        one_to_many = NULL,
        keep_untranslated = FALSE,
        translate_complexes = FALSE,
        uniprot_by_id_type = 'entrez'
) {

    column <- .nse_ensure_str(!!enquo(column))
    id_type <- .nse_ensure_str(!!enquo(id_type)) %>% if_null(column)
    uniprot <- id_type == 'uniprot'
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

    orthology_param <-
        list(
            source_organism,
            target_organism
        ) %>%
        set_names(`if`(
            resource == 'oma',
            c('organism_a', 'organism_b'),
            c('source', 'target')
        )) %>%
        {`if`(
            resource == 'homologene',
            c(
                .,
                `if`(
                    uniprot,
                    list(by = uniprot_by_id_type),
                    list(id_type = id_type)
                )
            ),
           .
        )}

    db_name <- `if`(
        resource == 'oma',
        'oma',
        `if`(uniprot, 'homologene_uniprot', 'homologene')
    )

    orthologous_pairs <-
        get_db(db_name, param = orthology_param) %>%
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
    join(
         orthologous_pairs,
         by = join_by(!!sym(column) == !!sym(ORTHO_SOURCE_COL))
    ) %>%
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
#' @importFrom magrittr %>% extract
#' @importFrom tibble tibble
#' @importFrom tidyr separate_rows unnest_longer
#' @importFrom dplyr filter mutate left_join n pull group_by distinct rowwise
#' @importFrom dplyr ungroup mutate_all distinct c_across summarise_at
#' @importFrom tidyselect everything
#' @importFrom rlang set_names !! := sym
#' @importFrom stringr str_detect str_replace
#' @noRd
complex_orthology <- function(
    orthology,
    identifiers,
    one_to_many = 1L
) {

    # NSE vs. R CMD check workaround
    components <- complexes <- vars <- NULL

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
