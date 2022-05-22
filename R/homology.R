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


#' Homology translation
#'
#' Translates identifiers between organisms using orthology data from Ensembl.
#'
#' @param d Data frame or character vector.
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
#' @param ensembl_orthology_types Character vector: use only this orthology
#'     relationship types. Possible values are "one2one", "one2many" and
#'     "many2many".
#' @param ensembl_min_orthology_confidence Integer: use only orthology
#'     relations with at least this level of confidence. In Ensembl the
#'     confidence can be either 0 or 1, so only these values make sense.
#'     If 0, all the orthology records will be used, if 1, only the ones
#'     with higher confidence.
#'
#' @return Data frame with the translated columns or character vector with
#'     translated identifiers.
#'
#' @examples
#' \dontrun{
#' # these proteins are ULK1, IFNG, EGFR, TGFB1, IL1R1
#' human_uniprots <- c("O75385", "P01579", "P00533", "P01137", "P14778")
#' homology_translate(human_uniprots)
#' }
#'
#' @importFrom magrittr %<>% %>% extract
#' @importFrom purrr map_chr
#' @importFrom rlang enquos !! enquo sym %||% :=
#' @importFrom dplyr rename_with inner_join rename filter select mutate pull
#' @importFrom tidyselect ends_with
#' @export
homology_translate <- function(
    d,
    ...,
    target = 10090,
    source = 9606,
    ensembl_orthology_types = c('one2one', 'one2many'),
    ensembl_min_orthology_confidence = 1L
){

    UNIPROT_DEFAULTS <- c('uniprot', 'uniprot', 'genesymbol')
    GENESYMBOL_DEFAULTS <- c('genesymbol')

    ensembl_orthology_types %<>% map_chr(~sprintf('ortholog_%s', .x))

    ids <-
        enquos(...) %>%
        map(.nse_ensure_str) %>%
        if_null_len0(
            `if`(
                'uniprot' %in% colnames(d) || is_uniprot(d),
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
    target_ens_org <- ensembl_name(target)
    target_up_col <- sym(sprintf('%s_uniprot', target_ens_org))
    target_col <- sym(sprintf('%s_%s', target_ens_org, target_id_type))

    target2_col <- id_cols[3]
    target2_id_type <- sym(id_types[3])
    target2_col_tmp <- sym(sprintf('%s__tmp', target2_col))


    source_col_str <- id_cols[1]
    source_col <- sym(source_col_str)
    source_col_pos <- d %>% colnames %>% {which(. == source_col_str)}

    # to handle single vectors as inputs, we convert them to data frames
    vector_input <- !is.data.frame(d)

    if(vector_input){

        d <- tibble(!!source_col := d)

    }

    # have to keep these to move the new column
    # to the position of the original column
    before_col <-
        d %>% colnames %>%
        extract(source_col_pos - 1L) %>%
        {`if`(length(.) == 0L, NULL, .)}
    after_col <-
        d %>% colnames %>%
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
                    is_uniprot(d$.),
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

    # homology table from Ensembl, filtered:
    ensembl_orthology(
        organism_a = source,
        organism_b = target
    ) %>%
    filter(
        b_orthology_type %in% ensembl_orthology_types &
        b_orthology_confidence >= ensembl_min_orthology_confidence
    ) %>%
    # translating source Ensembl peptide IDs to UniProts
    translate_ids(
        ensembl_peptide_id = ensp,
        uniprot,
        ensembl = TRUE
    ) %>%
    # removing the untranslated ones and
    # the ones without orthologous peptide
    filter(
        !is.na(b_ensembl_peptide) &
        !is.na(uniprot)
    ) %>%
    # adding an ugly label so we avoid any duplicate column names
    rename_with(~sprintf('%s__homology', .x)) %>%
    # joining the input data frame by source side uniprots
    inner_join(
        d,
        by = c(uniprot__homology = source_col_str)
    ) %>%
    # we want to keep the name from the original data frame
    # but that was on the left side, so we rename
    rename(!!source_col := uniprot__homology) %>%
    # translating the target Ensembl peptide IDs to UniProts
    translate_ids(
        b_ensembl_peptide__homology := ensp,
        !!target_up_col := uniprot,
        organism = target_ncbi,
        ensembl = TRUE
    ) %>%
    # if the final ID type is not UniProt, translating the target
    # UniProts to the desired ID type
    {`if`(
        target_id_type != 'uniprot',
        translate_ids(
            .,
            !!target_up_col := uniprot,
            !!target_col := !!sym(target_id_type),
            organism = target_ncbi
        ),
        .
    )} %>%
    # now we have the translation done
    # let's get rid of the redundant homology columns
    select(-ends_with('__homology')) %>%
    select(., -!!source_col) %>%
    # if a secondary identifier type is requested for the target
    # organism, let's translate it from UniProt:
    {`if`(
        is.na(target2_col),
        .,
        translate_ids(
            .,
            !!target_up_col := uniprot,
            !!target2_col_tmp := !!target2_id_type,
            organism = target_ncbi
        ) %>%
        mutate(
            !!sym(target2_col) := !!target2_col_tmp
        ) %>%
        select(-!!target2_col_tmp)
    )} %>%
    # if the target ID type is not UniProt,
    # we can remove the target UniProt column
    {`if`(
        target_id_type != 'uniprot',
        select(., -!!target_up_col),
        .
    )} %>%
    # remove the untranslated records
    filter(!is.na(!!target_col)) %>%
    # replace the original column with the translated one
    rename(!!source_col := !!target_col) %>%
    {`if`(
        vector_input,
        # output is a vector
        pull(., !!source_col) %>% unique,
        # move it to its original position
        {`if`(
            is.null(after_col),
            {`if`(
                is.null(before_col),
                # this happens to single column data frames only
                .,
                relocate(., !!source_col, .after = before_col)
            )},
            relocate(., !!source_col, .before = after_col)
        )}
    )}


}
