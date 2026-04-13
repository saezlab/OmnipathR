#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Alberto Valdeolivas
#                  Denes Turei (turei.denes@gmail.com)
#                  Attila Gabor
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Download one MACdb table
#'
#' @param file_type Character: one of `"metabolite"` or `"study"`.
#'
#' @return A tibble containing one MACdb table.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom rlang arg_match
#' @importFrom dplyr mutate across
#' @importFrom tidyselect where
#'
#' @noRd
macdb_download_table <- function(file_type = c('metabolite', 'study')) {

    file_type <- arg_match(file_type)

    generic_downloader(
        url_key = 'macdb_download_file',
        url_param = list(file_type),
        resource = sprintf('MACdb (%s)', file_type)
    ) %>%
    mutate(across(where(is.character), .macdb_normalize_missing)) %T>%
    load_success()

}


#' Normalize MACdb missing value markers
#'
#' @param x A vector.
#'
#' @return Character vector with common missing markers converted to `NA`.
#'
#' @noRd
.macdb_normalize_missing <- function(x) {

    if (!is.character(x)) {

        return(x)

    }

    x <- trimws(x)
    x[tolower(x) %in% c('', 'na', 'nan', 'n/a')] <- NA_character_

    x

}


#' Build a metabolite-centered MACdb metabolite-cancer association table
#'
#' Downloads MACdb metabolite and study tables and joins them by `Cohort_id`.
#' The output is ordered as metabolite identifiers, study information,
#' association statistics, and data-availability flags.
#'
#' @return A tibble with one row per metabolite association record.
#'
#' @examples
#' assoc <- macdb_metabolite_cancer_associations()
#' assoc
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of mutate
#' @importFrom logger log_error
#' @export
macdb_metabolite_cancer_associations <- function() {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    `case_control_p-value` <- log2FC <- case_concentration <-
        control_concentration <- NULL

    met <- macdb_download_table('metabolite')
    study <- macdb_download_table('study')

    met_cols <- names(met)
    duplicate_delta <- grepl('^Delta_concentration', met_cols)

    if (any(duplicate_delta)) {

        met <- met[, !duplicate_delta, drop = FALSE]

    }

    required_study <- c(
        'Cohort_id',
        'Trait_onto_ID',
        'Cancer_type',
        'Cancer_subtype',
        'Cancer_DOID',
        'Tissue',
        'pubmed_id'
    )

    required_met <- c(
        'Cohort_id',
        'original_metabolite_name',
        'pubchem_CID',
        'case_control_p-value',
        'log2FC'
    )

    missing_study <- setdiff(required_study, names(study))
    missing_met <- setdiff(required_met, names(met))

    if (length(missing_study) > 0L) {

        msg <- sprintf(
            'MACdb study table missing required columns: %s.',
            paste(missing_study, collapse = ', ')
        )
        log_error(msg)
        stop(msg)

    }

    if (length(missing_met) > 0L) {

        msg <- sprintf(
            'MACdb metabolite table missing required columns: %s.',
            paste(missing_met, collapse = ', ')
        )
        log_error(msg)
        stop(msg)

    }

    met_keep <- c(
        'Cohort_id',
        'original_metabolite_name',
        'pubchem_CID',
        'case_control_p-value',
        'log2FC',
        'case_concentration',
        'control_concentration'
    )

    study_keep <- c(
        'Cohort_id',
        'Trait_onto_ID',
        'Cancer_type',
        'Cancer_subtype',
        'Cancer_DOID',
        'Tissue',
        'pubmed_id'
    )

    met <- met %>%
    select(all_of(intersect(met_keep, names(met)))) %>%
    mutate(
        `case_control_p-value` = as.numeric(`case_control_p-value`),
        log2FC = as.numeric(log2FC),
        case_concentration = as.numeric(case_concentration),
        control_concentration = as.numeric(control_concentration)
    )

    study <- study %>%
    select(all_of(study_keep))

    duplicate_cohort <- unique(study$Cohort_id[duplicated(study$Cohort_id)])

    if (length(duplicate_cohort) > 0L) {

        msg <- sprintf(
            'MACdb study table has duplicated Cohort_id values (first 10): %s.',
            paste(head(duplicate_cohort, 10L), collapse = ', ')
        )
        log_error(msg)
        stop(msg)

    }

    study_idx <- match(met$Cohort_id, study$Cohort_id)
    study_side <- study[study_idx, setdiff(names(study), 'Cohort_id'), drop = FALSE]

    joined <- cbind(met, study_side)

    names(joined)[names(joined) == 'original_metabolite_name'] <- 'metabolite_name'
    names(joined)[names(joined) == 'pubchem_CID'] <- 'metabolite_pubchem_cid'
    names(joined)[names(joined) == 'case_control_p-value'] <- 'association_pvalue'
    names(joined)[names(joined) == 'log2FC'] <- 'association_log2fc'
    names(joined)[names(joined) == 'case_concentration'] <- 'association_case_concentration'
    names(joined)[names(joined) == 'control_concentration'] <- 'association_control_concentration'
    names(joined)[names(joined) == 'Trait_onto_ID'] <- 'study_trait_onto_id'
    names(joined)[names(joined) == 'Cancer_type'] <- 'study_cancer_type'
    names(joined)[names(joined) == 'Cancer_subtype'] <- 'study_cancer_subtype'
    names(joined)[names(joined) == 'Cancer_DOID'] <- 'study_cancer_doid'
    names(joined)[names(joined) == 'Tissue'] <- 'study_tissue'
    names(joined)[names(joined) == 'pubmed_id'] <- 'study_pubmed_id'

    nrows <- nrow(joined)

    joined$has_pValue <- if (
        'association_pvalue' %in% names(joined)
    ) {
        !is.na(joined$association_pvalue)
    } else {
        rep(FALSE, nrows)
    }

    joined$has_log2FC <- if (
        'association_log2fc' %in% names(joined)
    ) {
        !is.na(joined$association_log2fc)
    } else {
        rep(FALSE, nrows)
    }

    joined$has_caseConcentration <- if (
        'association_case_concentration' %in% names(joined)
    ) {
        !is.na(joined$association_case_concentration)
    } else {
        rep(FALSE, nrows)
    }

    joined$has_controlConcentration <- if (
        'association_control_concentration' %in% names(joined)
    ) {
        !is.na(joined$association_control_concentration)
    } else {
        rep(FALSE, nrows)
    }

    joined$has_bothConcentrations <- (
        joined$has_caseConcentration & joined$has_controlConcentration
    )

    ordered_cols <- c(
        'metabolite_pubchem_cid',
        'metabolite_name',
        'Cohort_id',
        'study_trait_onto_id',
        'study_cancer_type',
        'study_cancer_subtype',
        'study_cancer_doid',
        'study_tissue',
        'study_pubmed_id',
        'association_pvalue',
        'association_log2fc',
        'association_case_concentration',
        'association_control_concentration',
        'has_pValue',
        'has_log2FC',
        'has_caseConcentration',
        'has_controlConcentration',
        'has_bothConcentrations'
    )

    joined <- joined[, intersect(ordered_cols, names(joined)), drop = FALSE]

    joined

}
