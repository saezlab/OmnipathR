#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Diego Mananes
#                  Alberto Valdeolivas
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

#' List of GEMs in Metabolic Atlas
#'
#' Metabolic Atlas is a repository of genome scale models of metabolism (GEMs)
#' for various organisms, tissues and conditions. This function returns a list
#' of the available GEMs. Read more at \url{https://metabolicatlas.org}.
#'
#' @return A data frame (tibble) of GEMs.
#'
#' @examples
#' metabolic_atlas_models()
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom tidyr unnest_wider
#' @importFrom dplyr rename mutate select coalesce
#' @export
metabolic_atlas_list_models <- function() {

    # NSE vs. R CMD check workaround
    gemodelset <- sample <- description <- description2 <- NULL

    'metatlas_models' %>%
    download_to_cache() %>%
    safe_json() %>%
    tibble() %>%
    rename(description2 = description) %>%
    unnest_wider(c(gemodelset, sample)) %>%
    mutate(description = coalesce(description, description2)) %>%
    select(-description2)

}


#' Download and load one or more SBML models from Metabolic Atlas
#'
#' Downloads a genome-scale metabolic model (GEM) from Metabolic Atlas and
#' loads it using the SBMLR package (optional dependency).
#'
#' @param ... Arguments to filter models. Either:
#'     - A data frame of filtered models from metabolic_atlas_models()
#'     - Integer vector: model ID(s) from the `id` column
#'     - Expressions passed to dplyr::filter() (e.g., organism == 'Homo sapiens')
#' @param return_xml Logical: return an XML document object even if the SBMLR
#'     package is installed.
#'
#' @return A model object loaded by SBMLR, or an XML document object
#'     from xml2 if SBMLR is not available.
#'
#' @examples
#' # Load model by ID
#' model <- metabolic_atlas_model(1)
#'
#' # Load first human model
#' model <- metabolic_atlas_model(tissue = "Cervix")
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom rlang eval_tidy enexprs
#' @importFrom dplyr filter slice pull
#' @importFrom purrr map
#' @importFrom logger log_info log_error log_warn
#' @export
metabolic_atlas_models <- function(..., return_xml = FALSE) {

    # NSE vs. R CMD check workaround
    name <- files <- path <- NULL

    sbmlr_available <- requireNamespace("SBMLR", quietly = TRUE)

    if (!sbmlr_available) {
        msg <- "SBMLR package not available, install with: BiocManager::install('SBMLR')"
        log_warn(msg)
        warning(msg, call. = FALSE)
    }

    # Capture expressions first to prevent immediate evaluation
    filter_exprs <- enexprs(...)

    # Handle different argument types
    models <- if (length(filter_exprs) == 1L && is.data.frame(filter_exprs[[1L]])) {
        filter_exprs[[1L]]
    } else if (length(filter_exprs) == 1L && is.numeric(filter_exprs[[1L]])) {
        metabolic_atlas_list_models() %>% filter(id %in% filter_exprs[[1L]])
    } else {
        models <- metabolic_atlas_list_models()
        for (expr in filter_exprs) {
            models <- models %>% filter(eval_tidy(expr, data = models))
        }
        models
    }

    if (nrow(models) == 0L) {
        msg <- sprintf("No models found matching the criteria")
        log_error(msg)
        stop(msg, call. = FALSE)
    }

    models %>%
    {map(
        seq_len(nrow(.)),
        ~metabolic_atlas_model(slice(models, .x), return_xml = return_xml)
    )} %>%
    {`if`(length(.) == 1L, .[[1L]], .)}

}


#' Download and load an SBML model from Metabolic Atlas
#'
#' @importFrom xml2 read_xml
#' @noRd
metabolic_atlas_model <- function(model_row, return_xml = FALSE) {

    # NSE vs. R CMD check workaround
    id <- organism <- tissue <- cell_type <- condition <-
        reaction_count <- metabolite_count <- gene_count <- year <- NULL

    use_sbmlr <- requireNamespace('SBMLR', quietly = TRUE) & !return_xml

    sbml_path <-
        model_row$files[[1L]] %>%
        filter(format == 'SBML') %>%
        slice(1L) %>%
        pull(path)

    log_info(
        'Downloading model: %s %s',
        model_row$name,
        model_row %>%
            select(
                id, organism, tissue, cell_type, condition,
                reaction_count, metabolite_count, gene_count, year
            ) %>%
            as.list %>%
            compact_repr
    )

    method <- archive_extractor
    dl_args <- list(
        url_key = 'metatlas_model',
        url_param = list(sbml_path),
        reader = `if`(use_sbmlr, SBMLR::readSBML, xml2::read_xml),
        to_tempdir = use_sbmlr
    )

    exec(method, !!!dl_args)

}
