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
metabolic_atlas_models <- function() {

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


#' Download and load a Metabolic Atlas model
#'
#' Downloads a genome-scale metabolic model (GEM) from Metabolic Atlas and
#' loads it using the SBMLR package (optional dependency).
#'
#' @param ... Arguments to filter models. Either:
#'     - A data frame of filtered models from metabolic_atlas_models()
#'     - Integer vector: model ID(s) from the `id` column
#'     - Expressions passed to dplyr::filter() (e.g., organism == 'Homo sapiens')
#'
#' @return A model object loaded by SBMLR, or the path to the downloaded
#'     SBML file if SBMLR is not available.
#'
#' @examples
#' \dontrun{
#' # Load model by ID
#' model <- metabolic_atlas_model(1)
#'
#' # Load first human model
#' model <- metabolic_atlas_model(organism == "Homo sapiens")
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang eval_tidy enexprs
#' @importFrom dplyr filter slice
#' @importFrom logger log_info log_error log_warn
#' @export
metabolic_atlas_model <- function(...) {

    sbmlr_available <- requireNamespace("SBMLR", quietly = TRUE)

    if (!sbmlr_available) {
        msg <- "SBMLR package not available, install with: install.packages('SBMLR')"
        log_warn(msg)
        warning(msg, call. = FALSE)
    }

    args <- list(...)

    # Handle different argument types
    models <- if (length(args) == 1 && is.data.frame(args[[1]])) {
        args[[1]]
    } else if (length(args) == 1 && is.numeric(args[[1]])) {
        metabolic_atlas_models() %>% filter(id %in% args[[1]])
    } else {
        filter_exprs <- enexprs(...)
        models <- metabolic_atlas_models()
        for (expr in filter_exprs) {
            models <- models %>% filter(eval_tidy(expr, data = models))
        }
        models
    }

    # NSE vs. R CMD check workaround
    id <- name <- files <- NULL

    models <- models %>% slice(1)

    sbml_path <- models$files[[1]] %>%
        filter(format == "SBML") %>%
        slice(1) %>%
        pull(path)

    log_info("Downloading model: %s (ID: %s)", models$name[1], models$id[1])

    if (sbmlr_available) {
        archive_extractor(
            url_key = "metatlas_model",
            url_param = list(sbml_path),
            reader = SBMLR::readSBML
        )
    } else {
        download_to_cache(
            url_key = "metatlas_model",
            url_param = list(sbml_path)
        )
    }
}
