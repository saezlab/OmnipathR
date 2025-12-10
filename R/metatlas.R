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
