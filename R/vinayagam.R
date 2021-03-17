#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2021
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


#' Protein-protein interactions from Vinayagam 2011
#'
#' Retrieves the Supplementary Table S6 from Vinayagam et al. 2011.
#' Find out more at \url{https://doi.org/10.1126/scisignal.2001699}.
#'
#' @return A data frame (tibble) with interactions.
#'
#' @examples
#' vinayagam_interactions <- vinayagam_download()
#' vinayagam_interactions
#' # # A tibble: 34,814 x 5
#' #    `Input-node Gen. `Input-node Gen. `Output-node Ge. `Output-node Ge.
#' #    <chr>                       <dbl> <chr>                       <dbl>
#' #  1 C1orf103                    55791 MNAT1                        4331
#' #  2 MAST2                       23139 DYNLL1                       8655
#' #  3 RAB22A                      57403 APPL2                       55198
#' #  4 TRAP1                       10131 EXT2                         2132
#' #  5 STAT2                        6773 COPS4                       51138
#' # # . with 34,804 more rows, and 1 more variable:
#' # # `Edge direction score` <dbl>
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readxl read_xls
#' @importFrom utils unzip
#' @export
vinayagam_download <- function(){

    # NSE vs. R CMD check workaround
    from_cache <- vinayagam_url <- NULL

    top_env <- environment()

    'omnipath.vinayagam_url' %>%
    archive_downloader() %T>%
    {assign('from_cache', .$from_cache, envir = top_env)} %T>%
    {assign('vinayagam_url', .$url, envir = top_env)} %>%
    `$`('path') %>%
    unzip(
        files = '2001699_Tables_S1_S2_S6.xls',
        exdir = tempdir()
    ) %>%
    `[`(1) %>%
    read_xls(
        sheet = 'S6',
        progress = FALSE
    ) %>%
    origin_cache(from_cache) %>%
    source_attrs('Vinayagam et al. 2011', url = vinayagam_url) %T>%
    load_success()

}