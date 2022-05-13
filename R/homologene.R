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


#' Orthology data from NCBI HomoloGene
#'
#' Retrieves NCBI HomoloGene data without any processing. Processed tables
#' are more useful for most purposes, see below other functions that provide
#' those.
#'
#' @return A data frame as provided by NCBI HomoloGene.
#'
#' @examples
#' hg <- homologene_raw()
#' hg
#' # # A tibble: 275,237 × 6
#' #    hgroup ncbi_taxid  entrez genesymbol             gi refseqp
#' #     <dbl>      <dbl>   <dbl> <chr>               <dbl> <chr>
#' #  1      3       9606      34 ACADM             4557231 NP_000007.1
#' #  2      3       9598  469356 ACADM           160961497 NP_001104286.1
#' #  3      3       9544  705168 ACADM           109008502 XP_001101274.1
#' #  4      3       9615  490207 ACADM           545503811 XP_005622188.1
#' #  5      3       9913  505968 ACADM           115497690 NP_001068703.1
#' # # . with 275,232 more rows
#'
#' @importFrom magrittr %>%
#' @importFrom readr cols
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{}}}
#' }
homologene_raw <- function(){

    hdr <- c('hgroup', 'ncbi_taxid', 'entrez', 'genesymbol', 'gi', 'refseqp')

    'homologene' %>%
    generic_downloader(
        reader_param = list(
            col_names = hdr,
            col_types = cols()
        ),
        resource = 'NCBI HomoloGene'
    ) %T>%
    load_success()

}
