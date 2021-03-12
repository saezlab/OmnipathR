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


#' Downloads ligand-receptor interactions from Ramilowski et al. 2015
#'
#' Curated ligand-receptor pairs from Supplementary Table 2 of the article
#' "A draft network of ligand-receptor mediated multicellular signaling in
#' human" (\url{https://www.nature.com/articles/ncomms8866}).
#'
#' @return A data frame (tibble) with interactions.
#'
#' @examples
#' rami_interactions <- ramilowski_download()
#' rami_interactions
#' # # A tibble: 2,557 x 16
#' #    Pair.Name Ligand.Approved. Ligand.Name Receptor.Approv.
#' #    <chr>     <chr>            <chr>       <chr>
#' #  1 A2M_LRP1  A2M              alpha-2-ma. LRP1
#' #  2 AANAT_MT. AANAT            aralkylami. MTNR1A
#' #  3 AANAT_MT. AANAT            aralkylami. MTNR1B
#' #  4 ACE_AGTR2 ACE              angiotensi. AGTR2
#' #  5 ACE_BDKR. ACE              angiotensi. BDKRB2
#' # # . with 2,547 more rows, and 12 more variables: Receptor.Name <chr>,
#' # #   DLRP <chr>, HPMR <chr>, IUPHAR <chr>, HPRD <chr>,
#' # #   STRING.binding <chr>, STRING.experiment <chr>, HPMR.Ligand <chr>,
#' # #   HPMR.Receptor <chr>, PMID.Manual <chr>, Pair.Source <chr>,
#' # #   Pair.Evidence <chr>
#'
#' @export
#' @importFrom magrittr %T>%
ramilowski_download <- function(){

    xls_downloader(
        url_key = 'omnipath.ramilowski_url',
        sheet = 'All.Pairs',
        resource = 'Ramilowski et al. 2015'
    ) %T>%
    load_success()

}