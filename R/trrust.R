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


#' Downloads TF-target interactions from TRRUST
#'
#' TRRUST v2 (https://www.grnpedia.org/trrust/) is a database of literature
#' mined TF-target interactions for human and mouse.
#'
#' @param organism Character: either "human" or "mouse".
#'
#' @return A data frame of TF-target interactions.
#'
#' @examples
#' trrust_interactions <- trrust_download()
#' trrust_interactions
#' # # A tibble: 11,698 x 4
#' #    source_genesymbol target_genesymbol effect reference
#' #    <chr>             <chr>              <dbl> <chr>
#' #  1 AATF              BAX                   -1 22909821
#' #  2 AATF              CDKN1A                 0 17157788
#' #  3 AATF              KLK3                   0 23146908
#' #  4 AATF              MYC                    1 20549547
#' #  5 AATF              TP53                   0 17157788
#' #  6 ABL1              BAX                    1 11753601
#' #  7 ABL1              BCL2                  -1 11753601
#' # # . with 11,688 more rows
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom dplyr mutate recode
#' @importFrom tidyr separate_rows
trrust_download <- function(organism = 'human'){

    # NSE vs. R CMD check workaround
    effect <- reference <- NULL

    generic_downloader(
        url_key = 'omnipath.trrust_url',
        url_param = list(organism),
        reader_param = list(
            col_names = c(
                'source_genesymbol',
                'target_genesymbol',
                'effect',
                'reference'
            )
        ),
        resource = 'TRRUST'
    ) %>%
    mutate(
        effect = recode(
            effect,
            Repression = -1,
            Unknown    =  0,
            Activation =  1
        )
    ) %>%
    separate_rows(
        reference,
        sep = ';'
    ) %T>%
    load_success()

}