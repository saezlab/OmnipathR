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


#' Try how non standard evaluation works
#'
#' @noRd
nse_test <- function(a, ..., b = FALSE){

    print(a)
    print(b)

    print(enquos(...) %>% map(.nse_ensure_str))

}


#' A table with all UniProt IDs
#'
#' Retrieves a table from UniProt with all proteins for a certain organism.
#'
#' @param fields Character vector of fields as defined by UniProt. For
#'     possible values please refer to
#'     \url{https://www.uniprot.org/help/uniprotkb\%5Fcolumn\%5Fnames}
#' @param reviewed Retrieve only reviewed (`TRUE`), only unreviewed (`FALSE`)
#'     or both (`NULL`).
#' @param organism Integer, NCBI Taxonomy ID of the organism (by default
#'     9606 for human).
#'
#' @return Data frame (tibble) with the requested UniProt entries and fields.
#'
#' @importFrom rlang exec !!!
#' @export
#'
#' @examples
#' human_swissprot_ac <- all_uniprots(fields = 'entry name')
#' human_swissprot_ac
#' # # A tibble: 20,396 x 1
#' #    `Entry name`
#' #    <chr>
#' #  1 OR4K3_HUMAN
#' #  2 O52A1_HUMAN
#' #  3 O2AG1_HUMAN
#' #  4 O10S1_HUMAN
#' #  5 O11G2_HUMAN
#' # # . with 20,386 more rows
all_uniprots <- function(fields = 'id', reviewed = TRUE, organism = 9606){

    exec(.all_uniprots, !!!as.list(environment()))

}


#' R CMD check workaround, see details at \code{all_uniprots}
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom logger log_trace
#'
#' @noRd
.all_uniprots <- uniprot_domains %@% function(
    fields = 'id',
    reviewed = TRUE,
    organism = 9606,
    .subdomain = 'www'
){

    fields <- fields %>% paste(collapse = ',')

    log_trace(
        paste0(
            'Loading all UniProt records for organism %d ',
            '(only reviewed: %s); fields: %s'
        ),
        organism, reviewed, fields
    )

    reviewed <- `if`(
        is.null(reviewed),
        '',
        sprintf(' AND reviewed:%s', `if`(reviewed, 'yes', 'no'))
    )

    generic_downloader(
        url_key = 'all_uniprots',
        url_param = list(.subdomain, fields, organism, reviewed),
        reader_param = list(progress = FALSE),
        resource = 'UniProt'
    ) %T>%
    load_success()

}
