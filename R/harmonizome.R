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


#' Downloads a Harmonizome network dataset
#'
#' Downloads a single network dataset from Harmonizome
#' \url{https://maayanlab.cloud/Harmonizome}.
#'
#' @param dataset The dataset part of the URL. Please refer to the download
#'     section of the Harmonizome webpage.
#'
#' @return Data frame (tibble) with interactions.
#'
#' @examples
#' harmonizome_data <- harmonizome_download('phosphositeplus')
#' harmonizome_data
#' # # A tibble: 6,013 x 7
#' #    source   source_desc source_id target target_desc target_id weight
#' #    <chr>    <chr>           <dbl> <chr>  <chr>           <dbl>  <dbl>
#' #  1 TP53     na               7157 STK17A na               9263      1
#' #  2 TP53     na               7157 TP53RK na             112858      1
#' #  3 TP53     na               7157 SMG1   na              23049      1
#' #  4 UPF1     na               5976 SMG1   na              23049      1
#' # # . with 6,003 more rows
#'
#' @importFrom readr read_tsv read_lines
#' @importFrom magrittr %>% %T>%
#' @export
harmonizome_download <- function(dataset){

    path <-
        download_to_cache(
            url_key = 'omnipath.harmonizome_url',
            url_param = list(dataset)
        )

    path %>%
    gzfile() %>%
    read_lines(progress = FALSE) %>%
    `[`(-2) %>%
    read_tsv(col_types = cols(), progress = FALSE) %>%
    copy_source_attrs(path, resource = 'Harmonizome') %T>%
    load_success()

}