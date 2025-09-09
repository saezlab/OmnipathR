#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
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
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Kinase-PTM relationships from KinasePhos 3.0
#'
#' Predicted kinase-PTM interactions from https://awi.cuhk.edu.cn/KinasePhos/.
#'
#' @examples
#' kinasephos()
#'
#' @return A data frame (tibble) of kinase-PTM interactions.
#'
#' @importFrom magrittr %>% %T>%
#' @export
kinasephos <- function() {

    'kinasephos' %>%
    generic_downloader(
        resource = 'KinasePhos3.0',
        reader_param = list(delim = ',')
    ) %T>%
    load_success()

}
