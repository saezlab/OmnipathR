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


#' Makes sure we have a string even if the argument was passed by NSE
#'
#' @importsFrom magrittr %>%
#' @importsFrom rlang enquo quo_get_expr quo_text is_symbol
.nse_ensure_str <- function(arg){

    enquo(arg) %>%
    {`if`(
        is_symbol(quo_get_expr(.)),
        quo_text(.),
        quo_get_expr(.)
    )}

}


.ensure_dir <- function(path){

    dir_path <- dirname(path)
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

    if(dir.exists(dir_path)){
        logger::log_debug('Created directory `%s`', dir_path)
    }else{
        logger::log_warning('Failed to create directory `%s`', dir_path)
    }

}