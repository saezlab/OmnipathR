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

#' Ensures that `level` is a `loglevel` class object, converts if necessary
ensure_loglevel <- function(level){

    `if`(
        'loglevel' %in% class(level),
        level,
        get(toupper(level), envir = getNamespace('logger'))
    )

}


#' Returns the log level threshold for console messages, according to the
#' current settings
omnipath_console_loglevel <- function(){

    ensure_loglevel(options('omnipath.console_loglevel')[[1]])

}


#' Returns the path to the current OmnipathR log file
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \code{\link{omnipath_log}}
omnipath_logfile <- function(){

    logger:::get_logger_definitions(
        namespace = 'OmnipathR'
    )$default$appender %>%
    environment() %>%
    `$`('file') %>%
    normalizePath()

}


#' Browse the current OmnipathR log file
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \code{\link{omnipath_logfile}}
omnipath_log <- function(){

    omnipath_logfile() %>%
    file.show(title = 'OmnipathR log')

}