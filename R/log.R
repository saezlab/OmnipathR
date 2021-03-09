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
#'
#' @noRd
ensure_loglevel <- function(level){

    `if`(
        'loglevel' %in% class(level),
        level,
        get(toupper(level), envir = getNamespace('logger'))
    )

}


#' Returns the log level threshold for console messages, according to the
#' current settings
#'
#' @noRd
omnipath_console_loglevel <- function(){

    ensure_loglevel(options('omnipath.console_loglevel')[[1]])

}


#' Returns the path to the current OmnipathR log file
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \code{\link{omnipath_log}}
omnipath_logfile <- function(){

    # NSE vs. R CMD check workaround
    get_logger_definitions <- NULL

    (logger%:::%get_logger_definitions)(
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


#' Sets the log level for the package logger
#'
#' @param level Character or class `loglevel`. The desired log level.
#' @param target Character, either 'logfile' or 'console'
#'
#' @export
#' @importFrom magrittr %<>% %>%
omnipath_set_loglevel <- function(level, target = 'logfile'){

    # NSE vs. R CMD check workaround
    namespaces <- NULL

    level %<>% ensure_loglevel
    i_logger <- target %>% `==`(c('logfile', 'console')) %>% which

    omnipathr_loggers <- (logger%:::%namespaces)$OmnipathR
    omnipathr_loggers[[i_logger]]$threshold <- level

    assign(
        'OmnipathR',
        omnipathr_loggers,
        envir = logger%:::%namespaces
    )

    invisible(NULL)

}


#' Sets the log level for the console
#'
#' Use this method to change during a session which messages you want to be
#' printed on the console. Before loading the package, you can set it also by
#' the config file, with the omnipath.console_loglevel key.
#'
#' @param level Character or class `loglevel`. The desired log level.
#'
#' @export
#' @seealso \code{\link{omnipath_set_logfile_loglevel}}
omnipath_set_console_loglevel <- function(level){

    omnipath_set_loglevel(level = level, target = 'console')

}


#' Sets the log level for the logfile
#'
#' Use this method to change during a session which messages you want to be
#' written into the logfile. Before loading the package, you can set it also
#' by the config file, with the omnipath.loglevel key.
#'
#' @param level Character or class `loglevel`. The desired log level.
#'
#' @export
#' @seealso \code{\link{omnipath_set_console_loglevel}}
omnipath_set_logfile_loglevel <- function(level){

    omnipath_set_loglevel(level = level, target = 'logfile')

}