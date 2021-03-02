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
#' @importFrom magrittr %>%
#' @importFrom rlang enquo quo_get_expr quo_text is_symbol
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


#' Inserts an element to a list if the value is not null.
insert_if_not_null <- function(l, ...){

    elements <- list(...)

    for(name in names(elements)){

        if(!is.null(elements[[name]])){

            l[[name]] <- elements[[name]]

        }

    }

    return(l)

}


#' Makes sure value is a list
#'
#' This function is different from `as.list` because it returns a list with
#' a single NULL element if `value` is NULL
#'
#' @return A list
#'
ensure_list <- function(value){

    `if`(
        is.null(value),
        list(NULL),
        as.list(value)
    )

}


#' Returns NULL if `value` is a list with a single NULL element otherwise
#' the value itself
#'
list_null <- function(value){

    `if`(
        is.list(value) && length(value) == 1 && is.null(value[[1]]),
        NULL,
        value
    )

}


#' Returns the extension of a file name
#'
#' @importFrom magrittr %>%
file_extension <- function(name){

    name %>%
    strsplit('\\?') %>%
    `[[`(1) %>%
    `[`(1) %>%
    strsplit('\\.') %>%
    `[[`(1) %>%
    tail(
        `if`(
            length(.) > 1 && nchar(tail(., 1)) < 5,
            `if`(
                .[length(.) - 1] == 'tar',
                2,
                1
            ),
            0
        )
    ) %>%
    paste(collapse = '.')

}


#' Adds an extension
file_add_extension <- function(fname, ext){

    ifelse(
        is.character(ext) & ext != '' & !endsWith(fname, ext),
        paste(fname, ext, sep = '.'),
        fname
    )

}


#' Closes a connection when the parent call exits
close_on_exit <- function(con, envir = parent.frame()){

    assign('con', con, envir)
    do.call('on.exit', list(quote(close(con))), envir = envir)

}


#' Reads RDS from an URL
#'
#' Unfortunately the APIs of R built in functions and packages are very
#' diverse, for example, readRDS can not understand if an URL is passed to it.
#'
#' @importFrom magrittr %>%
url_rds <- function(URL){

    URL %>%
    url %>%
    readRDS

}