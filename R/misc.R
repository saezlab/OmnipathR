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
#'
#' @noRd
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
        logger::log_warn('Failed to create directory `%s`', dir_path)
    }

}


#' Inserts an element to a list if the value is not null.
#'
#' @noRd
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
#' @noRd
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
#' @importFrom utils tail
#'
#' @noRd
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
#'
#' @noRd
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


#' Adds an extension to a path or file name
#'
#' @noRd
file_add_extension <- function(fname, ext){

    ifelse(
        is.character(ext) & ext != '' & !endsWith(fname, ext),
        paste(fname, ext, sep = '.'),
        fname
    )

}


#' Converts any path to absolute path
#'
#' This function is used to convert relative file paths to absolute file paths
#' without checking if the file exists as \code{tools::file_path_as_absolute}
#' does. Modified from
#' https://github.com/sbfnk/RBi/blob/master/R/absolute_path.R
#'
#' @param filename name of a file, absolute or relative to a folder
#' @param dirname name of a folder where the file is supposed to be
#'
#' @return a character string containing the absolute path
#'
#' @noRd
absolute_path <- function(path, winslash = '\\'){

    if(
        (
            .Platform$OS.type == 'unix' &&
            substr(path, 1, 1) %in% c('/', '~')
        ) || (
            .Platform$OS.type == 'windows' &&
            grepl('^[A-z]:[\\\\/]', path)
        )
    ){
        # the path is already absolute
        result <- normalizePath(path, winslash, FALSE)
    } else {
        # prepending the current directory to make it absolute
        result <- file.path(normalizePath(getwd(), winslash), path)
    }

    return(result)

}


#' Closes a connection when the parent call exits
#'
#' @noRd
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
#'
#' @noRd
url_rds <- function(URL){

    URL %>%
    url %>%
    readRDS

}


#' Issues a log message about successful loading of a dataset
#'
#' @param data A data frame.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
load_success <- function(data){

    from_cache <- data %>% attr('origin') %>% {!is.null(.) && . == 'cache'}

    '%s: %sloaded %d records%s' %>%
    sprintf(
        attr(data, 'source'),
        `if`(from_cache, '', 'down'),
        data %>% {if_null(nrow(.), length(.))},
        `if`(from_cache, ' from cache', '')
    ) %>%
    logger::log_success()

}


#' Adds default parameters to a list of function arguments
#'
#' @param list The actual parameters typically provided by the user.
#' @param fun The function the parameters will be passed to.
#' @param defaults Named list with the default arguments which should be
#'     added to `param` if they are suitable for `fun` and if they haven't
#'     been overridden by `param`.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr walk2
#'
#' @noRd
add_defaults <- function(param, fun, defaults){

    this_env <- environment()

    defaults %>%
    walk2(
        names(.),
        function(value, key){
            if(
                !(key %in% names(param)) &&
                key %in% names(formals(fun))
            ){
                param[[key]] <- value
                assign('param', param, this_env)
            }
        }
    )

    return(param)

}


#' Copies attributes from one object to another
#'
#' @param to The object to copy the attribute to.
#' @param from The object to copy the attribute from.
#' @param names Character: the names of the attributes
#'
#' @importFrom magrittr %>%
#' @importFrom purrr reduce
#'
#' @noRd
copy_attrs <- function(to, from, names){

    names %>%
    reduce(
        function(target, name){
            target %>%
            `attr<-`(name, attr(from, name))
        },
        .init = to
    )

}


#' Returns `value1` if it's not NULL otherwise `value2`
#'
#' @importFrom magrittr %>%
#'
#' @noRd
if_null <- function(value1, value2){

    value1 %>%
    is.null %>%
    `if`(value2, value1)

}


#' Workaround against R CMD check notes about using `:::`
#'
#' @importFrom rlang enquo !!
#'
#' @noRd
`%:::%` <- function(pkg, fun){

    pkg <- .nse_ensure_str(!!enquo(pkg))
    fun <- .nse_ensure_str(!!enquo(fun))

    get(fun, envir = asNamespace(pkg), inherits = FALSE)

}