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
    grep('^[[:alnum:]]+$', ., value = TRUE) %>%
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
#' @importFrom rlang %||%
#'
#' @noRd
load_success <- function(data){

    from_cache <- data %>% attr('origin') %>% {!is.null(.) && . == 'cache'}

    '%s: %sloaded %d records%s' %>%
    sprintf(
        if_null_len0(attr(data, 'source'), 'Unknown source'),
        `if`(from_cache, '', 'down'),
        data %>% {if_null(nrow(.), length(.))} %||% 0,
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


#' Returns \code{NULL} if `value` is \code{NULL}, otherwise executes the
#' function `fun` on `value`.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
null_or_call <- function(value, fun, ...){

    value %>%
    is.null %>%
    `if`(NULL, fun(value, ...))

}


#' Returns `value1` if it's not zero length otherwise `value2`
#'
#' @importFrom magrittr %>%
#'
#' @noRd
if_empty <- function(value1, value2){

    value1 %>%
    length %>%
    `==`(0) %>%
    `if`(value2, value1)

}


#' Returns `value1` if it's not NULL or zero length otherwise `value2`
#'
#' @importFrom magrittr %>%
#'
#' @noRd
if_null_len0 <- function(value1, value2){

    value1 %>%
    {is.null(.) || length(.) == 0} %>%
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


#' Workaround against R CMD check `import not declared from` warnings
#'
#' @noRd
`%::%` <- function(pkg, name)
{
    pkg <- as.character(substitute(pkg))
    name <- as.character(substitute(name))
    getExportedValue(pkg, name)

}



#' Calls a function without throwing error on empty vectors
#'
#' Many functions can not tolerate if their input vector is empty. However,
#' especially in longer workflows, it can happen easily that we end up with
#' an empty vector. Wrapping the call into this method provides a control
#' over those cases and avoids error.
#'
#' @param v Vector, potentially empty or `NULL`.
#' @param fun Function to call.
#' @param .default The value to return if the vector is empty.
#' @param ... Arguments for `fun`.
#'
#' @noRd
empty_no_problem <- function(v, fun, .default = NULL, ...){

    `if`(
        is.null(v) || length(v) == 0,
        .default,
        fun(v, ...)
    )

}


#' Distinguish a list from a data frame
#'
#' \code{is.list} returns `TRUE` for data frames. If we want to check if
#' an object is a list but not a data frame we need to double check.
#'
#' @param obj The object to check, any kind of object.
#'
#' @noRd
just_a_list <- function(obj){

    is.list(obj) && !is.data.frame(obj)

}


#' Compare two lists for complete equality
#'
#' @param list1 A list.
#' @param list2 Another list.
#'
#' @return Logical.
#'
#' @importFrom magrittr %>%
#' @noRd
lists_identical <- function(list1, list2){

    keys1 <- list1 %>% names %>% sort
    keys2 <- list2 %>% names %>% sort
    keys <- union(keys1, keys2)

    if(is.null(keys1) && is.null(keys2)){
        return(identical(list1, list2))
    }

    if(any(keys1 != keys2)){
        return(FALSE)
    }

    for(key in keys){
        val1 <- list1[[key]]
        val2 <- list2[[key]]
        if(length(val1) != length(val2)){
            return(FALSE)
        }else if(is.list(val1) && is.list(val2)){
            if(!list_identical(val1, val2)){
                return(FALSE)
            }
        }else if(!is.atomic(val1) || !is.atomic(val2)){
            return(FALSE)
        }else if(any(sort(val1) != sort(val2))){
            return(FALSE)
        }
    }

    return(TRUE)

}