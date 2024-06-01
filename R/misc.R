#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
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


#' Makes sure we have a string even if the argument was passed by NSE
#'
#' @importFrom magrittr %>%
#' @importFrom rlang enquo quo_get_expr quo_text is_symbol
#' @importFrom stringr str_remove
#' @noRd
.nse_ensure_str <- function(arg){

    enquo(arg) %>%
    {`if`(
        is_symbol(quo_get_expr(.)),
        quo_text(.),
        quo_get_expr(.)
    )} %>%
    str_remove('`$') %>%
    str_remove('^`')

}


#' Makes sure the ellipsis is a named vector of strings
#'
#' @importFrom rlang %||% enquos
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @noRd
ellipsis_to_char <- function(...) {

    enquos(...) %>%
    map(.nse_ensure_str) %>%
    set_names(names(.) %||% unlist(.)) %>%
    set_names(ifelse(nchar(names(.)), names(.), unlist(.)))

}


#' @noRd
.ensure_dir <- function(path){

    dir_path <- dirname(path)
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

    if(dir.exists(dir_path)){
        logger::log_debug('Created directory `%s`', dir_path)
    }else{
        logger::log_warn('Failed to create directory `%s`', dir_path)
    }

}


#' @noRd
.path_writable <- function(path){

    path_prev <- ''

    while(is.character(path) && path != path_prev){

        if(file.access(path, mode = 2L) == 0L) {

            return(TRUE)

        }else if(file.exists(path)){

            return(FALSE)

        }else{

            prev_path <- path
            path <- dirname(path)

        }

    }

    return(FALSE)

}


#' @noRd
.ensure_safe_path <- function(path, directory = FALSE){

    if(!.path_writable(path)){

        fname <- `if`(directory, '', basename(path))
        fallback_path <- file.path(tempdir(check = TRUE), fname)
        warning(sprintf(
            'OmnipathR: Falling back to path `%s` instead of `%s`.',
            fallback_path,
            path
        ))
        path <- fallback_path

    }

    path

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


#' Makes sure value is a list
#'
#' If \code{value} is a list returns it unchanged, otherwise it wraps it
#' into a single element list.
#'
#' @return A list
#'
#' @importFrom magrittr %>% equals
#' @noRd
ensure_list_2 <- function(value){

    `if`(
        class(value)[1] %>% equals('list'),
        value,
        list(value)
    )

}


#' Merge two lists by name
#'
#' From the RCurl package. This is a method that merges the contents of one
#' list with another by adding the named elements in the second that are not
#' in the first. In other words, the first list is the target template, and
#' the second one adds any extra elements that it has.
#'
#' @param x The list to which elements will be added.
#' @param y The list which will supply additional elements to ‘x’ that
#'     are not already there by name.
#' @param ... Ignored.
#'
#' @return A named list whose name set is the union of the elements in names
#'    of x and y and whose values are those taken from y and then with
#'    those in x, overwriting if necessary.
#'
#' @author Duncan Temple Lang
#'
#' @noRd
merge_lists <- function (x, y, ...)
{
    if(length(x) == 0) return(y)

    if(length(y) == 0) return(x)

    i <- match(names(y), names(x))

    i <- is.na(i)

    if(any(i)){

        x[names(y)[which(i)]] <- y[which(i)]

    }

    x

}


#' Returns NULL if `value` is an empty list or a list with a single NULL
#' element otherwise the value itself
#'
#' @importFrom utils tail
#'
#' @noRd
list_null <- function(value){

    `if`(
        is.list(value) && (
            length(value) == 0 ||
            length(value) == 1 && is.null(value[[1]])
        ),
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
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang %||%
#'
#' @noRd
load_success <- function(data, resource = NULL, from_cache = NULL){

    from_cache %<>% if_null(data %>% is_from_cache)
    resource %<>% if_null(if_null_len0(attr(data, 'source'), 'Unknown source'))

    '%s: %sloaded %d records%s' %>%
    sprintf(
        resource,
        `if`(from_cache, '', 'down'),
        data %>% {if_null(nrow(.), length(.))} %||% 0,
        `if`(from_cache, ' from cache', '')
    ) %>%
    logger::log_success()

}


#' Is the object labelled as cache origin
#'
#' @param obj Any object, typically a data frame.
#'
#' @noRd
is_from_cache <- function(obj){

    obj %>% attr('origin') %>% {!is.null(.) && . == 'cache'}

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


#' Call a function if available
#'
#' @importFrom magrittr %>%
#' @importFrom rlang is_function
#' @noRd
maybe_call <- function(fun, ..., .return_args = FALSE) {

    fun %>%
    `if`(is_function(.), ., maybe_get(.)) %>%
    {`if`(
        is_function(.),
        .(...),
        `if`(
            is.null(.) || .return_args,
            list(...) %>% unlist_len1,
            .
        )
    )}

}


#' Try to get name, do nothing if it is not available
#'
#' @noRd
maybe_get <- function(name) {

    tryCatch(get(name), error = function(e) {return(name)})

}


#' Returns the first element for length one lists or vectors
#'
#' @importFrom dplyr first
#'
#' @noRd
unlist_len1 <- function(x) {

    `if`(length(x) <= 1L, first(x), x)

}


#' Returns `value1` if it's not zero length otherwise `value2`
#'
#' @importFrom magrittr %>%
#'
#' @noRd
if_empty <- function(value1, value2){

    value1 %>%
    is_empty %>%
    `if`(value2, value1)

}


#' Tells if value is NULL or a vector of length zero
#'
#' @importFrom magrittr %>% equals
#'
#' @noRd
is_empty <- function(value){

    value %>% length %>% equals(0L)

}



#' Tells if value is NULL or a vector of length zero or an empty string
#'
#' @importFrom magrittr %>% extract
#'
#' @noRd
is_empty_2 <- function(value){

    value %>%
    {
        is.null(.) ||
        length(.) == 0L ||
        !is_closure(.) &&
        !is.na(extract(., 1L)) &&
        extract(., 1L) == ''
    }

}


#' Tells if value is closure-like
#'
#' @importFrom magrittr %>% is_in
#' @noRd
is_closure <- function(value) {

    value %>%
    typeof %>%
    is_in(c('closure', 'special', 'builtin'))

}

#' Returns `value1` if it's not NULL or zero length otherwise `value2`
#'
#' @importFrom magrittr %>%
#'
#' @noRd
if_null_len0 <- function(value1, value2){

    value1 %>%
    is_empty_2 %>%
    `if`(value2, value1)

}


#' Same as `if_null_len0` but considers empty string empty too
#'
#' @importFrom magrittr %>%
#'
#' @noRd
if_null_len0_2 <- function(value1, value2){

    value1 %>%
    is_empty_2 %>%
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


#' Wrapper around loadNamespace to avoid R CMD check `call not declared from`
#' warnings.
#'
#' @param pkgname Character: name of the package.
#'
#' @return The package namespace loaded.
#'
#' @noRd
load_namespace <- function(pkgname){

    loadNamespace(pkgname)

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
            if(!lists_identical(val1, val2)){
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


#' Splits a vector into a list of chunks of a certain size
#'
#' @importFrom magrittr %>% equals
#'
#' @noRd
chunks <- function(v, size){

    v %>%
    split(seq_along(.) %>% `-`(1) %>% `%%`(size) %>% equals(0) %>% cumsum)

}


#' Pretty print a list of words
#'
#' @param words Character vector with words.
#' @param quotes Logical: wrap all words in backticks.
#'
#' @importFrom magrittr %>%
#' @importFrom utils head tail
#'
#' @noRd
pretty_list <- function(words, quotes = TRUE){

    words %>%
    sort %>%
    {`if`(quotes, sprintf('`%s`', .), .)} %>%
    {sprintf(
        '%s%s%s',
        paste(head(., -1), collapse = ', '),
        `if`(length(.) > 1, ' and ', ''),
        tail(., 1)
    )}

}


#' Plural form of a word depending on the count of the objects
#'
#' @param objects The objects to count.
#' @param word Character: the name of the objects.
#' @param irregular Character: irregular plural form (anything else than
#'     adding an "s").
#'
#' @importFrom rlang %||%
#'
#' @noRd
plural <- function(objects, word, irregular = NULL){

    `if`(
        length(objects) > 1,
        irregular %||% sprintf('%ss', word),
        word
    )

}


#' Convert an atomic variable to the type of another atomic variable
#'
#' @param var The variable to be converted.
#' @param type Convert to the type of this variable.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
as_type <- function(var, type){

    if(!is.atomic(var) || !is.atomic(type) || is.null(type)){

        var

    }else{

        (type %>% typeof %>% sprintf('as.%s', .) %>% get)(var)

    }

}


#' Retrieve the open connection(s) pointing to URI
#'
#' @param uri Character: path or URL the connection points to.
#'
#' @return A list of connection objects.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter pull
#' @importFrom purrr map
#' @noRd
get_connections <- function(uri){

    # NSE vs. R CMD check workaround
    description <- con_id <- NULL

    showConnections(all = TRUE) %>%
    as.data.frame %>%
    rownames_to_column('con_id') %>%
    filter(description == uri) %>%
    pull(con_id) %>%
    as.integer %>%
    map(getConnection)

}


#' Closes the open connection(s) pointing to URI
#'
#' @param uri Character: path or URL the connection points to.
#'
#' @return Invisible `NULL`.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr walk
#' @noRd
close_connection <- function(uri){

    uri %>%
    get_connections %>%
    walk(close)

    invisible(NULL)

}


#' Replace NULL elements by NAs in a list
#'
#' @param l A list.
#'
#' @return A list with NULL elements replaced by NAs.
#'
#' @importFrom purrr map
#' @noRd
null_to_na <- function(l){

    map(l, function(x){`if`(is.null(x), NA, x)})

}


#' Add indentation in front of lines
#'
#' @param lines Character vector.
#' @param depth Integer: number of spaces.
#'
#' @return Character vector of indented lines.
#'
#' @importFrom magrittr %>%
#' @noRd
indent <- function(lines, depth = 4L){

    ' ' %>%
    rep(depth) %>%
    paste(collapse = '') %>%
    paste0(lines)

}


#' Read a JSON
#'
#' @importFrom logger log_trace log_warn
#' @importFrom jsonlite validate fromJSON
#' @noRd
safe_json <- function(path, encoding = 'UTF-8', ...){

    lines <- suppressWarnings(readLines(con = path, encoding = encoding))

    log_trace('Reading JSON from `%s` (encoding: %s).', path, encoding)

    json_ok <- jsonlite::validate(txt = lines)

    log_trace('JSON validation successful: %s', json_ok)

    if(!json_ok){

        msg <- sprintf(
            'Failed to read JSON from `%s` (file exists: %s)',
            path,
            file.exists(path)
        )

        logger::log_warn(msg)
        warning(msg)

    }

    jsonlite::fromJSON(txt = `if`(json_ok, lines, '[]'), ...)

}


#' Is the package running on a build server?
#'
#' Build servers are for automated testing of packages. Bioconductor and
#' Saez Lab (developers & maintainers of this package) both run their own
#' build servers.
#'
#' @noRd
.on_buildserver <- function(){

    Sys.info()['user'] %in% c('biocbuild', 'omnipath')

}


#' Is the package running on a Bioconductor build server?'
#'
#' @noRd
.on_bioc_buildserver <- function(){

    Sys.info()['user'] == 'biocbuild'

}


#' Bypass slow doctests
#'
#' Makes the calling function return a value, bypassing further execution.
#' Used in functions that normally take a long time to run, and hence, due
#' to the 40 minutes timeout on Bioconductor build servers, we can not afford
#' to run them. We still run these doctests daily on our own build server,
#' so they fulfill their testing purpose.
#'
#' @param value The value to return.
#'
#' @return This function does not return anything, its parent returns
#'     `value`.
#'
#' @noRd
.slow_doctest <- function(value = NULL, env = parent.frame()){

    if(.on_buildserver()){

        log_trace('Bypassing call: `%s`.', deparse(sys.call(-1L)))
        log_trace(paste(
            'This behaviour is intended for running R CMD check',
            'within limited time and is triggered solely by the user name.',
            'Please report if you see this anywhere outside of a',
            'Bioconductor build server.'
        ))
        do.call(return, list(value), envir = env)

    }

}


#' Packages that seem to be missing from library
#'
#' @importFrom magrittr %>%
#' @importFrom utils installed.pacakges
#' @noRd
missing_packages <- function(pkgs) {

    installed.packages() %>%
    rownames %>%
    setdiff(pkgs, .)

}


#' Joins a list of words by comma and space
#'
#' @noRd
enum_format <- function(words) {

    paste(words, collapse = ', ')

}


#' Unlist, unique and sort a list or vector
#'
#' @param v List or vector.
#'
#' @return A vector of unique values sorted.
#'
#' @importFrom magrittr %>%
#' @noRd
unique_sorted <- function(v) {

    v %>% unlist %>% unique %>% sort

}


#' Joins a list of words by comma and space
#'
#' @noRd
enum_format <- function(words) {

    paste(words, collapse = ', ')

}


#' Uncompressed size of a connection
#'
#' @param con A connection.
#'
#' @return Integer: uncompressed size of the file the connection points to.
#'
#' @importFrom zip zip_list
#' @importFrom magrittr %>% extract extract2
#' @importFrom stringr str_split
#' @noRd
file_size <- function(con) {

    con %>%
    summary %>%
    extract2('description') %>%
    {`if`(
        class(con)[1L] == 'unz',

        str_split(., ':') %>%
        extract2(1L) %>%
        {extract(zip_list(.[1L]), 1L, 'uncompressed_size')},

        {`if`(
            class(con)[1L] == 'file',
            file.info(.) %>%
            extract(1L, 'size'),
            NA
        )}
    )}

}
