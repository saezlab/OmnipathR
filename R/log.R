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

    console_loglevel('OmnipathR')

}


#' Log level threshold for console messages from a certain package
#'
#' @importFrom magrittr %>%
#' @noRd
console_loglevel <- function(pkg = 'OmnipathR'){


    'console_loglevel' %>%
    pkg_getopt(pkg) %>%
    ensure_loglevel

}


#' Path to the current OmnipathR log file
#'
#' @return Character: path to the current logfile, or \code{NULL} if no
#'     logfile is available.
#'
#' @examples
#' omnipath_logfile()
#' # [1] "/home/denes/omnipathr/omnipathr-log/omnipathr-20210309-1642.log"
#'
#' @export
#' @seealso \code{\link{omnipath_log}}
omnipath_logfile <- function(){

    logfile('OmnipathR')

}


#' Path to the current logfile of a package
#'
#' @param pkg Character: name of a package.
#'
#' @importFrom magrittr %>%
#' @rdname omnipath_logfile
#' @export
logfile <- function(pkg = 'OmnipathR') {

    # NSE vs. R CMD check workaround
    get_logger_definitions <- NULL

    if(!has_logfile(pkg)) return()

    (logger%:::%get_logger_definitions)(
        namespace = pkg
    )[[2]]$appender %>%
    environment() %>%
    `$`('file') %>%
    normalizePath()

}


#' Browse the current OmnipathR log file
#'
#' @return Returns `NULL`.
#'
#' @examples
#' \dontrun{
#' omnipath_log()
#' # then you can browse the log file, and exit with `q`
#' }
#'
#' @export
#' @seealso \code{\link{omnipath_logfile}}
omnipath_log <- function(){

    read_log('OmnipathR')

}


#' Browse the latest log from a package
#'
#' @param pkg Character: name of a package.
#'
#' @importFrom magrittr %>%
#' @rdname omnipath_log
#' @export
read_log <- function(pkg = 'OmnipathR') {

    if(!has_logfile(pkg)) return()

    pkg %>%
    logfile() %>%
    file.show(title = sprintf('%s log', pkg))

}


#' Sets the log level for the package logger
#'
#' @param level Character or class `loglevel`. The desired log level.
#' @param target Character, either 'logfile' or 'console'
#'
#' @return Returns `NULL`.
#'
#' @examples
#' omnipath_set_loglevel(logger::FATAL, target = 'console')
#'
#' @export
omnipath_set_loglevel <- function(level, target = 'logfile'){

    set_loglevel(level, target = target, pkg = 'OmnipathR')

}


#' Sets the log level for any package
#'
#' @param pkg Character: name of the package.
#'
#' @importFrom magrittr %<>% %>% equals extract2
#' @export
#' @rdname omnipath_set_loglevel
set_loglevel <- function(level, target = 'logfile', pkg = 'OmnipathR') {

    # NSE vs. R CMD check workaround
    namespaces <- NULL

    level %<>% ensure_loglevel
    i_logger <- target %>% equals(c('console', 'logfile')) %>% which

    loggers <- (logger%:::%namespaces) %>% extract2(pkg)
    loggers[[i_logger]]$threshold <- level

    assign(pkg, loggers, envir = logger%:::%namespaces)

    invisible(NULL)

}


#' Sets the log level for the console
#'
#' Use this method to change during a session which messages you want to be
#' printed on the console. Before loading the package, you can set it also by
#' the config file, with the omnipathr.console_loglevel key.
#'
#' @param level Character or class `loglevel`. The desired log level.
#'
#' @return Returns `NULL`.
#'
#' @examples
#' omnipath_set_console_loglevel('warn')
#' # or:
#' omnipath_set_console_loglevel(logger::WARN)
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
#' by the config file, with the "omnipathr.loglevel" key.
#'
#' @param level Character or class `loglevel`. The desired log level.
#'
#' @return Returns `NULL`.
#'
#' @examples
#' omnipath_set_logfile_loglevel('info')
#' # or:
#' omnipath_set_logfile_loglevel(logger::INFO)
#'
#' @export
#' @seealso \code{\link{omnipath_set_console_loglevel}}
omnipath_set_logfile_loglevel <- function(level){

    omnipath_set_loglevel(level = level, target = 'logfile')

}


#' Dispatch a message to the OmnipathR logger
#'
#' Any package or script can easily send log messages and establish a logging
#' facility with the fantastic `logger` package. This function serves the
#' only purpose if you want to inject messages into the logger of OmnipathR.
#' Otherwise we recommend to use the `logger` package directly.
#'
#' @param level Character, numeric or class loglevel. A log level, if
#'     character one of the followings: "fatal", "error", "warn", "success",
#'     "info", "trace".
#' @param ... Arguments for string formatting, passed \code{sprintf} or
#'     \code{str_glue}.
#'
#' @return Returns `NULL`.
#'
#' @examples
#' omnipath_msg(
#'     level = 'success',
#'     'Talking to you in the name of OmnipathR, my favourite number is %d',
#'     round(runif(1, 1, 10))
#' )
#'
#' @importFrom magrittr %>%
#' @importFrom logger log_level
#' @export
omnipath_msg <- function(level, ...){

    level %>%
    ensure_loglevel %>%
    log_level(...)

}


#' Evaluates an expression with colorout off
#'
#' @param ex An expression.
#'
#' @return The value from the evaluated expression.
#'
#' @importFrom utils flush.console
#' @noRd
no_colorout <- function(ex){

    # to satisfy R CMD check:
    isColorOut <- noColorOut <- ColorOut <- NULL

    colorout_active <- 'colorout' %in% .packages() && isColorOut()

    if(colorout_active) noColorOut()
    result <- eval(ex)
    flush.console()
    if(colorout_active) ColorOut()

    invisible(result)

}


#' Patches logger::appender_console to not to interfere with colorout
#'
#' @importFrom crayon num_colors
#' @importFrom magrittr %<>%
#' @importFrom logger appender_console
#' @noRd
patch_logger_appender <- function(){

    ns <- loadNamespace('logger')

    original_appender_console <- logger::appender_console

    appender_patched <- function(lines){

        if(crayon::num_colors() > 1L){
            lines[1] %<>% {paste0('\033[0m', .)}
        }

        no_colorout(original_appender_console(lines))

    }

    patch_ns('appender_console', appender_patched, ns)
    patch_ns('appender_stderr', appender_patched, ns)

}


#' Parches logger::deparse_to_one_line
#'
#' @noRd
patch_logger_metavar <- function(){

    ns <- loadNamespace('logger')

    deparse_to_one_line <- function(x){

        return('')

#         gsub(
#             '\\s+',
#             ' ',
#             paste(deparse(x), collapse = ' '),
#             perl = TRUE
#         )

    }

    patch_ns('deparse_to_one_line', deparse_to_one_line, ns)

}


#' Is the logfile enabled?
#'
#' @noRd
omnipath_has_logfile <- function(){

    has_logfile('OmnipathR')

}


#' Is the logfile enabled for a certain package?
#'
#' @param pkg Character: name of a package
#'
#' @importFrom magrittr %>% equals not
#' @importFrom stringr str_to_lower
#' @noRd
has_logfile <- function(pkg = 'OmnipathR'){

    str_to_lower(pkg) %>%
    sprintf('%s.logfile', .) %>%
    getOption %>%
    if_null('') %>%
    str_to_lower %>%
    equals('none') %>%
    magrittr::not()

}


#' Shortcut for lowest log level (useful for development)
#'
#' @noRd
.optrace <- function() {

    omnipath_set_console_loglevel('trace')

}


#' Include a welcome message in the log with the package version
#'
#' @importFrom logger log_info
#' @noRd
log_welcome <- function(pkg) {

    log_info('Welcome to OmnipathR!')
    log_info('OmnipathR version: %s', packageVersion(pkgname))

}


#' Log appender that stores the messages for later use
#'
#' @importFrom magrittr %<>%
#' @noRd
appender_delay <- function(pkg) {

    force(pkg)

    structure(
        function(lines) {

            pkg_env(pkg)$delayed_log %<>% c(lines)

        },
        generator = deparse(match.call()),
        pkg = pkg
    )

}


#' Changes the logfile of a loaded package
#'
#' In `init_logger` we set up two log appanders, one for the  console and one
#' for a file output. In the beginning, we may not know the file yet, hence we
#' temporarily use an `appender_delay` instead. This function should be called
#' once the target file become known. Then we copy the delayed and stored
#' messages into the log file, create an appender for it, and continue logging
#' as normal.
#'
#' @noRd
switch_logfile <- function(pkg, path) {

}


#' Log the curl version
#'
#' @importFrom curl curl_version
#' @importFrom purrr map2_chr
#' @importFrom logger log_info
#' @importFrom magrittr %>%
#' @noRd
log_curl_version <- function() {

    curl_version() %>%
    map2_chr(names(.), ~sprintf('%s: %s', .y, paste(.x, collapse = ', '))) %>%
    paste0(collapse = '; ') %>%
    log_info('CURL: %s', .)

}


#' Log the session info
#'
#' @importFrom sessioninfo session_info
#' @importFrom logger log_info
#' @importFrom magrittr %>% extract2
#' @noRd
log_session_info <- function() {

    session_info() %>%
    extract2('platform') %>%
    compact_repr(limit = 999L, sep = '; ') %>%
    log_info('Session info: %s', .)

    session_info(info = 'external') %>%
    extract2('external') %>%
    compact_repr(limit = 999L, sep = '; ') %>%
    log_info('External libraries: %s', .)

}
