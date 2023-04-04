#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2022
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
#
#  This file is from https://github.com/klmr/decorator
#  Author: Kondrad Rudolph
#  License: Apache-2.0
#

#' @noRd
`%@%` <- function (decorator, f) UseMethod('%@%')


#' @noRd
`%@%.default` <- function (decorator, f)
    stop(deparse(substitute(decorator)), ' is not a decorator')


#' @noRd
`%@%.decorator` <- function (decorator, f) {
    pretty_decorators = as.list(match.call())[-1]
    # Patch delayed decorator.
    if (! is.null({pretty_patched = attr(decorator, 'calls')}))
        pretty_decorators = c(pretty_patched, pretty_decorators[-1])

    # Handle operator precedence so that we can chain decorators.
    if (inherits(f, 'decorator'))
        .delayed_decorate(decorator, f, pretty_decorators)
    else
        prettify(decorator(f), f, pretty_decorators[-length(pretty_decorators)])
}


#' @noRd
decorator <- function (f)
    structure(f, class = 'decorator')


decorator <- decorator(decorator)


#' @importFrom utils capture.output
#' @noRd
print.decorated <- function (x, useSource = TRUE, ...) {

    bare <- function (f) {
        bare = unclass(f)
        attr(bare, 'decorators') = NULL
        bare
    }

    fun_def <- capture.output(
        print.function(
            bare(x),
            useSource = useSource,
            ...
        )
    )
    for (decorator in attr(x, 'decorators'))
        cat(deparse(decorator), '%@%\n')
    cat(fun_def, sep = '\n')
    invisible(x)
}

attr(environment(print.decorated), 'S3') <- c(
    attr(environment(print.decorated), 'S3'),
    'print.decorated'
)
registerS3method(
    'print',
    'decorated',
    print.decorated,
    environment(print.decorated)
)


#' @noRd
prettify <- function (f, original, decorator_calls) {
    attr(f, 'srcref') = pretty_code(original)
    attr(f, 'decorators') = decorator_calls
    class(f) = c(class(f), 'decorated')
    f
}


#' @noRd
pretty_code <- function (f) {
    srcref = attr(f, 'srcref')
    if (is.null(srcref)) body(f) else srcref
}


#' @noRd
.delayed_decorate <- function (d1, d2, decorator_calls)
    structure(decorator(function (f) d1(d2(f))), calls = decorator_calls)
