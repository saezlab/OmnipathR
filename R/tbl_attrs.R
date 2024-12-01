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


#' @importFrom magrittr %>% %<>% extract
#' @noRd
copy_attrs <- function(to, from, .ignore = NULL) {

    attrs <-
        attributes(from) %>%
        extract(setdiff(names(.), .ignore)) %>%
        if_null(list())


    attributes(to) %<>% if_null(list()) %>% modifyList(attrs, .)

    return(to)

}


#' A patched data.frame (tibble) class
#'
#' For attribute preserving patching and dispatch of S3 methods
#'
#' @importFrom magrittr %<>%
#' @noRd
tbl_attrs <- function(data, original = NULL, .ignore = NULL) {
    data %<>% copy_attrs(original, .ignore = .ignore)
    class(data) %<>% c('tbl_attrs', .)
    return(data)
}


#' @importFrom magrittr %>%
#' @export
group_by.tbl_attrs <- function(.data, ...) {
    NextMethod() %>% tbl_attrs(.data)
}


#' @importFrom magrittr %>%
#' @export
ungroup.tbl_attrs <- function(.data, ...) {
    NextMethod() %>% tbl_attrs(.data, .ignore = 'groups')
}


#' @importFrom magrittr %>%
#' @export
summarise.tbl_attrs <- function(.data, ...) {
    NextMethod() %>% tbl_attrs(.data, .ignore = 'groups')
}


#' @importFrom magrittr %>%
#' @export
mutate.tbl_attrs <- function(.data, ...) {
    NextMethod() %>% tbl_attrs(.data)
}


#' @importFrom magrittr %>%
#' @export
rowwise.tbl_attrs <- function(data, ...) {
    class(data) %<>% setdiff('tbl_attrs')
    NextMethod() %>% tbl_attrs(data)
}


#' @importFrom magrittr %>%
#' @export
reframe.tbl_attrs <- function(.data, ...) {
    class(.data) %<>% setdiff('tbl_attrs')
    NextMethod() %>% tbl_attrs(.data, .ignore = 'groups')
}


ATTR_PRESERVING = list(
    list(pkg = 'tidyr', fn = 'unnest_longer'),
    list(pkg = 'tidyr', fn = 'unnest_wider')
)


#' Patch functions to preserve attributes
#'
#' As we see above most of the dplyr verbs which do not presesve attributes
#' by default can be patched using the custom tbl_attrs tibble subclass and S3
#' method dispatch. However, some other verbs are not S3 methods but simple
#' functions, and as such, we can make them attribute preserving by patching
#' them in the namespaces of the providing packages.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr walk
#' @noRd
patch_attr_preserving_all <- function() {

    ATTR_PRESERVING %>%
    walk(patch_attr_preserving)

}


#' @importFrom magrittr %<>% extract2
#' @importFrom rlang list2 fn_fmls_names
#' @noRd
patch_attr_preserving <- function(f) {

    PATCHED <- 'omnipathr_patched'

    if (requireNamespace(f$pkg, quietly = TRUE)) {

        ns <- asNamespace(f$pkg)
        original <- get(f$fn, ns)
        tbl_name <- original %>% fn_fmls_names %>% extract(1L)

        if (is.null(attr(original, PATCHED))) {

            patched <- function(data, ...) {
                call <- match.call()
                call[[1L]] <- original
                obj <- eval(call[[tbl_name]], parent.frame())
                result <- eval(call, parent.frame())
                result %<>% copy_attrs(obj, .ignore = 'groups')
                if (inherits(obj, 'tbl_attrs')) {
                    class(result) %<>% c('tbl_attrs', .)
                }
                result
            }

            formals(patched) <- formals(original)
            attr(patched, PATCHED) <- TRUE

            # patch_ns(f$fn, patched, ns)

        } else {

            patched <- original

        }

        patch_ns(f$fn, patched, asNamespace('OmnipathR'), add = TRUE)

    }

}
