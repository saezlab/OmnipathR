#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
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


#' Keep the curl handle in httr2 under the `handle` attribute of `request`
#'
#' @importFrom rlang !!! fn_fmls_syms new_function expr current_env
#' @noRd
patch_httr2_keep_handle <- function() {

    PATCHED <- 'omnipathr_patched'

    if (requireNamespace('httr2', quietly = TRUE)) {

        ns <- asNamespace('httr2')
        original <- get('req_perform1', ns)

        if (is.null(attr(original, PATCHED))) {

            patched <- new_function(
                args = formals(original),
                body = expr({
                    if (!is.null(handle)) {
                        attr(req, 'handle') <- handle
                    }
                    ORIGINAL(!!!fn_fmls_syms())
                }),
                env = current_env()
            )

            environment(patched)$ORIGINAL <- original

            attr(patched, PATCHED) <- TRUE

            patch_ns('req_perform1', patched, ns)

        } else {

            patched <- original

        }

        patch_ns('req_perform1', patched, asNamespace('OmnipathR'), add = TRUE)

    }

}
