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


#' Unlocks a binding in a namespace, changes its value and locks the
#' binding again
#'
#' @noRd
patch_ns <- function(name, patched, ns, add = FALSE) {

    ulb <- get('unlockBinding')
    lb <- get('lockBinding')

    if (add || name %in% names(ns)) {
        if (name %in% names(ns)) {
            ulb(name, as.environment(ns))
        }
        assign(name, patched, ns)
        lb(name, as.environment(ns))
    }

}
