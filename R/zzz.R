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


#' @importFrom rmarkdown pandoc_version
#' @importFrom logger log_trace log_info
#' @importFrom utils packageVersion
#' @noRd
.onLoad <- function(libname, pkgname){

    # Suppress common "incomplete final line" warnings globally
    .suppress_readlines_warnings()

    omnipath_init_config()
    patch_logger_metavar()
    patch_logger_appender()
    patch_attr_preserving_all()
    patch_httr2_keep_handle()
    omnipath_init_log()

    buildserver <- .on_buildserver()
    omnipathr.env$doctest_bypass <- TRUE

    if(buildserver){

        omnipath_set_console_loglevel('trace')

    }

    omnipath_init_cache()
    log_all_info(pkgname)

    if(buildserver) {

        logger::log_trace('Running on a build server, wiping cache.')
        cachedir <- omnipath_get_cachedir()
        logger::log_trace('Cache is at `%s`.', cachedir)
        logger::log_trace('Contains %i files.', length(list.files(cachedir)))
        logger::log_trace(
            'Cache is locked: %s.',
            'cache.lock' %in% list.files(cachedir)
        )
        tryCatch(
            omnipath_cache_wipe(),
            error = function(e){
                log_error_with_info('Failed to wipe cache: %s.', e)
                logger::log_trace('On a build server, unlocking cache db.')
            }
        )
        omnipath_unlock_cache_db()

        # report pandoc version for Bioc build server debugging
        logger::log_trace('Pandoc version: `%s`.', pandoc_version())

    }

    omnipath_init_db(pkgname)
    .load_magic_bytes(pkgname)
    .load_urls(pkgname)
    .load_id_types(pkgname)
    .load_organisms(pkgname)
    .load_entity_types(pkgname)

    if(buildserver){

        omnipath_unlock_cache_db()

        logger::log_trace(
            'Cache locked: %s',
            file.exists(omnipath_cache_lock_path())
        )

    }

}


#' Suppress common "incomplete final line" warnings globally
#'
#' This function sets up a global warning handler that specifically filters
#' out the harmless "incomplete final line found" warning that occurs when
#' reading files without trailing newlines, while preserving other important
#' warnings.
#'
#' @importFrom logger log_trace
#' @noRd
.suppress_readlines_warnings <- function() {

    # Store original warning handler
    original_warning_handler <- getOption("warning.expression")

    # Set up custom warning handler
    options(warning.expression = expression({
        # Call the original warning handler if it exists
        if (!is.null(original_warning_handler)) {
            eval(original_warning_handler)
        }

        # Filter out specific "incomplete final line" warnings
        # Pattern matches variations like "incomplete final line found on 'file'"
        if (grepl("incomplete final line", conditionMessage(condition), fixed = TRUE)) {
            # Log the suppression at trace level for debugging
            logger::log_trace("Suppressed readLines warning: %s", conditionMessage(condition))
            invokeRestart("muffleWarning")
        }
    }))

    logger::log_trace("Global readLines warning suppression enabled")

}
