#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Diego Mananes
#                  Alberto Valdeolivas
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


#' Read an SBML file
#'
#' Common function for reading SBML files, used by both Metabolic Atlas
#' integrated models and standard-GEM repositories. Uses SBMLR if available,
#' otherwise falls back to xml2 for basic XML parsing.
#'
#' @param path Character: path to the SBML file (local file or URL).
#' @param return_xml Logical: if TRUE, return an XML document object even if
#'     SBMLR is available.
#'
#' @return An SBML model object (if SBMLR is available) or an XML document
#'     object (from xml2).
#'
#' @importFrom xml2 read_xml
#' @importFrom logger log_warn log_trace
#' @noRd
sbml_read <- function(path, return_xml = FALSE) {

    sbmlr_available <- requireNamespace('SBMLR', quietly = TRUE)

    if (!sbmlr_available && !return_xml) {
        msg <- "SBMLR package not available, install with: BiocManager::install('SBMLR')"
        log_warn(msg)
        warning(msg, call. = FALSE)
    }

    use_sbmlr <- sbmlr_available && !return_xml

    if (use_sbmlr) {
        log_trace('Reading SBML with SBMLR: %s', path)
        # Suppress SBMLR's stdout output
        capture.output(
            result <- SBMLR::readSBML(path),
            type = 'output'
        )
        result
    } else {
        log_trace('Reading SBML with xml2: %s', path)
        xml2::read_xml(path)
    }

}


#' Callback for reading SBML from downloaded content
#'
#' Writes the content to a temp file and reads it with sbml_read.
#' This is needed because SBMLR requires a file path.
#'
#' @param content Character: the SBML content as a string.
#' @param return_xml Logical: passed to sbml_read.
#'
#' @return An SBML model object or XML document.
#'
#' @importFrom logger log_trace
#' @noRd
sbml_reader_callback <- function(return_xml = FALSE) {

    function(content) {
        temp_file <- tempfile(fileext = '.xml')
        on.exit(unlink(temp_file), add = TRUE)
        writeLines(content, temp_file)
        sbml_read(temp_file, return_xml = return_xml)
    }

}


#' Load an SBML model from a standard-GEM repository
#'
#' Downloads and parses the SBML (.xml) model file from a Metabolic Atlas
#' standard-GEM repository. Uses the SBMLR package if available, otherwise
#' falls back to xml2 for basic XML parsing.
#'
#' @param gem Character: name of the GEM (e.g., "Human-GEM") or full repo path.
#' @param ref Character: git reference (tag, branch, or commit SHA). If NULL,
#'     defaults to the latest version.
#' @param return_xml Logical: if TRUE, return an XML document object even if
#'     SBMLR is available. Useful for custom parsing.
#'
#' @return An SBML model object (if SBMLR is available) or an XML document
#'     object (from xml2), or NULL if the file is not found.
#'
#' @examples
#' human_gem <- metatlas_gem_sbml("Human-GEM")
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom logger log_info log_warn
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metatlas_gem_reactions}}}
#'     \item{\code{\link{metatlas_gem_metabolites}}}
#'     \item{\code{\link{metatlas_gem_genes}}}
#' }
metatlas_gem_sbml <- function(gem, ref = NULL, return_xml = FALSE) {

    .slow_doctest()

    # Resolve GEM info
    gem_info <- gem_repo_info(gem)

    if (is.null(gem_info)) {
        return(NULL)
    }

    # Default to latest version if ref not specified
    ref %<>% if_null(gem_info$latest_version)

    # Find the SBML file in the repo tree
    repo_tree <- gem_info$repo_tree[[1L]]
    file_path <- gem_find_file(repo_tree, '\\.xml$', 'model')

    if (is.null(file_path)) {
        # Try compressed version
        file_path <- gem_find_file(repo_tree, '\\.xml\\.gz$', 'model')
    }

    if (is.null(file_path)) {
        msg <- sprintf(
            'SBML file not found in GEM `%s` at ref `%s`.',
            gem_info$model_name, ref
        )
        log_warn(msg)
        warning(msg)
        return(NULL)
    }

    # Use first match if multiple found
    file_path <- file_path[1L]

    log_info(
        'Loading SBML model from %s (ref: %s)',
        gem_info$model_name, ref
    )

    gem_download_file(
        git_host = gem_info$git_host,
        git_repo = gem_info$git_repo,
        ref = ref,
        path = file_path,
        reader = sbml_reader_callback(return_xml = return_xml)
    )

}
