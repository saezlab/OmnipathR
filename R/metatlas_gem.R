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


#' Resolve GEM repository information
#'
#' Looks up a GEM by name in the Metabolic Atlas standard-GEM list and returns
#' repository metadata including host, path, latest version, and file tree.
#'
#' @param gem Character: name of the GEM (e.g., "Human-GEM") or full repo path
#'     (e.g., "SysBioChalmers/Human-GEM").
#'
#' @return A single-row tibble with columns: git_host, git_repo, model_name,
#'     latest_version, repo_tree, and other metadata from metabolic_atlas_list_gems().
#'     Returns NULL if the GEM is not found.
#'
#' @examples
#' \dontrun{
#' info <- gem_repo_info("Human-GEM")
#' info$latest_version
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom stringr str_detect fixed
#' @importFrom logger log_trace log_warn
#' @noRd
gem_repo_info <- function(gem) {

    # NSE vs. R CMD check workaround
    model_name <- git_repo <- NULL

    log_trace('Looking up GEM repository info for `%s`.', gem)

    gems <- get_db('metatlas_gems')

    # Try matching by model_name first, then by git_repo

    result <- gems %>% filter(model_name == gem)

    if (nrow(result) == 0L) {
        result <- gems %>% filter(git_repo == gem)
    }

    # Also try partial match on git_repo (e.g., "Human-GEM" matches "SysBioChalmers/Human-GEM")
    if (nrow(result) == 0L) {
        result <- gems %>% filter(str_detect(git_repo, fixed(gem)))
    }

    if (nrow(result) == 0L) {
        log_warn('GEM `%s` not found in Metabolic Atlas standard-GEM list.', gem)
        return(NULL)
    }

    if (nrow(result) > 1L) {
        log_warn(
            'Multiple GEMs match `%s`, using first match: %s',
            gem, result$model_name[1L]
        )
        result <- result[1L, ]
    }

    log_trace(
        'Found GEM `%s`: repo=%s, latest_version=%s',
        result$model_name, result$git_repo, result$latest_version
    )

    result

}


#' Download a raw file from a GEM repository
#'
#' Downloads a file from a GitHub or GitLab repository at a specific git
#' reference (tag, branch, or commit).
#'
#' @param git_host Character: git host ("github" or "gitlab").
#' @param git_repo Character: repository path in "owner/repo" format.
#' @param ref Character: git reference (tag, branch, or commit SHA).
#' @param path Character: path to the file within the repository.
#' @param reader Function: reader function to parse the downloaded content.
#'     Defaults to identity (returns raw text).
#' @param reader_param List: additional parameters passed to the reader function.
#' @param ... Additional arguments passed to generic_downloader.
#'
#' @return The file content, processed by the reader function.
#'
#' @examples
#' \dontrun{
#' # Download a TSV file
#' reactions <- gem_download_file(
#'     "github", "SysBioChalmers/Human-GEM", "v1.17.0",
#'     "model/reactions.tsv",
#'     reader = readr::read_tsv
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom logger log_trace log_warn
#' @importFrom utils URLencode
#' @noRd
gem_download_file <- function(
        git_host,
        git_repo,
        ref,
        path,
        reader = identity,
        reader_param = list(),
        ...
    ) {

    log_trace(
        'Downloading file from %s/%s at ref %s: %s',
        git_host, git_repo, ref, path
    )

    if (git_host == 'github') {

        url_key <- 'github_raw_file'
        url_param <- list(git_repo, ref, path)

    } else if (git_host == 'gitlab') {

        url_key <- 'gitlab_raw_file'
        # GitLab uses URL-encoded repo path in raw file URLs
        encoded_repo <- URLencode(git_repo, reserved = TRUE)
        url_param <- list(git_repo, ref, path)

    } else {

        msg <- sprintf('Unsupported git host: %s', git_host)
        log_warn(msg)
        warning(msg)
        return(NULL)

    }

    generic_downloader(
        url_key = url_key,
        reader = reader,
        url_param = url_param,
        reader_param = reader_param,
        use_httr = TRUE,
        ...
    )

}


#' Find a file in the GEM repository tree
#'
#' Searches the repository file tree for files matching a pattern.
#'
#' @param repo_tree Character vector: list of file paths in the repository.
#' @param pattern Character: regex pattern to match file paths.
#' @param directory Character: optional directory prefix to filter by.
#'
#' @return Character vector of matching file paths, or NULL if none found.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr keep
#' @importFrom stringr str_detect
#' @importFrom logger log_trace
#' @noRd
gem_find_file <- function(repo_tree, pattern, directory = NULL) {

    if (is.null(repo_tree)) {
        return(NULL)
    }

    matches <-
        repo_tree %>%
        {`if`(
            is.null(directory),
            .,
            keep(., ~str_detect(.x, sprintf('^%s', directory)))
        )} %>%
        keep(~str_detect(.x, pattern))

    if (length(matches) == 0L) {
        log_trace('No files matching pattern `%s` in repo tree.', pattern)
        return(NULL)
    }

    log_trace('Found %d files matching pattern `%s`.', length(matches), pattern)
    matches

}
