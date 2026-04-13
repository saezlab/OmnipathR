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

#' List of original GEMs in Metabolic Atlas
#'
#' Metabolic Atlas is a repository of genome scale models of metabolism (GEMs)
#' for various organisms, tissues and conditions. This function returns a list
#' of the available GEMs. Read more at \url{https://metabolicatlas.org}.
#'
#' @param integrated Logical: list the integrated (standard) GEMs. These are
#'     available by a separate API in Metabolic Atlas. There are 7 integrated
#'     GEMs available by the API, and 23 standard GEMs listed in Metabolic
#'     Atlas in total. In the repository a total of 360 models are available.
#'
#' @return A data frame (tibble) of GEMs.
#'
#' @examples
#' \dontrun{
#' metabolic_atlas_models()
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom tidyr unnest_wider
#' @importFrom dplyr rename mutate select coalesce
#' @export
metabolic_atlas_list_models <- function(integrated = FALSE) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    gemodelset <- sample <- description <- description2 <- NULL

    c('metatlas', `if`(integrated, 'integrated', NULL), 'models') %>%
    paste0(collapse = '_') %>%
    download_to_cache() %>%
    safe_json() %>%
    tibble() %>%
    {`if`(
        integrated,
        unnest_wider(., sample),
        rename(., description2 = description) %>%
        unnest_wider(c(gemodelset, sample)) %>%
        mutate(description = coalesce(description, description2)) %>%
        select(-description2)
    )}

}


#' Download and load one or more SBML models from Metabolic Atlas
#'
#' Downloads a genome-scale metabolic model (GEM) from Metabolic Atlas and
#' loads it using the SBMLR package (optional dependency).
#'
#' @param ... Arguments to filter models. Either:
#'     - A data frame of filtered models from metabolic_atlas_models()
#'     - Integer vector: model ID(s) from the `id` column
#'     - Expressions passed to dplyr::filter() (e.g., organism == 'Homo sapiens')
#' @param return_xml Logical: return an XML document object even if the SBMLR
#'     package is installed.
#'
#' @return A model object loaded by SBMLR, or an XML document object
#'     from xml2 if SBMLR is not available.
#'
#' @examples
#' \dontrun{
#' # Load model by ID
#' model <- metabolic_atlas_models(1)
#'
#' # Load first human model
#' model <- metabolic_atlas_models(tissue = "Cervix")
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom rlang eval_tidy enexprs
#' @importFrom dplyr filter slice pull
#' @importFrom purrr map
#' @importFrom logger log_info log_error log_warn
#' @export
metabolic_atlas_models <- function(..., return_xml = FALSE) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    name <- files <- path <- id <- NULL

    sbmlr_available <- requireNamespace("SBMLR", quietly = TRUE)

    if (!sbmlr_available && !return_xml) {
        msg <- "SBMLR package not available, install with: BiocManager::install('SBMLR')"
        log_warn(msg)
        warning(msg, call. = FALSE)
    }

    # Capture expressions first to prevent immediate evaluation
    filter_exprs <- enexprs(...)

    # Handle different argument types
    models <- if (length(filter_exprs) == 1L && is.data.frame(filter_exprs[[1L]])) {
        filter_exprs[[1L]]
    } else if (length(filter_exprs) == 1L && is.numeric(filter_exprs[[1L]])) {
        metabolic_atlas_list_models() %>% filter(id %in% filter_exprs[[1L]])
    } else {
        models <- metabolic_atlas_list_models()
        for (expr in filter_exprs) {
            models <- models %>% filter(eval_tidy(expr, data = models))
        }
        models
    }

    if (nrow(models) == 0L) {
        msg <- sprintf("No models found matching the criteria")
        log_error(msg)
        stop(msg, call. = FALSE)
    }

    models %>%
    {map(
        seq_len(nrow(.)),
        ~metabolic_atlas_model(slice(models, .x), return_xml = return_xml)
    )} %>%
    {`if`(length(.) == 1L, .[[1L]], .)}

}


#' List standard GEMs from Metabolic Atlas
#'
#' Retrieves information about standard GEMs from Metabolic Atlas.
#'
#' @return A data frame (tibble) with standard GEMs metadata and release information.
#'     Each model may have multiple rows, one for each release. The `latest` column
#'     indicates the most recent release for each model. The `repo_tree` column
#'     contains a character vector of file paths in the repository at the latest
#'     version, or NULL if the `gert` or `git2r` packages are not available.
#'
#' @examples
#' \dontrun{
#' metabolic_atlas_list_gems()
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom purrr map pmap
#' @importFrom tidyr unnest
#' @export
metabolic_atlas_list_gems <- function() {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    git_host <- git_repo <- gem_info <- latest_version <- NULL

    download_to_cache("metatlas_standard_gems_index") %>%
    safe_json() %>%
    tibble(git_host = names(.), git_repo = .) %>%
    unnest(git_repo) %>%
    mutate(
        gem_info = map(git_repo, metabolic_atlas_gem_info),
        git_url = sprintf('https://%s.com/%s', git_host, git_repo)
    ) %>%
    unnest(gem_info) %>%
    mutate(
        repo_tree = pmap(
            list(git_host, git_repo, latest_version),
            .git_repo_tree
        )
    )

}


#' Download and parse JSON metadata for one standard GEM from Metabolic Atlas
#'
#' @importFrom magrittr %>% extract extract2
#' @importFrom dplyr bind_rows mutate first
#' @importFrom purrr map map_chr
#' @importFrom stringr str_split
#' @importFrom logger log_trace
#' @importFrom tibble tibble
#' @importFrom rlang %||% set_names
#' @noRd
metabolic_atlas_gem_info <- function(repo_name) {

    gem_name <- repo_name %>% str_split('/') %>% unlist %>% last

    log_trace('Loading metadata about `%s` GEM.', gem_name)

    gem_json <-
        download_to_cache(
            "metatlas_standard_gem",
            url_param = list(gem_name)
        ) %>%
        safe_json(simplifyDataFrame = FALSE) %>%
        extract2(gem_name)

    metadata <- gem_json$metadata

    releases <-
        gem_json$releases %>%
        set_names(map_chr(., ~names(.x) %>% extract(1L))) %>%
        map(
            function(release) {
                release %>%
                extract2(1L) %>%
                extract2('standard-GEM') %>%
                extract2(2L) %>%
                extract2('test_results')
            }
        )

    tibble(
        model_name = gem_name,
        releases = list(releases),
        commits = metadata$commits %||% NA_integer_,
        contributors = metadata$contributors %||% NA_integer_,
        latest_commit_date = metadata$latest_commit_date %||% NA_character_,
        owner = metadata$owner %||% NA_character_,
        avatar = metadata$avatar %||% NA_character_
    ) %>%
    mutate(
        latest_version = map_chr(releases, ~names(.x) %>% first)
    )

}


#' Get the directory tree of a git repository at a specific reference
#'
#' Uses the GitHub or GitLab API to retrieve the repository tree at a given
#' reference.
#'
#' @param git_host Character: git host ("github" or "gitlab").
#' @param git_repo Character: repository path in "owner/repo" format.
#' @param ref Character: git reference (tag, branch, or commit SHA).
#'
#' @return A character vector of file paths in the repository, or NULL if
#'     an error occurred.
#'
#' @importFrom logger log_trace log_warn
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @noRd
.git_repo_tree <- function(git_host, git_repo, ref) {

    # Handle missing reference
    if (is.na(ref) || is.null(ref) || ref == '') {
        return(NULL)
    }

    log_trace('Getting repository tree for %s/%s at ref %s', git_host, git_repo, ref)

    tryCatch({

        if (git_host == 'github') {

            url_key <- 'github_repo_tree'
            url_param <- list(git_repo, ref)
            http_headers <- list(
                Accept = 'application/vnd.github.v3+json',
                `User-Agent` = 'OmnipathR'
            )

        } else if (git_host == 'gitlab') {

            url_key <- 'gitlab_repo_tree'
            # GitLab requires URL-encoded project path
            encoded_repo <- URLencode(git_repo, reserved = TRUE)
            url_param <- list(encoded_repo, ref)
            http_headers <- list(`User-Agent` = 'OmnipathR')

        } else {

            log_warn('Unsupported git host: %s', git_host)
            return(NULL)

        }

        tree_data <- generic_downloader(
            url_key = url_key,
            reader = jsonlite::fromJSON,
            url_param = url_param,
            reader_param = list(simplifyDataFrame = TRUE),
            http_headers = http_headers,
            use_httr = TRUE
        )

        # Extract paths - GitHub uses tree$path, GitLab returns path directly
        if (git_host == 'github') {
            if (!is.null(tree_data$tree) && nrow(tree_data$tree) > 0) {
                tree_data$tree$path
            } else {
                NULL
            }
        } else {
            if (is.data.frame(tree_data) && nrow(tree_data) > 0) {
                tree_data$path
            } else {
                NULL
            }
        }

    }, error = function(e) {
        log_warn('Failed to get repository tree for %s/%s: %s', git_host, git_repo, e$message)
        NULL
    })

}


#' Download and load an SBML model from Metabolic Atlas
#'
#' @importFrom xml2 read_xml
#' @noRd
metabolic_atlas_model <- function(model_row, return_xml = FALSE) {

    # NSE vs. R CMD check workaround
    id <- organism <- tissue <- cell_type <- condition <-
        reaction_count <- metabolite_count <- gene_count <- year <-
        format <- NULL

    use_sbmlr <- requireNamespace('SBMLR', quietly = TRUE) & !return_xml

    sbml_path <-
        model_row$files[[1L]] %>%
        filter(format == 'SBML') %>%
        slice(1L) %>%
        pull(path)

    if (length(sbml_path) == 0L) {
        log_warn(
            'No SBML file available for model `%s` (id=%s).',
            model_row$name,
            model_row$id
        )
        return(NULL)
    }

    log_info(
        'Downloading model: %s %s',
        model_row$name,
        model_row %>%
            select(
                id, organism, tissue, cell_type, condition,
                reaction_count, metabolite_count, gene_count, year
            ) %>%
            as.list %>%
            compact_repr
    )

    readSBML <- SBMLR %:::% readSBML

    silent_sbml <- function(x) {
        capture.output(
            result <- readSBML(x),
            type = 'output'
        )
        result
    }

    method <- archive_extractor
    dl_args <- list(
        url_key = 'metatlas_model',
        url_param = list(sbml_path),
        reader = `if`(use_sbmlr, silent_sbml, xml2::read_xml),
        to_tempdir = use_sbmlr
    )

    exec(method, !!!dl_args)

}
