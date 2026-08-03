#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2026
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): OmniPath Team (omnipathdb@gmail.com)
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


# Client for the `omnipath-metabo` web service (COSMOS PKN endpoints).
# See `R/metabo_networks.R` for the generic `/networks` endpoints, which
# share the internal engine defined below.


#' Base URL of the omnipath-metabo web service
#'
#' @param path Character: path to append to the omnipath-metabo base URL
#'     (\code{getOption('omnipathr.metabo_url')}).
#'
#' @return Character: the base URL joined with \code{path}.
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace
#'
#' @noRd
.metabo_url <- function(path){

    'omnipathr.metabo_url' %>%
    getOption %>%
    str_replace('/+$', '') %>%
    c(str_replace(path, '^/+', '')) %>%
    paste(collapse = '/')

}


#' Builds a query URL against the omnipath-metabo web service
#'
#' @param path Character: the endpoint path, e.g. \code{'cosmos/pkn'}.
#' @param ... Named query string parameters. \code{NULL} or \code{NA}
#'     values are omitted; vectors are collapsed with commas; mirrors
#'     \code{omnipath_url_add_param}, the same helper the main OmniPath
#'     query engine uses.
#'
#' @return Character: the complete, URL-encoded request URL.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr reduce
#' @importFrom utils URLencode
#'
#' @noRd
.metabo_build_url <- function(path, ...){

    params <- list(...)

    params %>%
    names %>%
    reduce(
        function(url, key){
            omnipath_url_add_param(url, key, params[[key]])
        },
        .init = .metabo_url(path)
    ) %>%
    URLencode()

}


#' Queries an endpoint of the omnipath-metabo web service
#'
#' Thin wrapper around \code{\link{generic_downloader}}: caching, retry and
#' logging all come from the existing OmnipathR download infrastructure, no
#' new machinery is introduced here.
#'
#' @param path Character: the endpoint path, e.g. \code{'cosmos/pkn'}.
#' @param ... Named query string parameters, passed to
#'     \code{.metabo_build_url}.
#' @param reader_param List: passed to \code{jsonlite::fromJSON} via
#'     \code{curl_read_json}.
#'
#' @return The parsed JSON response.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
.metabo_query <- function(path, ..., reader_param = list()){

    generic_downloader(
        url_key = .metabo_build_url(path, ...),
        reader = curl_read_json,
        reader_param = reader_param,
        resource = 'omnipath-metabo'
    )

}


#' Detects a 404 (not found) condition raised by the download layer
#'
#' The download machinery (\code{download_base}, via the \code{curl}
#' package) raises an R error for any non-2xx HTTP status before this
#' package ever sees a response body, so a 404 from the omnipath-metabo
#' service cannot be distinguished by inspecting parsed JSON -- it must be
#' caught here instead.
#'
#' @param cond A condition caught from a failed \code{.metabo_query} call.
#'
#' @return Logical.
#'
#' @noRd
.metabo_is_not_found <- function(cond){

    grepl('error: 404', conditionMessage(cond), fixed = TRUE)

}


#' Empty COSMOS PKN interaction table
#'
#' @return An empty tibble with the columns \code{\link{metabo_cosmos_pkn}}
#'     returns.
#'
#' @importFrom tibble tibble
#'
#' @noRd
.metabo_empty_pkn <- function(){

    tibble(
        source = character(),
        target = character(),
        source_type = character(),
        target_type = character(),
        id_type_a = character(),
        id_type_b = character(),
        interaction_type = character(),
        resource = character(),
        mor = integer(),
        locations = list(),
        attrs = list()
    )

}


#' COSMOS PKN categories available from omnipath-metabo
#'
#' Lists the PKN subset ("category") names currently served by the
#' omnipath-metabo web service (\url{https://metabo.omnipathdb.org}), e.g.
#' \code{'transporters'}, \code{'receptors'}. Use these names with the
#' \code{categories} argument of \code{\link{metabo_cosmos_pkn}}. The list
#' reflects whatever the service currently supports, not a fixed list
#' bundled with the package.
#'
#' @return Character vector of category names.
#'
#' @examples
#' \dontrun{
#' metabo_cosmos_categories()
#' }
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_cosmos_pkn}}}
#'     \item{\code{\link{metabo_cosmos_resources}}}
#'     \item{\code{\link{metabo_cosmos_organisms}}}
#' }
metabo_cosmos_categories <- function(){

    .slow_doctest()

    'cosmos/categories' %>%
    .metabo_query() %>%
    as.character

}


#' Organisms with a COSMOS PKN available from omnipath-metabo
#'
#' Lists the NCBI Taxonomy IDs currently pre-built and cached by the
#' omnipath-metabo web service.
#'
#' @return Integer vector of NCBI Taxonomy IDs.
#'
#' @examples
#' \dontrun{
#' metabo_cosmos_organisms()
#' }
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_cosmos_pkn}}}
#'     \item{\code{\link{metabo_cosmos_categories}}}
#' }
metabo_cosmos_organisms <- function(){

    .slow_doctest()

    'cosmos/organisms' %>%
    .metabo_query() %>%
    as.integer

}


#' Resources contributing to each COSMOS PKN category
#'
#' Lists, for each PKN category, the contributing resource names currently
#' available from the omnipath-metabo web service. Use the \code{resource}
#' column values with the \code{resources} argument of
#' \code{\link{metabo_cosmos_pkn}}.
#'
#' @return A tibble with columns \code{category} and \code{resource}.
#'
#' @examples
#' \dontrun{
#' metabo_cosmos_resources()
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom purrr imap_dfr
#' @importFrom tibble tibble
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_cosmos_pkn}}}
#'     \item{\code{\link{metabo_cosmos_categories}}}
#' }
metabo_cosmos_resources <- function(){

    .slow_doctest()

    'cosmos/resources' %>%
    .metabo_query() %>%
    imap_dfr(~ tibble(category = .y, resource = as.character(.x)))

}


#' Cache status of the omnipath-metabo COSMOS PKN build
#'
#' Reports which category/organism combinations are currently pre-built and
#' cached by the omnipath-metabo web service, and their sizes.
#'
#' @return A tibble with columns \code{category}, \code{organism},
#'     \code{size_mb} and \code{modified}, with the overall
#'     \code{total_size_mb} and \code{n_components} attached as attributes.
#'
#' @examples
#' \dontrun{
#' status <- metabo_cosmos_status()
#' attr(status, 'total_size_mb')
#' }
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr mutate
#' @importFrom tibble tibble as_tibble
#' @importFrom rlang .data
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_cosmos_pkn}}}
#'     \item{\code{\link{metabo_cosmos_categories}}}
#' }
metabo_cosmos_status <- function(){

    .slow_doctest()

    data <- .metabo_query('cosmos/status')

    components <- `if`(
        is.null(data$components) || length(data$components) == 0L,
        tibble(
            category = character(),
            organism = integer(),
            size_mb = double(),
            modified = double()
        ),
        data$components %>%
        as_tibble %>%
        mutate(organism = as.integer(.data$organism))
    )

    attr(components, 'total_size_mb') <- data$total_size_mb
    attr(components, 'n_components') <- data$n_components

    components

}


#' COSMOS prior-knowledge network from the omnipath-metabo web service
#'
#' Retrieves the COSMOS PKN, or a chosen subset of it, for one organism
#' from the pre-built, pre-cached omnipath-metabo web service
#' (\url{https://metabo.omnipathdb.org}). Unlike \code{\link{cosmos_pkn}},
#' which builds a PKN locally from STITCH, the Chalmers GEM and OmniPath,
#' this function is a thin, cached client for the separately maintained
#' Python \code{omnipath-metabo} package -- the two are unrelated
#' implementations and neither replaces the other.
#'
#' @param organism Character or integer: name, common name, latin name or
#'     NCBI Taxonomy ID of the organism.
#' @param categories Character vector: one or more PKN category names (see
#'     \code{\link{metabo_cosmos_categories}}), or \code{'all'} (default)
#'     for the complete PKN.
#' @param resources Character vector: restrict the result to interactions
#'     contributed by these resources (see
#'     \code{\link{metabo_cosmos_resources}}). By default all resources
#'     are included.
#' @param ... Ignored; reserved for future extensions of the web service.
#'
#' @return A tibble of PKN interactions with columns \code{source},
#'     \code{target}, \code{source_type}, \code{target_type},
#'     \code{id_type_a}, \code{id_type_b}, \code{interaction_type},
#'     \code{resource}, \code{mor}, \code{locations} and \code{attrs} (a
#'     nested list column preserving resource-specific metadata, such as
#'     \code{reaction_id}, for future use). Empty (zero-row, same columns)
#'     if the service has no data for the requested organism/category
#'     combination.
#'
#' @details The output is always a binary (flat source-target edge list)
#'     representation. A reaction-grouped ("hypergraph") representation,
#'     collapsing rows that share the same \code{attrs$reaction_id} into
#'     one multi-substrate/product hyperedge, is not yet implemented -- the
#'     data needed for it is already present in the \code{attrs} column,
#'     so this can be added later without a web service change.
#'
#' @examples
#' \dontrun{
#' transporters <- metabo_cosmos_pkn(
#'     organism = 9606,
#'     categories = 'transporters'
#' )
#' transporters
#' }
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble as_tibble
#' @importFrom logger log_warn
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_cosmos_categories}}}
#'     \item{\code{\link{metabo_cosmos_resources}}}
#'     \item{\code{\link{metabo_metabolic_signaling_pkn}}}
#'     \item{\code{\link{metabo_metabolite_protein_pkn}}}
#'     \item{\code{\link{cosmos_pkn}}}
#' }
metabo_cosmos_pkn <- function(
    organism = 9606,
    categories = 'all',
    resources = NULL,
    ...
){

    .slow_doctest()

    organism %<>% ncbi_taxid()

    data <- tryCatch(
        .metabo_query(
            'cosmos/pkn',
            organism = organism,
            categories = categories,
            resources = resources,
            format = 'json'
        ),
        error = function(e){

            if (!.metabo_is_not_found(e)) stop(e)

            log_warn(
                'omnipath-metabo: no cached PKN for organism = %s, categories = %s.',
                paste(organism, collapse = ','),
                paste(categories, collapse = ',')
            )
            NULL

        }
    )

    if (is.null(data) || !is.null(data$error)) {
        return(.metabo_empty_pkn())
    }

    if (is.null(data$network) || length(data$network) == 0L) {
        return(.metabo_empty_pkn())
    }

    as_tibble(data$network)

}


#' Metabolic-centric signalling PKN from omnipath-metabo
#'
#' Convenience wrapper around \code{\link{metabo_cosmos_pkn}} that
#' retrieves the complete COSMOS PKN (all categories) for one organism:
#' transport, enzyme catalysis, receptor and allosteric regulation, protein
#' signalling and gene regulation combined into one network connecting
#' metabolism to cell signalling.
#'
#' @param organism Character or integer: name, common name, latin name or
#'     NCBI Taxonomy ID of the organism.
#' @param ... Passed to \code{\link{metabo_cosmos_pkn}}.
#'
#' @return A tibble of PKN interactions; see \code{\link{metabo_cosmos_pkn}}.
#'
#' @examples
#' \dontrun{
#' pkn <- metabo_metabolic_signaling_pkn(organism = 9606)
#' pkn
#' }
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_cosmos_pkn}}}
#'     \item{\code{\link{metabo_metabolite_protein_pkn}}}
#' }
metabo_metabolic_signaling_pkn <- function(organism = 9606, ...){

    .slow_doctest()

    metabo_cosmos_pkn(organism = organism, categories = 'all', ...)

}


#' Metabolite-protein interaction PKN from omnipath-metabo
#'
#' Convenience wrapper around \code{\link{metabo_cosmos_pkn}} that
#' retrieves only the metabolite-protein interaction layers of the COSMOS
#' PKN (receptor and allosteric-regulation categories) for one organism,
#' excluding transport, metabolic-reaction and pure protein-protein
#' signalling edges. Comparable in scope to a metabolite-protein
#' interaction resource such as MetalinksDB, but derived from the COSMOS
#' PKN build rather than being MetalinksDB itself. For the
#' separately-curated MetalinksDB network, see
#' \code{\link{metabo_network_interactions}}.
#'
#' @param organism Character or integer: name, common name, latin name or
#'     NCBI Taxonomy ID of the organism.
#' @param ... Passed to \code{\link{metabo_cosmos_pkn}}.
#'
#' @return A tibble of PKN interactions; see \code{\link{metabo_cosmos_pkn}}.
#'
#' @examples
#' \dontrun{
#' pkn <- metabo_metabolite_protein_pkn(organism = 9606)
#' pkn
#' }
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_cosmos_pkn}}}
#'     \item{\code{\link{metabo_metabolic_signaling_pkn}}}
#'     \item{\code{\link{metabo_network_interactions}}}
#' }
metabo_metabolite_protein_pkn <- function(organism = 9606, ...){

    .slow_doctest()

    metabo_cosmos_pkn(
        organism = organism,
        categories = c('receptors', 'allosteric'),
        ...
    )

}
