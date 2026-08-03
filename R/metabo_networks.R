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


# Client for the generic `/networks` API of the `omnipath-metabo` web
# service: any specialized network registered on the server (e.g.
# MetalinksDB, LIANA) is reachable here by name, with no package update
# needed when new networks are added server-side. Shares the internal
# engine (`.metabo_url`, `.metabo_build_url`, `.metabo_query`) defined in
# `R/metabo.R`.


#' @noRd
.metabo_assert_network_name <- function(name){

    if (!is.character(name) || length(name) != 1L || is.na(name)) {

        msg <- '`name` must be a single, non-missing character string.'
        log_error(msg)
        stop(msg)

    }

    invisible(name)

}


#' @noRd
.metabo_assert_count <- function(x, arg_name){

    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0L) {

        msg <- sprintf('`%s` must be a single non-negative number.', arg_name)
        log_error(msg)
        stop(msg)

    }

    invisible(x)

}


#' Raises a clear error if a `/networks` response is an error body
#'
#' The server responds to an unknown network name with a JSON body
#' containing a `detail` field (e.g. `"Unknown network: foo"`); this
#' surfaces that as an R error pointing at \code{\link{metabo_networks}}
#' instead of returning the error body as if it were data.
#'
#' @importFrom logger log_error
#'
#' @noRd
.metabo_stop_if_error <- function(data){

    if (!is.null(data$detail)) {

        msg <- sprintf(
            '%s (use `metabo_networks()` to list valid network names).',
            data$detail
        )
        log_error(msg)
        stop(msg)

    }

    invisible(data)

}


#' Queries a `/networks/{name}/...` endpoint, translating a 404 into a
#' clear error
#'
#' The download machinery raises an R error for a 404 HTTP status before
#' this package ever sees a response body (see \code{.metabo_is_not_found}
#' in \code{R/metabo.R}), which is how the server responds to an unknown
#' network name -- so that condition is caught and re-raised here with a
#' message naming the problem and pointing at
#' \code{\link{metabo_networks}}, per the package's error-handling
#' convention (\code{log_error()} then \code{stop()}). Any other failure
#' (e.g. the service being unreachable) is left untouched.
#'
#' @importFrom logger log_error
#'
#' @noRd
.metabo_query_network <- function(path, name, ...){

    data <- tryCatch(
        .metabo_query(path, ...),
        error = function(e){

            if (!.metabo_is_not_found(e)) stop(e)

            msg <- sprintf(
                "Unknown network: '%s' (use `metabo_networks()` to list valid names).",
                name
            )
            log_error(msg)
            stop(msg)

        }
    )

    .metabo_stop_if_error(data)

}


#' @importFrom tibble as_tibble tibble
#' @noRd
.metabo_rows_to_tibble <- function(rows){

    `if`(
        is.null(rows) || length(rows) == 0L,
        tibble(),
        as_tibble(rows)
    )

}


#' Specialized networks registered on the omnipath-metabo web service
#'
#' Lists every specialized network currently registered on the
#' omnipath-metabo web service (\url{https://metabo.omnipathdb.org}) --
#' independently curated networks outside the COSMOS PKN family, such as
#' MetalinksDB (metabolite-protein interactions) or LIANA
#' (ligand-receptor). New networks appear here as soon as they are
#' registered server-side, without a package update. Pass a \code{name}
#' from the result to \code{\link{metabo_network_interactions}},
#' \code{\link{metabo_network_status}} or
#' \code{\link{metabo_network_resources}}.
#'
#' @return A tibble with columns \code{name}, \code{kind}, \code{schema_name},
#'     \code{combined_relation}, \code{included_sources} (list column) and
#'     \code{built_at}.
#'
#' @examples
#' \dontrun{
#' metabo_networks()
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_network_interactions}}}
#'     \item{\code{\link{metabo_network_status}}}
#'     \item{\code{\link{metabo_network_resources}}}
#' }
metabo_networks <- function(){

    .slow_doctest()

    'networks/' %>%
    .metabo_query() %>%
    as_tibble

}


#' Build status of a specialized network on omnipath-metabo
#'
#' @param name Character: a network name, as listed by
#'     \code{\link{metabo_networks}}.
#'
#' @return A named list with elements \code{name}, \code{present},
#'     \code{row_count} and \code{build_id}.
#'
#' @examples
#' \dontrun{
#' metabo_network_status('metalinksdb')
#' }
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_networks}}}
#'     \item{\code{\link{metabo_network_interactions}}}
#' }
metabo_network_status <- function(name){

    .slow_doctest()

    .metabo_assert_network_name(name)

    .metabo_query_network(sprintf('networks/%s/status', name), name)

}


#' Contributing resources of a specialized network on omnipath-metabo
#'
#' @param name Character: a network name, as listed by
#'     \code{\link{metabo_networks}}.
#'
#' @return A named list with elements \code{name}, \code{kind},
#'     \code{included_sources} and \code{combined_relation}.
#'
#' @examples
#' \dontrun{
#' metabo_network_resources('metalinksdb')
#' }
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_networks}}}
#'     \item{\code{\link{metabo_network_interactions}}}
#' }
metabo_network_resources <- function(name){

    .slow_doctest()

    .metabo_assert_network_name(name)

    .metabo_query_network(sprintf('networks/%s/resources', name), name)

}


#' Interactions of a specialized network on omnipath-metabo
#'
#' Retrieves interaction rows from any network registered on the
#' omnipath-metabo web service (see \code{\link{metabo_networks}} for
#' available names). The column set is not fixed by this function -- it
#' reflects whatever the requested network exposes, so new networks are
#' served correctly without a package update.
#'
#' @param name Character: a network name, as listed by
#'     \code{\link{metabo_networks}}.
#' @param source Character: restrict the result to interactions from this
#'     contributing source. \code{NULL} (default) returns all sources.
#' @param limit Integer: maximum rows per request (server default 1000,
#'     maximum 100000).
#' @param offset Integer: number of rows to skip (for manual paging).
#' @param paginate Logical: if \code{TRUE}, repeatedly request subsequent
#'     pages (starting after \code{offset}, one request per page, each
#'     individually cached) and combine them into a single
#'     table, until a page shorter than \code{limit} is returned. If
#'     \code{FALSE} (default), returns only the first page.
#' @param ... Ignored; reserved for future extensions of the web service.
#'
#' @return A tibble of interaction rows; the column set varies by network,
#'     as described above.
#'
#' @examples
#' \dontrun{
#' page <- metabo_network_interactions('metalinksdb', limit = 10)
#' all_rows <- metabo_network_interactions('metalinksdb', paginate = TRUE)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metabo_networks}}}
#'     \item{\code{\link{metabo_network_status}}}
#' }
metabo_network_interactions <- function(
    name,
    source = NULL,
    limit = 1000L,
    offset = 0L,
    paginate = FALSE,
    ...
){

    .slow_doctest()

    .metabo_assert_network_name(name)
    .metabo_assert_count(limit, 'limit')
    .metabo_assert_count(offset, 'offset')

    fetch_page <- function(offset){

        .metabo_query_network(
            sprintf('networks/%s/interactions', name),
            name,
            source = source,
            limit = limit,
            offset = offset,
            format = 'json'
        )

    }

    rows <- fetch_page(offset)$rows %>% .metabo_rows_to_tibble()

    if (paginate) {

        repeat {

            offset <- offset + limit
            page <- fetch_page(offset)$rows %>% .metabo_rows_to_tibble()

            if (nrow(page) == 0L) break

            rows <- bind_rows(rows, page)

            if (nrow(page) < limit) break

        }

    }

    rows

}
