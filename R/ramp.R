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


#' Download and open RaMP database SQLite
#'
#' @param version Character. The version of RaMP to download.
#'
#' @return SQLite connection.
#'
#' @examples
#' sqlite_con <- ramp_sqlite()
#'
#' @importFrom magrittr %>%
#' @importFrom R.utils gunzip
#' @importFrom RSQLite SQLite dbConnect
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{ramp_tables}}}
#' }
ramp_sqlite <- function(version = '2.5.4') {

    .slow_doctest()

    fake_url <- version %>% sprintf('RAMP_%s.sqlite', .)
    cache_record <- omnipath_cache_latest_or_new(url = fake_url)

    if (!cache_record$status == CACHE_STATUS$READY) {

        'ramp' %>%
        download_to_cache(url_param = list(version)) %>%
        gunzip(destname = cache_record$path, remove = FALSE)

        omnipath_cache_download_ready(cache_record)

    }

    dbConnect(SQLite(), cache_record$path)

}


#' List tables in RaMP database

#' @param version Character. The version of RaMP to download.
#'
#' @return Character vector of table names in the RaMP SQLite database.
#'
#' @examples
#' ramp_tables()
#'
#' @importFrom magrittr %>%
#' @importFrom RSQLite dbListTables
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{ramp_sqlite}}}
#' }
ramp_tables <- function(version = '2.5.4') {

    .slow_doctest()

    version %>% ramp_sqlite %>% dbListTables

}


#' Return table from RaMP database

#' @param version Character. The version of RaMP to download.
#'
#' @return Character vector of table names in the RaMP SQLite database.
#'
#' @examples
#' ramp_table('source')
#'
#' @importFrom magrittr %>%
#' @importFrom RSQLite dbReadTable
#' @importFrom tibble as_tibble
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{ramp_sqlite}}}
#'     \item{\code{\link{ramp_tables}}}
#' }
ramp_table <- function(name, version = '2.5.4') {

    .slow_doctest()

    version %>% ramp_sqlite %>% dbReadTable(name) %>% as_tibble()

}

