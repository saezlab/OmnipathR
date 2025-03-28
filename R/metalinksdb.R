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


#' Open MetalinksDB as an SQLite3 connection
#'
#' MetalinksDB is a database of metabolite-protein and small molecule
#' ligand-receptor interactions.
#'
#' @return An SQLite3 connection.
#'
#' @examples
#' con <- metalinksdb_sqlite()
#' con
#'
#' @importFrom magrittr %>%
#' @export
metalinksdb_sqlite <- function() {

    'metalinksdb' %>%
    sqlite_downloader()

}


#' List tables in MetalinksDB
#'
#' @return Character vector of table names in the MetalinksDB SQLite database.
#'
#' @examples
#' metalinksdb_tables()
#'
#' @importFrom magrittr %>%
#' @importFrom RSQLite dbListTables dbDisconnect
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metalinksdb_sqlite}}}
#' }
metalinksdb_tables <- function() {

    con <- metalinksdb_sqlite()
    on.exit(dbDisconnect(con))

    con %>% dbListTables

}


#' A table from MetalinksDB

#' @param name Character. The name of the MetalinksDB table to fetch.
#'
#' @return A data frame (tibble) of one table from the MetalinksDB SQLite database.
#'
#' @examples
#' metalinksdb_table('pathway')
#'
#' @importFrom magrittr %>%
#' @importFrom RSQLite dbReadTable dbDisconnect
#' @importFrom tibble as_tibble
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{metalinksdb_sqlite}}}
#'     \item{\code{\link{metalinksdb_tables}}}
#' }
metalinksdb_table <- function(name) {

    con <- metalinksdb_sqlite()
    on.exit(dbDisconnect(con))

    con %>% dbReadTable(name) %>% as_tibble()

}
