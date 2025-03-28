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


RAMP_METABOLITE_ID_TYPES <- c(
    "hmdb",
    "chebi",
    "chemspider",
    "kegg",
    "pubchem",
    "pubchem_cid",
    "CAS",
    "cas",
    "wikidata",
    "LIPIDMAPS",
    "lipidmaps",
    "lipidbank",
    "swisslipids",
    "plantfa",
    "kegg_glycan",
    "rhea-comp",
    "rhea",
    "polymer"
)

RAMP_GENE_ID_TYPES <- c(
    "hmdb",
    "wikidata",
    "ensg",
    "ensembl",
    "genesymbol",
    "gene_symbol",
    "uniprot",
    "entrez",
    "EN",
    "en",
    "ncbiprotein",
    "brenda"
)

RAMP_LATEST_VERSION <- '2.5.4'


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
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{ramp_tables}}}
#' }
ramp_sqlite <- function(version = RAMP_LATEST_VERSION) {

    .slow_doctest()

    sqlite_downloader('ramp', version = version, compr = 'gz')

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
ramp_tables <- function(version = RAMP_LATEST_VERSION) {

    .slow_doctest()

    version %>% ramp_sqlite %>% dbListTables

}


#' Return table from RaMP database

#' @param name Character. The name of the RaMP table to fetch.
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
ramp_table <- function(name, version = RAMP_LATEST_VERSION) {

    .slow_doctest()

    version %>% ramp_sqlite %>% dbReadTable(name) %>% as_tibble()

}


#' Pairwise ID translation table from RaMP database

#' @param from Character or Symbol. Name of an identifier type.
#' @param to Character or Symbol. Name of an identifier type.
#' @param version Character. The version of RaMP to download.
#'
#' @return Dataframe of pairs of identifiers.
#'
#' @examples
#' ramp_id_mapping_table('hmdb', 'kegg')
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr inner_join mutate filter select distinct
#' @importFrom tidyr separate_wider_delim
#' @importFrom stringr str_replace str_detect
#' @importFrom rlang !! enquo
#' @noRd
#' @seealso \itemize{
#'     \item{\code{\link{ramp_sqlite}}}
#'     \item{\code{\link{ramp_tables}}}
#'     \item{\code{\link{ramp_table}}}
#' }
ramp_id_mapping_table_impl <- function(from, to, version = RAMP_LATEST_VERSION) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    rampId <- sourceId <- resource <- ID <- NULL

    from <- .nse_ensure_str(!!enquo(from)) %>% ramp_id_type
    to <- .nse_ensure_str(!!enquo(to)) %>% ramp_id_type

    version %>%
    ramp_table('source', .) %>%
    mutate(
        sourceId = str_replace(sourceId, 'swisslipids:SLM', 'swisslipids')
    ) %>%
    separate_wider_delim(
        sourceId,
        names = c('resource', 'ID'),
        delim = ':',
        too_many = 'merge'
    ) %>%
    filter(!str_detect(ID, ':')) %>%
    mutate(ID = ifelse(resource == 'swisslipids', sprintf('SLM:%s', ID), ID)) %>%
    {inner_join(
        filter(., resource == from) %>% select(From = ID, rampId),
        filter(., resource == to) %>% select(To = ID, rampId),
        by = 'rampId'
    )} %>%
    select(-rampId) %>%
    distinct

}
