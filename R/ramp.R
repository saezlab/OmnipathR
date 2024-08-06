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
#' @importFrom R.utils gunzip
#' @importFrom RSQLite SQLite dbConnect
#' @noRd
ramp_sqlite <- function(version = '2.5.4') {

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

#' @importFrom RSQLite SQLite dbListTables
#' @noRd
ramp_tables <- function(version = '2.5.4') {

    version %>% ramp_sqlite %>% dbListTables

}