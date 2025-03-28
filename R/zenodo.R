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


#' Retrieves data from Zenodo
#'
#' Zenodo is a repository of large scientific datasets. Many projects and
#' publications make their datasets available at Zenodo. This function
#' downloads an archive from Zenodo and extracts the requested file.
#'
#' @param path Character: path to the file within the archive.
#' @param reader Optional, a function to read the connection.
#' @param reader_param List: arguments for the reader function.
#' @param url_key Character: name of the option containing the URL
#' @param zenodo_record The Zenodo record ID, either integer or character.
#' @param zenodo_fname The file name within the record.
#' @param url_param List: variables to insert into the URL string (which is
#' returned from the options).
#' @param url_key_param List: variables to insert into the `url_key`.
#' @param ... Passed to \code{archive_extractor}
#'
#' @examples
#' # an example from the OmnipathR::remap_tf_target_download function:
#' remap_dorothea <- zenodo_download(
#'     zenodo_record = 3713238,
#'     zenodo_fname = 'tf_target_sources.zip',
#'     path = (
#'         'tf_target_sources/chip_seq/remap/gene_tf_pairs_genesymbol.txt'
#'     ),
#'     reader = read_tsv,
#'     reader_param = list(
#'         col_names = c(
#'             'source_genesymbol',
#'             'target_genesymbol',
#'             'target_ensembl',
#'             'score'
#'         ),
#'         col_types = cols(),
#'         progress = FALSE
#'     ),
#'   resource = 'ReMap'
#' )
#'
#' @return A connection
#'
#' @importFrom magrittr %<>% %>%
#' @export
zenodo_download <- function(
    path,
    reader = NULL,
    reader_param = list(),
    url_key = NULL,
    zenodo_record = NULL,
    zenodo_fname = NULL,
    url_param = list(),
    url_key_param = list(),
    ...
){

    .slow_doctest()

    if(
        !is.null(zenodo_record) &&
        !is.null(zenodo_fname) &&
        is.null(url_key)
    ){

        url_key <- 'zenodo'
        url_param <- list(
            as.character(zenodo_record),
            zenodo_fname
        )

    }

    archive_extractor(
        path = path,
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param,
        reader = reader,
        reader_param = reader_param,
        ...
    )

}
