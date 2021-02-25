#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2021
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
#  Website: https://saezlab.github.io/omnipathr
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
    url_key_param = list()
){

    if(
        !is.null(zenodo_record) &&
        !is.null(zenodo_fname) &&
        is.null(url_key)
    ){

        url_key <- 'omnipath.zenodo_url'
        url_param <- list(
            as.character(zenodo_record),
            zenodo_fname
        )

    }

    zip_extractor(
        path = path,
        url_key = url_key,
        url_key_param = url_key_param,
        url_param = url_param,
        reader = reader,
        reader_param = reader_param
    )

}