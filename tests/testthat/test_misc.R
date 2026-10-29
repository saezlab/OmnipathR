#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2026
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


library(OmnipathR)
library(testthat)


test_that('UniProt downloads use HTTP/1.1, other hosts are unchanged', {

    expect_identical(
        OmnipathR:::host_curlopt(
            'https://rest.uniprot.org/uniprotkb/stream?format=tsv'
        ),
        list(http_version = 2L)
    )
    expect_identical(
        OmnipathR:::host_curlopt('https://omnipathdb.org/interactions'),
        list()
    )

})


test_that('On a build server, a failed UniProt download gives an empty table', {

    env <- OmnipathR:::omnipathr.env
    withr::defer(env$uniprot_down <- NULL)
    calls <- 0L

    local_mocked_bindings(
        .on_buildserver = function(...) TRUE,
        generic_downloader = function(...) {
            calls <<- calls + 1L
            stop('HTTP/2 stream 1 was not closed cleanly')
        },
        .package = 'OmnipathR'
    )

    result <- expect_warning(
        OmnipathR:::uniprot_download('new_uniprot', list(), 'accession,id'),
        'Failed to download from UniProt'
    )
    expect_identical(names(result), c('accession', 'id'))
    expect_identical(nrow(result), 0L)

    # the next query does not try UniProt again
    expect_identical(
        OmnipathR:::uniprot_download('new_uniprot', list(), 'accession'),
        tibble::tibble(accession = character(0L))
    )
    expect_identical(calls, 1L)

})


test_that('Outside build servers, a failed UniProt download is an error', {

    withr::local_options(omnipathr.console_loglevel = 'fatal')

    local_mocked_bindings(
        generic_downloader = function(...) stop('connection reset'),
        .on_buildserver = function(...) FALSE,
        .package = 'OmnipathR'
    )

    expect_error(
        OmnipathR:::uniprot_download('new_uniprot', list(), 'accession'),
        'connection reset'
    )

})
