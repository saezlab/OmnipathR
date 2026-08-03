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


library(OmnipathR)
library(logger)

# not to interfere with testthat's console display
omnipath_set_console_loglevel(logger::FATAL)


#' Can we reach the omnipath-metabo web service at all?
can_connect_metabo <- function(){

    tryCatch(
        as.logical(nrow(metabo_networks())),
        error = function(...) FALSE
    )

}


if (can_connect_metabo()) {

    test_that(
        'metabo_networks lists the known registered networks',
        {
            networks <- metabo_networks()
            expect_true(all(
                c('name', 'kind', 'included_sources') %in% names(networks)
            ))
            expect_true(all(
                c('metalinksdb', 'liana') %in% networks$name
            ))
        }
    )

    test_that(
        'metabo_network_status reports a present, non-empty network',
        {
            status <- metabo_network_status('metalinksdb')
            expect_true(status$present)
            expect_gt(status$row_count, 0L)
        }
    )

    test_that(
        'metabo_network_resources reports contributing sources',
        {
            resources <- metabo_network_resources('metalinksdb')
            expect_true(length(resources$included_sources) > 0L)
        }
    )

    test_that(
        'metabo_network_interactions respects the page size limit',
        {
            page <- metabo_network_interactions('metalinksdb', limit = 10L)
            expect_lte(nrow(page), 10L)
            expect_gt(nrow(page), 0L)
        }
    )

    test_that(
        'metabo_network_interactions raises a clear error for an unknown network',
        {
            expect_error(
                metabo_network_interactions('not-a-real-network'),
                regexp = 'Unknown network'
            )
        }
    )

} else {

    test_that(
        'Skipping omnipath-metabo networks tests, could not establish connection.',
        expect_true(TRUE)
    )

}
