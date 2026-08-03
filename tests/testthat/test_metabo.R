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
        as.logical(length(metabo_cosmos_categories())),
        error = function(...) FALSE
    )

}


if (can_connect_metabo()) {

    test_that(
        'metabo_cosmos_categories returns the known category names',
        {
            categories <- metabo_cosmos_categories()
            expect_true(is.character(categories))
            expect_true(all(
                c('transporters', 'receptors', 'enzyme_metabolite') %in%
                categories
            ))
        }
    )

    test_that(
        'metabo_cosmos_organisms returns human',
        {
            organisms <- metabo_cosmos_organisms()
            expect_true(is.integer(organisms))
            expect_true(9606L %in% organisms)
        }
    )

    test_that(
        'metabo_cosmos_resources returns a category/resource tibble',
        {
            resources <- metabo_cosmos_resources()
            expect_true(all(c('category', 'resource') %in% names(resources)))
            expect_gt(nrow(resources), 0L)
        }
    )

    test_that(
        'metabo_cosmos_pkn returns interactions for a single category',
        {
            transporters <- metabo_cosmos_pkn(
                organism = 9606,
                categories = 'transporters'
            )
            expect_gt(nrow(transporters), 0L)
            expect_true(all(
                c('source', 'target', 'resource', 'mor', 'attrs') %in%
                names(transporters)
            ))
            expect_true(!anyNA(transporters$interaction_type))
        }
    )

    test_that(
        'metabo_cosmos_pkn returns an empty result for an organism with no data',
        {
            empty <- metabo_cosmos_pkn(organism = 1L, categories = 'transporters')
            expect_equal(nrow(empty), 0L)
            expect_true('source' %in% names(empty))
        }
    )

    test_that(
        'metabo_metabolic_signaling_pkn equals the "all" category query',
        {
            signaling <- metabo_metabolic_signaling_pkn(organism = 9606)
            full <- metabo_cosmos_pkn(organism = 9606, categories = 'all')
            expect_equal(nrow(signaling), nrow(full))
        }
    )

    test_that(
        'metabo_metabolite_protein_pkn is restricted to metabolite-protein edges',
        {
            metalinks_like <- metabo_metabolite_protein_pkn(organism = 9606)
            expect_gt(nrow(metalinks_like), 0L)
            expect_true(all(
                !metalinks_like$interaction_type %in%
                c('transport', 'catalysis', 'signaling', 'gene_regulation')
            ))
        }
    )

} else {

    test_that(
        'Skipping omnipath-metabo tests, could not establish connection.',
        expect_true(TRUE)
    )

}
