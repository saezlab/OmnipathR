#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2020
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


library(OmnipathR)
library(dplyr)
library(magrittr)
library(tidyr)
library(jsonlite)
library(logger)

# not to interfere with testthat's console display
omnipath_set_console_loglevel(logger::FATAL)


#' Retrieves a list of resources for a certain query type and dataset.
#' We test each query with downloading only a few resources.
get_resources_test <- function(
        query_type,
        dataset = NULL,
        license = NULL,
        top = NULL
    ){

    license_levels <- list(
        academic = 1,
        commercial = 2
    )

    resources_url <- 'https://omnipathdb.org/resources'
    resources <- jsonlite::fromJSON(txt = resources_url)

    resources <- Filter(
        function(res_name){

            res <- resources[[res_name]]

            (
                query_type %in% names(res$queries) && (
                    is.null(dataset) ||
                    dataset %in% res$queries[[query_type]]$datasets
                ) && (
                    is.null(license) || (
                        license_levels[[res$license$purpose]] >=
                        license_levels[[license]]
                    )
                )
            )

        },
        names(resources)
    )

    resources <- sort(unique(resources))

    `if`(is.null(top), resources, head(resources, top))

}


#' We want to avoid that tests fail just because some generic, system
#' level problem with network or HTTP(S) connections.
can_we_download_anything_at_all <- function(){

    tryCatch(
        as.logical(
            nrow(
                read.csv('https://omnipathdb.org/queries/interactions')
            )
        ),
        error = function(...){FALSE}
    )

}


#' Try also with jsonlite
can_jsonlite_download <- function(){

    tryCatch(
        as.logical(
            length(
                jsonlite::fromJSON('https://omnipathdb.org/resources')
            )
        ),
        error = function(...){FALSE}
    )

}


#' Test if we can send a log message
test_log <- function(){

    test_that(
        'Testing logger',
        expect_error(logger::log_info('From tests with love!'), NA)
    )

}


test_log()

# Here we build a list with parameters that we later iterate through
# in order to test all the various download methods.
test_items <- list(
    complexes = list(qt = 'complexes'),
    enzyme_substrate = list(qt = 'enzsub'),
    intercell = list(qt = 'intercell'),
    annotations = list(qt = 'annotations'),
    omnipath = list(datasets = 'omnipath'),
    kinaseextra = list(datasets = 'kinaseextra'),
    pathwayextra = list(datasets = 'pathwayextra'),
    ligrecextra = list(datasets = 'ligrecextra'),
    tf_target = list(datasets = 'tf_target'),
    tf_mirna = list(datasets = 'tf_mirna'),
    mirna_target = list(datasets = 'mirnatarget')
)

if(can_we_download_anything_at_all() && can_jsonlite_download()){

    top_env <- environment()

    omnipath_cache_clean_db()

    for(item in names(test_items)){

        resource_col <- list(
            complexes = 'sources',
            interactions = 'sources',
            enzsub = 'sources',
            annotations = 'source',
            intercell = 'database'
        )

        method_name <- item
        method <- get(method_name)
        dataset <- test_items[[item]]$datasets
        query_type <-
            test_items[[item]]$qt %>%
            OmnipathR:::if_null('interactions')

        resources <- get_resources_test(
            query_type = query_type,
            dataset = dataset,
            license = 'academic',
            top = 5
        )

        test_that(
            sprintf('Testing `%s`', method_name),
            expect_error(
                assign(
                    'response',
                    method(resources = resources),
                    envir = top_env
                ),
                NA
            )
        )

        if(query_type %in% c('enzsub', 'interactions', 'complexes')){
            response %<>% separate_rows(sources, sep = ';')
        }

        resources_in_response <- response %>%
            pull(resource_col[[query_type]]) %>%
            unique()

        label <- sprintf(
            '%s%s',
            query_type,
            `if`(is.null(dataset), '', sprintf(' (%s)', dataset))
        )

        test_that(
            sprintf('%s: resources', label),
            expect_true(
                all(resources %in% resources_in_response)
            )
        )

        test_that(
            sprintf('%s: numof rows', label),
            expect_gt(nrow(response), 10)
        )

    }

}else{

    test_that(
        'Skipping tests, could not establish connection.',
        expect_true(TRUE)
    )

}
