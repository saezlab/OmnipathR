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


test_that('The Ensembl organisms table is read from species.json', {

    json <- withr::local_tempfile(fileext = '.json')
    writeLines(
        '{"last_updated": "2026-10-05", "species": {
            "Homo_sapiens": {"taxid": 9606, "species_taxonomy_id": 9606,
                "scientific_name": "Homo sapiens", "common_name": "Human",
                "assemblies": {}},
            "Escherichia_coli_str_K_12_substr_MG1655": {"taxid": 511145,
                "species_taxonomy_id": 562,
                "scientific_name": "Escherichia coli str. K-12 substr. MG1655",
                "common_name": null, "assemblies": {}}
        }}',
        json
    )

    local_mocked_bindings(
        download_to_cache = function(...) json,
        .package = 'OmnipathR'
    )

    result <- ensembl_organisms_raw()

    expect_identical(
        names(result),
        c('species', 'ncbi_tax_id', 'species_tax_id', 'latin_name', 'common_name')
    )
    expect_identical(result$ncbi_tax_id, c(9606L, 511145L))
    expect_identical(result$species_tax_id, c(9606L, 562L))
    expect_identical(result$common_name, c('Human', NA_character_))

})


test_that('Ensembl organism IDs are the BioMart dataset names', {

    local_mocked_bindings(
        ensembl_organisms_raw = function() tibble::tibble(
            species = c(
                'Homo_sapiens', 'Malus_domestica',
                'Faecalibacterium_prausnitzii_str_CNCM_I_4541',
                'Faecalibacterium_prausnitzii'
            ),
            ncbi_tax_id = c(9606L, 3750L, 853L, 853L),
            species_tax_id = c(9606L, 3750L, 853L, 853L),
            latin_name = c(
                'Homo sapiens', 'Malus domestica',
                'Faecalibacterium prausnitzii str. CNCM I 4541',
                'Faecalibacterium prausnitzii'
            ),
            common_name = c('Human', 'Apple', NA, NA)
        ),
        .package = 'OmnipathR'
    )

    result <- ensembl_organisms()

    expect_identical(
        names(result),
        c(
            'common_name', 'latin_name', 'ncbi_tax_id',
            'ensembl_id', 'common_name_alt'
        )
    )
    expect_false(anyDuplicated(result$ncbi_tax_id) > 0L)
    expect_identical(
        result$ensembl_id[result$ncbi_tax_id == 9606L],
        'hsapiens'
    )
    # the same short name in BioMart belongs to the opossum, not the apple
    expect_identical(
        result$ensembl_id[result$ncbi_tax_id == 3750L],
        NA_character_
    )
    # organisms only in BioMart are kept
    expect_identical(
        result$ensembl_id[result$ncbi_tax_id == 13616L],
        'mdomestica'
    )
    # of two entries with the same taxon, the species is kept, not the strain
    expect_identical(
        result$latin_name[result$ncbi_tax_id == 853L],
        'Faecalibacterium prausnitzii'
    )

})


test_that('A failed organism name source gives an empty table', {

    withr::local_options(omnipathr.console_loglevel = 'fatal')

    result <- expect_warning(
        OmnipathR:::organism_table_safe(
            function() stop('HTTP 504 Gateway Timeout'),
            ncbi_tax_id = integer(),
            latin_name = character()
        ),
        'HTTP 504'
    )
    expect_identical(
        result,
        tibble::tibble(ncbi_tax_id = integer(), latin_name = character())
    )

    expect_warning(
        OmnipathR:::organism_table_safe(
            function() tibble::tibble(taxon = 1L),
            ncbi_tax_id = integer()
        ),
        'missing columns'
    )

})
