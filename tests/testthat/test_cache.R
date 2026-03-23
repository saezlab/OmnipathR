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

omnipath_set_console_loglevel(logger::FATAL)


# Helper: set up a temporary cache dir, return the original cachedir for
# restoration. This avoids interfering with the user's real cache.
setup_test_cache <- function() {
    original <- omnipath_get_cachedir()
    test_dir <- file.path(tempdir(), paste0('omnipathr_test_cache_', Sys.getpid()))
    dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
    omnipath_set_cachedir(test_dir)
    original
}


# Helper: restore the original cache dir and clean up the temp dir
teardown_test_cache <- function(original) {
    test_dir <- omnipath_get_cachedir()
    omnipath_set_cachedir(original)
    unlink(test_dir, recursive = TRUE)
}


test_that('omnipath_cache_key generates deterministic SHA1 hashes', {

    key1 <- omnipath_cache_key(url = 'https://example.com/data.tsv')
    key2 <- omnipath_cache_key(url = 'https://example.com/data.tsv')
    key3 <- omnipath_cache_key(url = 'https://example.com/other.tsv')

    expect_equal(key1, key2)
    expect_false(key1 == key3)
    expect_equal(nchar(key1), 40L) # SHA1 hex digest

})


test_that('omnipath_cache_key differs with POST params', {

    key1 <- omnipath_cache_key(url = 'https://example.com/api')
    key2 <- omnipath_cache_key(
        url = 'https://example.com/api',
        post = list(query = 'test')
    )

    expect_false(key1 == key2)

})


test_that('omnipath_cache_get creates and retrieves records', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    url <- 'https://test.example.com/resource1.tsv'
    record <- omnipath_cache_get(url = url)

    expect_type(record, 'list')
    expect_equal(record$url, url)
    expect_true(nchar(record$key) == 40L)
    expect_equal(record$ext, 'tsv')
    expect_type(record$versions, 'list')

    # Retrieving the same URL should return the same record
    record2 <- omnipath_cache_get(url = url)
    expect_equal(record$key, record2$key)

})


test_that('omnipath_cache_get with create=FALSE returns NULL for missing', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    record <- omnipath_cache_get(
        url = 'https://test.example.com/nonexistent',
        create = FALSE
    )

    expect_null(record)

})


test_that('omnipath_cache_save and omnipath_cache_load round-trip', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    url <- 'https://test.example.com/roundtrip'
    test_data <- data.frame(a = 1:3, b = c('x', 'y', 'z'))

    omnipath_cache_save(test_data, url = url)

    loaded <- omnipath_cache_load(url = url)
    expect_true(all.equal(loaded, test_data, check.attributes = FALSE))
    expect_equal(attr(loaded, 'origin'), 'cache')

})


test_that('omnipath_cache_load returns NULL for missing entries', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    result <- omnipath_cache_load(url = 'https://test.example.com/missing')
    expect_null(result)

})


test_that('omnipath_cache_remove by key removes only the target entry', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    urls <- paste0('https://test.example.com/resource', 1:3)
    for (url in urls) {
        omnipath_cache_save(
            data.frame(src = url),
            url = url
        )
    }

    cache_before <- OmnipathR:::omnipathr.env$cache
    expect_equal(length(cache_before), 3L)

    # Remove the first entry by key
    target_key <- omnipath_cache_key(url = urls[1])
    omnipath_cache_remove(key = target_key)

    cache_after <- OmnipathR:::omnipathr.env$cache
    expect_equal(length(cache_after), 2L)
    expect_false(target_key %in% names(cache_after))

    # Verify the remaining entries are intact with their versions
    for (remaining_url in urls[2:3]) {
        remaining_key <- omnipath_cache_key(url = remaining_url)
        expect_true(remaining_key %in% names(cache_after))
        expect_gte(length(cache_after[[remaining_key]]$versions), 1L)
    }

})


test_that('omnipath_cache_remove by URL removes only the target entry', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    urls <- paste0('https://test.example.com/by_url_', 1:4)
    for (url in urls) {
        omnipath_cache_save(data.frame(x = 1), url = url)
    }

    expect_equal(length(OmnipathR:::omnipathr.env$cache), 4L)

    # Remove by URL
    omnipath_cache_remove(url = urls[2])

    cache_after <- OmnipathR:::omnipathr.env$cache
    expect_equal(length(cache_after), 3L)

    removed_key <- omnipath_cache_key(url = urls[2])
    expect_false(removed_key %in% names(cache_after))

    # All other entries still present with versions
    for (url in urls[c(1, 3, 4)]) {
        key <- omnipath_cache_key(url = url)
        expect_true(key %in% names(cache_after))
        expect_gte(length(cache_after[[key]]$versions), 1L)
    }

})


test_that('omnipath_cache_remove with no args does not wipe cache', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    urls <- paste0('https://test.example.com/noop_', 1:3)
    for (url in urls) {
        omnipath_cache_save(data.frame(x = 1), url = url)
    }

    n_before <- length(OmnipathR:::omnipathr.env$cache)
    expect_equal(n_before, 3L)

    # Calling remove with no arguments should be a no-op
    omnipath_cache_remove()

    n_after <- length(OmnipathR:::omnipathr.env$cache)
    expect_equal(n_after, n_before)

})


test_that('omnipath_cache_remove with wipe=TRUE clears everything', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    urls <- paste0('https://test.example.com/wipe_', 1:3)
    for (url in urls) {
        omnipath_cache_save(data.frame(x = 1), url = url)
    }

    expect_gte(length(OmnipathR:::omnipathr.env$cache), 1L)

    omnipath_cache_remove(wipe = TRUE)

    expect_equal(length(OmnipathR:::omnipathr.env$cache), 0L)

})


test_that('omnipath_cache_search finds records by URL pattern', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    omnipath_cache_save(data.frame(x = 1), url = 'https://example.com/alpha')
    omnipath_cache_save(data.frame(x = 2), url = 'https://example.com/beta')
    omnipath_cache_save(data.frame(x = 3), url = 'https://other.com/alpha')

    results <- omnipath_cache_search('example\\.com')
    expect_equal(length(results), 2L)

    results_alpha <- omnipath_cache_search('alpha')
    expect_equal(length(results_alpha), 2L)

    results_beta <- omnipath_cache_search('beta')
    expect_equal(length(results_beta), 1L)

})


test_that('omnipath_cache_clean_db removes entries for missing files', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    urls <- paste0('https://test.example.com/cleandb_', 1:3)
    for (url in urls) {
        omnipath_cache_save(data.frame(x = 1), url = url)
    }

    expect_equal(length(OmnipathR:::omnipathr.env$cache), 3L)

    # Delete the file for the second entry
    key2 <- omnipath_cache_key(url = urls[2])
    record2 <- OmnipathR:::omnipathr.env$cache[[key2]]
    file_path <- record2$versions[['1']]$path
    file.remove(file_path)

    omnipath_cache_clean_db()

    cache_after <- OmnipathR:::omnipathr.env$cache
    expect_equal(length(cache_after), 2L)
    expect_false(key2 %in% names(cache_after))

})


test_that('omnipath_cache_filter_versions filters by status', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    url <- 'https://test.example.com/filter_status'
    omnipath_cache_save(data.frame(x = 1), url = url)

    key <- omnipath_cache_key(url = url)
    record <- OmnipathR:::omnipathr.env$cache[[key]]

    ready_versions <- omnipath_cache_filter_versions(
        record,
        status = 'ready'
    )
    expect_equal(length(ready_versions), 1L)

    started_versions <- omnipath_cache_filter_versions(
        record,
        status = 'started'
    )
    expect_equal(length(started_versions), 0L)

})


test_that('omnipath_cache_latest_version returns the most recent version', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    url <- 'https://test.example.com/latest_version'

    # Save twice to create 2 versions
    omnipath_cache_save(data.frame(x = 1), url = url)
    omnipath_cache_save(data.frame(x = 2), url = url)

    key <- omnipath_cache_key(url = url)
    record <- OmnipathR:::omnipathr.env$cache[[key]]
    latest <- omnipath_cache_latest_version(record)

    expect_type(latest, 'character')

})


test_that('omnipath_cache_update_status changes version status', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    url <- 'https://test.example.com/update_status'
    omnipath_cache_save(data.frame(x = 1), url = url)

    key <- omnipath_cache_key(url = url)
    record <- OmnipathR:::omnipathr.env$cache[[key]]
    expect_equal(record$versions[['1']]$status, 'ready')

    omnipath_cache_update_status(
        version = '1',
        key = key,
        status = 'failed'
    )

    record_after <- OmnipathR:::omnipathr.env$cache[[key]]
    expect_equal(record_after$versions[['1']]$status, 'failed')

})


test_that('omnipath_set_cachedir changes and respects cache directory', {

    original <- omnipath_get_cachedir()

    test_dir <- file.path(tempdir(), 'omnipathr_cachedir_test')
    dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
    on.exit({
        omnipath_set_cachedir(original)
        unlink(test_dir, recursive = TRUE)
    })

    omnipath_set_cachedir(test_dir)
    expect_equal(omnipath_get_cachedir(), test_dir)

    omnipath_set_cachedir(original)
    expect_equal(omnipath_get_cachedir(), original)

})


test_that('cache records preserve URL and extension', {

    original <- setup_test_cache()
    on.exit(teardown_test_cache(original))

    url_tsv <- 'https://test.example.com/data.tsv'
    url_csv <- 'https://test.example.com/data.csv'

    record_tsv <- omnipath_cache_get(url = url_tsv)
    record_csv <- omnipath_cache_get(url = url_csv)

    expect_equal(record_tsv$ext, 'tsv')
    expect_equal(record_csv$ext, 'csv')
    expect_equal(record_tsv$url, url_tsv)
    expect_equal(record_csv$url, url_csv)

})
