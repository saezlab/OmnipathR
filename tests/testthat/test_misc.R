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


test_that('The environment variable marks a build server', {

    withr::local_envvar(OMNIPATHR_BUILDSERVER = 'true')
    expect_true(OmnipathR:::.on_buildserver(user = 'someone'))

    withr::local_envvar(OMNIPATHR_BUILDSERVER = 'TRUE')
    expect_true(OmnipathR:::.on_buildserver(user = 'someone'))

})


test_that('Without the variable, only the build users are build servers', {

    withr::local_envvar(OMNIPATHR_BUILDSERVER = NA)
    expect_false(OmnipathR:::.on_buildserver(user = 'someone'))
    expect_true(OmnipathR:::.on_buildserver(user = 'biocbuild'))
    expect_true(OmnipathR:::.on_buildserver(user = 'omnipath'))

    withr::local_envvar(OMNIPATHR_BUILDSERVER = 'false')
    expect_false(OmnipathR:::.on_buildserver(user = 'someone'))

})
