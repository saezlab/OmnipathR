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

library(OmnipathR)
library(tibble)
library(dplyr)
library(magrittr)
library(tidyselect)


omnipath_set_console_loglevel('fatal')

input0 <- tibble(
    a = c('A', 'A', 'A', 'B', 'C', 'D', 'E', 'F'),
    b = c('a', 'b', 'c', 'd', 'a', 'e', 'e', NA)
)

input1 <-
    input0 %>%
    group_by(a) %>%
    mutate(b = list(b)) %>%
    summarize_all(~extract(., 1L))

output_expanded <- tibble(
    a = c('A', 'A', 'A', 'B', 'C', 'D', 'E', 'F'),
    b = c('a', 'b', 'c', 'd', 'a', 'e', 'e', NA),
    a_b_to_ambiguity = c(3L, 3L, 3L, 1L, 1L, 1L, 1L, 0L),
    a_b_from_ambiguity = c(2L, 1L, 1L, 1L, 2L, 2L, 2L, 1L),
    a_b_ambiguity = c(
        'many-to-many',
        'one-to-many',
        'one-to-many',
        'one-to-one',
        'many-to-one',
        'many-to-one',
        'many-to-one',
        'one-to-none'
    )
)

output_collapsed <-
    output_expanded %>%
    group_by(a) %>%
    reframe(across(everything(), ~list(.x))) %>%
    ungroup %>%
    rowwise %>%
    mutate(across(ends_with('ambiguity'), ~list(set_names(.x, replace_na(b, ''))))) %>%
    ungroup


test_that(
    'Translation ambiguity [expand-expanded]',
    expect_equal(output_expanded, ambiguity(input0, a, b))
)

test_that(
    'Translation ambiguity [expand-collapsed]',
    expect_equal(output_expanded, ambiguity(input1, a, b, expand = TRUE))
)

test_that(
    'Translation ambiguity [collapse-expanded]',
    expect_equal(output_collapsed, ambiguity(input0, a, b, expand = FALSE))
)

test_that(
    'Translation ambiguity [collapse-collapsed]',
    expect_equal(output_collapsed, ambiguity(input1, a, b))
)
