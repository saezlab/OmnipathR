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


#' Generic OBO parser
#'
#' Reads the contents of an OBO file and processes it into a list based
#' data structure.
#'
#' @param path Path to the OBO file.
#' @param relations Character vector: process only these relations.
#' @param shorten_namespace Logical: shorten the namespace to a single
#'     code (as usual for Gene Ontology, e.g. cellular_component = "C").
#' @param tables Logical: return data frames (tibbles) instead of nested
#'     lists.
#'
#' @return A list with the following elements: 1) "names" a list with
#'     terms as names and names as values; 2) "namespaces" a list with
#'     terms as names and namespaces as values; 3) "relations" a list with
#'     relations between terms: terms are keys, values are lists with
#'     relations as names and character vectors of related terms as
#'     values; 4) "subsets" a list with terms as keys and character
#'     vectors of subset names as values (or \code{NULL} if the term
#'     does not belong to any subset); 5) "obsolete" character vector
#'     with all the terms labeled as obsolete. If the \code{tables}
#'     parameter is \code{TRUE},
#'
#' @examples
#' goslim_url <-
#'     "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"
#' path <- tempfile()
#' download.file(goslim_url, destfile = path, quiet = TRUE)
#' obo <- obo_reader(path, tables = FALSE)
#' unlink(path)
#' names(obo)
#' # [1] "names"      "namespaces" "relations"  "subsets"    "obsolete"
#' head(obo$relations, n = 2)
#' # $`GO:0000001`
#' # $`GO:0000001`$is_a
#' # [1] "GO:0048308" "GO:0048311"
#' #
#' # $`GO:0000002`
#' # $`GO:0000002`$is_a
#' # [1] "GO:0007005"
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom stringr str_split_fixed str_trim str_sub
#' @importFrom stringr str_to_upper str_split
#' @importFrom dplyr mutate filter select pull summarize group_by
#' @importFrom tidyr separate chop
#' @importFrom purrr map_chr map2
#' @importFrom logger log_trace
#' @export
obo_parser <- function(
    path,
    relations = c(
        'is_a', 'part_of', 'occurs_in', 'regulates',
        'positively_regulates', 'negatively_regulates'
    ),
    shorten_namespace = TRUE,
    tables = TRUE
){

    short_namespace <- function(ns){

        ns %>%
        str_split('_') %>%
        map_chr(last) %>%
        str_sub(1, 1) %>%
        str_to_upper

    }

    log_trace('Processing OBO file `%s`.', path)

    proc_namespace <- `if`(shorten_namespace, short_namespace, identity)

    con <- file(path, 'r')
    on.exit(close(con))

    term_name <- list()
    term_relations <- list()
    term_namespace <- list()
    obsolete_terms <- NULL
    term_subset <- list()

    raw <-
        readLines(con) %>%
        str_split_fixed(':', n = 2) %>%
        `colnames<-`(c('key', 'value')) %>%
        as_tibble() %>%
        mutate(
            value = str_trim(value),
            row = 1:n()
        )

    terms <-
        raw %>%
        filter(key == 'id')

    raw %<>%
        mutate(
            term = cut(
                row,
                c(terms$row - 1, Inf),
                terms$value,
                right = FALSE
            )
        ) %>%
        select(-row) %>%
        filter(!is.na(term))

    externals <-
        raw %>%
        filter(key == 'namespace' & value == 'external') %>%
        pull(term)

    raw %<>% filter(!term %in% externals)

    term_name <-
        raw %>%
        filter(key == 'name') %>%
        {`if`(
            tables,
            select(., term, name = value),
            term_value_list(.)
        )}

    term_namespace <-
        raw %>%
        filter(key == 'namespace') %>%
        mutate(value = proc_namespace(value)) %>%
        {`if`(
            tables,
            select(., term, namespace = value),
            term_value_list(.)
        )}


    obsolete_terms <-
        raw %>%
        filter(key == 'is_obsolete' & value == 'true') %>%
        pull(term) %>%
        as.character

    term_subset <-
        raw %>%
        filter(key == 'subset') %>%
        {`if`(
            tables,
            select(., term, subset = value),
            group_by(., term) %>%
            summarize(value = list(value)) %>%
            term_value_list(.)
        )}

    term_relations <-
        raw %>%
        filter(key %in% relations | key == 'relationship') %>%
        mutate(
            value =
                str_split(value, '!', n = 2) %>%
                map_chr(first) %>%
                str_trim
        ) %>%
        {bind_rows(
            filter(., key %in% relations),
            filter(., key == 'relationship') %>%
            separate(value, c('key', 'value'), sep = ' ')
        )} %>%
        group_by(term, key) %>%
        summarize(value = list(value), .groups = 'drop') %>%
        select(., term, relation = key, parents = value) %>%
        {`if`(
            tables,
            .,
            relations_table_to_list(.)
        )}

    log_trace('Finished processing OBO file `%s`.', path)

    list(
        names = term_name,
        namespaces = term_namespace,
        relations = term_relations,
        subsets = term_subset,
        obsolete = obsolete_terms
    )

}


#' Transform `term` and `value` columns of a tibble to a named list
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @noRd
term_value_list <- function(d){

    d %>%
    {setNames(
        as.list(pull(., value)),
        pull(., term)
    )}

}