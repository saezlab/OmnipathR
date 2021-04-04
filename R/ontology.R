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


#' Nested list from a table of ontology relations
#'
#' @param relations A data frame of ontology relations (the "relations"
#'     element of the list returned by \code{\link{obo_parser}} in case
#'     its argument `tables` is \code{TRUE}).
#'
#' @return The relations converted to a nested list.
#'
#' @examples
#' goslim_url <-
#'     "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"
#' path <- tempfile()
#' download.file(goslim_url, destfile = path, quiet = TRUE)
#' obo <- obo_reader(path, tables = TRUE)
#' unlink(path)
#' rel_list <- relations_table_to_list(obo$relations)
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !! sym
#' @importFrom tidyr chop replace_na
#' @importFrom purrr map2
#' @importFrom dplyr mutate
#' @importFrom logger log_trace
#'
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{relations_list_to_table}}}
#'     \item{\code{\link{swap_relations}}}
#'     \item{\code{\link{obo_parser}}}
#' }
relations_table_to_list <- function(relations){

    direction <- c('parents', 'children')
    to_str <-
        intersect(names(relations), direction)[1] %>%
        replace_na('side2')
    to_sym <- sym(to_str)
    direction <- `if`(
        to_str == 'parents',
        direction,
        `if`(
            to_str == 'children',
            rev(direction),
            NULL
        )
    )

    log_trace('Converting ontology relations from table to list.')

    relations %>%
    chop(c('relation', to_str)) %>%
    mutate(value = map2(!!to_sym, relation, setNames)) %>%
    term_value_list %>%
    `attr<-`('direction', direction)

}


#' Table from a nested list of ontology relations
#'
#' Converting the nested list to a table is a more costly operation, it takes
#' a few seconds. Best to do it only once, or pass \code{tables = TRUE} to
#' \code{\link{obo_parser}}, and convert the data frame to list, if you
#' also need it in list format.
#'
#' @param relations A nested list of ontology relations (the "relations"
#'     element of the list returned by \code{\link{obo_parser}} in case
#'     its argument `tables` is \code{FALSE}).
#' @param direction Override the direction (i.e. child -> parents or parent
#'     -> children). The nested lists produced by functions in the current
#'     package add an attribute "direction" thus no need to pass this value.
#'     If the attribute and the argument are both missing, the column will
#'     be named simply "side2" and it won't be clear whether the relations
#'     point from "term" to "side2" or the other way around. The direction
#'     should be a character vector of length 2 with the values "parents"
#'     and "children".
#'
#' @return The relations converted to a data frame (tibble).
#'
#' @examples
#' goslim_url <-
#'     "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"
#' path <- tempfile()
#' download.file(goslim_url, destfile = path, quiet = TRUE)
#' obo <- obo_reader(path, tables = FALSE)
#' unlink(path)
#' rel_tbl <- relations_list_to_table(obo$relations)
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !! sym %||%
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest_longer
#' @importFrom logger log_trace
#'
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{swap_relations}}}
#'     \item{\code{\link{relations_table_to_list}}}
#'     \item{\code{\link{obo_parser}}}
#' }
relations_list_to_table <- function(relations, direction = NULL){

    direction <-
        direction %||%
        attr(relations, 'direction') %||%
        c('side1', 'side2')
    to_str <- direction[2]
    to_sym <- sym(to_str)

    log_trace('Converting ontology relations from list to table.')

    relations %>%
    tibble(rel = .) %>%
    mutate(term = names(rel)) %>%
    unnest_longer(
        rel,
        values_to = to_str,
        indices_to = 'relation'
    ) %>%
    select(term, relation, !!to_sym)

}


#' Reverse the direction of ontology relations
#'
#' @param relations The `relations` component of the data returned by
#'     \code{\link{obo_parser}} or any `...ontology_download` function
#'     such as \code{\link{go_ontology_download}}. Depending on the
#'     \code{tables} argument of those functions the `relations` can be
#'     a data frame or a nested list.
#'
#' @return Same type as the input, but the relations swapped: if in the input
#'     these pointed from each child to the parents, in the output they
#'     point from each parent to their children, and vice versa.
#'
#' @examples
#' goslim_url <-
#'     "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"
#' path <- tempfile()
#' download.file(goslim_url, destfile = path, quiet = TRUE)
#' obo <- obo_reader(path)
#' unlink(path)
#' rel_swapped <- swap_relations(obo$relations)
#'
#' @importFrom magrittr %>%
#' @importFrom rlang sym !! :=
#' @importFrom dplyr select summarize group_by
#' @importFrom tidyr unnest
#' @importFrom logger log_fatal log_trace
#'
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{relations_list_to_table}}}
#'     \item{\code{\link{relations_table_to_list}}}
#'     \item{\code{\link{obo_parser}}}
#' }
swap_relations <- function(relations){

    dfclass <- inherits(relations, 'data.frame')

    if(dfclass){
        if('parents' %in% names(relations)){
            c_in <- sym('parents')
            c_out <- sym('children')
        }else if('children' %in% names(relations)){
            c_in <- sym('children')
            c_out <- sym('parents')
        }else if('side2' %in% names(relations)){
            c_in <- sym('side2')
            c_out <- sym('side2')
        }else{
            msg <- paste0(
                'swap_relations: the input data frame must have one of the ',
                'following columns: "parents", "children" or "side2".'
            )
            log_fatal(msg)
            stop(msg)
        }
    }

    log_trace('Swapping direction of ontology relations.')

    `if`(
        dfclass,
        relations %>%
            unnest(!!c_in) %>%
            group_by(!!c_in, relation) %>%
            summarize(term = list(term), .groups = 'drop') %>%
            select(term = !!c_in, relation, !!c_out := term),
        relations %>%
            relations_list_to_table() %>%
            swap_relations() %>%
            relations_table_to_list()
    )

}


descendants <- function(
    descendants_of,
    onto
){}