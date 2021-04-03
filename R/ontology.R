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
#' @importFrom dplyr mutate select summarize group_by
#' @importFrom tidyr chop unnest unnest_longer
#' @importFrom purrr map2
#' @export
swap_relations <- function(relations){

    dfclass <- inherits(relations, 'data.frame')

    if(dfclass){
        if('parents' %in% names(relations)){
            c_in <- sym('parents')
            c_out <- sym('children')
        }else{
            c_in <- sym('children')
            c_out <- sym('parents')
        }
    }

    `if`(
        dfclass,
        relations %>%
            unnest(!!c_in) %>%
            group_by(!!c_in, relation) %>%
            summarize(term = list(term), .groups = 'drop') %>%
            select(term = !!c_in, relation, !!c_out := term),
        relations %>%
            tibble(rel = .) %>%
            mutate(term = names(rel)) %>%
            unnest_longer(
                rel,
                values_to = 'parents',
                indices_to = 'relation'
            ) %>%
            swap_relations() %>%
            chop(c('relation', 'children')) %>%
            mutate(value = map2(children, relation, setNames)) %>%
            term_value_list
    )

}


descendants <- function(
    descendants_of,
    onto
){}