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


#' Creates an igraph object which helps to control transformations of
#' ontology relation data structures
#'
#' @importFrom igraph graph_from_literal as.directed V E get.edgelist
#' @importFrom purrr pmap_chr
#' @importFrom stringr str_sub
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#'
#' @noRd
get_ontology_db_variants_graph <- function(){

    g <-
        graph_from_literal(
            rel_tbl_c2p-rel_tbl_p2c,
            rel_tbl_c2p-rel_lst_c2p,
            rel_tbl_p2c-rel_lst_p2c
        ) %>%
        as.directed

    V(g)$tbl <- grepl('tbl', V(g)$name)
    V(g)$c2p <- grepl('c2p', V(g)$name)
    E(g)$fun <-
        g %>%
        get.edgelist %>%
        `colnames<-`(c('from', 'to')) %>%
        as_tibble %>%
        pmap_chr(
            function(from, to){
                if(str_sub(from, -3) != str_sub(to, -3)){
                    'swap_relations'
                }else if(str_sub(from, 5, 7) == 'tbl'){
                    'relations_table_to_list'
                }else{
                    'relations_list_to_table'
                }
            }
        )

    E(g)$weight <-
        ifelse(E(g)$fun == 'swap_relations', 1, 2)

    return(g)

}


.ontology_db_variants_graph <- get_ontology_db_variants_graph()


#' Finds the most efficient way to transform ontology relationships into the
#' desired format
#'
#' @param db An ontology database (as produced by \code{\link{obo_parser}}
#'     and accessed by \code{\link{get_db}}.
#' @param tbl Logical: the data structure should be a data frame
#'     (\code{TRUE}) or a list (\code{FALSE}).
#' @param c2p Logical: the data structure should contain child-to-parents
#'     (\code{TRUE}) or parent-to-children (\code{FALSE}) relations.
#'
#'
#' @return A list with the following elements: "operations" character
#'     vector with the function names; "start" name of the starting data
#'     structure; "end" name of the target data structure. Note: if the
#'     target data structure already exists "operations" a zero length
#'     vector.
#'
#' @importFrom magrittr %>%
#' @importFrom igraph V shortest_paths
#' @importFrom dplyr first last
#' @importFrom purrr map_dbl
#'
#' @noRd
ontology_db_transformations <- function(db, tbl, c2p){

    g <- .ontology_db_variants_graph

    to <- (V(g)$tbl == tbl & V(g)$c2p == c2p) %>% which
    from <- V(g)$name %in% names(db) %>% which

    paths <- shortest_paths(g, to, from, mode = 'in', output = 'both')
    idx <- paths$epath %>% map_dbl(function(e){sum(e$weight)}) %>% which.min
    operations <- paths$epath[[idx]]$fun
    start <- paths$vpath[[idx]] %>% last %>% `$`('name')
    end <- paths$vpath[[idx]] %>% first %>% `$`('name')

    list(
        operations = operations,
        start = start,
        end = end
    )

}


#' Access an ontology database
#'
#' Retrieves an ontology database with relations in the desired data
#' structure. The database is automatically loaded and the requested data
#' structure is constructed if necessary. The databases stay loaded up to a
#' certain time period (see the option \code{omnipath.db_lifetime}). Hence
#' the first one of repeated calls to this function might take long and the
#' subsequent ones should be really quick.
#'
#' @param key Character: key of the ontology database. For the available keys
#'     see \code{\link{omnipath_show_db}}.
#' @param rel_tbl Logical: wheter the ontology relations data structure
#'     should be a data frame or a list.
#' @param child_parents Logical: whether the ontology relations should point
#'     from child to parents (\code{TRUE}) or from parent to children
#'     (\code{FALSE}).
#'
#' @examples
#' go <- get_ontology_db('go_basic', child_parents = FALSE)
#'
#' @importFrom magrittr %<>%
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{omnipath_show_db}}}
#'     \item{\code{\link{get_db}}}
#' }
get_ontology_db <- function(key, rel_tbl = TRUE, child_parents = TRUE){

    db <- get_db(key)

    transf <- ontology_db_transformations(db, rel_tbl, child_parents)

    relations <- db[[transf$start]]

    for(op in transf$operations){

        relations %<>% (get(op))

    }

    db[[transf$end]] <- relations
    omnipath.env$db[[key]]$db <- db

    get_db(key)

}


descendants <- function(
    terms,
    db_key = 'go_basic',
    db_param = list(tables = TRUE),
    reload = FALSE
){

    db <- get_db(key = db_key, param = db_param, reload = reload)



}