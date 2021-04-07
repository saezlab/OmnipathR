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


#' Graph from a table of ontology relations
#'
#' @param relations A data frame of ontology relations (the "relations"
#'     element of the list returned by \code{\link{obo_parser}} in case
#'     its argument `tables` is \code{TRUE}).
#'
#' @return The relations converted to an igraph graph object.
#'
#' @details
#' By default the relations point from child to parents, the edges in the
#' graph will be of the same direction. Use \code{\link{swap_relations}}
#' on the data frame to reverse the direction.
#'
#' @examples
#' go <- get_db('go_basic')
#' go_graph <- relations_table_to_graph(go$relations)
#'
#' @importFrom logger log_trace
#' @importFrom magrittr %>%
#' @importFrom dplyr rename
#' @importFrom tidyr unnest
#' @importFrom igraph graph_from_data_frame
#' @export
relations_table_to_graph <- function(relations){

    log_trace('Converting ontology relations from table to graph.')

    relations %>%
    rename(term2 = 3) %>%
    unnest(term2) %>%
    select(term, term2, relation) %>%
    graph_from_data_frame

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
#' @importFrom igraph graph_from_literal as.directed V E V<- E<-
#' @importFrom igraph get.edgelist delete_edges
#' @importFrom purrr pmap_chr
#' @importFrom stringr str_sub
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble as_tibble
#'
#' @noRd
get_ontology_db_variants_graph <- function(){

    g <-
        graph_from_literal(
            rel_tbl_c2p-rel_tbl_p2c,
            rel_tbl_c2p-rel_lst_c2p,
            rel_tbl_p2c-rel_lst_p2c,
            rel_tbl_c2p-rel_gra_c2p,
            rel_tbl_p2c-rel_gra_p2c
        ) %>%
        as.directed

    V(g)$fmt <- str_sub(V(g)$name, 5, 7)
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
                }else if(
                    str_sub(from, 5, 7) == 'tbl' &&
                    str_sub(to, 5, 7) == 'lst'
                ){
                    'relations_table_to_list'
                }else if(
                    str_sub(from, 5, 7) == 'lst' &&
                    str_sub(to, 5, 7) == 'tbl'
                ){
                    'relations_list_to_table'
                }else if(
                    str_sub(from, 5, 7) == 'tbl' &&
                    str_sub(to, 5, 7) == 'gra'
                ){
                    'relations_table_to_graph'
                }else(
                    NA
                )
            }
        )

    g %<>% delete_edges(E(g)$fun %>% is.na %>% which)

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
#' @param fmt Character: the data structure should be 1) "tbl" a data frame,
#'     2) "lst" a list or 3) "gra" a graph.
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
ontology_db_transformations <- function(db, fmt, c2p){

    g <- .ontology_db_variants_graph

    to <- (V(g)$fmt == fmt & V(g)$c2p == c2p) %>% which
    from <- V(g)$name %in% names(db) %>% which

    paths <- shortest_paths(g, to, from, mode = 'in', output = 'both')
    idx <- paths$epath %>% map_dbl(function(e){sum(e$weight)}) %>% which.min

    operations <- paths$epath[[idx]]$fun %>% rev
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
#' @param rel_fmt Character: the data structure of the ontology relations.
#'     Posible values are 1) "tbl" a data frame, 2) "lst" a list or 3) "gra"
#'     a graph.
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
get_ontology_db <- function(key, rel_fmt = 'tbl', child_parents = TRUE){

    db <- get_db(key)

    transf <- ontology_db_transformations(db, rel_fmt, child_parents)

    relations <- db[[transf$start]]

    for(op in transf$operations){

        relations %<>% (get(op))

    }

    db[[transf$end]] <- relations
    omnipath.env$db[[key]]$db <- db

    get_db(key)

}


#' All nodes of a subtree starting from the selected nodes
#'
#' Starting from the selected nodes, recursively walks the ontology tree
#' until it reaches either the root or leaf nodes. Collects all visited
#' nodes.
#'
#' @param terms Character vector of ontology term IDs or names. A mixture of
#'     IDs and names can be provided.
#' @param ancestors Logical: if \code{FALSE} the ontology tree is traversed
#'     towards the leaf nodes; if \code{TRUE}, the tree is traversed until
#'     the root. The former returns the ancestors (parents), the latter
#'     the descendants (children).
#' @param db_key Character: key to identify the ontology database. For the
#'     available keys see \code{\link{omnipath_show_db}}.
#' @param ids Logical: whether to return IDs or term names.
#' @param relations Character vector of ontology relation types. Only these
#'     relations will be used.
#'
#' @return Character vector of ontology IDs. If the input terms are all
#'     leaves or roots \code{NULL} is returned. The starting nodes won't
#'     be included in the result unless they fall onto the traversal path
#'     from other nodes.
#'
#' @details
#' Note: this function relies on the database manager, the first call might
#' take long because of the database load process. Subsequent calls within
#' a short period should be faster. See \code{\link{get_ontology_db}}.
#'
#' @examples
#' walk_ontology_tree(c('GO:0006241', 'GO:0044211'))
#' # [1] "GO:0006139" "GO:0006220" "GO:0006221" "GO:0006241" "GO:0006725"
#' # [6] "GO:0006753" "GO:0006793" "GO:0006796" "GO:0006807" "GO:0008150"
#' # ... (truncated)
#' walk_ontology_tree(c('GO:0006241', 'GO:0044211'), ancestors = FALSE)
#' # [1] "GO:0044210" "GO:0044211"
#' walk_ontology_tree(
#'     c('GO:0006241', 'GO:0044211'),
#'     ancestors = FALSE,
#'     ids = FALSE
#' )
#' # [1] "'de novo' CTP biosynthetic process" "CTP salvage"
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom igraph delete_edges
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{omnipath_show_db}}}
#'     \item{\code{\link{get_ontology_db}}}
#' }
walk_ontology_tree <- function(
    terms,
    ancestors = TRUE,
    db_key = 'go_basic',
    ids = TRUE,
    method = 'gra',
    relations = c(
        'is_a', 'part_of', 'occurs_in', 'regulates',
        'positively_regulates', 'negatively_regulates'
    )
){

    db <- get_ontology_db(
        key = db_key,
        rel_fmt = method,
        child_parents = ancestors
    )

    rel_key <- `if`(ancestors, 'rel_gra_c2p', 'rel_gra_p2c')

    rel <- db[[rel_key]]

    if(method == 'gra'){

        rel %<>%
            delete_edges(
                E(rel)$relation %in% relations %>% `!` %>% which
            )

    }

    terms %<>% ontology_ensure_id(db_key = db_key)

    fun <- `if`(
        method == 'gra',
        .walk_ontology_tree_graph,
        .walk_ontology_tree
    )

    fun(terms, rel, relations) %>%
    ontology_name_id(ids = ids, db_key = db_key)

}


#' See \code{link{walk_ontology_tree}}
#'
#' @importFrom igraph V ego
#' @noRd
.walk_ontology_tree_graph <- function(terms, rel, relations){

    V(rel)[terms, na_ok = TRUE] %>%
    na.omit %>%
    {ego(rel, order = 16, nodes = ., mode = 'out')} %>%
    unlist %>%
    names %>%
    unique

}


#' See \code{link{walk_ontology_tree}}
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr map
#' @noRd
.walk_ontology_tree <- function(
    terms,
    rel,
    relations
){

    if(length(terms) == 1L){

        if(!terms %in% names(rel)){

            return(NULL)

        }else{

            related_terms <- rel[[terms]][relations] %>% unlist
            c(
                related_terms,
                rel %>% names %>% intersect(related_terms) %>%
                .walk_ontology_tree(rel = rel, relations = relations)
            ) %>%
            unname %>%
            unique

        }

    }else{

        terms %>% map(
            .walk_ontology_tree,
            rel = rel,
            relations = relations
        ) %>%
        unlist %>%
        unique

    }

}


#' All descendants in the ontology tree
#'
#' Starting from the selected nodes, recursively walks the ontology tree
#' until it reaches the leaf nodes. Collects all visited nodes, which are
#' the descendants (children) of the starting nodes.
#'
#' @param terms Character vector of ontology term IDs or names. A mixture of
#'     IDs and names can be provided.
#' @param db_key Character: key to identify the ontology database. For the
#'     available keys see \code{\link{omnipath_show_db}}.
#' @param ids Logical: whether to return IDs or term names.
#' @param relations Character vector of ontology relation types. Only these
#'     relations will be used.
#'
#' @return Character vector of ontology IDs. If the input terms are all
#'     leaves \code{NULL} is returned. The starting nodes won't be included
#'     in the result unless some of them are descendants of other starting
#'     nodes.
#'
#' @details
#' Note: this function relies on the database manager, the first call might
#' take long because of the database load process. Subsequent calls within
#' a short period should be faster. See \code{\link{get_ontology_db}}.
#'
#' @examples
#' descendants('GO:0005035', ids = FALSE)
#' # [1] "tumor necrosis factor-activated receptor activity"
#' # [2] "TRAIL receptor activity"
#' # [3] "TNFSF11 receptor activity"
#'
#' @importFrom rlang exec !!!
#' @export
descendants <- function(
    terms,
    db_key = 'go_basic',
    ids = TRUE,
    relations = c(
        'is_a', 'part_of', 'occurs_in', 'regulates',
        'positively_regulates', 'negatively_regulates'
    )
){

    exec(walk_ontology_tree, ancestors = FALSE, !!!as.list(environment()))

}


#' All ancestors in the ontology tree
#'
#' Starting from the selected nodes, recursively walks the ontology tree
#' until it reaches the root. Collects all visited nodes, which are the
#' ancestors (parents) of the starting nodes.
#'
#' @param terms Character vector of ontology term IDs or names. A mixture of
#'     IDs and names can be provided.
#' @param db_key Character: key to identify the ontology database. For the
#'     available keys see \code{\link{omnipath_show_db}}.
#' @param ids Logical: whether to return IDs or term names.
#' @param relations Character vector of ontology relation types. Only these
#'     relations will be used.
#'
#' @return Character vector of ontology IDs. If the input terms are all
#'     root nodes, \code{NULL} is returned. The starting nodes won't be
#'     included in the result unless some of them are ancestors of other
#'     starting nodes.
#'
#' @details
#' Note: this function relies on the database manager, the first call might
#' take long because of the database load process. Subsequent calls within
#' a short period should be faster. See \code{\link{get_ontology_db}}.
#'
#' @examples
#' ancestors('GO:0005035', ids = FALSE)
#' # [1] "molecular_function"
#' # [2] "transmembrane signaling receptor activity"
#' # [3] "signaling receptor activity"
#' # [4] "molecular transducer activity"
#'
#' @importFrom rlang exec !!!
#' @export
ancestors <- function(
    terms,
    db_key = 'go_basic',
    ids = TRUE,
    relations = c(
        'is_a', 'part_of', 'occurs_in', 'regulates',
        'positively_regulates', 'negatively_regulates'
    )
){

    exec(walk_ontology_tree, ancestors = TRUE, !!!as.list(environment()))

}


#' Translate between ontology IDs and names
#'
#' Makes sure that the output contains only valid IDs or term names. The
#' input can be a mixture of IDs and names. The order of the input won't
#' be preserved in the output.
#'
#' @param terms Character: ontology IDs or term names.
#' @param ids Logical: the output should contain IDs or term names.
#' @param db_key Character: key to identify the ontology database. For the
#'     available keys see \code{\link{omnipath_show_db}}.
#'
#' @return Character vector of ontology IDs or term names.
#'
#' @examples
#' ontology_name_id(c('mitochondrion inheritance', 'reproduction'))
#' # [1] "GO:0000001" "GO:0000003"
#' ontology_name_id(c('GO:0000001', 'reproduction'), ids = FALSE)
#' # [1] "mitochondrion inheritance" "reproduction"
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull
#' @importFrom rlang !! sym
#' @importFrom tibble tibble
#' @export
ontology_name_id <- function(terms, ids = TRUE, db_key = 'go_basic'){

    db <- get_db(db_key)
    to_col <- `if`(ids, 'term', 'name')
    from_col <- `if`(ids, 'name', 'term')

    id_name <- `if`(
        inherits(db$names, 'data.frame'),
        db$names,
        db$names %>% {tibble(term = names(.), name = unlist(.))}
    )

    id_name %>%
    filter(!!sym(from_col) %in% terms) %>%
    pull(!!sym(to_col)) %>%
    union(
        id_name %>%
        pull(!!sym(to_col)) %>%
        intersect(terms)
    )

}


#' Only ontology IDs
#'
#' Converts a mixture of ontology IDs and names to only IDs. If an element
#' of the input is missing from the chosen ontology it will be dropped.
#' This can happen if the ontology is a subset (slim) version, but also if
#' the input is not a valid ID or name.
#'
#' @param terms Character: ontology IDs or term names.
#' @param db_key Character: key to identify the ontology database. For the
#'     available keys see \code{\link{omnipath_show_db}}.
#'
#' @return Character vector of ontology IDs.
#'
#' @examples
#' ontology_ensure_id(c('mitochondrion inheritance', 'GO:0001754'))
#' # [1] "GO:0000001" "GO:0001754"
#'
#' @export
ontology_ensure_id <- function(terms, db_key = 'go_basic'){

    ontology_name_id(terms, db_key = db_key)

}


#' Only ontology term names
#'
#' Converts a mixture of ontology IDs and names to only names. If an element
#' of the input is missing from the chosen ontology it will be dropped.
#' This can happen if the ontology is a subset (slim) version, but also if
#' the input is not a valid ID or name.
#'
#' @param terms Character: ontology IDs or term names.
#' @param db_key Character: key to identify the ontology database. For the
#'     available keys see \code{\link{omnipath_show_db}}.
#'
#' @return Character vector of ontology term names.
#'
#' @examples
#' ontology_ensure_name(c('reproduction', 'GO:0001754', 'foo bar'))
#' # [1] "eye photoreceptor cell differentiation" "reproduction"
#'
#' @export
ontology_ensure_name <- function(terms, db_key = 'go_basic'){

    ontology_name_id(terms, ids = FALSE, db_key = db_key)

}


#' Looks like an ontology ID
#'
#' Tells if the input has the typical format of ontology IDs, i.e. a code
#' of capital letters, a colon, followed by a numeric code.
#'
#' @param terms Character vector with strings to check.
#'
#' @return A logical vector with the same length as the input.
#'
#' @examples
#' is_ontology_id(c('GO:0000001', 'reproduction'))
#' # [1]  TRUE FALSE
#'
#' @export
is_ontology_id <- function(terms){

    grepl('^[A-Z]+:[0-9]+$', terms)

}