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
#
#  Bio Model Analyzer export: converts path to a BMA motif
#  Author: Ben Hall
#


#' Ends a function where something has gone wrong, printing information about the error
#' @param a string with information about why the error occurred
#' @noRd
#' @importFrom logger log_warn
wrong_input <- function(reason){
    logger::log_warn(reason)
    return(NULL)
}


#' Returns a formatted string describing a BMA interaction between variables
#'
#' @param a unique id, variable ids describing the source and targets, and
#' the edge description
#' @noRd
bma_relationship <- function(id, from, to, type){
    rel <- sprintf(
        '{"Id":%d, "FromVariable":%d, "ToVariable":%d, "Type":"%s"}',
        id, from, to, type
    )
    return(rel)
}


#' Returns a formatted string describing the model parameters of a BMA
#' variable
#'
#' @param a unique id, human readable name (e.g. JAG1), unique variable id,
#' granularity (number of levels) and the formula
#' @noRd
bma_variable_model <- function(id, name, granularity, formula = ""){
    var <- sprintf(
        '{"Name":"%s", "Id":%d, "RangeFrom":0, "RangeTo":%d, "Formula":"%s"}',
        name, id, granularity, formula
    )
    return(var)
}


#' Returns a formatted string describing the layout parameters of a BMA
#' variable
#'
#' @param a unique id, human readable name (e.g. JAG1), granularity (number
#' of levels) and the update formula
#' @noRd
bma_variable_layout <- function(id, name, x, y, description = "") {
    var <- sprintf(
        paste0(
            '{"Id":%d, "Name":"%s", "Type":"Constant", "ContainerId":0, ',
            '"PositionX":%f, "PositionY":%f, "CellX":0, "CellY":0, ',
            '"Angle":0, "Description":"%s"}'
        ),
        id, name, x, y, description
    )
    return(var)
}

#' Returns a string containing the target function of a variable
#'
#' @return Returns either empty string (interpreted as default function), or
#' ranularity - activity of upstream inhibitor
#' @param bool stating whether the interaciton is an inhibition, granularity
#' of variables (number of levels), and source of interaction
#' @noRd
bma_formula <- function(inhibitor, granularity, upstream){
    f <- ifelse(inhibitor, sprintf("%d-var(%s)", granularity, upstream), "")
    return(f)
}


#' Returns a string describing the evidence behind an interaction
#'
#' Contains all interaction types with a simple descriptor and PMIDs
#' @param takes an edge from omnipath "e", and optionally the name of the
#' upstream variable ("incoming")
#' @noRd
bma_description <- function(e, incoming = ""){
    sign <- ifelse(e$is_stimulation == 1,
            ifelse(e$is_inhibition == 1, "Mixed", "Activator"),
            ifelse(e$is_inhibition == 1, "Inhibitor", "Unknown"))
    refs <- paste(unlist(e$references), sep = '', collapse = ', ')
    incoming <- ifelse(incoming == "", "", paste("", incoming, "", sep = " "))
    return(sprintf("%s%s-PMID:%s.", incoming, sign, refs))
    }


#' Prints a BMA motif to the screen from a sequence of edges, which can be
#' copy/pasted into the BMA canvas
#'
#' Intended to parallel print_path_es
#' @param takes an sequence of edges, a graph, and a granularity
#' @export
#' @importFrom igraph tail_of head_of
bma_motif_es <- function(edge_seq, G, granularity = 2){
    if(length(edge_seq) == 0) {
        wrong_input("BMA motif: empty path")
    }
    #Process-
    ## Create list of variables
    ## Create layout of variables
    ## Create list of links
    ## Print format as follows (x is a string, xN is an integer, xF is a float)
    ### {"Model":     {"Name": "Omnipath motif",
    ###              "Variables": [{"Name":"x", "Id":xN, "RangeFrom" = 0,
    ###                             "RangeTo" = granularity}, ...]
    ###             "Relationships": [{"Id":xN, "FromVariable":xN,
    ###                            "ToVariable":xN, "Type":"Activator"}, ...]
    ###         }
    ###  "Layout":     {"Variables": [{"Id":xN, "Name":"x", "Type":"Constant",
    ###                                "ContainerId":0, "PositionX":xF,
    ###                                "PositionY":xF, "CellX":0, "CellY":0,
    ###                                "Angle":0, "Description":""}, ...]
    ###            "Containers":[]
    ###            }
    ### }

    #Code for identifying sign
    #signs <- ifelse(edge_seq$is_stimulation == 1,
    #    ifelse(edge_seq$is_inhibition == 1, "( + /-)", "( + )"),
    #    ifelse(edge_seq$is_inhibition == 1, "( - )", "( ? )"))
    #interaction <- paste0(" == ", signs, " == >")

    #relationships <- ifelse(edge_seq$is_stimulation == 1,
    #    ifelse(edge_seq$is_inhibition == 1, "Activator", "Activator"),
    #    ifelse(edge_seq$is_inhibition == 1, "Inhibitor",
    #    return(wrongInput("\nUnsigned input graph\n")))
    sources <- tail_of(G, edge_seq)$name
    variable_names <- c(sources, head_of(G, edge_seq)$name[length(edge_seq)])
    var_num <- length(variable_names)

    positions <- vector('list', var_num)
    variables <- vector('list', var_num)
    relationships <- vector('list', var_num - 1)


    formula = ""
    description = ""
    x = 125
    y = 140
    for (i in seq_along(variable_names))
        {
    v <- bma_variable_model(i, variable_names[i], granularity, formula)
    p <- bma_variable_layout(i, variable_names[i], x, y, description)
    if (i < var_num){
        # Simplified sign- if inhibition, inhibitor,
        # else (activator/mixed/unknown) activation
        r <- bma_relationship(
            i + var_num, i, i + 1,
            ifelse(edge_seq[i]$is_inhibition == 1, "Inhibitor", "Activator")
        )
        relationships[[i]] <- r
        formula <- bma_formula(
            (edge_seq[i]$is_inhibition == 1),
            granularity,
            variable_names[i]
        )
        description <- bma_description(edge_seq[i])
    }
    positions[[i]] <- p
    variables[[i]] <- v

    x <- x + 86.6025404
    ymod <- ifelse(i %% 2 == 0, 50, -50)
    y <- y + ymod
    }
    result <- sprintf(
        paste0(
            '{"Model": {"Name": "Omnipath motif", "Variables":[%s], ',
            '"Relationships":[%s]}, "Layout":{"Variables":[%s], ',
            '"Containers":[]}}\n'
        ),
        paste(variables, sep = '', collapse = ', '),
        paste(relationships, sep = '', collapse = ', '),
        paste(positions, sep = '', collapse = ', ')
    )
    cat(result)

}


#' Prints a BMA motif to the screen from a sequence of nodes, which can be
#' copy/pasted into the BMA canvas
#'
#' Intended to parallel print_path_vs
#' @param takes an sequence of nodes, and a granularity
#' @export
#' @importFrom logger log_warn
#' @importFrom igraph E
bma_motif_vs <- function(node_seq, G){

    if(length(node_seq) == 0){
        logger::log_warn("BMA motif: empty path")
        return(invisible(NULL))
    }
    node_seq_names <- unique_node_seq(node_seq)
    for(i in seq(node_seq_names)){
        print(
            paste0(
                "pathway ", i, ": ",
                paste(node_seq_names[[i]], collapse = " -> ")
            )
        )
        edge_set <- c()
        for(j in 2:length(node_seq_names[[i]])){
            edge_set <- c(edge_set, E(G)[node_seq_names[[i]][[j-1]]  %->%
                node_seq_names[[i]][[j]]])
        }
        bma_motif_es(E(G)[edge_set], G)
    }
}