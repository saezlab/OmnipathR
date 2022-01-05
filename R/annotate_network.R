#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2022
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


#' Network interactions with annotations
#'
#' Annotations are often useful in a network context, e.g. one might want to
#' label the interacting partners by their pathway membership. This function
#' takes a network data frame and joins an annotation data frame from both
#' the left and the right side, so both the source and target molecular
#' entities will be labeled by their annotations. If one entity has many
#' annotations these will yield many rows, hence the interacting pairs won't
#' be unique across the data frame any more. Also if one entity has really
#' many annotations the resulting data frame might be huge, we recommend to
#' be careful with that. Finally, if you want to do the same but with
#' intercell annotations, there is the \code{\link{import_intercell_network}}
#' function.
#'
#' @param network Behaviour depends on type: if list, will be passed as
#'     arguments to \code{\link{import_omnipath_interactions}} to obtain a
#'     network data frame; if a data frame or tibble, it will be used as a
#'     network data frame; if a character vector, will be assumed to be a
#'     set of resource names and interactions will be queried from these
#'     resources.
#' @param annot Either the name of an annotation resource (for a list of
#'     available resources call \code{\link{get_annotation_resources}}), or
#'     an annotation data frame. If the data frame contains more than one
#'     resources, only the first one will be used.
#' @param ... Column names selected from the annotation data frame (passed
#'     to \code{dplyr::select}, if empty all columns will be selected.)
#'
#' @return A data frame of interactions with annotations for both interacting
#'     entities.
#'
#' @examples
#' signalink_with_pathways <-
#'     annotated_network('SignaLink3', 'SignaLink_pathway')
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom checkmate assert_data_frame
#' @importFrom rlang enexprs !!!
#' @importFrom dplyr select left_join
#' @export
annotated_network <- function(
    network = NULL,
    annot = NULL,
    ...
){

    # NSE vs. R CMD check workaround
    entity_type <- uniprot <- genesymbol <- NULL

    annot_sel <- enexprs(...)

    network %<>%
    {`if`(is.character(.), list(resources = .), .)} %>%
    {`if`(
        just_a_list(.),
        do.call(import_omnipath_interactions, .),
        .
    )} %>%
    assert_data_frame

    annot %<>%
    {`if`(
        is.character(.),
        import_omnipath_annotations(resources = .),
        .
    )} %>%
    {`if`('record_id' %in% names(.), pivot_annotations(.), .)} %>%
    {`if`(just_a_list(.), .[[1]], .)} %>%
    assert_data_frame %>%
    {`if`(
        length(annot_sel) > 0,
        select(., uniprot, !!!annot_sel),
        select(., -genesymbol, -entity_type)
    )}

    network %>%
    left_join(annot, by = c('source' = 'uniprot')) %>%
    left_join(
        annot,
        by = c('target' = 'uniprot'),
        suffix = c('_source', '_target')
    )

}