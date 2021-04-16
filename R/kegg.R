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


#' Download the KEGG Pathways database
#'
#' Downloads all pathway diagrams in the KEGG Pathways database in KGML
#' format and processes the XML to extract the interactions.
#'
#' @return A data frame (tibble) of interactions.
#'
#' @export
kegg_pathways_download <- function(){

    pathways <- kegg_pathway_list()



}


#' List of KEGG pathways
#'
#' Retrieves a list of available KEGG pathways.
#'
#' @return Data frame of pathway names and identifiers.
#'
#' @examples
#' kegg_pws <- kegg_pathway_list()
#' kegg_pws
#' # # A tibble: 521 x 2
#' #    id       name
#' #    <chr>    <chr>
#' #  1 map01100 Metabolic pathways
#' #  2 map01110 Biosynthesis of secondary metabolites
#' #  3 map01120 Microbial metabolism in diverse environments
#' #  4 map01200 Carbon metabolism
#' #  5 map01210 2-Oxocarboxylic acid metabolism
#' #  6 map01212 Fatty acid metabolism
#' #  7 map01230 Biosynthesis of amino acids
#' # # . with 514 more rows
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom purrr map keep
#' @importFrom xml2 read_html xml_find_all xml_attr xml_text
#' @importFrom stringr str_extract
#' @importFrom tidyr unnest_wider
#' @export
kegg_pathway_list <- function(){

    # NSE vs. R CMD check workaround
    pathway <- NULL

    path <-
        download_to_cache(
            url_key = 'omnipath.kegg_list_url',
            ext = 'html'
        )

    path %>%
    read_html %>%
    xml_find_all('.//a') %>%
    map(
        function(a){
            list(
                id =
                    a %>%
                    xml_attr('href') %>%
                    str_extract('((?:hsa|map)[0-9]+)'),
                name = xml_text(a)
            )
        }
    ) %>%
    keep(
        function(a){!is.na(a$id)}
    ) %>%
    tibble(pathway = .) %>%
    unnest_wider(pathway) %>%
    copy_source_attrs(path, resource = 'KEGG') %T>%
    load_success()

}


#' Download one KEGG pathway
#'
#' Downloads one pathway diagram from the KEGG Pathways database in KGML
#' format and processes the XML to extract the interactions.
#'
#' @param pathway_id Character: a KEGG pathway identifier, for example
#'     "hsa04350".
#' @param process Logical: process the data or return it in raw format.
#'     processing means joining the entries and relations into a single
#'     data frame and adding UniProt IDs.
#'
#' @return A data frame (tibble) of interactions if \code{process} is
#'     \code{TRUE}, otherwise a list with two data frames: "entries" is
#'     a raw table of the entries while "relations" is a table of relations
#'     extracted from the KGML file.
#'
#' @importFrom httr add_headers
#' @importFrom rlang exec !!!
#' @importFrom xml2 read_xml xml_find_all xml_attr xml_child
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom tidyr unnest_wider separate_rows
#' @importFrom dplyr left_join mutate n
#' @export
kegg_pathway_download <- function(pathway_id, process = TRUE){

    # NSE vs. R CMD check workaround
    genesymbol <- NULL

    req_hdrs <- list(
        Referer = paste0(
            'http://www.genome.jp/kegg-bin/show_pathway',
            '?map=hsa04710&show_description=show'
        )
    )

    path <-
        download_to_cache(
            url_key = 'omnipath.kegg_kgml_url',
            url_param = list(pathway_id),
            http_param = exec(add_headers, !!!req_hdrs),
            ext = 'xml'
        )

    kgml <- path %>% read_xml

    entries <-
        kgml %>%
        xml_find_all('.//entry') %>%
        map(
            function(e){
                list(
                    kgml_id = xml_attr(e, 'id'),
                    kegg_id = xml_attr(e, 'name'),
                    genesymbol = e %>% xml_child %>% xml_attr('name')
                )
            }
        ) %>%
        tibble(entries = .) %>%
        unnest_wider(entries) %>%
        separate_rows(genesymbol, sep = '[, ]+')

    relations <-
        kgml %>%
        xml_find_all('.//relation') %>%
        map(
            function(r){
                subtype <- xml_child(r)
                list(
                    source = xml_attr(r, 'entry1'),
                    target = xml_attr(r, 'entry2'),
                    type = xml_attr(r, 'type'),
                    effect = xml_attr(subtype, 'name'),
                    arrow = xml_attr(subtype, 'value')
                )
            }
        ) %>%
        tibble(relations = .) %>%
        unnest_wider(relations)

    if(process){

        relations %>%
        mutate(relation_id = sprintf('%s:%d', pathway_id, n())) %>%
        left_join(
            entries,
            by = c('source' = 'kgml_id')
        ) %>%
        left_join(
            entries,
            by = c('target' = 'kgml_id'),
            suffix = c('_source', '_target')
        )

    }else{

        list(
            entries = entries,
            relations = relations
        )

    }

}