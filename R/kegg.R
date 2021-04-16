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
#' @param max_expansion Numeric: the maximum number of relations
#'     derived from a single relation record. As one entry might represent
#'     more than one molecular entities, one relation might yield a large
#'     number of relations in the processing. This happens in a combinatorial
#'     way, e.g. if the two entries represent 3 and 4 entities, that results
#'     12 relations. If \code{NULL}, all relations will be expanded.
#'
#' @return A data frame (tibble) of interactions if \code{process} is
#'     \code{TRUE}, otherwise a list with two data frames: "entries" is
#'     a raw table of the entries while "relations" is a table of relations
#'     extracted from the KGML file.
#'
#' @examples
#' tgf_pathway <- kegg_pathway_download('hsa04350')
#' tgf_pathway
#' # # A tibble: 50 x 12
#' #    source target type  effect arrow relation_id kegg_id_source
#' #    <chr>  <chr>  <chr> <chr>  <chr> <chr>       <chr>
#' #  1 51     49     PPrel activ. -->   hsa04350:1  hsa:7040 hsa:.
#' #  2 57     55     PPrel activ. -->   hsa04350:2  hsa:151449 hs.
#' #  3 34     32     PPrel activ. -->   hsa04350:3  hsa:3624 hsa:.
#' #  4 20     17     PPrel activ. -->   hsa04350:4  hsa:4838
#' #  5 60     46     PPrel activ. -->   hsa04350:5  hsa:4086 hsa:.
#' # # . with 45 more rows, and 5 more variables: genesymbol_source <chr>,
#' # #   uniprot_source <chr>, kegg_id_target <chr>,
#' # #   genesymbol_target <chr>, uniprot_target <chr>
#'
#' @importFrom httr add_headers
#' @importFrom rlang exec !!!
#' @importFrom xml2 read_xml xml_find_all xml_attr xml_child
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom tidyr unnest_wider separate_rows
#' @importFrom dplyr mutate row_number
#' @export
kegg_pathway_download <- function(
    pathway_id,
    process = TRUE,
    max_expansion = NULL
){

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
        unnest_wider(entries)

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
        unnest_wider(relations) %>%
        mutate(
            relation_id = sprintf('%s:%d', pathway_id, row_number())
        )

    if(process){

        kegg_process(entries, relations, max_expansion)

    }else{

        list(
            entries = entries,
            relations = relations
        )

    }

}


#' Interactions from KGML
#'
#' Processes KEGG Pathways data extracted from a KGML file. Joins the entries
#' and relations into a single data frame and translates the Gene Symbols to
#' UniProt IDs.
#'
#' @param entries A data frames with entries extracted from a KGML
#'     file by \code{\link{kegg_pathway_download}}.
#' @param relations A data frames with relations extracted from a KGML
#'     file by \code{\link{kegg_pathway_download}}.
#' @param max_expansion Numeric: the maximum number of relations
#'     derived from a single relation record. As one entry might represent
#'     more than one molecular entities, one relation might yield a large
#'     number of relations in the processing. This happens in a combinatorial
#'     way, e.g. if the two entries represent 3 and 4 entities, that results
#'     12 relations. If \code{NULL}, all relations will be expanded.
#'
#' @return A data frame (tibble) of interactions.
#'
#' @examples
#' hsa04350 <- kegg_pathway_download('hsa04350', process = FALSE)
#' tgf_pathway <- kegg_process(hsa04350$entries, hsa04350$relatiions)
#' tgf_pathway
#' # # A tibble: 50 x 12
#' #    source target type  effect arrow relation_id kegg_id_source
#' #    <chr>  <chr>  <chr> <chr>  <chr> <chr>       <chr>
#' #  1 51     49     PPrel activ. -->   hsa04350:1  hsa:7040 hsa:.
#' #  2 57     55     PPrel activ. -->   hsa04350:2  hsa:151449 hs.
#' #  3 34     32     PPrel activ. -->   hsa04350:3  hsa:3624 hsa:.
#' #  4 20     17     PPrel activ. -->   hsa04350:4  hsa:4838
#' #  5 60     46     PPrel activ. -->   hsa04350:5  hsa:4086 hsa:.
#' # # . with 45 more rows, and 5 more variables: genesymbol_source <chr>,
#' # #   uniprot_source <chr>, kegg_id_target <chr>,
#' # #   genesymbol_target <chr>, uniprot_target <chr>
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr separate_rows
#' @importFrom dplyr left_join mutate n rename group_by filter ungroup
#' @export
kegg_process <- function(entries, relations, max_expansion = NULL){

    uniprots <- get_db('up_gs_human')

    entries %<>%
        separate_rows(genesymbol, sep = '[, ]+') %>%
        left_join(
            uniprots %>% rename(uniprot = From),
            by = c('genesymbol' = 'To')
        )

    relations %>%
    left_join(
        entries,
        by = c('source' = 'kgml_id')
    ) %>%
    left_join(
        entries,
        by = c('target' = 'kgml_id'),
        suffix = c('_source', '_target')
    ) %>%
    {`if`(
        is.null(max_expansion),
        .,
        group_by(., relation_id) %>%
        filter(n() <= max_expansion) %>%
        ungroup
    )} %>%
    filter(
        !is.na(uniprot_source) &
        !is.na(uniprot_target)
    )

}