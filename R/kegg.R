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


#' Download the KEGG Pathways database
#'
#' Downloads all pathway diagrams in the KEGG Pathways database in KGML
#' format and processes the XML to extract the interactions.
#'
#' @param max_expansion Numeric: the maximum number of relations
#'     derived from a single relation record. As one entry might represent
#'     more than one molecular entities, one relation might yield a large
#'     number of relations in the processing. This happens in a combinatorial
#'     way, e.g. if the two entries represent 3 and 4 entities, that results
#'     12 relations. If \code{NULL}, all relations will be expanded.
#' @param simplify Logical: remove KEGG's internal identifiers and the
#'     pathway annotations, keep only unique interactions with direction
#'     and effect sign.
#'
#' @return A data frame (tibble) of interactions.
#'
#' @examples
#' \dontrun{
#' kegg_pw <- kegg_pathways_download(simplify = TRUE)
#' kegg_pw
#' # # A tibble: 6,765 x 6
#' #    uniprot_source uniprot_target type  effect genesymbol_source
#' #    <chr>          <chr>          <chr> <chr>  <chr>
#' #  1 Q03113         Q15283         PPrel activ. GNA12
#' #  2 Q9Y4G8         P62070         PPrel activ. RAPGEF2
#' #  3 Q13972         P62070         PPrel activ. RASGRF1
#' #  4 O95267         P62070         PPrel activ. RASGRP1
#' #  5 P62834         P15056         PPrel activ. RAP1A
#' # # . with 6,760 more rows, and 1 more variable: genesymbol_target <chr>
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom progress progress_bar
#' @importFrom dplyr filter mutate bind_rows
#' @importFrom purrr pmap
#' @importFrom stringr str_pad str_trunc
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{kegg_pathway_list}}}
#'     \item{\code{\link{kegg_process}}}
#'     \item{\code{\link{kegg_pathway_download}}}
#' }
kegg_pathways_download <- function(max_expansion = NULL, simplify = FALSE){

    # NSE vs. R CMD check workaround
    id <- NULL

    pathways <-
        kegg_pathway_list() %>%
        filter(
            !startsWith(id, 'map')
        )

    pb <- progress_bar$new(
        total = pathways %>% nrow,
        format = paste0(
            '  KEGG Pathways: :pw ',
            '[:bar] :current/:total, eta: :eta'
        )
    )

    pathways  %>%
    pmap(
        function(id, name){
            display_name <- name %>% str_trunc(22) %>% str_pad(22, 'right')
            pb$tick(tokens = list(pw = display_name))
            kegg_pathway_download(
                id,
                max_expansion = max_expansion,
                simplify = simplify
            ) %>%
            null_or_call(mutate, pathway = name, pathway_id = id)
        }
    ) %>%
    bind_rows %>%
    kegg_simplify(simplify = simplify)

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
#' @importFrom tibble tibble
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{kegg_process}}}
#'     \item{\code{\link{kegg_pathway_download}}}
#'     \item{\code{\link{kegg_pathways_download}}}
#'     \item{\code{\link{kegg_open}}}
#'     \item{\code{\link{kegg_picture}}}
#'     \item{\code{\link{kegg_info}}}
#' }
kegg_pathway_list <- function(){

    # NSE vs. R CMD check workaround
    pathway <- NULL

    path <-
        download_to_cache(
            url_key = 'kegg_list',
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
#' @param simplify Logical: remove KEGG's internal identifiers and the
#'     pathway annotations, keep only unique interactions with direction
#'     and effect sign.
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
#' @importFrom rlang exec !!! %||%
#' @importFrom xml2 read_xml xml_find_all xml_attr xml_child xml_children
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom tidyr unnest_wider separate_rows
#' @importFrom dplyr mutate row_number
#' @importFrom tibble tibble
#' @importFrom logger log_error log_trace
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{kegg_process}}}
#'     \item{\code{\link{kegg_pathways_download}}}
#'     \item{\code{\link{kegg_pathway_list}}}
#' }
kegg_pathway_download <- function(
    pathway_id,
    process = TRUE,
    max_expansion = NULL,
    simplify = FALSE
){

    # NSE vs. R CMD check workaround
    genesymbol <- NULL

    req_hdrs <- list(
        Referer = paste0(
            'http://www.genome.jp/kegg-bin/show_pathway',
            '?map=hsa04710&show_description=show'
        )
    )

    kgml <- tryCatch(

        {

            download_to_cache(
                url_key = 'kegg_kgml',
                url_param = list(pathway_id),
                http_param = exec(add_headers, !!!req_hdrs),
                ext = 'xml'
            ) %>%
            read_xml

        },

        error = function(err){

            log_error(
                paste0(
                    'Failed to download KEGG pathway `%s` (URL: %s). ',
                    'Unless you provided an invalid pathway_id or there ',
                    'is a problem with your network connection, this is a ',
                    'known issue on Windows and we will try to fix it as ',
                    'soon as possible. The original error message was: %s'
                ),
                pathway_id,
                get_url(key = 'kegg_kgml', param = list(pathway_id)),
                err
            )

            return(NULL)

        }

    )

    if(is.null(kgml)){

        log_trace('Exiting `kegg_pathway_download` due to previous error.')

        result <- `if`(
            process,
            tibble(),
            list(entries = tibble(), relations = tibble())
        )

        return(result)

    }

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
                subtype <- xml_children(r)[1]
                list(
                    source = xml_attr(r, 'entry1'),
                    target = xml_attr(r, 'entry2'),
                    type = xml_attr(r, 'type'),
                    effect =
                        null_or_call(subtype, xml_attr, 'name') %||%
                        'unknown',
                    arrow =
                        null_or_call(subtype, xml_attr, 'value') %||%
                        '---'
                )
            }
        ) %>%
        tibble(relations = .) %>%
        unnest_wider(relations) %>%
        mutate(
            relation_id = sprintf('%s:%d', pathway_id, row_number())
        )

    if(process){

        kegg_process(entries, relations, max_expansion, simplify)

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
#' @param simplify Logical: remove KEGG's internal identifiers and the
#'     pathway annotations, keep only unique interactions with direction
#'     and effect sign.
#'
#' @return A data frame (tibble) of interactions. In rare cases when a
#'     pathway doesn't contain any relation, returns \code{NULL}.
#'
#' @examples
#' hsa04350 <- kegg_pathway_download('hsa04350', process = FALSE)
#' tgf_pathway <- kegg_process(hsa04350$entries, hsa04350$relations)
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
#' @seealso \itemize{
#'     \item{\code{\link{kegg_pathway_download}}}
#'     \item{\code{\link{kegg_pathways_download}}}
#'     \item{\code{\link{kegg_pathway_list}}}
#' }
kegg_process <- function(
    entries,
    relations,
    max_expansion = NULL,
    simplify = FALSE
){

    # NSE vs. R CMD check workaround
    From <- uniprot_source <- uniprot_target <- relation_id <- NULL

    if(!nrow(relations)){

        return(NULL)

    }

    uniprots <- get_db('up_gs', organism = 9606)

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
    filter(
        !is.na(uniprot_source) &
        !is.na(uniprot_target)
    ) %>%
    {`if`(
        is.null(max_expansion),
        .,
        group_by(., relation_id) %>%
        filter(n() <= max_expansion) %>%
        ungroup
    )} %>%
    kegg_simplify(simplify = simplify)

}


#' Keep only the most important columns of KEGG data
#'
#' Removes KEGG's internal identifiers and the pathway annotations, keeps
#' only unique interactions with direction and effect sign.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct
#' @noRd
kegg_simplify <- function(tbl, simplify = TRUE){

    # NSE vs. R CMD check workaround
    uniprot_source <- uniprot_target <- type <- effect <-
    genesymbol_source <- genesymbol_target <- NULL

    `if`(
        simplify,
        tbl %>%
            select(
                uniprot_source,
                uniprot_target,
                type,
                effect,
                genesymbol_source,
                genesymbol_target
            ) %>%
            distinct,
        tbl
    )

}


#' Open a KEGG Pathway diagram in the browser
#'
#' @param pathway_id Character: a KEGG Pathway identifier, e.g. "hsa04710".
#'     For a complete list of IDs see \code{\link{kegg_pathway_list}}.
#'
#' @return Returns \code{NULL}.
#'
#' @details
#' To open URLs in the web browser the "browser" option must to be set to a
#' a valid executable. You can check the value of this option by
#' \code{getOption("browser")}. If your browser is firefox and the executable
#' is located in the system path, you can set the option to point to it:
#' \code{options(browser = "firefox")}. To make it a permanent setting, you
#' can also include this in your \code{.Rprofile} file.
#'
#' @examples
#' if(any(getOption('browser') != '')) kegg_open('hsa04710')
#'
#' @importFrom magrittr %>%
#' @importFrom utils browseURL
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{kegg_pathway_list}}}
#'     \item{\code{\link{kegg_picture}}}
#'     \item{\code{\link{kegg_info}}}
#' }
kegg_open <- function(pathway_id){

    'kegg_pathway' %>%
    url_parser(url_param = list(pathway_id)) %>%
    browseURL

}


#' Download a pathway diagram as a picture
#'
#' Downloads a KEGG Pathway diagram as a PNG image.
#'
#' @param pathway_id Character: a KEGG Pathway identifier, e.g. "hsa04710".
#'     For a complete list of IDs see \code{\link{kegg_pathway_list}}.
#' @param path Character: save the image to this path. If \code{NULL}, the
#'     image will be saved in the current directory under the name
#'     \code{<pathway_id>.png}.
#'
#' @return Invisibly returns the path to the downloaded file.
#'
#' @examples
#' kegg_picture('hsa04710')
#' kegg_picture('hsa04710', path = 'foo/bar')
#' kegg_picture('hsa04710', path = 'foo/bar/circadian.png')
#'
#' @importFrom magrittr %<>% %>%
#' @export
#' @seealso \code{\link{kegg_pathway_list}}
#' @seealso \itemize{
#'     \item{\code{\link{kegg_pathway_list}}}
#'     \item{\code{\link{kegg_open}}}
#'     \item{\code{\link{kegg_info}}}
#' }
kegg_picture <- function(pathway_id, path = NULL){

    imgname <- sprintf('%s.png', pathway_id)

    path %<>%
        if_null('.') %>%
        {`if`(endsWith(., '.png'), ., file.path(., imgname))}

    path %>% dirname %>% dir.create(recursive = TRUE, showWarnings = FALSE)

    'kegg_pw_png' %>%
    url_parser(url_param = list(substr(pathway_id, 1, 3), pathway_id)) %>%
    download_base(fun = NULL, path = path) %>%
    invisible

}


#' Information about a KEGG Pathway
#'
#' @param pathway_id Character: a KEGG Pathway identifier, e.g. "hsa04710".
#'     For a complete list of IDs see \code{\link{kegg_pathway_list}}.
#'
#' @return List with the pathway information.
#'
#' @examples
#' kegg_info('map00563')
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr map discard keep walk
#' @importFrom xml2 read_html xml_find_all xml_children
#' @importFrom xml2 xml_text xml_find_first
#' @importFrom stringr str_replace
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{kegg_pathway_list}}}
#'     \item{\code{\link{kegg_picture}}}
#'     \item{\code{\link{kegg_open}}}
#' }
kegg_info <- function(pathway_id){

    info <- list(
        id = pathway_id,
        name = NULL,
        desc = NULL,
        pubmed = NULL,
        diseases = NULL,
        rel_pathways = NULL,
        module = NULL
    )

    fields <- list(
        Description = 'desc',
        Disease = 'diseases',
        Reference = 'pubmed',
        Name = 'name',
        Module = 'module',
        Relatedpathway = 'rel_pathways'
    )

    top_env <- environment()

    path <-
        'kegg_pw_info' %>%
        download_to_cache(url_param = list(pathway_id), ext = 'html')

    path %>%
    read_html %>%
    xml_find_all('.//tr') %>%
    map(xml_children) %>%
    discard(function(ch){length(ch) < 2}) %>%
    map(
        function(ch){
            list(
                ch[1] %>% xml_text,
                ch[2]
            )
        }
    ) %>%
    keep(function(item){item[[1]] %in% names(fields)}) %>%
    walk(
        function(item){
            tbl <- item[[2]] %>% xml_find_first('.//table')
            key <- fields[[item[[1]]]]
            content <-
                item[[2]] %>%
                {`if`(
                    is.na(tbl),
                    xml_text(.),
                    kegg_info_tbl(.)
                )} %>%
                str_replace('^PMID:', '')
            info[[key]] %<>% c(content)
            assign('info', info, envir = top_env)
        }
    )

    return(info)

}


#' Process a table from a KEGG pathway info sheet
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map discard map_chr
#' @importFrom xml2 xml_find_all xml_children xml_text
#' @noRd
kegg_info_tbl <- function(tbl){

    tbl %>%
    xml_find_all('.//tr') %>%
    map(xml_children) %>%
    discard(function(ch){length(ch) < 2}) %>%
    map_chr(function(ch){ch[2] %>% xml_text})

}


#' Protein pathway annotations
#'
#' Downloads all KEGG pathways and creates a table of protein-pathway
#' annotations.
#'
#' @param pathways A table of KEGG pathways as produced by \code{
#'      \link{kegg_pathways_download}}.
#'
#' @return A data frame (tibble) with UniProt IDs and pathway names.
#'
#' @examples
#' \dontrun{
#' kegg_pw_annot <- kegg_pathway_annotations()
#' kegg_pw_annot
#' # # A tibble: 7,341 x 4
#' #    uniprot genesymbol pathway                pathway_id
#' #    <chr>   <chr>      <chr>                  <chr>
#' #  1 Q03113  GNA12      MAPK signaling pathway hsa04010
#' #  2 Q9Y4G8  RAPGEF2    MAPK signaling pathway hsa04010
#' #  3 Q13972  RASGRF1    MAPK signaling pathway hsa04010
#' #  4 O95267  RASGRP1    MAPK signaling pathway hsa04010
#' #  5 P62834  RAP1A      MAPK signaling pathway hsa04010
#' # # . with 7,336 more rows
#' }
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom rlang %||%
#' @importFrom dplyr select bind_rows distinct
#' @export
#' @seealso \code{\link{kegg_pathways_download}}
kegg_pathway_annotations <- function(pathways = NULL){

    # NSE vs. R CMD check workaround
    uniprot_source <- uniprot_target <- pathway <- pathway_id <-
    genesymbol_source <- genesymbol_target <- NULL

    pathways %<>% {. %||% kegg_pathways_download()}

    pathways %>%
    select(
        uniprot = uniprot_source,
        genesymbol = genesymbol_source,
        pathway,
        pathway_id
    ) %>%
    bind_rows(
        pathways %>%
        select(
            uniprot = uniprot_target,
            genesymbol = genesymbol_target,
            pathway,
            pathway_id
        )
    ) %>%
    distinct

}
