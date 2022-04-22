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


#' Table of Ensembl organisms
#'
#' A table with various taxon IDs and metadata about related Ensembl database
#' contents, as shown at https://www.ensembl.org/info/about/species.html.
#' The "Taxon ID" column contains the NCBI Taxonomy identifiers.
#'
#' @return The table described above as a data frame.
#'
#' @examples
#' ens_org <- ensembl_organisms_raw()
#' ens_org
#'
#' @importFrom magrittr %>%
#' @importFrom rvest read_html html_element html_table
#' @export
ensembl_organisms_raw <- function(){

    'ensembl_organisms' %>%
    download_to_cache %>%
    read_html %>%
    html_element('table') %>%
    html_table()

}


#' Organism names and identifiers from Ensembl
#'
#' A table with various taxon names and identifiers: English common names,
#' latin (scientific) names, Ensembl organism IDs and NCBI taxonomy IDs.
#'
#' @return A data frame with the above mentioned columns.
#'
#' @examples
#' ens_org <- ensembl_organisms()
#' ens_org
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate recode
#' @importFrom purrr pmap_chr map_chr
#' @importFrom stringr str_extract_all str_extract str_sub
#' @importFrom rlang !!!
#' @export
ensembl_organisms <- function(){

    SPECIALS <- list(
        bixbtaurus = 'bihybrid'
    )

    ensembl_organisms_raw() %>%
    select(
        common_name = `Common name`,
        latin_name = `Scientific name`,
        ncbi_tax_id = `Taxon ID`
    ) %>%
    mutate(
        ensembl_id =
            pmap_chr(
                list(
                    str_extract_all(latin_name, '\\b[A-z]') %>%
                    map_chr(paste0, collapse = '') %>%
                    str_to_lower %>%
                    str_sub(end = -2L),
                    str_extract(latin_name, '\\w+$')
                ),
                paste0
            ),
        ensembl_id = recode(ensembl_id, !!!SPECIALS)
    )

}


#' Query the Ensembl BioMart web service
#'
#' @param attrs Character vector: one or more Ensembl attribute names.
#' @param filters Character vector: one or more Ensembl filter names.
#' @param transcript Logical: include Ensembl transcript IDs in the result.
#' @param peptide Logical: include Ensembl peptide IDs in the result.
#' @param gene Logical: include Ensembl gene IDs in the result.
#' @param dataset Character: An Ensembl dataset name.
#'
#' @return Data frame with the query result
#'
#' @examples
#' cel_genes <- biomart_query(
#'     attrs = c("external_gene_name", "start_position", "end_position"),
#'     gene = TRUE,
#'     dataset = "celegans_gene_ensembl"
#' )
#' cel
#' # # A tibble: 46,934 × 4
#' #   ensembl_gene_id external_gene_name start_position end_position
#' #   <chr>           <chr>                       <dbl>        <dbl>
#' # 1 WBGene00000001  aap-1                     5107843      5110183
#' # 2 WBGene00000002  aat-1                     9599178      9601695
#' # 3 WBGene00000003  aat-2                     9244402      9246360
#' # 4 WBGene00000004  aat-3                     2552260      2557736
#' # 5 WBGene00000005  aat-4                     6272529      6275721
#' # # . with 46,924 more rows
#'
#' @importFrom magrittr %<>% %>% %T>%
#' @importFrom purrr map_chr
#' @importFrom readr cols col_character read_tsv type_convert
#' @importFrom dplyr slice_tail slice_head
#' @importFrom logger log_warn
#' @export
biomart_query <- function(
    attrs,
    filters = NULL,
    transcript = FALSE,
    peptide = FALSE,
    gene = FALSE,
    dataset = 'hsapiens_gene_ensembl'
){

    TEMPLATE <- biomart_xml_template()
    FILTER_TEMPLATE <- '<Filter name="%s" excluded="0"/>' %>% indent(8L)
    ATTR_TEMPLATE <- '<Attribute name="%s"/>' %>% indent(8L)

    local_env <- environment()

    attrs %<>%
        {`if`(peptide, c('ensembl_peptide_id',. ), .)} %>%
        {`if`(transcript, c('ensembl_transcript_id', .), .)} %>%
        {`if`(gene, c('ensembl_gene_id', .), .)} %>%
        assign('col_names', ., envir = local_env) %>%
        map_chr(~sprintf(ATTR_TEMPLATE, .x)) %>%
        paste0(collapse = '\n')

    filters %<>%
        map_chr(~sprintf(FILTER_TEMPLATE, .x)) %>%
        paste0(collapse = '\n') %>%
        {`if`(nchar(.), sprintf('\n%s', .), .)}

    query <-
        TEMPLATE %>%
        sprintf(dataset, attrs, filters)

    'biomart' %>%
    generic_downloader(
        url_param = list(query),
        reader = function(...){suppressWarnings(read_tsv(...))},
        reader_param = list(
            col_names = col_names,
            col_types = cols(.default = col_character()),
            progress = FALSE
        )
    ) %>%
    {`if`(
        slice_tail(., n = 1L) %>% first %>% equals('[success]'),
        slice_head(., n = -1L),
        {
            log_warn(
                'BioMart: missing success flag, data might be incomplete!'
            )
            .
        }
    )} %>%
    type_convert(col_types = cols()) %T>%
    load_success()

}


#' BioMart query XML template
#'
#' @return The XML as a single character string.
#'
#' @importFrom utils packageName
#' @importFrom magrittr %>%
#' @importFrom readr read_file
#' @noRd
biomart_xml_template <- function(){

    system.file(
        'internal',
        'biomart.xml',
        package = utils::packageName(),
        mustWork = TRUE
    ) %>%
    read_file

}
