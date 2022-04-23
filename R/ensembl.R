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
#' @importFrom logger log_warn log_trace
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

    log_trace('BioMart query: %s', query)

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
                'BioMart: missing success flag, data might ',
                'be incomplete or contain error message!'
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


#' Identifier translation table from Ensembl
#'
#' @param to Character or symbol: target ID type. See Details for possible
#'     values.
#' @param from Character or symbol: source ID type. See Details for possible
#'     values.
#' @param organism Character or integer: NCBI Taxonomy ID or name of the
#'     organism (by default 9606 for human).
#'
#' @return A data frame (tibble) with columns `From` and `To`.
#'
#' @details The arguments \code{to} and \code{from} can be provided either
#' as character or as symbol (NSE). Their possible values are either Ensembl
#' attribute names or synonyms listed at \code{\link{translate_ids}}.
#'
#' @examples
#' ensp_up <- ensembl_id_mapping_table("ensp")
#' ensp_up
#' # # A tibble: 119,129 × 2
#' #    From   To
#' #    <chr>  <chr>
#' #  1 P03886 ENSP00000354687
#' #  2 P03891 ENSP00000355046
#' #  3 P00395 ENSP00000354499
#' #  4 P00403 ENSP00000354876
#' #  5 P03928 ENSP00000355265
#' # # . with 119,124 more rows
#'
#' @importFrom dplyr recode
#' @importFrom magrittr %>%
#' @importFrom rlang enquo !! !!! set_names
#' @importFrom logger log_warn log_trace
#' @importFrom tibble tibble
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{translate_ids}}}
#'     \item{\code{\link{uniprot_full_id_mapping_table}}}
#'     \item{\code{\link{uniprot_id_mapping_table}}}
#' }
ensembl_id_mapping_table <- function(
    to,
    from = 'uniprot',
    organism = 9606
){

    get_field_name <- function(label){

        label %>% recode(!!!omnipath.env$id_types$ensembl)

    }

    to <-
        .nse_ensure_str(!!enquo(to)) %>%
        get_field_name()
    from <-
        .nse_ensure_str(!!enquo(from)) %>%
        get_field_name()

    organism %<>% ensembl_name_warn

    if(is.na(organism)){

        log_warn('Returning empty table.')

        return(tibble(From = character(0), To = character(0)))

    }

    log_trace(
        paste0(
            'Creating ID mapping table from `%s` to `%s`, ',
            'for organism %s'
        ),
        from, to, organism
    )

    dataset <- ensembl_dataset(organism)

    biomart_query(
        attrs = c(from, to),
        dataset = dataset
    ) %>%
    set_names(c('From', 'To'))

}

#' Ensembl dataset name from organism
#'
#' @param organism Character or integer: an organism (taxon) name or
#'     identifier. If an Ensembl dataset name is provided
#'
#' @return Character: name of an ensembl dataset.
#'
#' @examples
#' ensembl_dataset(10090)
#' # [1] "mmusculus_gene_ensembl"
#'
#' @importFrom magrittr %>% %T>%
#' @export
ensembl_dataset <- function(organism) {

    organism %>%
    as.character %>%
    {`if`(
        endsWith(., '_gene_ensembl'),
        .,
        ensembl_name_warn(.) %>%
        sprintf('%s_gene_ensembl', .)
    )}

}


#' Same as ensembl_name but warns upon failure
#'
#' @importFrom logger log_warn
#' @noRd
ensembl_name_warn <- function(name){

    ensname <- ensembl_name(name)

    if(is.na(ensname)){

        log_warn(
            'Could not find Ensembl name for organism `%s`.',
            as.character(name)
        )

    }

    ensname

}


#' Orthologous gene pairs from Ensembl
#'
#' @param organism_a Character or integer: organism name or identifier for
#'     the left side organism. We query the Ensembl dataset of this organism
#'     and add the orthologues of the other organism to it. Ideally this is
#'     the organism you translate from.
#' @param organism_b Character or integer: organism name or identifier for
#'     the right side organism. We add orthology information of this organism
#'     to the gene records of the left side organism.
#' @param attrs_a Further attributes about organism_a genes. Will be simply
#'     added to the attributes list.
#' @param attrs_b Further attributes about organism_b genes (orthologues).
#'     The available attributes are: "associated_gene_name", "chromosome",
#'     "chrom_start", "chrom_end",  "wga_coverage", "goc_score", "perc_id_r1",
#'     "perc_id", "subtype". Attributes included by default: "ensembl_gene",
#'     "ensembl_peptide", "canonical_transcript_protein",
#'     "orthology_confidence" and "orthology_type".
#' @param colrename Logical: replace prefixes from organism_b attribute
#'     column names, so the returned table always have the same column
#'     names, no matter the organism. E.g. for mouse these columns all
#'     have the prefix "mmusculus_homolog_", which this option changes
#'     to "b_".
#'
#' @details Only the records with orthology information are returned. The
#'     order of columns is the following: defaults of organism_a, extra
#'     attributes of organism_b, defaults of organism_b, extra attributes
#'     of organism_b.
#'
#' @examples
#' \dontrun{
#' sffish <- ensembl_homology(
#'     organism_b = 'Siamese fighting fish',
#'     attrs_a = 'external_gene_name',
#'     attrs_b = 'associated_gene_name'
#' )
#' sffish
#' # # A tibble: 175,608 × 10
#' #    ensembl_gene_id ensembl_transcript_id ensembl_peptide. external_gene_n.
#' #    <chr>           <chr>                 <chr>            <chr>
#' #  1 ENSG00000277196 ENST00000621424       ENSP00000481127  NA
#' #  2 ENSG00000277196 ENST00000615165       ENSP00000482462  NA
#' #  3 ENSG00000278817 ENST00000613204       ENSP00000482514  NA
#' #  4 ENSG00000274847 ENST00000400754       ENSP00000478910  MAFIP
#' #  5 ENSG00000273748 ENST00000612919       ENSP00000479921  NA
#' # # . with 175,603 more rows, and 6 more variables:
#' # #   b_ensembl_peptide <chr>, b_ensembl_gene <chr>,
#' # #   b_orthology_type <chr>, b_orthology_confidence <dbl>,
#' # #   b_canonical_transcript_protein <chr>, b_associated_gene_name <chr>
#' #
#' }
#'
#' @importFrom magrittr %<>% %>%
#' @importFrom dplyr rename_with
#' @importFrom stringr str_replace
#' @importFrom logger log_trace
#' @export
ensembl_homology <- function(
    organism_a = 9606,
    organism_b = 10090,
    attrs_a = NULL,
    attrs_b = NULL,
    colrename = TRUE
){

    HOMOLOGY_ATTRS_DEFAULT <- c(
        'homolog_ensembl_peptide',
        'homolog_ensembl_gene',
        'homolog_orthology_type',
        'homolog_orthology_confidence',
        'homolog_canonical_transcript_protein'
    )

    organism_a %<>% ensembl_name_warn
    organism_b %<>% ensembl_name_warn
    dataset <- organism_a %>% ensembl_dataset
    filters <- organism_b %>% sprintf('with_%s_homolog', .)

    attrs_b %<>% map_chr(~sprintf('homolog_%s', .x))

    attrs <-
        HOMOLOGY_ATTRS_DEFAULT %>%
        c(attrs_b) %>%
        map_chr(~sprintf('%s_%s', organism_b, .x)) %>%
        c(attrs_a, .)

    log_trace(
        'Querying orthology data from Ensembl, attributes: %s.',
        paste(attrs, collapse = ', ')
    )

    biomart_query(
        attrs = attrs,
        filters = filters,
        transcript = TRUE,
        peptide = TRUE,
        gene = TRUE,
        dataset = dataset
    ) %>%
    {`if`(
        colrename,
        rename_with(
            .,
            str_replace,
            pattern = sprintf('%s_homolog_', organism_b),
            replacement = 'b_'
        ),
        .
    )}

}
