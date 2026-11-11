#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
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
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Table of Ensembl organisms
#'
#' A table of the organisms in Ensembl, as listed in the `species.json` file
#' at https://ftp.ebi.ac.uk/pub/ensemblorganisms/. The table contains one row
#' for each organism (species or strain) with at least one genome in Ensembl.
#' If the EBI FTP is not available, a copy at https://rescued.omnipathdb.org/
#' is used. On build servers, the copy is the first choice.
#'
#' @return A data frame with the following columns:
#'     \describe{
#'         \item{species}{Character: the name of the organism in Ensembl,
#'             e.g. "Homo_sapiens".}
#'         \item{ncbi_tax_id}{Integer: NCBI Taxonomy ID of the organism.}
#'         \item{species_tax_id}{Integer: NCBI Taxonomy ID of the species.
#'             For strains, it is different from \code{ncbi_tax_id}.}
#'         \item{latin_name}{Character: scientific name.}
#'         \item{common_name}{Character: English common name, \code{NA} for
#'             most microorganisms.}
#'     }
#'
#' @examples
#' ens_org <- ensembl_organisms_raw()
#' ens_org
#'
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom jsonlite read_json
#' @importFrom purrr map_int map_chr
#' @importFrom tibble tibble
#' @export
ensembl_organisms_raw <- function(){

    # the rescued copy contains only the fields below, it is faster to
    # download than the original from the EBI FTP
    url_keys <- c('ensembl_organisms', 'ensembl_organisms_rescued')
    url_keys %<>% `if`(.on_buildserver(), rev(.), .)

    organisms <-
    tryCatch(
        url_keys[1L] %>%
        download_to_cache(http_headers = user_agent()),
        error = function(err) {

            log_error_with_info(
                'Failed to download Ensembl organisms table.',
                err
            )

            log_trace('Accessing Ensembl organisms table from `%s`.', url_keys[2L])
            url_keys[2L] %>%
            download_to_cache(http_headers = user_agent())

        }
    ) %>%
    read_json %>%
    extract2('species')

    tibble(
        species = names(organisms),
        ncbi_tax_id = organisms %>% map_int('taxid') %>% unname,
        species_tax_id = organisms %>% map_int('species_taxonomy_id') %>% unname,
        latin_name = organisms %>% map_chr('scientific_name') %>% unname,
        common_name =
            organisms %>%
            map_chr('common_name', .default = NA_character_) %>%
            unname
    )

}


#' Organism names and identifiers from Ensembl
#'
#' A table with various taxon names and identifiers: English common names,
#' latin (scientific) names, Ensembl organism IDs and NCBI taxonomy IDs.
#' The names and NCBI Taxonomy IDs are from \code{ensembl_organisms_raw}.
#' The Ensembl organism IDs (e.g. "hsapiens") are the names of the Ensembl
#' BioMart datasets, available only for the about 200 organisms in BioMart.
#'
#' @return A data frame with the above mentioned columns, one row for each
#'     NCBI Taxonomy ID. For the organisms in BioMart, the column
#'     `common_name` contains the name used in BioMart (e.g. "Mouse"), and
#'     `common_name_alt` the name in the Ensembl organisms table (e.g.
#'     "House mouse").
#'
#' @examples
#' ens_org <- ensembl_organisms()
#' ens_org
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate arrange distinct full_join coalesce
#' @export
ensembl_organisms <- function(){

    # NSE vs. R CMD check workaround
    latin_name <- ncbi_tax_id <- common_name <- ensembl_id <-
    latin_name_biomart <- common_name_biomart <- common_name_alt <- NULL

    ensembl_organisms_raw() %>%
    # one row per taxon: the shortest name, i.e. the species, not a strain
    arrange(nchar(latin_name)) %>%
    distinct(ncbi_tax_id, .keep_all = TRUE) %>%
    select(common_name, latin_name, ncbi_tax_id) %>%
    full_join(
        ensembl_biomart_organisms(),
        by = 'ncbi_tax_id',
        suffix = c('', '_biomart')
    ) %>%
    mutate(
        latin_name = coalesce(latin_name, latin_name_biomart),
        common_name_alt = common_name,
        # the BioMart names are the short, familiar ones (e.g. "Mouse")
        common_name = coalesce(common_name_biomart, common_name)
    ) %>%
    select(common_name, latin_name, ncbi_tax_id, ensembl_id, common_name_alt)

}


#' Organisms in Ensembl BioMart
#'
#' Ensembl BioMart is not developed any more, the latest release is at
#' https://jun2026.archive.ensembl.org/biomart/. The organism IDs in this
#' table are the dataset names in BioMart (without the "_gene_ensembl"
#' suffix). The table is part of the package, because the organisms in
#' BioMart do not change any more.
#'
#' @return A data frame with the columns `ensembl_id`, `ncbi_tax_id`,
#'     `latin_name` and `common_name`.
#'
#' @importFrom readr read_tsv cols col_character col_integer
#' @noRd
ensembl_biomart_organisms <- function(){

    system.file(
        'internal',
        'ensembl_biomart.tsv',
        package = 'OmnipathR',
        mustWork = TRUE
    ) %>%
    read_tsv(
        col_types = cols(
            ensembl_id = col_character(),
            ncbi_tax_id = col_integer(),
            latin_name = col_character(),
            common_name = col_character()
        ),
        progress = FALSE
    )

}


# Ensembl BioMart servers, in the order we try them
BIOMART_MIRRORS <- c('biomart', 'biomart_useast', 'biomart_asia')


#' Query the Ensembl BioMart web service
#'
#' The query goes to the main Ensembl site first. If it fails, the US East
#' and the Asia mirrors are tried, then the next round of attempts starts
#' with the main site again. On build servers, if all of them fail,
#' a warning is issued and an empty data frame is returned, so a temporary
#' outage of Ensembl does not fail the package check. After such a failure,
#' the further queries in the same R session return an empty data frame
#' without a download attempt.
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
#' cel_genes
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
#' @importFrom magrittr %<>% %>% %T>% extract2 equals
#' @importFrom purrr map_chr map
#' @importFrom readr cols col_character read_tsv type_convert
#' @importFrom dplyr slice_tail slice_head
#' @importFrom logger log_warn log_trace
#' @importFrom tibble tibble
#' @importFrom rlang set_names !!!
#' @export
biomart_query <- function(
    attrs = NULL,
    filters = NULL,
    transcript = FALSE,
    peptide = FALSE,
    gene = FALSE,
    dataset = 'hsapiens_gene_ensembl'
){

    # NSE vs. R CMD check workaround
    col_names <- NULL

    TEMPLATE <- biomart_xml_template()
    FILTER_TEMPLATE <- '<Filter name="%s" excluded="0"/>' %>% indent(8L)
    ATTR_TEMPLATE <- '<Attribute name="%s"/>' %>% indent(8L)
    NCBI_TAX_ID <-
        dataset %>%
        str_split('_', simplify = TRUE) %>%
        first %>%
        ncbi_taxid

    local_env <- environment()

    attrs %<>%
        {`if`(peptide, c('ensembl_peptide_id',. ), .)} %>%
        {`if`(transcript, c('ensembl_transcript_id', .), .)} %>%
        {`if`(gene, c('ensembl_gene_id', .), .)} %>%
        assign('col_names', ., envir = local_env) %>%
        map_chr(~sprintf(ATTR_TEMPLATE, .x)) %>%
        paste0(collapse = '\n')

    if(!nchar(attrs)){

        msg <- 'BioMart: the query must contain at least one attribute.'
        omnipath_error(msg)

    }

    filters %<>%
        map_chr(~sprintf(FILTER_TEMPLATE, .x)) %>%
        paste0(collapse = '\n') %>%
        {`if`(nchar(.), sprintf('\n%s', .), .)}

    query <-
        TEMPLATE %>%
        sprintf(dataset, attrs, filters)

    log_trace('BioMart query: %s', query)

    empty <-
        col_names %>%
        set_names(., .) %>%
        as.list %>%
        map(~character(0L)) %>%
        tibble(!!!.)

    if (.on_buildserver() && isTRUE(omnipathr.env$biomart_down)) {

        log_warn(
            'BioMart failed earlier in this session, ',
            'skipping query on build server.'
        )
        return(empty)

    }

    result <- tryCatch(
        BIOMART_MIRRORS %>%
        generic_downloader(
            url_param = list(query),
            reader = function(...){read_tsv(...)},
            reader_param = list(
                col_names = col_names,
                col_types = cols(.default = col_character()),
                progress = FALSE
            ),
            use_httr = TRUE,
            http_headers = user_agent()
        ),
        error = identity
    )

    if (inherits(result, 'condition')) {

        msg <- sprintf(
            paste0(
                'Failed BioMart query: all Ensembl mirrors (%s) failed. ',
                'This is likely a temporary issue on the Ensembl server ',
                'side. Last error: %s'
            ),
            paste(BIOMART_MIRRORS, collapse = ', '),
            conditionMessage(result)
        )

        if (.on_buildserver()) {

            # on build servers time is limited, if Ensembl does not respond
            # now, we don't try again in this session
            omnipathr.env$biomart_down <- TRUE
            log_warn(msg)
            warning(msg)
            return(empty)

        }

        log_error_with_info(msg)
        stop(result)

    }

    biomart_url <-
        BIOMART_MIRRORS[1L] %>%
        url_parser(url_param = list(query))

    result %<>%
        {`if`(
            slice_tail(., n = 1L) %>% extract2(1L) %>% equals('[success]'),
            slice_head(., n = -1L),
            {
                log_warn(
                    'BioMart: missing success flag, data might ',
                    'be incomplete or contain error message!'
                )
                log_warn(
                    .[[1]]
                )
                .
            }
        )}

    if(
        ncol(result) < length(col_names) ||
        any(grepl('Query ERROR', result[[1L]], fixed = TRUE))
    ){
        omnipath_cache_remove(url = biomart_url)
        msg <- sprintf(
            paste0(
                'Failed BioMart query: the Ensembl BioMart ',
                'returned an error instead of data. This is ',
                'likely a temporary issue on the Ensembl server ',
                'side. Response: %s'
            ),
            paste(result[[1L]], collapse = ' ')
        )
        if(.on_buildserver()){
            log_warn(msg)
            warning(msg)
            return(empty)
        }
        omnipath_error(msg)
    }

    result %>%
    type_convert(col_types = cols()) %>%
    `attr<-`('oraganism', NCBI_TAX_ID) %T>%
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
#' @return A data frame of orthologous gene pairs with gene, transcript
#'     and peptide identifiers and confidence values.
#'
#' @details Only the records with orthology information are returned. The
#'     order of columns is the following: defaults of organism_a, extra
#'     attributes of organism_b, defaults of organism_b, extra attributes
#'     of organism_b.
#'
#' @examples
#' \dontrun{
#' sffish <- ensembl_orthology(
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
ensembl_orthology <- function(
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
