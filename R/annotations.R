#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
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


#' Imports annotations from OmniPath
#'
#' Imports protein annotations about function, localization, expression,
#' structure and other properties of proteins from OmniPath
#' \url{https://omnipathdb.org/annotations}.
#' Note: there might be also a few miRNAs annotated; a vast majority of
#' protein complex annotations are inferred from the annotations of the
#' members: if all members carry the same annotation the complex inherits.
#'
#' @details Downloading the full annotations
#' dataset is disabled by default because the size of this data is
#' around 1GB. We recommend to retrieve the annotations for a set of proteins
#' or only from a few resources, depending on your interest. You can always
#' download the full database from
#' \url{https://archive.omnipathdb.org/omnipath_webservice_annotations__recent.tsv}
#' using any standard R or \code{readr} method.
#'
#' @return A data frame containing different gene and complex annotations.
#'
#' @param proteins Vector containing the genes or proteins for whom
#'     annotations will be retrieved (UniProt IDs or HGNC Gene Symbols or
#'     miRBase IDs). It is also possible to donwload annotations for protein
#'     complexes. To do so, write 'COMPLEX:' right before the genesymbols of
#'     the genes integrating the complex. Check the vignette for examples.
#' @param resources Load the annotations only from these databases.
#'     See \code{\link{get_annotation_resources}} for possible values.
#' @param wide Convert the annotation table to wide format, which
#'     corresponds more or less to the original resource. If the data comes
#'     from more than one resource a list of wide tables will be returned.
#'     See examples at \code{\link{pivot_annotations}}.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides intercell annotations
#'     only for human. For other organisms, the data will be translated
#'     by orthologous gene pairs from human.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... Additional arguments.
#'
#' @examples
#' annotations <- import_omnipath_annotations(
#'     proteins = c('TP53', 'LMNA'),
#'     resources = c('HPA_subcellular')
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_annotation_databases}}}
#'     \item{\code{\link{pivot_annotations}}}
#' }
#'
#' @importFrom magrittr %<>%
#' @importFrom rlang exec !!!
#' @importFrom logger log_fatal log_success
#' @export
import_omnipath_annotations <- function(
    proteins = NULL,
    resources = NULL,
    wide = FALSE,
    organism = 'human',
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(query_type = 'annotations')
    wide <- pop(args, wide)
    proteins <- pop(args, proteins)

    if(
        is.null(proteins) &&
        is.null(resources)
    ){

        msg <- paste(
            'Downloading the entire annotations database is not allowed',
            'by default because of its huge size (>1GB). If you really',
            'want to do that, you find static files at',
            'https://archive.omnipathdb.org/.',
            'However we recommend to query a set of proteins or a few',
            'resources, depending on your interest.'
        )

        log_fatal(msg)
        stop(msg)

    }

    if(length(proteins) < 600){

        result <- exec(import_omnipath, proteins = proteins, !!!args)

        if(!is.null(proteins)){
            result <- result[
                which(
                    result$uniprot %in% proteins |
                    result$genesymbol %in% proteins
                ),
            ]
        }

    }else{

        parts <- list()

        proteins_chunks <-
            split(
                proteins,
                cut(
                    seq_along(proteins),
                    breaks = length(proteins) / 500L,
                    labels = FALSE
                )
            )

        cat('Downloading ')

        for(proteins_chunk in proteins_chunks){

            cat('.')

            parts <- append(
                parts,
                list(
                    exec(
                        import_omnipath,
                        proteins = proteins_chunk,
                        recursive_call = TRUE,
                        silent = TRUE,
                        !!!args
                    )
                )
            )

        }

        cat(' ready.\n')

        result <- do.call(rbind, parts)

        log_success(
            'Downloaded %d annotation records.',
            nrow(result)
        )

    }

    result %<>% {`if`(wide, pivot_annotations(.), .)}

    return(result)

}


#' Retrieves a list of available resources in the annotations database
#' of OmniPath
#'
#' Get the names of the resources from \url{https://omnipath.org/annotations}.
#'
#' @param dataset ignored for this query type
#' @param ... optional additional arguments
#'
#' @return character vector with the names of the annotation resources
#'
#' @examples
#' get_annotation_resources()
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_omnipath_annotations}}}
#' }
#'
#' @export
get_annotation_resources <- function(dataset = NULL, ...){

    return(get_resources(query_type = 'annotations', datasets = dataset))

}


#' Annotation categories and resources
#'
#' A full list of annotation resources, keys and values.
#'
#' @return A data frame with resource names, annotation key labels and
#'     for each key all possible values.
#'
#' @examples
#' annot_cat <- annotation_categories()
#' annot_cat
#' # # A tibble: 46,307 x 3
#' #    source           label    value
#' #    <chr>            <chr>    <chr>
#' #  1 connectomeDB2020 role     ligand
#' #  2 connectomeDB2020 role     receptor
#' #  3 connectomeDB2020 location ECM
#' #  4 connectomeDB2020 location plasma membrane
#' #  5 connectomeDB2020 location secreted
#' #  6 KEGG-PC          pathway  Alanine, aspartate and glutamate metabolism
#' #  7 KEGG-PC          pathway  Amino sugar and nucleotide sugar metabolism
#' #  8 KEGG-PC          pathway  Aminoacyl-tRNA biosynthesis
#' #  9 KEGG-PC          pathway  Arachidonic acid metabolism
#' # 10 KEGG-PC          pathway  Arginine and proline metabolism
# . with 46,297 more rows
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr separate_rows
#' @export
annotation_categories <- function(){

    # NSE vs. R CMD check workaround
    value <- NULL

    import_omnipath('annotations_summary', license = NA) %>%
    separate_rows(value, sep = '#')

}


#' Converts annotation tables to a wide format
#'
#' Use this method to reconstitute the annotation tables into the format of
#' the original resources. With the `wide=TRUE` option
#' \code{\link{import_omnipath_annotations}} applies this function to the
#' downloaded data.
#'
#' @param annotations A data frame of annotations downloaded from the
#'    OmniPath web service by \code{\link{import_omnipath_annotations}}.
#'
#' @return A wide format data frame (tibble) if the provided data contains
#' annotations from one resource, otherwise a list of wide format tibbles.
#'
#' @examples
#' # single resource: the result is a data frame
#' disgenet <- import_omnipath_annotations(resources = 'DisGeNet')
#' disgenet <- pivot_annotations(disgenet)
#' disgenet
#' # # A tibble: 126,588 × 11
#' #    uniprot genesymbol entity_type disease      type  score   dsi   dpi
#' #    <chr>   <chr>      <chr>       <chr>        <chr> <dbl> <dbl> <dbl>
#' #  1 P04217  A1BG       protein     Schizophren. dise.  0.3  0.7   0.538
#' #  2 P04217  A1BG       protein     Hepatomegaly phen.  0.3  0.7   0.538
#' #  3 P01023  A2M        protein     Fibrosis, L. dise.  0.3  0.529 0.769
#' #  4 P01023  A2M        protein     Acute kidne. dise.  0.3  0.529 0.769
#' #  5 P01023  A2M        protein     Mental Depr. dise.  0.3  0.529 0.769
#' # # . with 126,583 more rows, and 3 more variables: nof_pmids <dbl>,
#' # #   nof_snps <dbl>, source <chr>
#'
#' # multiple resources: the result is a list
#' annotations <- import_omnipath_annotations(
#'     resources = c('DisGeNet', 'SignaLink_function', 'DGIdb', 'kinase.com')
#' )
#' annotations <- pivot_annotations(annotations)
#' names(annotations)
#' # [1] "DGIdb"              "DisGeNet"           "kinase.com"
#' # [4] "SignaLink_function"
#' annotations$kinase.com
#' # # A tibble: 825 x 6
#' #    uniprot genesymbol entity_type group family subfamily
#' #    <chr>   <chr>      <chr>       <chr> <chr>  <chr>
#' #  1 P31749  AKT1       protein     AGC   Akt    NA
#' #  2 P31751  AKT2       protein     AGC   Akt    NA
#' #  3 Q9Y243  AKT3       protein     AGC   Akt    NA
#' #  4 O14578  CIT        protein     AGC   DMPK   CRIK
#' #  5 Q09013  DMPK       protein     AGC   DMPK   GEK
#' # # . with 815 more rows
#'
#' @export
#' @importFrom magrittr %>% %<>% is_greater_than
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select pull group_split group_by group_keys
#' @importFrom purrr map
#' @importFrom rlang set_names
#' @importFrom readr type_convert cols
#'
#' @seealso \code{\link{import_omnipath_annotations}}
pivot_annotations <- function(annotations){

    # NSE vs. R CMD check workaround
    record_id <- entity_type <- NULL

    more_than_one_resources <-
        annotations %>%
        pull(source) %>%
        unique %>%
        length %>%
        is_greater_than(1L)

    if(more_than_one_resources) {

        annotations %>%
        group_by(source) %>%
        {set_names(
            group_split(.),
            group_keys(.)$source
        )} %>%
        map(pivot_annotations)

    } else {

        id_cols <- c(
            'record_id',
            'uniprot',
            'genesymbol'
        )
        entity_type_field <- 'entity_type' %in% annotations$label
        id_cols %<>% {`if`(entity_type_field, ., c(., 'entity_type'))}

        (
            annotations %>%
            {`if`(
                entity_type_field,
                select(., -entity_type),
                .
            )} %>%
            tidyr::pivot_wider(
                id_cols = id_cols,
                names_from = 'label',
                values_from = 'value'
            ) %>%
            select(-record_id) %>%
            type_convert(col_types = cols())
        )

    }

}
