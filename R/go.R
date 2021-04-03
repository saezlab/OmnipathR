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


#' Downloads gene annotations from Gene Ontology
#'
#' Gene Ontology is an ontology of gene subcellular localizations, molecular
#' functions and involvement in biological processes. Gene products across
#' many organisms are annotated with the ontology terms. This function
#' downloads the gene-ontology term associations for certain model organisms
#' or all organisms. For a description of the columns see
#' \url{http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/}.
#'
#' @param organism Character: either "chicken", "cow", "dog", "human", "pig"
#'     or "uniprot_all".
#' @param aspects Character vector with some of the following elements:
#'     "C" (cellular component), "F" (molecular function) and "P" (biological
#'     process). Gene Ontology is three separate ontologies called as three
#'     aspects. By this parameter you can control which aspects to include
#'     in the output.
#'
#' @return A tibble (data frame) of annotations as it is provided by the
#' database
#'
#' @examples
#' goa_data <- go_annot_download()
#' goa_data
#' # # A tibble: 606,840 x 17
#' #    db       db_object_id db_object_symbol qualifier go_id   db_ref
#' #    <fct>    <chr>        <chr>            <fct>     <chr>   <chr>
#' #  1 UniProt. A0A024RBG1   NUDT4B           NA        GO:000. GO_REF:00.
#' #  2 UniProt. A0A024RBG1   NUDT4B           NA        GO:000. GO_REF:00.
#' #  3 UniProt. A0A024RBG1   NUDT4B           NA        GO:004. GO_REF:00.
#' #  4 UniProt. A0A024RBG1   NUDT4B           NA        GO:005. GO_REF:00.
#' #  5 UniProt. A0A024RBG1   NUDT4B           NA        GO:005. GO_REF:00.
#' # # . with 606,830 more rows, and 11 more variables:
#' # #   evidence_code <fct>, with_or_from <chr>, aspect <fct>,
#' # #   db_object_name <chr>, db_object_synonym <chr>,
#' # #   db_object_type <fct>, taxon <fct>, date <date>,
#' # #   assigned_by <fct>, annotation_extension <chr>,
#' # #   gene_product_from_id <chr>
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom readr cols col_factor col_character col_date
#' @importFrom dplyr filter
go_annot_download <- function(
    organism = 'human',
    aspects = c('C', 'F', 'P')
){

    'omnipath.go_annot_url' %>%
    generic_downloader(
        url_param = list(organism),
        reader_param = list(
            skip = 40,
            col_names = c(
                'db',
                'db_object_id',
                'db_object_symbol',
                'qualifier',
                'go_id',
                'db_ref',
                'evidence_code',
                'with_or_from',
                'aspect',
                'db_object_name',
                'db_object_synonym',
                'db_object_type',
                'taxon',
                'date',
                'assigned_by',
                'annotation_extension',
                'gene_product_from_id'
            ),
            col_types = cols(
                db = col_factor(),
                qualifier = col_factor(),
                evidence_code = col_factor(),
                aspect = col_factor(),
                db_object_type = col_factor(),
                taxon = col_factor(),
                date = col_date('%Y%m%d'),
                assigned_by = col_factor(),
                annotation_extension = col_character(),
                gene_product_from_id = col_character()
            )
        ),
        resource = 'Gene Ontology'
    ) %>%
    filter(aspect %in% aspects) %T>%
    load_success()

}