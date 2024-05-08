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

PROTEIN_FIELDS <- c(
    'accession',
    'name',
    'version',
    'general_function',
    'specific_function',
    'protein_type'
)

PROTEIN_ID_FIELDS <- c(
    'gene_name',
    'genebank_protein_id',
    'uniprot_id',
    'uniprot_name',
    'hgnc_id',
    'genecard_id',
    'geneatlas_id'
)

METABOLITE_FIELDS <- c(
    'accession',
    'name',
    'version',
    'status',
    'description',
    'chemical_formula',
    'average_molecular_weight',
    'monisotopic_molecular_weight',
    'iupac_name',
    'traditional_iupac',
    'cas_registry_number',
    'smiles',
    'inchi',
    'inchikey',
    'state',
    'synthesis_reference'
)

METABOLITE_ID_FIELDS <- c(
    'chemspider_id',
    'drugbank_id',
    'foodb_id',
    'pubchem_compound_id',
    'pdb_id',
    'chebi_id',
    'phenol_explorer_compound_id',
    'knapsack_id',
    'kegg_id',
    'biocyc_id',
    'bigg_id',
    'wikipedia_id',
    'metlin_id',
    'vmh_id',
    'fbonto_id'
)

MULTI_FIELDS <- c(
    'secondary_accessions',
    'synonyms'
)


#' Field names for the HMDB metabolite dataset
#'
#' @return Character vector of field names.
#'
#' @examples
#' hmdb_metabolite_fields()
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{hmdb_table}}}
#'     \item{\code{\link{hmdb_protein_fields}}}
#' }
hmdb_metabolite_fields <- function() {

    return(c(
        METABOLITE_FIELDS,
        METABOLITE_ID_FIELDS,
        MULTI_FIELDS
    ))

}


#' Field names for the HMDB proteins dataset
#'
#' @return Character vector of field names.
#'
#' @examples
#' hmdb_protein_fields()
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{hmdb_table}}}
#'     \item{\code{\link{hmdb_metabolite_fields}}}
#' }
hmdb_protein_fields <- function() {

    return(c(
        PROTEIN_FIELDS,
        PROTEIN_ID_FIELDS,
        MULTI_FIELDS
    ))

}


#' Simple xml2 parser for HMDB
#'
#' @importFrom magrittr %>% extract equals
#' @importFrom purrr map map_chr map_int
#' @importFrom xml2 xml_find_all xml_find_first xml_children xml_text xml_ns
#' @noRd
hmdb_xml2_parse <- function(dataset, fields) {

    dataset %<>% the_record

    log_info('Creating xml2 parser for record `%s`.', dataset)

    function(xml) {

        xml %<>% read_xml
        ns <- xml %>% xml_ns %>% names %>% extract(1L)

        xml %>%
        xml_find_all(sprintf('//%s:%s', ns, dataset)) %>%
        {map(
            fields,
            \(field, nodeset = .) {
                xml_find_first(nodeset, sprintf('./%s:%s', ns, field))
            }
        )} %>%
        set_names(fields) %>%
        map(
            ~map(
                .x,
                ~`if`(
                    length(xml_children(.x)) == 0L,
                    xml_text(.x),
                    map_chr(xml_children(.x), xml_text)
                )
            ) %>%
            {`if`(
                map_int(., length) %>% equals(1L) %>% all,
                unlist(.),
                .
            )}
        )

    }

}


#' Download a HMDB XML file and process it into a table
#'
#' @param fields Character: fields to extract from the XML. This is a very
#'     minimal parser that is able to extract the text content of simple fields
#'     and multiple value fields which contain a list of leaves within one
#'     container tag under the record tag. A full list of fields available in
#'     HMDB is available by the \code{\link{hmdb_protein_fields}} and \code{
#'     \link{hmdb_metabolite_fields}} functions.
#' @param dataset Character: name of an HMDB XML dataset, such as
#'     "metabolites", "proteins", "urine", "serum", "csf", "saliva", "feces",
#'     "sweat".
#'
#' @return A data frame (tibble) with each column corresponding to a field.
#'
#' @examples
#' hmdb_table(fields())
#'
#' @importFrom magrittr %<>% %>% extract2
#' @importFrom XML xmlEventParse
#' @importFrom tibble as_tibble
#' @importFrom tidyr unnest_wider
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{hmdb_protein_fields}}}
#'     \item{\code{\link{hmdb_metabolite_fields}}}
#' }
hmdb_table <- function(
        fields = hmdb_metabolite_fields(),
        dataset = 'metabolites'
    ) {

    dataset %>%
    {`if`(
        . %in% c('metabolites', 'proteins'),
        sprintf('hmdb_%s', .),
        sprintf('%s_metabolites', .)
    )} %>%
    {archive_extractor(
        url_key = 'hmdb',
        path = sprintf('%s.xml', .),
        url_param = list(.)
    )} %>%
    parse_in_chunks(
        record = dataset %>% the_record,
        header_lines = 2L,
        parser = hmdb_xml2_parse(dataset, fields)
    ) %>%
    as_tibble

}


#' Record tag from dataset name
#'
#' @importFrom magrittr %>%
#' @noRd
the_record <- function(dataset) {

    dataset %>% {`if`(. == 'proteins', 'protein', 'metabolite')}

}
