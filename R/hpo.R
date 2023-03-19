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


#' Downloads protein annotations from Human Phenotype Ontology
#'
#' Human Phenotype Ontology (HPO) provides a standardized vocabulary of
#' phenotypic abnormalities encountered in human disease. Each term in the
#' HPO describes a phenotypic abnormality. HPO currently contains over 13,000
#' terms and over 156,000 annotations to hereditary diseases. See more at
#' \url{https://hpo.jax.org/app/}.
#'
#' @return A tibble (data frame) of annotations as it is provided by the
#' database
#'
#' @examples
#' hpo_data <- hpo_download()
#' hpo_data
#' # # A tibble: 231,738 x 9
#' #    entrez_gene_id entrez_gene_symb. hpo_term_id hpo_term_name
#' #             <dbl> <chr>             <chr>       <chr>
#' #  1           8192 CLPP              HP:0000013  Hypoplasia of the ute.
#' #  2           8192 CLPP              HP:0004322  Short stature
#' #  3           8192 CLPP              HP:0000786  Primary amenorrhea
#' #  4           8192 CLPP              HP:0000007  Autosomal recessive i.
#' #  5           8192 CLPP              HP:0000815  Hypergonadotropic hyp.
#' # # . with 231,733 more rows, and 5 more variables:
#' # #   frequency_raw <chr>, frequency_hpo <chr>, info_gd_source <chr>,
#' # #   gd_source <chr>, disease_id <chr>
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom dplyr mutate across na_if
#' @importFrom tidyselect vars_select_helpers
hpo_download <- function(){

    'hpo' %>%
    generic_downloader(
        reader_param = list(
            skip = 1,
            col_names = c(
                'entrez_gene_id',
                'entrez_gene_symbol',
                'hpo_term_id',
                'hpo_term_name',
                'frequency_raw',
                'frequency_hpo',
                'info_gd_source',
                'gd_source',
                'disease_id'
            )
        ),
        resource = 'Human Phenotype Ontology'
    ) %>%
    mutate(across(vars_select_helpers$where(is.character), na_if, '-')) %T>%
    load_success()

}
