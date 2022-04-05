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


#' Downloads interactions from the Guide to Pharmacology database
#'
#' Downloads ligand-receptor interactions from the Guide to Pharmacology
#' (IUPHAR/BPS) database (\url{https://www.guidetopharmacology.org/}).
#'
#' @return A tibble (data frame) of interactions as it is provided by the
#' database
#'
#' @examples
#' g2p_data <- guide2pharma_download()
#' g2p_data
#' # # A tibble: 21,586 x 38
#' #    target target_id target_gene_sym. target_uniprot target_ensembl_.
#' #    <chr>      <dbl> <chr>            <chr>          <chr>
#' #  1 12S-L.      1387 ALOX12           P18054         ENSG00000108839
#' #  2 15-LO.      1388 ALOX15           P16050         ENSG00000161905
#' #  3 15-LO.      1388 ALOX15           P16050         ENSG00000161905
#' #  4 15-LO.      1388 ALOX15           P16050         ENSG00000161905
#' # # . with 21,576 more rows, and 33 more variables: target_ligand <chr>,
#' # #   target_ligand_id <chr>, target_ligand_gene_symbol <chr>,
#' # ... (truncated)
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom readr cols col_character col_number
#' @importFrom tidyr separate_rows
#' @importFrom tibble tibble
#' @importFrom logger log_error
guide2pharma_download <- function(){

    # NSE vs. R CMD check workaround
    ligand_gene_symbol <- target_gene_symbol <- NULL

    tryCatch(

        {

            'guide2pharma' %>%
            generic_downloader(
                reader_param = list(
                    comment = '##',
                    col_types = cols(
                        ligand_context = col_character(),
                        receptor_site = col_character(),
                        target_ligand = col_character(),
                        target_ligand_id = col_character(),
                        target_ligand_pubchem_sid = col_number(),
                        target_ligand_gene_symbol = col_character(),
                        target_ligand_uniprot = col_character(),
                        target_ligand_ensembl_gene_id = col_character()
                    )
                ),
                resource = 'Guide to Pharmacology (IUPHAR/BPS)'
            ) %>%
            separate_rows(ligand_gene_symbol, sep = '\\|') %>%
            separate_rows(target_gene_symbol, sep = '\\|') %T>%
            load_success()

        },

        error = function(err){

            'guide2pharma' %>%
            get_url %>%
            close_connection

            log_error(
                paste0(
                    'Failed to download data from Guide to Pharmacology ',
                    '(guidetopharmacology.org). Most likely it is due to ',
                    'a temporary issue with the server`s SSL certificate. ',
                    'Returning an empty data frame. The original error ',
                    'message was: %s'
                ),
                err
            )

            tibble()

        }

    )

}
