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
#  Website: https://r.omnipathdb.org/
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
#' @importFrom dplyr rename_with
#' @importFrom stringr str_to_lower str_replace_all
#' @importFrom tidyr separate_rows
#' @importFrom tibble tibble
#' @importFrom logger log_error
guide2pharma_download <- function(){

    .slow_doctest()

    # NSE vs. R CMD check workaround
    ligand_gene_symbol <- target_gene_symbol <- NULL

    tryCatch(

        {

            'guide2pharma' %>%
            generic_downloader(
                reader_param = list(
                    comment = '"#',
                    col_types = cols(
                        `Target ID` = col_character(),
                        `Target Subunit IDs` = col_character(),
                        `Target Ligand ID` = col_character(),
                        `Target Ligand Subunit IDs` = col_character(),
                        `Ligand Context` = col_character(),
                        `Receptor Site` = col_character(),
                        `Target Ligand` = col_character(),
                        `Target Ligand ID` = col_character(),
                        `Target Ligand PubChem SID` = col_character(),
                        `Ligand PubChem SID` = col_character(),
                        `Target Ligand Gene Symbol` = col_character(),
                        `Target Ligand UniProt ID` = col_character(),
                        `Target Ligand Ensembl Gene ID` = col_character()
                    )
                ),
                resource = 'Guide to Pharmacology (IUPHAR/BPS)'
            ) %>%
            rename_with(~str_replace_all(str_to_lower(.), ' ', '_')) %>%
            {copy_attrs(
                separate_rows(., ligand_gene_symbol, sep = '\\|'),
                .,
                c('origin', 'source')
            )} %>%
            {copy_attrs(
                separate_rows(., target_gene_symbol, sep = '\\|'),
                .,
                c('origin', 'source')
            )} %T>%
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
                    'a temporary issue with the server`s SSL certificate, ',
                    'or an update in their data format. ',
                    'Returning an empty data frame. The original error ',
                    'message was: %s'
                ),
                err
            )

            tibble()

        }

    )

}
