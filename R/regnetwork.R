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


#' Interactions from RegNetwork
#'
#' Downloads transcriptional and post-transcriptional regulatory interactions
#' from the RegNetwork database (http://www.regnetworkweb.org/).
#'
#' @param organism Character: either human or mouse.
#'
#' @return Data frame with interactions
#'
#' @examples
#' \donttest{
#' regn_interactions <- regnetwork_download()
#' }
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom dplyr left_join mutate select
#' @importFrom readr read_tsv cols col_character
#' @importFrom utils download.file
regnetwork_download <- function(organism = 'human'){

    # NSE vs. R CMD check workaround
    source_entrez <- target_entrez <- NULL

    url <-
        'omnipath.regnetwork_url' %>%
        url_parser(url_param = list(organism))

    version <- omnipath_cache_latest_or_new(url = url)

    from_cache <- version$status == CACHE_STATUS$READY

    if(!from_cache){


        url %>% download.file(destfile = version$path, quiet = TRUE)
        omnipath_cache_download_ready(version)

    }

    version$path %>%
    unz(sprintf('%s.source', organism), open = 'rb') %T>%
    close_on_exit(envir = parent.frame()) %>%
    read_tsv(
        col_names = c(
            'source_genesymbol',
            'source_entrez',
            'target_genesymbol',
            'target_entrez'
        ),
        col_types = cols(
            .default = col_character()
        ),
        progress = FALSE
    ) %>%
    left_join(
        regnetwork_directions(organism = organism) %>%
        select(-source_entrez, -target_entrez),
        by = c('source_genesymbol', 'target_genesymbol')
    ) %>%
    mutate(
        source_type = mirna_or_protein(source_entrez),
        target_type = mirna_or_protein(target_entrez)
    ) %>%
    origin_cache(from_cache) %>%
    source_attrs('RegNetwork', url) %T>%
    load_success()

}


#' Transcription factor effects from RegNetwork
#'
#' @param organism Character: either human or mouse.
#'
#' @examples
#' \donttest{
#' regn_dir <- regnetwork_directions()
#' }
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_delim cols col_character
#' @importFrom dplyr mutate recode
#' @importFrom utils download.file
#'
#' @export
regnetwork_directions <- function(organism = 'human'){

    # NSE vs. R CMD check workaround
    effect <- NULL

    url <-
        'omnipath.regnetwork_url' %>%
        url_parser(url_param = list('RegulatoryDirections'))

    version <- omnipath_cache_latest_or_new(url = url)

    if(version$status != CACHE_STATUS$READY){


        url %>% download.file(destfile = version$path, quiet = TRUE)
        omnipath_cache_download_ready(version)

    }

    version$path %>%
    unz(sprintf('new_kegg.%s.reg.direction.txt', organism), open = 'rb') %T>%
    close_on_exit(envir = parent.frame()) %>%
    read_delim(
        delim = ' ',
        skip = 1,
        col_names = c(
            'source_genesymbol',
            'source_entrez',
            'target_genesymbol',
            'target_entrez',
            'effect'
        ),
        col_types = cols(
            .default = col_character()
        ),
        progress = FALSE
    ) %>%
    mutate(
        effect = recode(
            effect,
            `--|`    = -1,
            `-->`    =  1,
            .default =  0
        )
    )

}


#' For a vector of miRBase IDs and various other IDs mixed, returns a vector
#' of molecule types "mirna" and "protein"
#'
#' @noRd
mirna_or_protein <- function(ids){

    ifelse(
        grepl('^MI(?:MAT)?[0-9]', ids),
        'mirna',
        'protein'
    )

}