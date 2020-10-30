#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2020
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


#' Build NicheNet prior knowledge
#'
#' Builds all prior knowledge data required by NicheNet. For this it calls
#' a multitude of methods to download and combine data from various
#' databases according to the settings. The content of the prior knowledge
#' data is highly customizable, see more in Details.
#'
#' @param signaling_network A list of parameters for building the signaling
#'     network, passed to \code{\link{nichenet_signaling_network}}
#' @export
#'
#' @seealso \code{\link{nichenet_signaling_network}}
nichenet_prior_knowledge <- function(
    signaling_network = list()
){

    list(
        signaling_network = do.call(
            nichenet_signaling_network,
            signaling_network
        )
    )

}


#' Build a NicheNet signaling network
#'
#' Builds signaling network prior knowledge for NicheNet using multiple
#' resources.
#'
#' @param omnipath List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_omnipath}}
#' @param pathwaycommons List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_pathwaycommons}}
#' @param harmonizome List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_harmonizome}}
#' @export
#'
#' @seealso \code{\link{nichenet_signaling_network_omnipath},
#'     \link{nichenet_signaling_network_pathwaycommons},
#'     \link{nichenet_signaling_network_harmonizome}}
nichenet_signaling_network <- function(
    omnipath = list(),
    pathwaycommons = NULL,
    harmonizome = NULL,
    vinayagam = NULL,
    cpdb = NULL,
    evex = NULL,
    inweb = NULL
){



}


#' Builds signaling network prior knowledge for NicheNet using OmniPath
#'
#' @param min_curation_effort Lower threshold for curation effort
#' @param ... Passed to \code{\link{import_post_translational_interactions}}
#' @importsFrom dplyr %>% mutate select
#' @export
#'
#' @seealso
nichenet_signaling_network_omnipath <- function(
    min_curation_effort = 0,
    ...
){

    import_post_translational_interactions(entity_types = 'protein', ...) %>%
    select(from = source_genesymbol, to = target_genesymbol, is_directed) %>%
    mutate(
        source = ifelse(
            is_directed,
            'omnipath_directed',
            'omnipath_undirected'
        ),
        database = 'omnipath'
    ) %>%
    select(-is_directed)

}


#' NicheNet signaling network from PathwayCommons
#'
#' Builds signaling network prior knowledge for NicheNet using PathwayCommons.
#'
#' @param interaction_types Character vector with PathwayCommons interaction
#'     types. Please refer to the default value and the PathwayCommons
#'     webpage.
#' importsFrom dplyr %>% mutate
#' @importsFrom readr read_tsv
#' @export
nichenet_signaling_network_pathwaycommons <- function(
    interaction_types = c(
        'catalysis-precedes',
        'controls-phosphorylation-of',
        'controls-state-change-of',
        'controls-transport-of',
        'in-complex-with',
        'interacts-with'
    ),
    ...
){

    options('omnipath.pathwaycommons_url') %>%
    read_tsv(
        col_names = c('from', 'source', 'to'),
        col_types = cols()
    ) %>%
    filter(
        source %in% interaction_types
    ) %>%
    mutate(
        source = sprintf(
            'pathwaycommons_%s',
            gsub('-', '_', source, fixed = TRUE)
        ),
        database = 'pathwaycommons_signaling'
    )

}


#' NicheNet signaling network from Harmonizome
#'
#' Builds signaling network prior knowledge for NicheNet using Harmonizome
#'
#' @param datasets The datasets to use. For possible values please refer to
#'     default value and the Harmonizome webpage.
#' @importsFrom dplyr %>% mutate bind_rows
#' @export
nichenet_signaling_network_harmonizome <- function(
    datasets = c(
        'phosphosite',
        'kea',
        'depod',
    ),
    ...
){

    dataset_names <- list(
        phosphosite = 'PhosphoSite',
        kea = 'KEA',
        depod = 'DEPOD'
    )

    do.call(
        bind_rows,
        datasets %>% lapply(harmonizome_download)
    ) %>%
    mutate(
        source = sprintf('harmonizome_%s', dataset_names[source]),
        database = 'harmonizome'
    )

}


#' NicheNet signaling network from Vinayagam
#'
#' Builds signaling network prior knowledge for NicheNet using Vinayagam 2011
#' Supplementary Table S6
#'
#' Find out more at https://doi.org/10.1126/scisignal.2001699
#'
#' @importsFrom dplyr %>% select mutate
#' @importsFrom readxl read_xls
#' @export
nichenet_signaling_network_vinayagam <- function(...){

    tmp_zip <- tempfile(fileext = '.zip')
    tmp_xls <- tempfile(fileext = '.xls')
    xls_con <- file(tmp_xls, open = 'wb')

    on.exit(unlink(tmp_xls))

    'omnipath.vinayagam_url' %>%
    options() %>%
    as.character() %>%
    download.file(destfile = tmp_zip, quiet = TRUE)

    xls_unz_con <- unz(tmp_zip, '2001699_Tables_S1_S2_S6.xls', open = 'rb')

    xls_unz_con %>%
    readBin(raw(), n = 6000000) %>%
    writeBin(xls_con, useBytes = TRUE)

    close(xls_unz_con)
    close(xls_con)
    unlink(tmp_zip)

    tmp_xls %>%
    read_xls(sheet = 'S6', progress = FALSE) %>%
    select(
        from = `Input-node Gene Symbol`,
        to = `Output-node Gene Symbol`
    ) %>%
    mutate(
        source = 'vinayagam_ppi',
        database = 'vinayagam'
    )

}


#' Builds signaling network prior knowledge for NicheNet using ConsensusPathDB
#' (CPDB)
#'
#' @importsFrom dplyr select mutate distinct
#' @importsFrom magrittr %>%
#' @export
nichenet_signaling_network_cpdb <- function(...){

    consensuspathdb() %>%
    select(from = genesymbol_a, to = genesymbol_b) %>%
    distinct() %>%
    mutate(
        source = sprintf(
            'cpdb_%s',
            ifelse(in_complex, 'complex', 'interaction')
        ),
        database = 'cpdb'
    ) %>%
    select(-in_complex)

}


#' NicheNet signaling network from EVEX
#'
#' Builds signaling network prior knowledge for NicheNet from the EVEX
#' database.
#'
#' @importsFrom magrittr %>% `n'est pas`
#' @importsFrom dplyt select mutate filter
#' @export
#'
#' @seealso \code{\link{evex}}
nichenet_signaling_network_evex <- function(...){

    categories <- list(
        Binding = 'binding',
        `Regulation of binding` = 'regulation_binding',
        Regulation_of_phosphorylation = 'phosphorylation'
    )

    evex() %>%
    select(
        from = source_genesymbol,
        to = target_genesymbol,
        coarse_type,
        refined_type
    ) %>%
    filter(
        coarse_type != 'Indirect_regulation' &
        `n'est pas`(refined_type %in% c(
            # these belong to transcriptional regulation
            'Regulation of expression',
            'Regulation of transcription',
            'Catalysis of DNA methylation'
        ))
    ) %>%
    mutate(
        source = sprintf(
            'evex_%s',
            ifelse(
                refined_type %in% names(categories),
                categories[refined_type],
                ifelse(
                    startsWith(refined_type, 'Catalysis'),
                    'catalysis',
                    ifelse(
                        coarse_type == 'Regulation',
                        'regulation_other',
                        'binding'
                    )
                )
            )
        ),
        database = 'evex_signaling'
    ) %>%
    select(-coarse_type, -refined_type)


}