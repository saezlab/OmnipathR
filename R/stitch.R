#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Diego Mananes
#'                 Alberto Valdeolivas
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


#' Download STITCH link data frame from \url{http://stitch.embl.de/}
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return STITCH links data frame
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv
#'
#' @noRd
stitch_proteins <- function(organism) {

    .slow_doctest()

    'stitch_proteins' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = list(organism),
        reader_param = list(trim_ws = TRUE),
        resource = NULL,
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' Download STITCH actions data frame from \url{http://stitch.embl.de/}
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return STITCH actions data frame
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv
#'
#' @noRd
stitch_actions <- function(organism) {

    .slow_doctest()

    'stitch_actions' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = list(organism),
        reader_param = list(trim_ws = TRUE),
        resource = NULL,
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' Processing chemical-protein interactions from STITCH
#'
#' @param stitch.actions STITCH actions data frame obtained from the
#'   \pkg{OmnipathR} R package.
#' @param stitch.links STITCH links data frame obtained from the
#'   \pkg{OmnipathR} R package.
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#' @param omnipath.PKN Protein-protein interactions obtained using
#'   \code{omnipath_for_cosmos}.
#' @param mapping.biomart BioMart ontology mapping data frame. If \code{NULL},
#'   this info is obtained using the \ckg{bioMaRt} R package.
#' @param threshold Confidence cutoff used for STITCH connections
#'   (700 by default).
#'
#' @return List containing PKN with COSMOS and OCEAN format.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter case_when mutate select
#'
#' @noRd
stitch_format_gem <- function(
        stitch.actions,
        stitch.links,
        organism,
        omnipath.PKN,
        mapping.biomart = NULL,
        threshold = 700
) {
    dataset.biomart <- switch(
        as.character(organism),
        '9606' = 'hsapiens_gene_ensembl',
        '10090' = 'mmusculus_gene_ensembl',
        '10116' = 'rnorvegicus_gene_ensembl',
        '7955' = 'drerio_gene_ensembl',
        '7227' = 'dmelanogaster_gene_ensembl',
        '6239' = 'celegants_gene_ensembl'
    )
    if (is.null(dataset.biomart))
        stop('Chosen organism is not recognizable')

    log_info('Reading provided STITCH files.')

    links.detail <- as.data.frame(stitch.links) %>% filter(
        combined_score >= threshold,
        experimental >= threshold | database >= threshold
    ) %>% mutate(
        ID = paste(chemical, protein, sep = '_'),
        ID_reverse = paste(protein, chemical, sep = '_')
    )
    STITCH <- as.data.frame(stitch.actions) %>%
        filter(mode == 'activation' | mode == 'inhibition', a_is_acting) %>%
        mutate(ID = paste(item_id_a, item_id_b, sep = '_')) %>%
        filter(ID %in% links.detail$ID | ID %in% links.detail$ID_reverse)
    STITCH <- STITCH[,-7]
    ## df of proteins in STICH
    prots <- unique(c(STITCH$item_id_a, STITCH$item_id_b))
    prots <- prots[grepl(paste0(organism, '\\.'), prots)]
    prots <- as.data.frame(cbind(prots, gsub(paste0(organism, '\\.'), '', prots)))
    colnames(prots) <- c('original', 'ensembl_prots')
    ## getting info from Biomart
    log_info('Using information from BiomaRt')

    if (is.null(mapping.biomart)) {
        ensembl.link <- useEnsembl(biomart = 'ensembl', dataset = dataset.biomart)
        ensembl.df <- getBM(
            filters = 'ensembl_peptide_id',
            attributes = c(
                'ensembl_peptide_id','ensembl_gene_id', 'external_gene_name'# , 'entrezgene_id', 'description'
            ),
            values = prots[[2]],
            mart = ensembl.link
        )
        colnames(ensembl.df)[1] <- 'ensembl_prots'
    } else {
        ensembl.df <- mapping.biomart %>% filter(
            ensembl_peptide_id %in% prots[[2]]
        )
        colnames(ensembl.df)[1] <- 'ensembl_prots'
    }

    prots <- merge(prots, ensembl.df, by = 'ensembl_prots')
    ## external_gene_name for mouse, Idk in other cases, check this
    prots <- prots[prots$external_gene_name != '',]
    prots.vec <- prots$external_gene_name
    names(prots.vec) <- prots$original

    if (verbose) message('\n\t>>> Generating PKN network')

    STITCH <- STITCH %>% mutate(
        item_id_a = case_when(
            grepl('\\.', item_id_a) & (item_id_a %in% names(prots.vec)) ~
                prots.vec[item_id_a],
            grepl('^CID', item_id_a) ~ gsub('CID[a-z]0*', 'Metab__', item_id_a),
            TRUE ~ item_id_a
        ),
        item_id_b = case_when(
            grepl('\\.', item_id_b) & (item_id_b %in% names(prots.vec)) ~
                prots.vec[item_id_b],
            grepl('^CID', item_id_b) ~ gsub('CID[a-z]0*', 'Metab__', item_id_b),
            TRUE ~ item_id_b
        ),
        sign = case_when(action == 'inhibition' ~ -1, TRUE ~ 1)
    ) %>% select(1, 2, 7)
    colnames(STITCH) <- c('source', 'target', 'sign')
    CIDs <- unique(as.character(unlist(STITCH[,c(1,3)])))
    CIDs <- CIDs[grepl('Metab__', CIDs)] %>% gsub('Metab__', '', .)
    ## Convert CID to HMDB Id when available
    metabolitesMapping.mod <-
        metaboliteIDMapping::metabolitesMapping %>%
        filter(CID %in% CIDs, !is.na(HMDB)) %>%
        mutate(HMDB = paste0('Metab__', HMDB))
    metabolitesMapping.vec <- metabolitesMapping.mod$HMDB
    names(metabolitesMapping.vec) <- paste0('Metab__', metabolitesMapping.mod$CID)
    ## metabolites with no HMDB are kept
    STITCH <- STITCH %>% mutate(
        source = case_when(
            grepl('Metab__', source) & source %in% names(metabolitesMapping.vec) ~
                metabolitesMapping.vec[source],
            TRUE ~ source
        ),
        target = case_when(
            grepl('Metab__', target) & target %in% names(metabolitesMapping.vec) ~
                metabolitesMapping.vec[target],
            TRUE ~ target
        )
    )
    # TODO: at this point, STITCH contains metabolites in both columns of the dataframe
    ## this should be checked
    omn.prots <- unique(as.character(unlist(omnipath.PKN[,c(1, 2)])))
    STITCH <- unique(STITCH[which(STITCH$target %in% omn.prots),])

    STITCH$source <- paste(STITCH$source, '_c', sep = '')

    return(STITCH)
}
