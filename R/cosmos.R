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

#' Prior knowledge network (PKN) for COSMOS
#'
#' This function generates the prior knowledge network (PKN) needed to run
#' COSMOS using information from different resources through the \pkg{OmnipathR}
#' R package. Particularly, \code{cosmos_pkn} will obtain:
#' \itemize{ \item Genome-scale metabolic
#' model (GEM) of the required organism from Wang et al., 2021.
#' \item Interaction network of
#' chemicals and proteins from STITCH (\url{http://stitch.embl.de/}) for the
#' required organism. \item Protein-protein interactions from Omnipath (Türei
#' et al., 2021) for the required organism} With these three pieces of
#' information, the function will generate the required causal network for
#' COSMOS to run.
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'       Supported taxons are 9606 (\emph{Homo sapiens}), 10090
#'       (\emph{Mus musculus}), and 10116 (\emph{Rattus norvegicus}).
#' @param translate.genes Whether translating genes from ENSEMBL into SYMBOL.
#'     Only required when \code{organism == 9606} (\code{FALSE} by default).
#' @param biomart.use.omnipath Whether using BioMart information from OmnipathR
#'     (\code{TRUE} by default).
#' @param gem_reactions.map.col Column of reaction IDs in the GEM
#'     (\code{'rxns'} by default).
#' @param gem_metabolites.map.col Column of reaction IDs in the GEM
#'     (\code{'mets'} by default).
#' @param gem_list.params List containing the name of the slots where the
#'     information to construct the PKN is located in the gem_ If a matlab object
#'     is provided, this list parameter should not be modified.
#' @param gem_degree.mets.threshold Degree cutoff used to filter out
#'     metabolites (400 by default). The objective is to remove cofactors and
#'     over-promiscuous metabolites.
#' @param stitch.threshold Confidence cutoff used for STITCH connections
#'     (700 by default).
#'
#' @return List of 4 elements containing the necessary information for COSMOS to
#'     run: causal PKN, mapping data frame for metabolites from GEM,
#'     reaction-to-gene data frame from GEM, and mapping data frame for reactions
#'     from gem_
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M,
#'        Cholley PE, et al. Genome-scale metabolic network reconstruction of model
#'        animals as a platform for translational research. Proceedings of the
#'        National Academy of Sciences. 2021 Jul 27;118(30):e2102344118.
#'
#'        Türei D, Valdeolivas A, Gul L, Palacio‐Escat N, Klein M, Ivanova O, et al.
#'         Integrated intra‐ and intercellular signaling knowledge for multicellular
#'         omics analysis. Molecular Systems Biology. 2021 Mar;17(3):e9923.
#'
#' @examples
#' \dontrun{
#'     human.PKN.COSMOS <- cosmos_pkn(organism = 9606)
#' }
#'
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom rlang exec !!!
#'
#' @export
cosmos_pkn <- function(
    organism,
    gene_id_type = NULL,
    chalmers_gem_metab_max_degree = 400L,
    stitch_score = 700L
){

    .slow_doctest()

    organism %<>% ncbi_taxid()

    args <- environment() %>% as.list

    ## check dependencies (Suggests in DESCRIPTION)
    c('R.matlab', 'metaboliteIDMapping') %>%
    missing_packages %>%
    paste(collapse = ', ') %>%
    {`if`(nchar(.), sprintf('Missing packages: %s', .) %T>% log_error %>% stop)}

    cache_pseudo_url <- 'COSMOS_PKN_%s' %>% sprintf(organism)

    in_cache <-
        omnipath_cache_get(
            url = cache_pseudo_url,
            post = args,
            create = FALSE
        ) %>%
        omnipath_cache_latest_version

    if (is.null(in_cache)) {

        log_success(
            paste0(
                'Building COSMOS PKN (organism: %s). ',
                'This will take 10-30 min at the first ',
                'time, and will be saved in the cache for later use.'
            ),
            organism
        )

        result <-
            exec(.cosmos_pkn, !!!cache_pseudo_post) %>%
            omnipath_cache_save(url = cache_pseudo_url, post = args)

    } else {

        result <- omnipath_cache_load(url = cache_pseudo_url, post = args)

    }

    return(result)

}



#' Build prior knowledge network for COSMOS
#'
#' This function generates the prior knowledge network (PKN) needed to run
#' COSMOS using information from different resources. It will download the
#' required information through the \pkg{OmnipathR} R package. Particularly,
#' \code{cosmos_pkn} will obtain: \itemize{ \item Genome-scale metabolic
#' model (GEM) of the required organism from Wang et al., 2021.
#' \item Interaction network of
#' chemical and proteins from STITCH (\url{http://stitch.embl.de/}) for the
#' required organism. \item Protein-protein interactions from Omnipath (Türei
#' et al., 2021)} With these three pieces of information, the function will
#' generate the required causal network for COSMOS to run.
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#' @param gene_id_type Whether translating genes from ENSEMBL into SYMBOL.
#'     \code{FALSE} by default.
#' @param gem_reactions.map.col Column of reaction IDs in the GEM
#'     (\code{'rxns'} by default).
#' @param gem_metabolites.map.col Column of reaction IDs in the GEM
#'     (\code{'mets'} by default).
#' @param gem_list.params List containing the name of the slots where the
#'     information to construct the PKN is located in the gem_ If a matlab object
#'     is provided, this list parameter should not be modified.
#' @param met_max_degree Degree cutoff used to filter out
#'     metabolites (400 by default). The objective is to remove cofactors and
#'     other metabolites with many connections.
#' @param stitch_score Confidence cutoff used for STITCH connections
#'     (700 by default).
#'
#' @return List of 4 elements containing the necessary information for COSMOS to
#'     run: causal PKN, mapping data frame for metabolites from GEM,
#'     reaction-to-gene data frame from GEM, and mapping data frame for reactions
#'     from gem_
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M,
#'        Cholley PE, et al. Genome-scale metabolic network reconstruction of model
#'        animals as a platform for translational research. Proceedings of the
#'        National Academy of Sciences. 2021 Jul 27;118(30):e2102344118.
#'
#'        Türei D, Valdeolivas A, Gul L, Palacio‐Escat N, Klein M, Ivanova O, et al.
#'         Integrated intra‐ and intercellular signaling knowledge for multicellular
#'         omics analysis. Molecular Systems Biology. 2021 Mar;17(3):e9923.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
.cosmos_pkn <- function(
        organism,
        gene_id_type = NULL,
        chalmers_gem_metab_max_degree = 400L,
        stitch_score = 700L
) {

    organism %<>% organism_for('cosmos')

    ## download STITCH using OmnipathR (if already done, it will be taken from cache)
    stitch <- stitch_gem(organism = organism, stitch_score = stitch_score)
    chalmers <- chalmers_gem_network(
        organism = organism,
        metab_max_degree = chalmers_gem_metab_max_degree
    )
    ominpath <- omnipath_for_cosmos(organism)

    cosmos_combine_networks(
        chalmers = chalmers,
        omnipath = omnipath,
        stitch = stitch
    ) %T>%
    {log_success('COSMOS PKN ready, %i interactions.', nrow(.))}

}


#' OmniPath PPI for the COSMOS PKN
#'
#' @param organism Character or integer: name or NCBI Taxonomy ID of the
#'       organism.
#' @param resources Character: names of one or more resources. Correct spelling
#'       is important.
#' @param datasets Character: one or more network datasets in OmniPath.
#' @param interaction_types Character: one or more interaction type
#' @param ... Further parameters to \code{\link{import_omnipath_interactions}}.
#'
#' @return Data frame with the columns source, target and sign.
#'
#' @importFrom magrittr %<>% %>% %T>%
#' @importFrom logger log_info
#' @importFrom dplyr mutate filter select bind_rows select
#' @noRd
omnipath_for_cosmos <- function(
        organism = 9606L,
        resources = NULL,
        datasets = NULL,
        interaction_types = NULL,
        ...
    ) {

    # NSE vs. R CMD check workaround
    consensus_stimulation <- consensus_inhibition <- NULL

    organism %<>% organism_supported('omnipath')
    paste0(
        'OmniPath network for COSMOS PKN; datasets: %s; ',
        'resources: %s; interaction types: %s; organism: %s.'
    ) %>%
    log_info(
        if_null_len0_2(enum_format(datasets), 'omnipath'),
        if_null_len0_2(enum_format(resources), 'all'),
        if_null_len0_2(enum_format(interaction_types), 'post-translational (PPI)'),
        common_name(organism)
    )

    import_omnipath_interactions(
        organism = organism,
        resources = resources,
        datasets = datasets,
        types = interaction_types,
        ...
    ) %>%
    filter(
        !is.na(references) &
        (
            consensus_stimulation == 1 |
            consensus_inhibition == 1
        )
    ) %>%
    mutate(sign = consensus_stimulation - consensus_inhibition) %>%
    select(
        source = source_genesymbol,
        target = target_genesymbol,
        sign
    ) %>%
    bind_rows(
        mutate(., sign = ifelse(sign, sign, 1)),
        filter(., sign == 0) %>% mutate(sign = -1)
    ) %T>%
    {log_info('OmniPath PPI for COSMOS PKN ready: %i interactions.', nrow(.))}

}


#' Translating metabolites from the Metabolic Atlas ID to the HMDB and KEGG IDs
#'
#' It transforms the Metabolic Atlas IDs into the HMDBs IDs. If a metabolite has
#' no HMDB ID, then the KEGG IDs are used. In case there is not entry for HMDB
#' or KEGG, the metabolite keeps the Metabolic Atlas ID.
#'
#' @param list.network List obtained using \code{gem_basal_pkn}.
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'     translated into the desired ontology, gene-to-reactions data frame,
#'     metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr case_when mutate
#' @importFrom stringr str_sub
#'
#' @noRd
translate_to_hmdb <- function(list.network) {
    metab.map <- list.network[[2]]
    list.network[[1]] <- list.network[[1]] %>% mutate(
        source = case_when(
            grepl('Metab__', source) ~ case_when(
                !is.na(
                    metab.map[gsub('Metab__', replacement = '', x = source), 'metHMDBID']
                ) ~ paste0(
                    'Metab__',
                    metab.map[gsub('Metab__', replacement = '', x = source), 'metHMDBID'],
                    '_', str_sub(source, start = nchar(source))
                ),
                !is.na(
                    metab.map[gsub('Metab__', replacement = '', x = source), 'metKEGGID']
                ) ~ paste0(
                    'Metab__',
                    metab.map[gsub('Metab__', replacement = '', x = source), 'metKEGGID'],
                    '_', str_sub(source, start = nchar(source))
                ), TRUE ~ source
            ), TRUE ~ source
        ),
        target = case_when(
            grepl('Metab__', target) ~ case_when(
                !is.na(
                    metab.map[gsub('Metab__', replacement = '', x = target), 'metHMDBID']
                ) ~ paste0(
                    'Metab__',
                    metab.map[gsub('Metab__', replacement = '', x = target), 'metHMDBID'],
                    '_', str_sub(target, start = nchar(target))
                ),
                !is.na(
                    metab.map[gsub('Metab__', replacement = '', x = target), 'metKEGGID']
                ) ~ paste0(
                    'Metab__',
                    metab.map[gsub('Metab__', replacement = '', x = target), 'metKEGGID'],
                    '_', str_sub(target, start = nchar(target))
                ), TRUE ~ target
            ), TRUE ~ target
        )
    )

    return(
        list(
            gem_pkn = list.network[[1]],
            mets.map = metab.map,
            reac.to.gene = list.network[[3]],
            reac.map = list.network[[4]]
        )
    )
}


#' Translating gene names from one ontology to another
#'
#' @param list.network List obtained using \code{gem_basal_pkn}.
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#' @param mapping.biomart BioMart ontology mapping data frame. If \code{NULL},
#'     this info is obtained using the \ckg{bioMaRt} R package.
#' @param ont.from Ontology to translate genes from
#'     (\code{'ensembl_gene_id' by default}).
#' @param ont.to Ontology to translate genes into
#'     (\code{'external_gene_name' by default}).
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'     translated into the desired ontology, gene-to-reactions data frame,
#'     metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr case_when mutate select
#' @importFrom stringr str_sub
#'
#' @noRd
translate_to_genesymbol <- function(
        list.network,
        organism,
        mapping.biomart = NULL,
        ont.from = 'ensembl_gene_id',
        ont.to = 'external_gene_name'
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

    regex <- '(Gene\\d+__)|(_reverse)'
    ## getting biomart info
    genes.GEM <- grep('Gene\\d+__', unlist(list.network[[1]]), value = TRUE) %>%
        gsub(regex, '', .) %>%
        ifelse(grepl('_', .), sapply(strsplit(., split = '_'), \(x) x), .) %>%
        unlist()
    if (is.null(mapping.biomart)) {
        ensembl.link <- useEnsembl(biomart = 'ensembl', dataset = dataset.biomart)
        ensembl.df <- getBM(
            filters = ont.from,
            attributes = c('ensembl_gene_id', 'external_gene_name'),
            values = genes.GEM,
            mart = ensembl.link
        ) %>% unique()
        rownames(ensembl.df) <- ensembl.df$ensembl_gene_id
    } else {
        ensembl.df <- mapping.biomart %>% select(-ensembl_peptide_id) %>%
            unique() %>% filter(
                !is.na(.data[[ont.from]]), !is.na(.data[[ont.to]]),
            )
        rownames(ensembl.df) <- ensembl.df[[ont.from]]
    }
    ## translating genes when possible (not found genes are not removed)
    ## when complexes are present (several genes concatenated), this code does not work
    list.network[[1]] <- list.network[[1]] %>% mutate(
        source = case_when(
            ## cases with a single gene
            grepl('Gene\\d+__', source) ~ case_when(
                !is.na(
                    ensembl.df[gsub(regex, replacement = '', x = source), ont.to]
                ) ~ paste0(
                    'Gene',
                    gsub('\\D', '', sapply(strsplit(x = source, split = '__'), \(x) x[1])),
                    '__',
                    ensembl.df[gsub(regex, replacement = '', x = source), ont.to]
                ),
                ## cases with complexes: more than 1 gene
                grepl('[0-9]_[E]', source) ~
                    paste0(
                        'Gene',
                        gsub('\\D', '', sapply(strsplit(x = target, split = '__'), \(x) x[1])),
                        '__',
                        unlist(
                            strsplit(
                                gsub(
                                    pattern = 'reverse', replacement = '',
                                    grep('[0-9]_[E]', source, value = T)[1]
                                ),
                                split = '_'
                            )
                        )[-c(1:2)] %>% ensembl.df[., ont.to] %>% paste(collapse = '_')
                    ),
                TRUE ~ source
            ), TRUE ~ source
        ),
        target = case_when(
            ## cases with a single gene
            grepl('Gene\\d+__', target) ~ case_when(
                !is.na(
                    ensembl.df[gsub(regex, replacement = '', x = target), ont.to]
                ) ~ paste0(
                    'Gene',
                    gsub('\\D', '', sapply(strsplit(x = target, split = '__'), \(x) x[1])),
                    '__',
                    ensembl.df[gsub(regex, replacement = '', x = target), ont.to]
                ),
                ## cases with complexes: more than 1 gene
                grepl('[0-9]_[E]', target) ~
                    paste0(
                        'Gene',
                        gsub('\\D', '', sapply(strsplit(x = target, split = '__'), \(x) x[1])),
                        '__',
                        unlist(
                            strsplit(
                                gsub(
                                    pattern = 'reverse', replacement = '',
                                    grep('[0-9]_[E]', target, value = T)[1]
                                ),
                                split = '_'
                            )
                        )[-c(1:2)] %>% ensembl.df[., ont.to] %>% paste(collapse = '_')
                    ),
                TRUE ~ target
            ), TRUE ~ target
        )
    )
    list.network[[3]] <- list.network[[3]] %>% mutate(
        Gene = case_when(
            !is.na(
                ensembl.df[gsub(regex, replacement = '', x = Gene), ont.to]
            ) ~ paste0(
                'Gene',
                gsub('\\D', '', sapply(strsplit(x = Gene, split = '__'), \(x) x[1])),
                '__',
                ensembl.df[gsub(regex, replacement = '', x = Gene), ont.to]
            ),
            TRUE ~ Gene
        )
    )

    return(
        list(
            gem_pkn = list.network[[1]],
            mets.map = list.network[[2]],
            reac.to.gene = list.network[[3]],
            reac.map = list.network[[4]]
        )
    )
}


#' Connecting PKN derived from GEM with Omnipath protein-protein interactions
#'
#' @param gem_pkn List obtained using \code{gem_basal_pkn}.
#' @param omnipath_pkn Protein-protein interactions obtained using
#'     \code{omnipath_for_cosmos}.
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'     translated into the desired ontology, gene-to-reactions data frame,
#'     metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
cosmos_connect_gem_omnipath <- function(gem_pkn, omnipath_pkn) {
    elements <- unique(as.character(unlist(gem_pkn)))
    elements <- elements[grepl('^Gene\\d+__', elements)]
    elements <- gsub('(.*__)|(_TRANSPORTER[0-9]+)|(_reverse$)', '', elements)
    ## this function can be vectorized
    connectors.df <- sapply(
        X = elements,
        FUN = function(ele) {
            if (grepl('_', ele)) {
                genes.sep <- str_split(string = ele, pattern = '_')[[1]]
                if(length(genes.sep) < 10) {
                    genes_connector_list <- sapply(
                        X = genes.sep,
                        FUN = function(gene) {
                            return(c(gene, ele))
                        }
                    )
                    return(t(genes_connector_list))
                }
            } else {
                return(c(ele, ele))
            }
        }
    ) %>% do.call(rbind, .) %>% as.data.frame()
    names(connectors.df) <- c('source', 'target')
    connectors.df <- connectors.df[which(
        connectors.df$source %in% omnipath_pkn$source |
            connectors.df$source %in% omnipath_pkn$target
    ),]
    network.df.new <- as.data.frame(rbind(gem_pkn, connectors.df))

    return(network.df.new)
}


#' Combine GEM-, Omnipath- and STITCH- derived PKNs
#'
#' @param gem_pkn Metabolic PKN obtained using \code{gem_basal_pkn}.
#' @param omnipath_pkn Protein-protein interactions obtained using
#'     \code{omnipath_for_cosmos}.
#' @param stitch_pkn Chemical-protein PKN obtained using \code{stitch_format_gem}.
#'
#' @return List containing PKN with COSMOS and OCEAN format.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
cosmos_combine_networks <- function(gem_pkn, omnipath_pkn, stitch_pkn) {
    ## connecting Omnipath and GEM
    gem_pkn <- cosmos_connect_gem_omnipath(
        gem_pkn = gem_pkn, omnipath_pkn = omnipath_pkn
    )
    gem_pkn$sign <- 1
    meta.PKN <- as.data.frame(
        rbind(omnipath_pkn, stitch_pkn, gem_pkn)
    ) %>% unique()
    meta.PKN <- meta.PKN[, c(1, 3, 2)]
    names(meta.PKN) <- c('source', 'interaction', 'target')
    #TODO: manual correction: difficult to generalize for different organisms, shall I remove it?
    #probably erroneous interaction (WHY?? this only works for human / mouse)
    # meta.network <- meta.network[-which(
    #       meta.network$source == 'Prkca' & meta.network$target == 'Src'
    # ),]
    # meta.PKN <- meta.PKN[-which(
    #       meta.PKN$source == 'Prkca' & meta.PKN$target == 'Src'
    # ),]
    #probably erroneous interaction
    # meta_network <- meta_network[-which(meta.network$source == 'Ltc45'),]
    #I don't know where this interaction comes from, the sources are wrong (https://www.nature.com/articles/onc2008228)
    # meta_network <- meta.network[!(grepl('Cad_reverse', meta.network$source) | grepl('Cad_reverse', meta.network$target)) ,]
    #redHuman confirms that the reaction is actually not reversible: NOT FOUND IN MOUSE EITHER

    return(meta.PKN)
}

