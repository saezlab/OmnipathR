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

SUPPORTED_ORGANISMS <- list(
    gem = c('Human', 'Mouse', 'Rat', 'Zebrafish', 'Fruitfly', 'Worm'),
)

#' Organism-specific prior knowledge networks for COSMOS
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
#'   Supported taxons are 9606 (\emph{Homo sapiens}), 10090
#'   (\emph{Mus musculus}), and 10116 (\emph{Rattus norvegicus}).
#' @param translate.genes Whether translating genes from ENSEMBL into SYMBOL.
#'   Only required when \code{organism == 9606} (\code{FALSE} by default).
#' @param biomart.use.omnipath Whether using BioMart information from OmnipathR
#'   (\code{TRUE} by default).
#' @param gem_reactions.map.col Column of reaction IDs in the GEM
#'   (\code{'rxns'} by default).
#' @param gem_metabolites.map.col Column of reaction IDs in the GEM
#'   (\code{'mets'} by default).
#' @param gem_list.params List containing the name of the slots where the
#'   information to construct the PKN is located in the gem_ If a matlab object
#'   is provided, this list parameter should not be modified.
#' @param gem_degree.mets.threshold Degree cutoff used to filter out
#'   metabolites (400 by default). The objective is to remove cofactors and
#'   over-promiscuous metabolites.
#' @param stitch.threshold Confidence cutoff used for STITCH connections
#'   (700 by default).
#'
#' @return List of 4 elements containing the necessary information for COSMOS to
#'   run: causal PKN, mapping data frame for metabolites from GEM,
#'   reaction-to-gene data frame from GEM, and mapping data frame for reactions
#'   from gem_
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M,
#'      Cholley PE, et al. Genome-scale metabolic network reconstruction of model
#'      animals as a platform for translational research. Proceedings of the
#'      National Academy of Sciences. 2021 Jul 27;118(30):e2102344118.
#'
#'      Türei D, Valdeolivas A, Gul L, Palacio‐Escat N, Klein M, Ivanova O, et al.
#'       Integrated intra‐ and intercellular signaling knowledge for multicellular
#'       omics analysis. Molecular Systems Biology. 2021 Mar;17(3):e9923.
#'
#' @examples
#' \dontrun{
#'   human.PKN.COSMOS <- cosmos_pkn(organism = 9606)
#' }
#'
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom rlang exec !!!
#'
#' @export
cosmos_pkn <- function(
    organism,
    translate.genes = FALSE,
    biomart.use.omnipath = TRUE,
    gem_reactions.map.col = 'rxns',
    gem_metabolites.map.col = 'mets',
    gem_list.params = list(
        stoich.name = 'S',
        reaction.name = 'grRules',
        lb.name = 'lb',
        ub.name = 'ub',
        rev.name = 'rev',
        reaction.ID.name = 'rxns',
        metabolites.ID.name = 'mets',
        metabolites.names.name = 'metNames',
        metabolites.fomulas.name = 'metFormulas',
        metabolites.inchi.name = 'inchis'
    ),
    gem_degree.mets.threshold = 400,
    stitch.threshold = 700
){

    organism %<>%
    c(
        original = .,
        ensembl = ensembl_name(.)
    ) %>%
    {`if`(
        is.na(extract2(., 'ensembl')),
        {
            extract2(., 'original') %>%
            sprintf('Could not recognize organism `%s`.', .) %T>%
            log_error
            stop
        },
        extract2(., 'ensembl')
    )}

    gene_ensembl <- sprintf('%s_gene_ensembl', ens_organism)

    ## check dependencies (Suggests in DESCRIPTION)
    c('R.matlab', 'metaboliteIDMapping') %>%
    missing_packages %>%
    paste(collapse = ', ') %>%
    {`if`(nchar(.), sprintf('Missing packages: %s', .) %T>% log_error %>% stop)}

    .slow_doctest()

    cache_pseudo_url <- 'PKN_COSMOS_%s' %>% sprintf(organism)
    cache_pseudo_post <- list(
        organism = organism,
        translate.genes = translate.genes,
        biomart.use.omnipath = biomart.use.omnipath,
        gem_reactions.map.col = gem_reactions.map.col,
        gem_metabolites.map.col = gem_metabolites.map.col,
        gem_list.params = gem_list.params,
        gem_degree.mets.threshold = gem_degree.mets.threshold,
        stitch.threshold = stitch.threshold
    )

    in_cache <- omnipath_cache_get(
            url = cache_pseudo_url,
            post = cache_pseudo_post,
            create = FALSE
        ) %>% omnipath_cache_latest_version

    if (is.null(in_cache)) {
        log_success(
            paste0(
                'Building COSMOS PKN (organism: %s). ',
                'This will take 10-30 min at the first ',
                'time, and will be saved in the cache for later use.'
            ),
            organism
        )

        res <- exec(.cosmos_pkn, !!!cache_pseudo_post) %>%
            omnipath_cache_save(
                url = cache_pseudo_url,
                post = cache_pseudo_post
            )
    } else {
        res <- omnipath_cache_load(
            url = cache_pseudo_url,
            post = cache_pseudo_post
        )
    }

    return(res)
}



#' Generating COSMOS' PKN for different organisms
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
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#' @param translate.genes Whether translating genes from ENSEMBL into SYMBOL.
#'   \code{FALSE} by default.
#' @param biomart.use.omnipath Whether using BiomaRt information from OmnipathR
#'   (\code{TRUE} by default).
#' @param gem_reactions.map.col Column of reaction IDs in the GEM
#'   (\code{'rxns'} by default).
#' @param gem_metabolites.map.col Column of reaction IDs in the GEM
#'   (\code{'mets'} by default).
#' @param gem_list.params List containing the name of the slots where the
#'   information to construct the PKN is located in the gem_ If a matlab object
#'   is provided, this list parameter should not be modified.
#' @param gem_degree.mets.threshold Degree cutoff used to filter out
#'   metabolites (400 by default). The objective is to remove cofactors and
#'   other metabolites with many connections.
#' @param stitch.threshold Confidence cutoff used for STITCH connections
#'   (700 by default).
#'
#' @return List of 4 elements containing the necessary information for COSMOS to
#'   run: causal PKN, mapping data frame for metabolites from GEM,
#'   reaction-to-gene data frame from GEM, and mapping data frame for reactions
#'   from gem_
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M,
#'      Cholley PE, et al. Genome-scale metabolic network reconstruction of model
#'      animals as a platform for translational research. Proceedings of the
#'      National Academy of Sciences. 2021 Jul 27;118(30):e2102344118.
#'
#'      Türei D, Valdeolivas A, Gul L, Palacio‐Escat N, Klein M, Ivanova O, et al.
#'       Integrated intra‐ and intercellular signaling knowledge for multicellular
#'       omics analysis. Molecular Systems Biology. 2021 Mar;17(3):e9923.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
.cosmos_pkn <- function(
        organism,
        translate.genes = FALSE,
        biomart.use.omnipath = TRUE,
        gem_reactions.map.col = 'rxns',
        gem_metabolites.map.col = 'mets',
        gem_list.params = list(
            stoich.name = 'S',
            reaction.name = 'grRules',
            lb.name = 'lb',
            ub.name = 'ub',
            rev.name = 'rev',
            reaction.ID.name = 'rxns',
            metabolites.ID.name = 'mets',
            metabolites.names.name = 'metNames',
            metabolites.fomulas.name = 'metFormulas',
            metabolites.inchi.name = 'inchis'
        ),
        gem_degree.mets.threshold = 400,
        stitch.threshold = 700
) {
    ## check organisms
    dataset.biomart <- switch(
        as.character(organism),
        '9606' = 'hsapiens_gene_ensembl',
        '10090' = 'mmusculus_gene_ensembl',
        '10116' = 'rnorvegicus_gene_ensembl'
    )
    if (is.null(dataset.biomart))
        stop(
            'Chosen organism is not recognizable Available options are: ',
            paste(c(9606, 10090, 10116, 7955, 7227, 6239), collapse = ', ')
        )

    log_info('Loading GEM obtained from Wang et al., 2021...')
    ## download GEM using OmnipathR (if already done, it will be taken from cache)
    gem_raw <- OmnipathR:::gem_raw(organism = organism) %>% as.data.frame()
    gem_metabolites <- OmnipathR:::gem_metabolites(organism = organism) %>% as.data.frame()
    gem_reactions <- OmnipathR:::gem_reactions(organism = organism) %>% as.data.frame()

    log_info('Loading protein-chemical interactions from STITCH...')
    ## download STITCH using OmnipathR (if already done, it will be taken from cache)
    stitch.actions <- OmnipathR:::stitch_actions(organism = organism) %>%
        as.data.frame()
    stitch.prot.details <- OmnipathR:::stitch_proteins(
        organism = organism
    ) %>% as.data.frame()

    if (biomart.use.omnipath == TRUE) {
        log_info('Using the OmnipathR to retrieve BioMart information')
        ## get info from BiomartR using OmnipathR
        mapping.biomart <- OmnipathR::biomart_query(
            attrs = c(
                'ensembl_peptide_id','ensembl_gene_id', 'external_gene_name'
            ),
            dataset = dataset.biomart
        ) %>% as.data.frame()
    } else {
        log_info('Using the BiomaRt R package when needed')

        mapping.biomart <- NULL
    }
    ## Omnipath data
    log_info('Loading protein-protein interactions from Omnipath...')
    omnipath.PKN <- omnipath_for_cosmos(organism)
    ## Getting GEM PKN
    log_info('Processing gem_..')
    gem_PKN.list <- gem_basal_pkn(
        matlab.object = gem_raw,
        reactions.map = gem_reactions,
        reactions.map.col = gem_reactions.map.col,
        metabolites.map = gem_metabolites,
        metabolites.map.col = gem_metabolites.map.col,
        list.params.GEM = gem_list.params,
        degree.mets.cutoff = gem_degree.mets.threshold
    )
    log_info('Formatting GEM PKN for COSMOS...')
    gem_PKN.list <- translate_to_hmdb(gem_PKN.list)

    if (as.character(organism) == '9606') translate.genes <- TRUE

    if (translate.genes){
        gem_PKN.list <- translate_to_genesymbol(
            gem_PKN.list, organism = organism,
            mapping.biomart = mapping.biomart
        )
    }
    gem_PKN.list <- cosmos_format_gem(gem_PKN.list)

    log_info('Getting STITCH PKN...')

    stitch.PKN <- stitch_format(
        stitch.actions = stitch.actions,
        stitch.links = stitch.prot.details,
        organism = organism,
        omnipath.PKN = omnipath.PKN,
        mapping.biomart = mapping.biomart,
        threshold = stitch.threshold
    )

    output.final <- cosmos_combine_networks(
        gem_network = gem_PKN.list[[1]],
        omnipath.PKN = omnipath.PKN,
        stitch.PKN = stitch.PKN
    )

    log_success('COSMOS PKN ready.')

    return(
        list(
            COSMOS.PKN = output.final,
            gem_mets.map = gem_PKN.list[[2]],
            gem_reac.to.gene = gem_PKN.list[[3]],
            reac.map = gem_PKN.list[[4]]
        )
    )
}


#' OmniPath PPI for the COSMOS PKN
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Data frame with the columns source, target and sign.
#'
#' @noRd
omnipath_for_cosmos <- function(
        organism = 9606
) {
    full_pkn_mm <- as.data.frame(import_omnipath_interactions(organism = organism))
    full_pkn_mm <- full_pkn_mm[!is.na(full_pkn_mm$references),]
    clean_PKN_mm <- full_pkn_mm[
        full_pkn_mm$consensus_stimulation == 1 |
            full_pkn_mm$consensus_inhibition == 1,
    ]
    clean_PKN_mm$sign <- clean_PKN_mm$consensus_stimulation -
        clean_PKN_mm$consensus_inhibition
    clean_PKN_mm <- clean_PKN_mm[, c(3 ,4, 16)]
    clean_PKN_supp_mm <- clean_PKN_mm[clean_PKN_mm$sign == 0,]
    clean_PKN_supp_mm$sign <- -1
    clean_PKN_mm[clean_PKN_mm$sign == 0, 'sign'] <- 1
    clean_PKN_mm <- as.data.frame(rbind(clean_PKN_mm, clean_PKN_supp_mm))
    names(clean_PKN_mm) <- c('source', 'target', 'sign')

    return(clean_PKN_mm)
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
#'   translated into the desired ontology, gene-to-reactions data frame,
#'   metabolite-mapping data frame, and reactions-mapping data frame.
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
            gem_PKN = list.network[[1]],
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
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#' @param mapping.biomart BioMart ontology mapping data frame. If \code{NULL},
#'   this info is obtained using the \ckg{bioMaRt} R package.
#' @param ont.from Ontology to translate genes from
#'   (\code{'ensembl_gene_id' by default}).
#' @param ont.to Ontology to translate genes into
#'   (\code{'external_gene_name' by default}).
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'   translated into the desired ontology, gene-to-reactions data frame,
#'   metabolite-mapping data frame, and reactions-mapping data frame.
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
            gem_PKN = list.network[[1]],
            mets.map = list.network[[2]],
            reac.to.gene = list.network[[3]],
            reac.map = list.network[[4]]
        )
    )
}


#' Connecting PKN derived from GEM with Omnipath protein-protein interactions
#'
#' @param gem_PKN List obtained using \code{gem_basal_pkn}.
#' @param omnipath.PKN Protein-protein interactions obtained using
#'   \code{omnipath_for_cosmos}.
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'   translated into the desired ontology, gene-to-reactions data frame,
#'   metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
cosmos_connect_gem_omnipath <- function(gem_PKN, omnipath.PKN) {
    elements <- unique(as.character(unlist(gem_PKN)))
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
        connectors.df$source %in% omnipath.PKN$source |
            connectors.df$source %in% omnipath.PKN$target
    ),]
    network.df.new <- as.data.frame(rbind(gem_PKN, connectors.df))

    return(network.df.new)
}


#' Combine GEM-, Omnipath- and STITCH- derived PKNs
#'
#' @param gem_network Metabolic PKN obtained using \code{gem_basal_pkn}.
#' @param omnipath.PKN Protein-protein interactions obtained using
#'   \code{omnipath_for_cosmos}.
#' @param stitch.PKN Chemical-protein PKN obtained using \code{stitch_format_gem}.
#'
#' @return List containing PKN with COSMOS and OCEAN format.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
cosmos_combine_networks <- function(gem_network, omnipath.PKN, stitch.PKN) {
    ## connecting Omnipath and GEM
    gem_network <- cosmos_connect_gem_omnipath(
        gem_PKN = gem_network, omnipath.PKN = omnipath.PKN
    )
    gem_network$sign <- 1
    meta.PKN <- as.data.frame(
        rbind(omnipath.PKN, stitch.PKN, gem_network)
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

