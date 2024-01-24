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

#' COSMOS PKN from the genome scale metabolic model in Wang et al. 2021
#'
#' Processing GEMs from Wang et al., 2021
#' (\url{https://github.com/SysBioChalmers}) to generate PKN for COSMOS
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicus), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#' @param reaction_id_type Column from \code{reactions.map} used as ID in
#'     \code{matlab.object}. This parameter should not be modified.
#' @param metabolite_id_type Column from \code{metabolites.map} used as ID in
#'     \code{matlab.object}. This parameter should not be modified.
#' @param list.params.GEM List of parameters to correctly get information from
#'     the matlab object. This parameter should not be modified.
#' @param degree.mets.cutoff Degree cutoff used to prune metabolites with high
#'     degree assuming they are cofactors (400 by default).
#'
#' @return List containing PKN with COSMOS and OCEAN format, gene-to-reactions
#'     data frame, metabolite-mapping data frame, and reactions-mapping data
#'     frame.
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'     PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'     Genome-scale metabolic network reconstruction of model animals as a
#'     platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'     27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% extract extract2 equals
#' @importFrom dplyr pull
#' @importFrom tibble as_tibble
#' @importFrom purrr map
#'
#' @export
chalmers_gem <- function(
        organism = 'Human',
        reaction_id_type = 'rxns',
        metabolite_id_type = 'mets',
        list.params.GEM = list(
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
        metab_max_degree = 400L
    ) {

    log_info('Processing Chalmers SysBio GEM.')

    raw <- gem_raw(organism = organism)
    metabolites <- gem_metabolites(organism = organism)
    reactions <- gem_reactions(organism = organism)

    ## check parameters
    err = NULL
    if (metab_max_degree < 1L) {
        err = '`metab_max_degree` cannot be less than 1.'
    }

    if (!is.null(err)) {
        log_error(err)
        stop(err)
    }


    raw %>%
    as_tibble %>%
    pull(1L) %>%
    extract(,,1L) %>%
    map(
        function(x) {
            if (x %>% dim %>% equals(1L) %>% all) x %>% extract(1L,1L) else x
        }
    ) %>%
    extract(c('rxns', 'S', 'lb', 'ub', 'rev', 'grRules')) %>%
    map(extract, ,1L) %>%
    as_tibble %>%
    mutate(
        rxns = unlist(rxns),
        direction = ifelse(ub + lb >= 0L, 'forward', 'backward'),
        rev = as.logical(rev),
        grRules = (
            unlist(grRules, recursive = FALSE) %>%
            map_chr(extract, 1L) %>%
            str_replace_all(' and ', '_') %>%
            str_replace_all('\\[\\(\\)\\]|_AT\\d+', '') %>%
            str_split(' or ') %>%
            map(str_split, ' and ') %>%
            map(discard, is_empty_2)
        ),
        orphan = map_lgl(
            grRules,
            ~any(map_lgl(.x, ~any(str_detect(.x, '^\\d+$'))))
        )
    )


    ##############################################################################
    ## metabolites
    metabolites.IDs <- unlist(
        matlab.object[[which(attribs.mat == list.params.GEM$metabolites.ID.name)]]
    )
    metabolites.names <- .metab_info(
        matlab.object = matlab.object, attribs.mat = attribs.mat,
        name = list.params.GEM$metabolites.names.name
    )
    ## check if IDs are the same and show number of lost metabolites
    # metabolites.map[[metabolites.map.col]]
    inter.metab <- intersect(
        metabolites.map[[metabolites.map.col]], metabolites.IDs
    )
    rownames(metabolites.map) <- metabolites.map[[metabolites.map.col]]
    metabolites.map <- metabolites.map[metabolites.IDs, ]
    ## adding additional information
    metabolites.formulas <- .metab_info(
        matlab.object = matlab.object, attribs.mat = attribs.mat,
        name = list.params.GEM$metabolites.fomulas.name
    )
    metabolites.inchi <- .metab_info(
        matlab.object = matlab.object, attribs.mat = attribs.mat,
        name = list.params.GEM$metabolites.inchi.name
    )
    metabolites.map <- cbind(
        metabolites.map,
        Metabolite.Name = metabolites.IDs,
        Metabolite.Formula = metabolites.formulas,
        Metabolite.Inchi = metabolites.inchi
    )
    metabolites.map[metabolites.map == ''] <- NA
    ##############################################################################
    ## SIF file: PKN
    log_info('Generating PKN')

    reaction.to.genes.df.reac <- reaction.to.genes.df
    reaction.list <- list()
    for (reac.idx in seq(ncol(s.matrix))) {
        reaction <- s.matrix[, reac.idx]
        #modify gene name so reactions that are catalised by same enzyme stay separated
        reaction.to.genes.df.reac[reaction.to.genes.df$Reaction == reac.idx, 1] <- paste(
            paste0('Gene', reac.idx),
            reaction.to.genes.df[reaction.to.genes.df$Reaction == reac.idx, 1],
            sep = '__'
        )
        # get the enzymes associated with reaction
        genes <- reaction.to.genes.df.reac[reaction.to.genes.df.reac$Reaction == reac.idx, 1]
        if (as.vector(lbs[reac.idx, 4] == 'forward')) {
            reactants <- metabolites.IDs[reaction == -1]
            products <- metabolites.IDs[reaction == 1]
        } else {
            reactants <- metabolites.IDs[reaction == 1]
            products <- metabolites.IDs[reaction == -1]
        }
        reactants <- paste0('Metab__', reactants)
        products <- paste0('Metab__', products)
        number_of_interations <- length(reactants) + length(products)
        # now for each enzyme, we create a two column dataframe recapitulating the
        # interactions between the metabolites and this enzyme
        reaction.df <- lapply(
            X = as.list(genes),
            FUN = \(gene) {
                gene.df <- data.frame(
                    # reactants followed by the enzyme (the enzyme is repeated as many time as they are products)
                    source = c(reactants, rep(gene, number_of_interations - length(reactants))),
                    # enzyme(repeated as many time as they are reactants) followed by products
                    target = c(rep(gene, number_of_interations - length(products)), products)
                )
                if (reversible[reac.idx]) {
                    gene.df.reverse <- data.frame(
                        source = c(
                            rep(
                                paste(gene, '_reverse', sep = ''),
                                number_of_interations - length(products)
                            ),
                            products
                        ),
                        target = c(
                            reactants,
                            rep(
                                paste(gene, '_reverse', sep = ''),
                                number_of_interations - length(reactants)
                            )
                        )
                    )
                    gene.df <- rbind(gene.df, gene.df.reverse)
                }
                return(gene.df)
            }
        ) %>% do.call(rbind, .)
        reaction.list[[reac.idx]] <- reaction.df
    }
    reaction.df.all <- do.call(rbind, reaction.list)
    ## removing those reactions with no metab <--> gene
    reaction.df.all <- reaction.df.all[reaction.df.all$source != 'Metab__' &
                                                                             reaction.df.all$target != 'Metab__',]
    ## only complete cases
    reaction.df.all <- reaction.df.all[complete.cases(reaction.df.all),]
    ##############################################################################
    ## removing cofactors (metabolites with a high degree)
    metabs.degree <- sort(
        table(
            grep(
                '^Metab__', c(reaction.df.all$source, reaction.df.all$target),
                value = TRUE
            )
        ),
        decreasing = TRUE
    )
    log_info(
        'Number of metabolites removed after degree >',
        degree.mets.cutoff,  ': ', sum(metabs.degree >= degree.mets.cutoff)
    )

    metabs.degree.f <- metabs.degree[metabs.degree < degree.mets.cutoff]
    reactions.df.no.cofac <- reaction.df.all[
        reaction.df.all$source %in% names(metabs.degree.f) |
            reaction.df.all$target %in% names(metabs.degree.f),
    ]
    mets <-
        grep(
            pattern = 'Metab__',
            x = unique(c(reactions.df.no.cofac[[1]], reactions.df.no.cofac[[1]])),
            value = TRUE
        ) %>%
        gsub('Metab__', '', .)

    metabolites.map <- metabolites.map[mets, ]
    log_info('Final number of connections: ', nrow(reactions.df.no.cofac))

    return(
        list(
            gem_pkn = reactions.df.no.cofac,
            mets.map = metabolites.map,
            reac.to.gene = reaction.to.genes.df.reac,
            reac.map = reactions.map
        )
    )
}


#' Formatting PKN derived from GEM for COSMOS
#'
#' It determines and marks transporters and reverse reactions.
#'
#' @param list.network List obtained using \code{gem_basal_pkn}.
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'   translated into the desired ontology, gene-to-reactions data frame,
#'   metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @noRd
cosmos_format_gem <- function(list.network) {
    reaction.network <- list.network[[1]]
    enzyme_reacs <- unique(c(reaction.network$source, reaction.network$target))
    enzyme_reacs <- enzyme_reacs[grepl('^Gene', enzyme_reacs)]
    enzyme_reacs_reverse <- enzyme_reacs[grepl('_reverse',enzyme_reacs)]
    enzyme_reacs <- enzyme_reacs[!grepl('_reverse',enzyme_reacs)]

    log_info('Step 1: Defining transporters')

    new_df_list <- sapply(
        X = enzyme_reacs,
        FUN = function(enzyme_reac, reaction.network) {
            df <- reaction.network[which(
                reaction.network$source == enzyme_reac |
                    reaction.network$target == enzyme_reac
            ),]
            if (dim(df)[1] < 2) {
                return(NA)
            } else {
                if (dim(df)[1] < 3) {
                    return(df)
                } else {
                    for(i in 1:dim(df)[1]) {
                        if(grepl('Metab__', df[i, 1])) {
                            counterpart <- which(
                                gsub('_[a-z]$','',df[,2]) == gsub('_[a-z]$','',df[i,1])
                            )
                            if(length(counterpart) > 0) {
                                df[i, 2] <- paste0(df[i, 2], paste0('_TRANSPORTER', i))
                                df[counterpart, 1] <- paste0(
                                    df[counterpart, 1], paste0('_TRANSPORTER', i)
                                )
                            }
                        }
                    }
                    return(df)
                }
            }
        },
        reaction.network = reaction.network
    )
    new_df <- as.data.frame(do.call(rbind, new_df_list))

    log_info('Step 2: Defining reverse reactions')

    new_df_reverse <- sapply(
        X = enzyme_reacs_reverse,
        FUN = function(enzyme_reac_reverse, reaction.network) {
            df <- reaction.network[which(
                reaction.network$source == enzyme_reac_reverse |
                    reaction.network$target == enzyme_reac_reverse
            ),]
            if(dim(df)[1] < 2) {
                return(NA)
            } else {
                if(dim(df)[1] < 3) {
                    return(df)
                } else {
                    for(i in 1:dim(df)[1]) {
                        if(grepl('Metab__',df[i,1])) {
                            counterpart <- which(
                                gsub('_[a-z]$','',df[,2]) == gsub('_[a-z]$','',df[i,1])
                            )
                            if(length(counterpart) > 0) {
                                transporter <- gsub('_reverse', '', df[i, 2])
                                transporter <- paste0(
                                    transporter, paste0(paste0('_TRANSPORTER', i), '_reverse')
                                )
                                df[i, 2] <- transporter
                                df[counterpart, 1] <- transporter
                            }
                        }
                    }
                    return(df)
                }
            }
        }, reaction.network = reaction.network
    )
    new_df_reverse <- as.data.frame(do.call(rbind, new_df_list))
    reaction.network.new <- as.data.frame(rbind(new_df, new_df_reverse))
    reaction.network.new <- reaction.network.new[complete.cases(reaction.network.new),]
    ## filter metabolites in mapping mets
    metabs <- c(
        grep('Metab__', reaction.network.new[[1]], value = TRUE),
        grep('Metab__', reaction.network.new[[2]], value = TRUE)
    ) %>% unique() %>% gsub('(Metab__)|(_[a-z])', '', .)
    list.network[[2]] <- list.network[[2]] %>%
        filter(metHMDBID %in% metabs | metBiGGID %in% metabs | mets %in% metabs)

    return(
        list(
            gem_pkn = reaction.network.new,
            mets.map = list.network[[2]],
            reac.to.gene = list.network[[3]],
            reac.map = list.network[[4]]
        )
    )
}


#' Reactions from the Chalmers SysBio GEM (Wang et al., 2021)
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Data frame of reaction identifiers.
#'
#' @examples
#' chalmers_gem_reactions()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'   PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'   Genome-scale metabolic network reconstruction of model animals as a
#'   platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'   27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_metabolites}}}
#'     \item{\code{\link{chalmers_gem_raw}}}
#'     \item{\code{\link{chalmers_gem_cosmos}}}
#' }
chalmers_gem_reactions <- function(organism = 'Human') {

    .slow_doctest()

    organism %<>% organism_for('chalmers-gem')

    'chalmers_gem' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = organism %>% list('reactions.tsv'),
        reader_param = list(),
        resource = 'Chalmers GEM',
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' Metabolites from the Chalmers SysBio GEM (Wang et al., 2021)
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Data frame of metabolite identifiers.
#'
#' @examples
#' chalmers_gem_metabolites()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'   PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'   Genome-scale metabolic network reconstruction of model animals as a
#'   platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'   27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom readr read_tsv
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_reactions}}}
#'     \item{\code{\link{chalmers_gem_raw}}}
#'     \item{\code{\link{chalmers_gem_cosmos}}}
#' }
chalmers_gem_metabolites <- function(organism = 'Human') {

    .slow_doctest()

    organism %<>% organism_for('chalmers-gem')

    suppressWarnings(generic_downloader(
        'chalmers_gem',
        reader = read_tsv,
        url_key_param = list(),
        url_param = organism %>% list('metabolites.tsv'),
        reader_param = list(),
        resource = 'Chalmers GEM',
        post = NULL,
        use_httr = FALSE
    )) %T>%
    load_success()

}


#' GEM matlab file from Chalmers Sysbio (Wang et al., 2021)
#'
#' Downloads and imports the matlab file containing the genome scale metabolic
#' models created by Chalmers SysBio.
#'
#' @param organism Character or integer: name or identifier of the organism.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Matlab object containing the GEM.
#'
#' @examples
#' chalmers_gem_raw()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'   PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'   Genome-scale metabolic network reconstruction of model animals as a
#'   platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'   27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>% %<>%
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_reactions}}}
#'     \item{\code{\link{chalmers_gem_metabolites}}}
#'     \item{\code{\link{chalmers_gem_cosmos}}}
#' }
chalmers_gem_raw <- function(organism = 'Human') {

    .slow_doctest()

    organism %<>%
        organism_for('chalmers-gem') %T>%
        log_info('Downloading GEM from Chalmers Sysbio for organism `%s`.', .)

    'chalmers_gem' %>%
    generic_downloader(
        reader = R.matlab::readMat,
        url_key_param = list(),
        url_param = organism %>% list(sprintf('%s-GEM.mat', .)),
        reader_param = list(),
        resource = 'Chalmers GEM',
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' Keep entries without elements in GEM processing
#'
#' @param gem_raw GEM from a Matlab object
#' @param attribs_mat Atribute from the Matlab object to be parsed
#' @param name Vector of elements to be checked in the Matlab object
#'
#' @return Vector with NAs in those entries with no element
#'
#' @noRd
metab_info <- function(gem_raw, attribs_mat, name) {
    unlist(
        sapply(
            X = gem_raw[[which(attribs_mat == name)]],
            FUN = \(elem) {
                if (length(unlist(elem) != 0)) {
                    return(elem)
                } else {
                    return(NA)
                }
            }
        )
    )
}
