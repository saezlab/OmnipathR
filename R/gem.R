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
#' @param matlab.object Matlab object containing a gem_
#' @param reactions.map Data frame with information to map reaction names using
#'   different ontologies.
#' @param metabolites.map Data frame with information to map metabolite names
#'   using different ontologies.
#' @param reactions.map.col Column from \code{reactions.map} used as ID in
#'   \code{matlab.object}. This parameter should not be modified.
#' @param metabolites.map.col Column from \code{metabolites.map} used as ID in
#'   \code{matlab.object}. This parameter should not be modified.
#' @param list.params.GEM List of parameters to correctly get information from
#'   the matlab object. This parameter should not be modified.
#' @param degree.mets.cutoff Degree cutoff used to prune metabolites with high
#'   degree assuming they are cofactors (400 by default).
#'
#' @return List containing PKN with COSMOS and OCEAN format, gene-to-reactions
#'   data frame, metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'   PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'   Genome-scale metabolic network reconstruction of model animals as a
#'   platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'   27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
gem_basal_pkn <- function(
        matlab.object,
        reactions.map,
        metabolites.map,
        reactions.map.col = 'rxns',
        metabolites.map.col = 'mets',
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
        degree.mets.cutoff = 400
) {
    ## check parameters
    if (!reactions.map.col %in% colnames(reactions.map)) {
        stop('reactions.map.col cannot be found in reactions.map data.frame')
    } else if (!metabolites.map.col %in% colnames(metabolites.map)) {
        stop('metabolites.map.col cannot be found in metabolites.map data.frame')
    } else if (degree.mets.cutoff < 1) {
        stop('degree.mets.cutoff cannot be less than 1')
    }

    attribs.mat <- rownames(matlab.object)
    matlab.object <- matlab.object[[1]]
    ## check elements are in the object
    invisible(
        sapply(
            names(list.params.GEM), \(idx) {
                if (!list.params.GEM[[idx]] %in% attribs.mat) {
                    stop(
                        paste0(
                            idx, 'element in list.params.GEM (',
                            list.params.GEM[[idx]], ') is not in matlab object'
                        )
                    )
                }
            }
        )
    )
    ## obtaining data
    s.matrix <- matlab.object[[which(attribs.mat == list.params.GEM$stoich.name)]]
    reaction.list <- matlab.object[[which(attribs.mat == list.params.GEM$reaction.name)]]
    ##############################################################################
    ## reactions
    # direction reactions
    lbs <- as.data.frame(
        cbind(
            matlab.object[[which(attribs.mat == list.params.GEM$lb.name)]],
            matlab.object[[which(attribs.mat == list.params.GEM$ub.name)]],
            matlab.object[[which(attribs.mat == list.params.GEM$rev.name)]]
        )
    )
    ## this could be done with mutate
    lbs$direction <- ifelse(
        (matlab.object[[which(attribs.mat == list.params.GEM$ub.name)]] +
             matlab.object[[which(attribs.mat == list.params.GEM$lb.name)]]) >= 0,
        'forward', 'backward'
    )
    reversible <- ifelse(
        matlab.object[[which(attribs.mat == list.params.GEM$rev.name)]] == 1,
        TRUE, FALSE
    )
    reaction.ids <- unlist(
        matlab.object[[which(attribs.mat == list.params.GEM$reaction.ID.name)]]
    )
    ## reaction to genes df
    reaction.to.genes.df <- lapply(
        seq_along(reaction.list),
        \(idx) {
            genes.reac <- unlist(reaction.list[[idx]], recursive = FALSE)
            if (length(genes.reac) != 0) {
                genes <- unique(
                    gsub(
                        ' and ', '_',
                        gsub(
                            '[()]','',
                            gsub('_AT[0-9]+','', strsplit(genes.reac, split = ' or ')[[1]])
                        )
                    )
                )
                return(
                    data.frame(
                        Gene = genes, Reaction = rep(idx, length(genes)),
                        Reaction.ID = rep(reaction.ids[idx], length(genes))
                    )
                )
            } else {
                return(
                    data.frame(
                        Gene = idx, Reaction = idx, Reaction.ID = reaction.ids[idx]
                    )
                )
            }
        }
    ) %>% do.call(rbind, .)
    orphan.reacts <- grepl(pattern = '^\\d+$', reaction.to.genes.df$Gene)
    reaction.to.genes.df[orphan.reacts, 'Reaction.ID'] <- paste0(
        'orphanReac.', reaction.to.genes.df[orphan.reacts, 'Reaction.ID']
    )
    reaction.to.genes.df[orphan.reacts, 'Gene'] <- paste0(
        reaction.to.genes.df[orphan.reacts, 'Gene'], '.',
        reaction.to.genes.df[orphan.reacts, 'Reaction.ID']
    )
    reaction.to.genes.df <- unique(reaction.to.genes.df)


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
            gem_PKN = reactions.df.no.cofac,
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
#' @param list.network List obtained using \code{.create_gem_basal_PKN}.
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'   translated into the desired ontology, gene-to-reactions data frame,
#'   metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @noRd
.format_gem_cosmos <- function(list.network) {
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
            gem_PKN = reaction.network.new,
            mets.map = list.network[[2]],
            reac.to.gene = list.network[[3]],
            reac.map = list.network[[4]]
        )
    )
}


#' Download GEM reactions file from Wang et al., 2021
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Mapping reactions data frame for GEM processing.
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
#' @noRd
gem_reacts <- function(organism) {

    dataset.github <- switch(
        as.character(organism),
        '9606' = 'Human',
        '10090' = 'Mouse',
        '10116' = 'Rat',
        '7955' = 'Zebrafish',
        '7227' = 'Fruitfly',
        '6239' = 'Worm'
    )

    .slow_doctest()

    'gem_github' %>%
        generic_downloader(
            reader = read_tsv,
            url_key_param = list(),
            url_param = list(dataset.github, 'reactions.tsv'),
            reader_param = list(),
            resource = NULL,
            post = NULL,
            use_httr = FALSE
        ) %T>% load_success()
}


#' Download GEM metabolites file from Wang et al., 2021
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Mapping metabolites data frame for GEM processing.
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
#' @noRd
gem_metabs <- function(organism) {

    dataset.github <- switch(
        as.character(organism),
        '9606' = 'Human',
        '10090' = 'Mouse',
        '10116' = 'Rat',
        '7955' = 'Zebrafish',
        '7227' = 'Fruitfly',
        '6239' = 'Worm'
    )

    .slow_doctest()

    'gem_github' %>%
        generic_downloader(
            reader = read_tsv,
            url_key_param = list(),
            url_param = list(dataset.github, 'metabolites.tsv'),
            reader_param = list(),
            resource = NULL,
            post = NULL,
            use_httr = FALSE
        ) %T>% load_success()
}


#' Download GEM matlab file from Wang et al., 2021
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Matlab object containing GEM information.
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'   PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'   Genome-scale metabolic network reconstruction of model animals as a
#'   platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'   27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom stringr str_to_title
#'
#' @noRd
gem_matlab <- function(organism) {

    organism %<>% common_name %>% str_to_title

    .slow_doctest()

    'gem_github' %>%
        generic_downloader(
            reader = R.matlab::readMat,
            url_key_param = list(),
            url_param = list(dataset.github, paste0(dataset.github, '-gem_mat')),
            reader_param = list(),
            resource = NULL,
            post = NULL,
            use_httr = FALSE
        ) %T>% load_success()

}


#' Keep entries without elements in GEM processing
#'
#' @param matlab.object GEM from a Matlab object
#' @param attribs.mat Atribute from the Matlab object to be parsed
#' @param name Vector of elements to be checked in the Metlab object
#'
#' @return Vector with NAs in those entries with no element
#'
#' @noRd
metab_info <- function(matlab.object, attribs.mat, name) {
    unlist(
        sapply(
            X = matlab.object[[which(attribs.mat == name)]],
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
