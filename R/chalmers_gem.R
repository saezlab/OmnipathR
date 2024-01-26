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
#' Process the GEMs from Wang et al., 2021
#' (\url{https://github.com/SysBioChalmers}) into convenient tables.
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicus), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return List containing the following elements: \itemize{
#'     \item{reactions: tibble of reaction data;}
#'     \item{metabolites: tibble of metabolite data;}
#'     \item{reaction_ids: translation table of reaction identifiers;}
#'     \item{metabolite_ids: translation table of metabolite identifiers;}
#'     \item{S: Stoichiometric matrix (sparse).}
#' }
#'
#' @examples
#' gem <- chalmers_gem()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'     PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'     Genome-scale metabolic network reconstruction of model animals as a
#'     platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'     27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>% extract extract2 equals
#' @importFrom dplyr pull across na_if first mutate slice n
#' @importFrom tibble as_tibble
#' @importFrom purrr map map_chr map_lgl map2 pmap_int
#' @importFrom logger log_info
#' @importFrom tidyselect everything
#' @importFrom stringr str_split str_replace_all str_detect
#'
#' @export
chalmers_gem <- function(organism = 'Human') {

    log_info('Processing Chalmers SysBio GEM.')

    raw <- chalmers_gem_raw(organism = organism)
    metabolite_ids <- chalmers_gem_metabolites(organism = organism)
    reaction_ids <- chalmers_gem_reactions(organism = organism)

    S <- raw %>% extract2(1L) %>% extract(,,1L) %>% extract2('S')

    metabolites <-
        raw %>%
        chalmers_gem_matlab_tibble(mets, metNames, metFormulas, inchis) %>%
        mutate(
            across(
                everything(),
                ~na_if(map_chr(unlist(.x, recursive = FALSE), first), '')
            )
        ) %T>%
        {log_info('Chalmers GEM: %i metabolites.', nrow(.))}

    reactions <-
        raw %>%
        chalmers_gem_matlab_tibble(rxns, lb, ub, reversible, grRules) %>%
        mutate(
            rxns = unlist(rxns),
            direction = ifelse(ub + lb >= 0L, 1L, -1L),
            reversible = as.logical(reversible),
            grRules = (
                unlist(grRules, recursive = FALSE) %>%
                map_chr(extract, 1L) %>%
                str_replace_all('\\[\\(\\)\\]|_AT\\d+', '') %>%
                str_split(' or ') %>%
                map(str_split, ' and ') %>%
                map(discard, is_empty_2)
            ),
            orphan = map_lgl(
                grRules,
                ~any(map_lgl(.x, ~any(str_detect(.x, '^\\d+$'))))
            ),
            reactants = map2(
                1L:n(),
                direction,
                ~.x %>% extract(S,,.) %>% equals(-.y) %>% which %>%
                     slice(metabolites, .) %>% pull(mets)
            ),
            products = map2(
                1L:n(),
                direction,
                ~.x %>% extract(S,,.) %>% equals(.y) %>% which %>%
                    slice(metabolites, .) %>% pull(mets)
            ),
            numof_interactions = pmap_int(
                list(reactants, products),
                ~length(c(..1, ..2))
            )
        ) %T>%
        {log_info(
            'Chalmers GEM: %i reactions, %i orphan.',
            nrow(.),
            pull(., orphan) %>% sum(na.rm = TRUE)
        )}

    return(
        list(
            reactions = reactions,
            metabolites = metabolites,
            reaction_ids = reaction_ids,
            metabolite_ids = metabolite_ids,
            S = S
        )
    )

}


#' Chalmers SysBio GEM in the form of gene-metabolite interactions
#'
#' Processing GEMs from Wang et al., 2021
#' (\url{https://github.com/SysBioChalmers}) to generate PKN for COSMOS
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicus), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#' @param metab_max_degree Degree cutoff used to prune metabolites with high
#'     degree assuming they are cofactors (400 by default).
#'
#' @return Data frame (tibble) of gene-metabolite interactions.
#'
#' @examples
#' gem <- chalmers_gem_network()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'     PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'     Genome-scale metabolic network reconstruction of model animals as a
#'     platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'     27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>% extract2
#' @importFrom purrr map map_chr
#' @importFrom dplyr mutate across select bind_rows filter
#' @importFrom tidyr unnest_longer
#' @importFrom logger log_info log_error
#' @export
chalmers_gem_network <- function(
        organism = 'Human',
        metab_max_degree = 400L
    ) {

    .slow_doctest()

    if (metab_max_degree < 1L) {
        '`metab_max_degree` cannot be less than 1.' %T>% log_error %>% stop
    }

    organism %>%
    chalmers_gem %T>%
    {log_info('Chalmers GEM: compiling gene-metabolite interactions.')} %>%
    extract2('reactions') %>%
    # once I see better what's the purpose of this labeling,
    # I might move it a bit further down the pipeline - denes
    mutate(
        genes = map(
            grRules,
            ~map(.x, ~map_chr(~sprintf('Gene%i__%s', ri, .x)))
        ) #,
        # # this labeling is not necessary at all
        # # at least I hope so -denes
        # across(
        #     c(reactants, products),
        #     ~map(.x, ~paste0('Metab__', .x))
        # )
    ) %>%
    select(ri, genes, reactants, products, rev) %>%
    unnest_longer(grRules, simplify = FALSE) %>%
    unnest_longer(reactants) %>%
    unnest_longer(products) %>%
    mutate(
        # note: this represents the AND relationship between genes,
        # i.e. both genes together are required for the reaction
        grRules = map_chr(~paste(.x, collapse = '_')),
        reverse = FALSE
    ) %>%
    {bind_rows(
        select(., ri, reactants, grRules, reverse),
        select(., ri, grRules, products, reverse),
        filter(., rev) %>% select(ri, grRules, reactants) %>% mutate(reverse = TRUE),
        filter(., rev) %>% select(ri, products, grRules) %>% mutate(reverse = TRUE)
    )} %T>%
    {log_info(
        'Chalmers GEM: %i records in gene-metabolite interactions table.',
        nrow(.)
    )} %>%
    list(
        lo_degree =
            c(
                pull(., mets = reactants),
                pull(., mets = products)
            ) %>%
            table %>%
            keep(~.x > metab_max_degree) %>%
            names
    ) %>%
    {filter(
        extract2(., 1L),
        reactants %in% .$lo_degree &
        products %in% .$lo_degree
    )} %T>%
    {log_info(
        paste0(
            'Chalmers GEM: %i gene-metabolite interactions after removing ',
            'metabolites with more than %i interactions.'
        ),
        nrow(.),
        metab_max_degree
    )} %T>%
    # here was removal of records with  missing metabolite or gene
    # I'm not aware of any reason such records should exist
    # so I removed the filter for now - denes
    {log_info('Chalmers GEM: gene-metabolite network is ready.')}

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


#' Tibble from Chalmers GEM Matlab object
#'
#' @param matlab Chalmers GEM in an R object loaded from the Matlab
#'     dump.
#' @param ... Variable names: should contain either only reaction or
#'     metabolite variables, otherwise num of rows won't be uniform.
#'
#' @return Tibble with the requested variables.
#'
#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom dplyr pull rename
#' @importFrom magrittr %>% extract extract2 equals
#' @importFrom rlang enquos
#'
#' @noRd
chalmers_gem_matlab_tibble <- function(matlab, ...) {

    cols <-
        enquos(...) %>%
        map_chr(.nse_ensure_str)

    matlab %>%
    as_tibble %>%
    pull(1L) %>%
    extract(,,1L) %>%
    map(
        function(x) {
            if (x %>% dim %>% equals(1L) %>% all) x %>% extract(1L,1L) else x
        }
    ) %>%
    extract(cols) %>%
    map(extract,,1L) %>%
    as_tibble %>%
    rename(reversible = rev)

}
