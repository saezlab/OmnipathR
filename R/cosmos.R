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
#'   (\code{"rxns"} by default).
#' @param gem_metabolites.map.col Column of reaction IDs in the GEM
#'   (\code{"mets"} by default).
#' @param gem_list.params List containing the name of the slots where the
#'   information to construct the PKN is located in the gem_ If a matlab object
#'   is provided, this list parameter should not be modified.
#' @param gem_degree.mets.threshold Degree cutoff used to filter out
#'   metabolites (400 by default). The objective is to remove cofactors and
#'   over-promiscuous metabolites.
#' @param stitch.threshold Confidence cutoff used for STITCH connections
#'   (700 by default).
#' @param verbose Whether showing messages during execution (\code{TRUE} by
#'   default).
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
    gem_reactions.map.col = "rxns",
    gem_metabolites.map.col = "mets",
    gem_list.params = list(
        stoich.name = "S",
        reaction.name = "grRules",
        lb.name = "lb",
        ub.name = "ub",
        rev.name = "rev",
        reaction.ID.name = "rxns",
        metabolites.ID.name = "mets",
        metabolites.names.name = "metNames",
        metabolites.fomulas.name = "metFormulas",
        metabolites.inchi.name = "inchis"
    ),
    gem_degree.mets.threshold = 400,
    stitch.threshold = 700,
    verbose = TRUE
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
        stitch.threshold = stitch.threshold,
        verbose = verbose
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
#'   (\code{"rxns"} by default).
#' @param gem_metabolites.map.col Column of reaction IDs in the GEM
#'   (\code{"mets"} by default).
#' @param gem_list.params List containing the name of the slots where the
#'   information to construct the PKN is located in the gem_ If a matlab object
#'   is provided, this list parameter should not be modified.
#' @param gem_degree.mets.threshold Degree cutoff used to filter out
#'   metabolites (400 by default). The objective is to remove cofactors and
#'   other metabolites with many connections.
#' @param stitch.threshold Confidence cutoff used for STITCH connections
#'   (700 by default).
#' @param verbose Whether showing messages during execution (\code{TRUE} by
#'   default).
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
        gem_reactions.map.col = "rxns",
        gem_metabolites.map.col = "mets",
        gem_list.params = list(
            stoich.name = "S",
            reaction.name = "grRules",
            lb.name = "lb",
            ub.name = "ub",
            rev.name = "rev",
            reaction.ID.name = "rxns",
            metabolites.ID.name = "mets",
            metabolites.names.name = "metNames",
            metabolites.fomulas.name = "metFormulas",
            metabolites.inchi.name = "inchis"
        ),
        gem_degree.mets.threshold = 400,
        stitch.threshold = 700,
        verbose = TRUE
) {
    ## check organisms
    dataset.biomart <- switch(
        as.character(organism),
        "9606" = "hsapiens_gene_ensembl",
        "10090" = "mmusculus_gene_ensembl",
        "10116" = "rnorvegicus_gene_ensembl"
    )
    if (is.null(dataset.biomart))
        stop(
            "Chosen organism is not recognizable Available options are: ",
            paste(c(9606, 10090, 10116, 7955, 7227, 6239), collapse = ", ")
        )

    if (verbose) message(">>> Loading GEM obtained from Wang et al., 2021...")
    ## download GEM using OmnipathR (if already done, it will be taken from cache)
    gem_matlab <- OmnipathR:::gem_matlab(organism = organism) %>% as.data.frame()
    gem_metabs <- OmnipathR:::gem_metabs(organism = organism) %>% as.data.frame()
    gem_reacts <- OmnipathR:::gem_reacts(organism = organism) %>% as.data.frame()

    if (verbose) message("\n>>> Loading protein-chemical interactions from STITCH...")
    ## download STITCH using OmnipathR (if already done, it will be taken from cache)
    stitch.actions <- OmnipathR:::stitch_actions(organism = organism) %>%
        as.data.frame()
    stitch.prot.details <- OmnipathR:::stitch_prot_details(
        organism = organism
    ) %>% as.data.frame()

    if (biomart.use.omnipath == TRUE) {
        if (verbose) message("\n>>> Using the OmnipathR to retrieve BioMart information")
        ## get info from BiomartR using OmnipathR
        mapping.biomart <- OmnipathR::biomart_query(
            attrs = c(
                "ensembl_peptide_id",'ensembl_gene_id', 'external_gene_name'
            ),
            dataset = dataset.biomart
        ) %>% as.data.frame()
    } else {
        if (verbose) message("\n>>> Using the BiomaRt R package when needed")

        mapping.biomart <- NULL
    }
    ## Omnipath data
    if (verbose) message("\n>>> Loading protein-protein interactions from Omnipath...")
    omnipath.PKN <- .retrievingOmnipath(organism)
    ## Getting GEM PKN
    if (verbose) message("\n>>> Processing gem_..")
    gem_PKN.list <- .create_gem_basal_PKN(
        matlab.object = gem_matlab,
        reactions.map = gem_reacts,
        reactions.map.col = gem_reactions.map.col,
        metabolites.map = gem_metabs,
        metabolites.map.col = gem_metabolites.map.col,
        list.params.GEM = gem_list.params,
        degree.mets.cutoff = gem_degree.mets.threshold,
        verbose = verbose
    )
    if (verbose) message("\n>>> Formatting GEM PKN for COSMOS...")
    gem_PKN.list <- .mets_to_HMDB(gem_PKN.list)

    if (as.character(organism) == "9606") translate.genes <- TRUE

    if (translate.genes){
        gem_PKN.list <- .genes_to_symbol(
            gem_PKN.list, organism = organism,
            mapping.biomart = mapping.biomart
        )
    }
    gem_PKN.list <- .format_gem_COSMOS(gem_PKN.list, verbose = TRUE)

    if (verbose) message("\n>>> Getting STITCH PKN...")
    stitch.PKN <- .formatSTITCH(
        stitch.actions = stitch.actions,
        stitch.links = stitch.prot.details,
        organism = organism,
        omnipath.PKN = omnipath.PKN,
        mapping.biomart = mapping.biomart,
        threshold = stitch.threshold,
        verbose = verbose
    )
    output.final <- .mixing_resources(
        gem_network = gem_PKN.list[[1]],
        omnipath.PKN = omnipath.PKN,
        stitch.PKN = stitch.PKN
    )
    # if (dim(output.final)[[1]][1] == 0) {
    #       stop(
    #           paste(
    #               "Output incorrectly generated. This might happen when the used",
    #               "gene ontology by GEM and OmnipathR/STITCH is not the same. Check",
    #               "translate.genes parameter"
    #           )
    #       )
    # }
    if (verbose) message("\nDONE")

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
.omnipath_for_cosmos <- function(
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
    clean_PKN_mm[clean_PKN_mm$sign == 0, "sign"] <- 1
    clean_PKN_mm <- as.data.frame(rbind(clean_PKN_mm, clean_PKN_supp_mm))
    names(clean_PKN_mm) <- c("source", "target", "sign")

    return(clean_PKN_mm)
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
.metab_info <- function(matlab.object, attribs.mat, name) {
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
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
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
.gem_basal_pkn <- function(
        matlab.object,
        reactions.map,
        metabolites.map,
        reactions.map.col = "rxns",
        metabolites.map.col = "mets",
        list.params.GEM = list(
            stoich.name = "S",
            reaction.name = "grRules",
            lb.name = "lb",
            ub.name = "ub",
            rev.name = "rev",
            reaction.ID.name = "rxns",
            metabolites.ID.name = "mets",
            metabolites.names.name = "metNames",
            metabolites.fomulas.name = "metFormulas",
            metabolites.inchi.name = "inchis"
        ),
        degree.mets.cutoff = 400,
        verbose = TRUE
) {
    ## check parameters
    if (!reactions.map.col %in% colnames(reactions.map)) {
        stop("reactions.map.col cannot be found in reactions.map data.frame")
    } else if (!metabolites.map.col %in% colnames(metabolites.map)) {
        stop("metabolites.map.col cannot be found in metabolites.map data.frame")
    } else if (degree.mets.cutoff < 1) {
        stop("degree.mets.cutoff cannot be less than 1")
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
                            idx, "element in list.params.GEM (",
                            list.params.GEM[[idx]], ") is not in matlab object"
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
        "forward", "backward"
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
                        " and ", "_",
                        gsub(
                            "[()]","",
                            gsub("_AT[0-9]+","", strsplit(genes.reac, split = " or ")[[1]])
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
    orphan.reacts <- grepl(pattern = "^\\d+$", reaction.to.genes.df$Gene)
    reaction.to.genes.df[orphan.reacts, "Reaction.ID"] <- paste0(
        "orphanReac.", reaction.to.genes.df[orphan.reacts, "Reaction.ID"]
    )
    reaction.to.genes.df[orphan.reacts, "Gene"] <- paste0(
        reaction.to.genes.df[orphan.reacts, "Gene"], ".",
        reaction.to.genes.df[orphan.reacts, "Reaction.ID"]
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
    metabolites.map[metabolites.map == ""] <- NA
    ##############################################################################
    ## SIF file: PKN
    if (verbose) message("\n>>> Generating PKN")

    reaction.to.genes.df.reac <- reaction.to.genes.df
    reaction.list <- list()
    for (reac.idx in seq(ncol(s.matrix))) {
        reaction <- s.matrix[, reac.idx]
        #modify gene name so reactions that are catalised by same enzyme stay separated
        reaction.to.genes.df.reac[reaction.to.genes.df$Reaction == reac.idx, 1] <- paste(
            paste0("Gene", reac.idx),
            reaction.to.genes.df[reaction.to.genes.df$Reaction == reac.idx, 1],
            sep = "__"
        )
        # get the enzymes associated with reaction
        genes <- reaction.to.genes.df.reac[reaction.to.genes.df.reac$Reaction == reac.idx, 1]
        if (as.vector(lbs[reac.idx, 4] == "forward")) {
            reactants <- metabolites.IDs[reaction == -1]
            products <- metabolites.IDs[reaction == 1]
        } else {
            reactants <- metabolites.IDs[reaction == 1]
            products <- metabolites.IDs[reaction == -1]
        }
        reactants <- paste0("Metab__", reactants)
        products <- paste0("Metab__", products)
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
                                paste(gene, "_reverse", sep = ""),
                                number_of_interations - length(products)
                            ),
                            products
                        ),
                        target = c(
                            reactants,
                            rep(
                                paste(gene, "_reverse", sep = ""),
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
    reaction.df.all <- reaction.df.all[reaction.df.all$source != "Metab__" &
                                                                             reaction.df.all$target != "Metab__",]
    ## only complete cases
    reaction.df.all <- reaction.df.all[complete.cases(reaction.df.all),]
    ##############################################################################
    ## removing cofactors (metabolites with a high degree)
    metabs.degree <- sort(
        table(
            grep(
                "^Metab__", c(reaction.df.all$source, reaction.df.all$target),
                value = TRUE
            )
        ),
        decreasing = TRUE
    )
    if (verbose)
        message(
            "\t>>> Number of metabolites removed after degree >",
            degree.mets.cutoff,  ": ", sum(metabs.degree >= degree.mets.cutoff)
        )
    metabs.degree.f <- metabs.degree[metabs.degree < degree.mets.cutoff]
    reactions.df.no.cofac <- reaction.df.all[
        reaction.df.all$source %in% names(metabs.degree.f) |
            reaction.df.all$target %in% names(metabs.degree.f),
    ]
    mets <- grep(
        pattern = "Metab__",
        x = unique(c(reactions.df.no.cofac[[1]], reactions.df.no.cofac[[1]])),
        value = TRUE
    ) %>% gsub("Metab__", "", .)
    metabolites.map <- metabolites.map[mets, ]
    if (verbose) {
        message(
            "\t>>> Final number of connections: ", nrow(reactions.df.no.cofac)
        )
    }

    return(
        list(
            gem_PKN = reactions.df.no.cofac,
            mets.map = metabolites.map,
            reac.to.gene = reaction.to.genes.df.reac,
            reac.map = reactions.map
        )
    )
}


#' Translating metabolites from the Metabolic Atlas ID to the HMDB and KEGG IDs
#'
#' It transforms the Metabolic Atlas IDs into the HMDBs IDs. If a metabolite has
#' no HMDB ID, then the KEGG IDs are used. In case there is not entry for HMDB
#' or KEGG, the metabolite keeps the Metabolic Atlas ID.
#'
#' @param list.network List obtained using \code{.create_gem_basal_PKN}.
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
.mets_to_HMDB <- function(list.network) {
    metab.map <- list.network[[2]]
    list.network[[1]] <- list.network[[1]] %>% mutate(
        source = case_when(
            grepl("Metab__", source) ~ case_when(
                !is.na(
                    metab.map[gsub("Metab__", replacement = "", x = source), "metHMDBID"]
                ) ~ paste0(
                    "Metab__",
                    metab.map[gsub("Metab__", replacement = "", x = source), "metHMDBID"],
                    "_", str_sub(source, start = nchar(source))
                ),
                !is.na(
                    metab.map[gsub("Metab__", replacement = "", x = source), "metKEGGID"]
                ) ~ paste0(
                    "Metab__",
                    metab.map[gsub("Metab__", replacement = "", x = source), "metKEGGID"],
                    "_", str_sub(source, start = nchar(source))
                ), TRUE ~ source
            ), TRUE ~ source
        ),
        target = case_when(
            grepl("Metab__", target) ~ case_when(
                !is.na(
                    metab.map[gsub("Metab__", replacement = "", x = target), "metHMDBID"]
                ) ~ paste0(
                    "Metab__",
                    metab.map[gsub("Metab__", replacement = "", x = target), "metHMDBID"],
                    "_", str_sub(target, start = nchar(target))
                ),
                !is.na(
                    metab.map[gsub("Metab__", replacement = "", x = target), "metKEGGID"]
                ) ~ paste0(
                    "Metab__",
                    metab.map[gsub("Metab__", replacement = "", x = target), "metKEGGID"],
                    "_", str_sub(target, start = nchar(target))
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
#' @param list.network List obtained using \code{.create_gem_basal_PKN}.
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#' @param mapping.biomart BioMart ontology mapping data frame. If \code{NULL},
#'   this info is obtained using the \ckg{bioMaRt} R package.
#' @param ont.from Ontology to translate genes from
#'   (\code{"ensembl_gene_id" by default}).
#' @param ont.to Ontology to translate genes into
#'   (\code{"external_gene_name" by default}).
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
.genes_to_symbol <- function(
        list.network,
        organism,
        mapping.biomart = NULL,
        ont.from = "ensembl_gene_id",
        ont.to = "external_gene_name"
) {
    dataset.biomart <- switch(
        as.character(organism),
        "9606" = "hsapiens_gene_ensembl",
        "10090" = "mmusculus_gene_ensembl",
        "10116" = "rnorvegicus_gene_ensembl",
        "7955" = "drerio_gene_ensembl",
        "7227" = "dmelanogaster_gene_ensembl",
        "6239" = "celegants_gene_ensembl"
    )
    if (is.null(dataset.biomart))
        stop("Chosen organism is not recognizable")

    regex <- "(Gene\\d+__)|(_reverse)"
    ## getting biomart info
    genes.GEM <- grep("Gene\\d+__", unlist(list.network[[1]]), value = TRUE) %>%
        gsub(regex, "", .) %>%
        ifelse(grepl("_", .), sapply(strsplit(., split = "_"), \(x) x), .) %>%
        unlist()
    if (is.null(mapping.biomart)) {
        ensembl.link <- useEnsembl(biomart = "ensembl", dataset = dataset.biomart)
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
            grepl("Gene\\d+__", source) ~ case_when(
                !is.na(
                    ensembl.df[gsub(regex, replacement = "", x = source), ont.to]
                ) ~ paste0(
                    "Gene",
                    gsub("\\D", "", sapply(strsplit(x = source, split = "__"), \(x) x[1])),
                    "__",
                    ensembl.df[gsub(regex, replacement = "", x = source), ont.to]
                ),
                ## cases with complexes: more than 1 gene
                grepl("[0-9]_[E]", source) ~
                    paste0(
                        "Gene",
                        gsub("\\D", "", sapply(strsplit(x = target, split = "__"), \(x) x[1])),
                        "__",
                        unlist(
                            strsplit(
                                gsub(
                                    pattern = "reverse", replacement = "",
                                    grep("[0-9]_[E]", source, value = T)[1]
                                ),
                                split = "_"
                            )
                        )[-c(1:2)] %>% ensembl.df[., ont.to] %>% paste(collapse = "_")
                    ),
                TRUE ~ source
            ), TRUE ~ source
        ),
        target = case_when(
            ## cases with a single gene
            grepl("Gene\\d+__", target) ~ case_when(
                !is.na(
                    ensembl.df[gsub(regex, replacement = "", x = target), ont.to]
                ) ~ paste0(
                    "Gene",
                    gsub("\\D", "", sapply(strsplit(x = target, split = "__"), \(x) x[1])),
                    "__",
                    ensembl.df[gsub(regex, replacement = "", x = target), ont.to]
                ),
                ## cases with complexes: more than 1 gene
                grepl("[0-9]_[E]", target) ~
                    paste0(
                        "Gene",
                        gsub("\\D", "", sapply(strsplit(x = target, split = "__"), \(x) x[1])),
                        "__",
                        unlist(
                            strsplit(
                                gsub(
                                    pattern = "reverse", replacement = "",
                                    grep("[0-9]_[E]", target, value = T)[1]
                                ),
                                split = "_"
                            )
                        )[-c(1:2)] %>% ensembl.df[., ont.to] %>% paste(collapse = "_")
                    ),
                TRUE ~ target
            ), TRUE ~ target
        )
    )
    list.network[[3]] <- list.network[[3]] %>% mutate(
        Gene = case_when(
            !is.na(
                ensembl.df[gsub(regex, replacement = "", x = Gene), ont.to]
            ) ~ paste0(
                "Gene",
                gsub("\\D", "", sapply(strsplit(x = Gene, split = "__"), \(x) x[1])),
                "__",
                ensembl.df[gsub(regex, replacement = "", x = Gene), ont.to]
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


#' Formatting PKN derived from GEM for COSMOS
#'
#' It determines and marks transporters and reverse reactions.
#'
#' @param list.network List obtained using \code{.create_gem_basal_PKN}.
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'   translated into the desired ontology, gene-to-reactions data frame,
#'   metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @noRd
.format_gem_cosmos <- function(list.network, verbose = TRUE) {
    reaction.network <- list.network[[1]]
    enzyme_reacs <- unique(c(reaction.network$source, reaction.network$target))
    enzyme_reacs <- enzyme_reacs[grepl("^Gene", enzyme_reacs)]
    enzyme_reacs_reverse <- enzyme_reacs[grepl("_reverse",enzyme_reacs)]
    enzyme_reacs <- enzyme_reacs[!grepl("_reverse",enzyme_reacs)]

    if (verbose) message("\t>>> Step 1: Defining transporters")

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
                        if(grepl("Metab__", df[i, 1])) {
                            counterpart <- which(
                                gsub("_[a-z]$","",df[,2]) == gsub("_[a-z]$","",df[i,1])
                            )
                            if(length(counterpart) > 0) {
                                df[i, 2] <- paste0(df[i, 2], paste0("_TRANSPORTER", i))
                                df[counterpart, 1] <- paste0(
                                    df[counterpart, 1], paste0("_TRANSPORTER", i)
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

    if (verbose) message("\t>>> Step 2: Defining reverse reactions")

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
                        if(grepl("Metab__",df[i,1])) {
                            counterpart <- which(
                                gsub("_[a-z]$","",df[,2]) == gsub("_[a-z]$","",df[i,1])
                            )
                            if(length(counterpart) > 0) {
                                transporter <- gsub("_reverse", "", df[i, 2])
                                transporter <- paste0(
                                    transporter, paste0(paste0("_TRANSPORTER", i), "_reverse")
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
        grep("Metab__", reaction.network.new[[1]], value = TRUE),
        grep("Metab__", reaction.network.new[[2]], value = TRUE)
    ) %>% unique() %>% gsub("(Metab__)|(_[a-z])", "", .)
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


#' Connecting PKN derived from GEM with Omnipath protein-protein interactions
#'
#' @param gem_PKN List obtained using \code{.create_gem_basal_PKN}.
#' @param omnipath.PKN Protein-protein interactions obtained using
#'   \code{.retrievingOmnipath}.
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
#'
#' @return List containing PKN with COSMOS and OCEAN format with genes
#'   translated into the desired ontology, gene-to-reactions data frame,
#'   metabolite-mapping data frame, and reactions-mapping data frame.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
.connect_gem_omnipath <- function(
        gem_PKN,
        omnipath.PKN,
        verbose = TRUE
) {
    elements <- unique(as.character(unlist(gem_PKN)))
    elements <- elements[grepl("^Gene\\d+__", elements)]
    elements <- gsub("(.*__)|(_TRANSPORTER[0-9]+)|(_reverse$)", "", elements)
    ## this function can be vectorized
    connectors.df <- sapply(
        X = elements,
        FUN = function(ele) {
            if (grepl("_", ele)) {
                genes.sep <- str_split(string = ele, pattern = "_")[[1]]
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
    names(connectors.df) <- c("source", "target")
    connectors.df <- connectors.df[which(
        connectors.df$source %in% omnipath.PKN$source |
            connectors.df$source %in% omnipath.PKN$target
    ),]
    network.df.new <- as.data.frame(rbind(gem_PKN, connectors.df))

    return(network.df.new)
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
#'   \code{.retrievingOmnipath}.
#' @param mapping.biomart BioMart ontology mapping data frame. If \code{NULL},
#'   this info is obtained using the \ckg{bioMaRt} R package.
#' @param threshold Confidence cutoff used for STITCH connections
#'   (700 by default).
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
#'
#' @return List containing PKN with COSMOS and OCEAN format.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter case_when mutate select
#'
#' @noRd
.format_gem_stitch <- function(
        stitch.actions,
        stitch.links,
        organism,
        omnipath.PKN,
        mapping.biomart = NULL,
        threshold = 700,
        verbose = TRUE
) {
    dataset.biomart <- switch(
        as.character(organism),
        "9606" = "hsapiens_gene_ensembl",
        "10090" = "mmusculus_gene_ensembl",
        "10116" = "rnorvegicus_gene_ensembl",
        "7955" = "drerio_gene_ensembl",
        "7227" = "dmelanogaster_gene_ensembl",
        "6239" = "celegants_gene_ensembl"
    )
    if (is.null(dataset.biomart))
        stop("Chosen organism is not recognizable")
    if (verbose) {
        message("\t>>> Reading provided STITCH files\n")
    }
    links.detail <- as.data.frame(stitch.links) %>% filter(
        combined_score >= threshold,
        experimental >= threshold | database >= threshold
    ) %>% mutate(
        ID = paste(chemical, protein, sep = "_"),
        ID_reverse = paste(protein, chemical, sep = "_")
    )
    STITCH <- as.data.frame(stitch.actions) %>%
        filter(mode == "activation" | mode == "inhibition", a_is_acting) %>%
        mutate(ID = paste(item_id_a, item_id_b, sep = "_")) %>%
        filter(ID %in% links.detail$ID | ID %in% links.detail$ID_reverse)
    STITCH <- STITCH[,-7]
    ## df of proteins in STICH
    prots <- unique(c(STITCH$item_id_a, STITCH$item_id_b))
    prots <- prots[grepl(paste0(organism, "\\."), prots)]
    prots <- as.data.frame(cbind(prots, gsub(paste0(organism, "\\."), "", prots)))
    colnames(prots) <- c("original", "ensembl_prots")
    ## getting info from Biomart
    if (verbose) message("\t>>> Using information from BiomaRt")

    if (is.null(mapping.biomart)) {
        ensembl.link <- useEnsembl(biomart = "ensembl", dataset = dataset.biomart)
        ensembl.df <- getBM(
            filters = "ensembl_peptide_id",
            attributes = c(
                "ensembl_peptide_id",'ensembl_gene_id', 'external_gene_name'# , 'entrezgene_id', "description"
            ),
            values = prots[[2]],
            mart = ensembl.link
        )
        colnames(ensembl.df)[1] <- "ensembl_prots"
    } else {
        ensembl.df <- mapping.biomart %>% filter(
            ensembl_peptide_id %in% prots[[2]]
        )
        colnames(ensembl.df)[1] <- "ensembl_prots"
    }

    prots <- merge(prots, ensembl.df, by = "ensembl_prots")
    ## external_gene_name for mouse, Idk in other cases, check this
    prots <- prots[prots$external_gene_name != "",]
    prots.vec <- prots$external_gene_name
    names(prots.vec) <- prots$original

    if (verbose) message("\n\t>>> Generating PKN network")

    STITCH <- STITCH %>% mutate(
        item_id_a = case_when(
            grepl("\\.", item_id_a) & (item_id_a %in% names(prots.vec)) ~
                prots.vec[item_id_a],
            grepl("^CID", item_id_a) ~ gsub("CID[a-z]0*", "Metab__", item_id_a),
            TRUE ~ item_id_a
        ),
        item_id_b = case_when(
            grepl("\\.", item_id_b) & (item_id_b %in% names(prots.vec)) ~
                prots.vec[item_id_b],
            grepl("^CID", item_id_b) ~ gsub("CID[a-z]0*", "Metab__", item_id_b),
            TRUE ~ item_id_b
        ),
        sign = case_when(action == "inhibition" ~ -1, TRUE ~ 1)
    ) %>% select(1, 2, 7)
    colnames(STITCH) <- c("source", "target", "sign")
    CIDs <- unique(as.character(unlist(STITCH[,c(1,3)])))
    CIDs <- CIDs[grepl("Metab__", CIDs)] %>% gsub("Metab__", "", .)
    ## Convert CID to HMDB Id when available
    metabolitesMapping.mod <-
        metaboliteIDMapping::metabolitesMapping %>%
        filter(CID %in% CIDs, !is.na(HMDB)) %>%
        mutate(HMDB = paste0("Metab__", HMDB))
    metabolitesMapping.vec <- metabolitesMapping.mod$HMDB
    names(metabolitesMapping.vec) <- paste0("Metab__", metabolitesMapping.mod$CID)
    ## metabolites with no HMDB are kept
    STITCH <- STITCH %>% mutate(
        source = case_when(
            grepl("Metab__", source) & source %in% names(metabolitesMapping.vec) ~
                metabolitesMapping.vec[source],
            TRUE ~ source
        ),
        target = case_when(
            grepl("Metab__", target) & target %in% names(metabolitesMapping.vec) ~
                metabolitesMapping.vec[target],
            TRUE ~ target
        )
    )
    # TODO: at this point, STITCH contains metabolites in both columns of the dataframe
    ## this should be checked
    omn.prots <- unique(as.character(unlist(omnipath.PKN[,c(1, 2)])))
    STITCH <- unique(STITCH[which(STITCH$target %in% omn.prots),])

    STITCH$source <- paste(STITCH$source, "_c", sep = "")

    return(STITCH)
}


#' Combine GEM-, Omnipath- and STITCH- derived PKNs
#'
#' @param gem_network Metabolic PKN obtained using \code{.create_gem_basal_PKN}.
#' @param omnipath.PKN Protein-protein interactions obtained using
#'   \code{.retrievingOmnipath}.
#' @param stitch.PKN Chemical-protein PKN obtained using \code{.formatSTITCH}.
#'
#' @return List containing PKN with COSMOS and OCEAN format.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
.combine_resources <- function(gem_network, omnipath.PKN, stitch.PKN) {
    ## connecting Omnipath and GEM
    gem_network <- .connecting_gem_omnipath(
        gem_PKN = gem_network, omnipath.PKN = omnipath.PKN
    )
    gem_network$sign <- 1
    meta.PKN <- as.data.frame(
        rbind(omnipath.PKN, stitch.PKN, gem_network)
    ) %>% unique()
    meta.PKN <- meta.PKN[, c(1, 3, 2)]
    names(meta.PKN) <- c("source", "interaction", "target")
    #TODO: manual correction: difficult to generalize for different organisms, shall I remove it?
    #probably erroneous interaction (WHY?? this only works for human / mouse)
    # meta.network <- meta.network[-which(
    #       meta.network$source == "Prkca" & meta.network$target == "Src"
    # ),]
    # meta.PKN <- meta.PKN[-which(
    #       meta.PKN$source == "Prkca" & meta.PKN$target == "Src"
    # ),]
    #probably erroneous interaction
    # meta_network <- meta_network[-which(meta.network$source == "Ltc45"),]
    #I don't know where this interaction comes from, the sources are wrong (https://www.nature.com/articles/onc2008228)
    # meta_network <- meta.network[!(grepl("Cad_reverse", meta.network$source) | grepl("Cad_reverse", meta.network$target)) ,]
    #redHuman confirms that the reaction is actually not reversible: NOT FOUND IN MOUSE EITHER

    return(meta.PKN)
}


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
stitch_prot_details <- function(organism) {

    .slow_doctest()

    'stitch_prot_details' %>%
        generic_downloader(
            reader = read_tsv,
            url_key_param = list(),
            url_param = list(organism),
            reader_param = list(trim_ws = TRUE),
            resource = NULL,
            post = NULL,
            use_httr = FALSE
        ) %T>% load_success()
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
        ) %T>% load_success()
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
        "9606" = "Human",
        "10090" = "Mouse",
        "10116" = "Rat",
        "7955" = "Zebrafish",
        "7227" = "Fruitfly",
        "6239" = "Worm"
    )

    .slow_doctest()

    'gem_github' %>%
        generic_downloader(
            reader = read_tsv,
            url_key_param = list(),
            url_param = list(dataset.github, "reactions.tsv"),
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
        "9606" = "Human",
        "10090" = "Mouse",
        "10116" = "Rat",
        "7955" = "Zebrafish",
        "7227" = "Fruitfly",
        "6239" = "Worm"
    )

    .slow_doctest()

    'gem_github' %>%
        generic_downloader(
            reader = read_tsv,
            url_key_param = list(),
            url_param = list(dataset.github, "metabolites.tsv"),
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
#'
#' @noRd
gem_matlab <- function(organism) {

    dataset.github <- switch(
        as.character(organism),
        "9606" = "Human",
        "10090" = "Mouse",
        "10116" = "Rat",
        "7955" = "Zebrafish",
        "7227" = "Fruitfly",
        "6239" = "Worm"
    )

    .slow_doctest()

    'gem_github' %>%
        generic_downloader(
            reader = R.matlab::readMat,
            url_key_param = list(),
            url_param = list(dataset.github, paste0(dataset.github, "-gem_mat")),
            reader_param = list(),
            resource = NULL,
            post = NULL,
            use_httr = FALSE
        ) %T>% load_success()

}
