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
#' @param protein_ids Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the "a"
#'     and "b" sides of the interaction, respectively. The default ID type for
#'     proteins is Esembl Gene ID, and by default UniProt IDs and Gene
#'     Symbols are included.
#' @param metabolite_ids Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the "a"
#'     and "b" sides of the interaction, respectively. The default ID type for
#'     metabolites is Metabolic Atlas ID, and HMDB IDs and KEGG IDs are included.
#' @param chalmers_gem_metab_max_degree Degree cutoff used to filter out
#'     metabolites (400 by default). The objective is to remove cofactors and
#'     over-promiscuous metabolites.
#' @param stitch_score Confidence cutoff used for STITCH connections (700 by
#'     default).
#' @param ... Further parameters to \code{\link{import_omnipath_interactions}}.
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
#'     human_cosmos <- cosmos_pkn(organism = 9606)
#' }
#'
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom rlang exec !!!
#'
#' @export
cosmos_pkn <- function(
    organism = 'human',
    protein_ids = c('uniprot', 'genesymbol'),
    metabolite_ids = c('hmdb', 'kegg'),
    chalmers_gem_metab_max_degree = 400L,
    stitch_score = 700L,
    ...
){

    .slow_doctest()

    organism %<>% ncbi_taxid()

    args <- environment() %>% as.list %>% c(list(...))

    ## check dependencies (Suggests in DESCRIPTION)
    'R.matlab' %>%
    missing_packages %>%
    paste(collapse = ', ') %>%
    {`if`(nchar(.), sprintf('Missing packages: %s', .) %T>% log_error %>% stop)}

    with_cache(
        name = 'COSMOS_PKN_%s' %>% sprintf(organism),
        args = args,
        callback = .cosmos_pkn
    )

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
#' @param ... Further parameters to \code{\link{import_omnipath_interactions}}.
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
    organism = 'human',
    protein_ids = c('uniprot', 'genesymbol'),
    metabolite_ids = c('hmdb', 'kegg'),
    chalmers_gem_metab_max_degree = 400L,
    stitch_score = 700L,
    ...
) {

    organism %<>% organism_for('cosmos')

    log_success(
        paste0(
            'Building COSMOS PKN (organism: %s). ',
            'This will take 10-30 min at the first ',
            'time, and will be saved in the cache for later use.'
        ),
        organism
    )

    stitch <- stitch_gem(
        organism = organism,
        min_score = stitch_score,
        protein_ids = protein_ids,
        metabolite_ids = metabolite_ids,
        cosmos = TRUE
    )

    chalmers <- chalmers_gem_network(
        organism = organism,
        protein_ids = protein_ids,
        metabolite_ids = metabolite_ids,
        metab_max_degree = chalmers_gem_metab_max_degree
    )

    omnipath <- omnipath_for_cosmos(organism, id_types = protein_ids, ...)

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
#' @param id_types Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the
#'     "source" and "target" sides of the interaction, respectively. The
#'     default ID type for proteins is Esembl Gene ID, and by default UniProt
#'     IDs and Gene Symbols are included. The UniProt IDs returned by the web
#'     service are left intact, while the Gene Symbols are queried from
#'     Ensembl. These Gene Symbols are different from the ones returned from
#'     the web service, and match the Ensembl Gene Symbols used by other
#'     components of the COSMOS PKN.
#' @param ... Further parameters to \code{\link{import_omnipath_interactions}}.
#'
#' @return Data frame with the columns source, target and sign.
#'
#' @examples
#' op_cosmos <- omnipath_for_cosmos()
#' op_cosmos
#'
#' @importFrom magrittr %<>% %>% %T>%
#' @importFrom logger log_info
#' @importFrom dplyr mutate filter select bind_rows select row_number
#' @importFrom rlang !!! syms
#' @export
#' @seealso \code{\link{cosmos_pkn}}
omnipath_for_cosmos <- function(
        organism = 9606L,
        resources = NULL,
        datasets = NULL,
        interaction_types = NULL,
        id_types = c('uniprot', 'genesymbol'),
        ...
    ) {

    # NSE vs. R CMD check workaround
    consensus_stimulation <- consensus_inhibition <- NULL

    organism %<>% organism_for('omnipath')
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
    mutate(
        sign = consensus_stimulation - consensus_inhibition,
        record_id = row_number()
    ) %>%
    select(source, target, sign, record_id) %>%
    {`if`(
        'uniprot' %in% id_types,
        mutate(., uniprot_source = source, uniprot_target = target),
        .
    )} %>%
    translate_ids_multi(
        source = uniprot,
        target,
        !!!syms(setdiff(id_types, 'uniprot')),
        ensembl = TRUE,
        organism = organism
    ) %>%
    bind_rows(
        mutate(., sign = ifelse(sign, sign, 1)),
        filter(., sign == 0) %>% mutate(sign = -1)
    ) %T>%
    {log_info('OmniPath PPI for COSMOS PKN ready: %i interactions.', nrow(.))}

}


#' Combine components of COSMOS PKN
#'
#' @param chalmers Data frame: the Chalmers Sysbiol GEM PKN as produced by
#'     \{code{\link{chalmers_gem_network}}.
#' @param omnipath Data frame: the OmniPath PKN as produced by
#'     \{code{\link{omnipath_for_cosmos}}.
#' @param stitch Data frame: the STITCH PKN as produced by
#'     \{code{\link{stitch_gem}}.
#'
#' @return Data frame: the combined PKN suitable for COSMOS and OCEAN.
#'
#' @importFrom magrittr %<>% %>% %T>%
#' @importFrom dplyr rename rename_with bind_rows mutate relocate select across
#' @importFrom tidyselect starts_with
#' @importFrom stringr str_replace str_detect
#' @importFrom logger log_trace log_info
#' @noRd
cosmos_combine_networks <- function(
    chalmers,
    omnipath,
    stitch
){

    .slow_doctest()

    # NSE vs. R CMD check workaround
    item_id_a <- item_id_b <- score <- ri <- ci <- .x <- comp <- reverse <-
    transporter <- met_to_gene <- target <- resource <- entity_type_source <-
    entity_type_target <- score_stitch <- ci_chalmers <- comp_chalmers <-
    reverse_chalmers <- transporter_chalmers <- record_id <- NULL

    stitch_etype <- function(x) {
        ifelse(str_detect(x, '^ENS'), 'protein', 'metabolite')
    }

    log_info('Combining components of COSMOS PKN.')
    log_trace('Preparing STITCH')

    stitch %<>%
        rename(
            source = item_id_a,
            target = item_id_b,
            score_stitch = score
        ) %>%
        rename_with(~str_replace(.x, '_a$', '_source')) %>%
        rename_with(~str_replace(.x, '_b$', '_target')) %>%
        mutate(
            entity_type_source = stitch_etype(source),
            entity_type_target = stitch_etype(target),
            resource = 'STITCH'
        )

    log_trace('Preparing Chalmers Sysbio GEM')

    chalmers %<>%
        rename(
            record_id = ri,
            ci_chalmers = ci,
            comp_chalmers = comp,
            reverse_chalmers = reverse,
            transporter_chalmers = transporter
        ) %>%
        mutate(
            entity_type_source = ifelse(met_to_gene, 'metabolite', 'protein'),
            entity_type_target = ifelse(met_to_gene, 'protein', 'metabolite'),
            resource = 'Chalmers Sysbiol GEM'
        ) %>%
        select(-met_to_gene)

    log_trace('Preparing OmniPath')

    omnipath %<>%
        mutate(
            entity_type_source = 'protein',
            entity_type_target = 'protein',
            resource = 'OmniPath'
        )

    log_trace('Combining data frames')

    bind_rows(chalmers, stitch, omnipath) %>%
    relocate(
        sign,
        record_id,
        resource,
        entity_type_source,
        entity_type_target,
        score_stitch,
        ci_chalmers,
        comp_chalmers,
        reverse_chalmers,
        transporter_chalmers,
        .after = target
    ) %>%
    mutate(across(starts_with('entity_type_'), as.factor)) %T>%
    {log_info('COSMOS combined PKN ready')}

}
