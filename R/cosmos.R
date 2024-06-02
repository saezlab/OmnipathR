#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Diego Mananes
#                 Alberto Valdeolivas
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
#' The prior knowledge network (PKN) used by COSMOS is a network of
#' heterogenous causal interactions: it contains protein-protein,
#' reactant-enzyme and enzyme-product interactions. It is a combination of
#' multiple resources:
#' \itemize{
#'     \item Genome-scale metabolic model (GEM) from Chalmers Sysbio (Wang et
#'     al., 2021.)
#'     \item Network of chemical-protein interactions from STITCH
#'     (\url{https://stitch.embl.de/})
#'     \item Protein-protein interactions from Omnipath (Türei
#'     et al., 2021)
#' }
#' This function downloads, processes and combines the resources above. With
#' all downloads and processing the build might take 30-40 minutes. Data is
#' cached at various levels of processing, shortening processing times. With
#' all data downloaded and HMDB ID translation data preprocessed, the build
#' takes 3-4 minutes; the complete PKN is also saved in the cache, if this is
#' available, loading it takes only a few seconds.
#'
#' @param organism Character or integer: name or NCBI Taxonomy ID of an
#'      organism. Supported organisms vary by resource: the Chalmers GEM is
#'      available only for human, mouse, rat, fish, fly and worm. OmniPath can
#'      be translated by orthology, but for non-vertebrate or less researched
#'      taxa very few orthologues are available. STITCH is available for a
#'      large number of organisms, please refer to their web page:
#'      \url{https://stitch.embl.de/}.
#' @param protein_ids Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the
#'     "source" and "target" sides of the interaction, respectively. The
#'     default ID type for proteins depends on the resource, hence the "source"
#'     and "target" columns are heterogenous. By default UniProt IDs and Gene
#'     Symbols are included. The Gene Symbols used in the COSMOS PKN are
#'     provided by Ensembl, and do not completely agree with the ones provided
#'     by UniProt and used in OmniPath data by default.
#' @param metabolite_ids Character: translate the metabolite identifiers to
#'     these ID types. Each ID type results two extra columns in the output,
#'     for the "source" and "target" sides of the interaction, respectively.
#'     The default ID type for metabolites depends on the resource, hence the
#'     "source" and "target" columns are heterogenous. By default HMDB IDs and
#'     KEGG IDs are included.
#' @param chalmers_gem_metab_max_degree Numeric: remove metabolites from the
#'     Chalmers GEM network with defgrees larger than this. Useful to remove
#'     cofactors and over-promiscuous metabolites.
#' @param stitch_score Include interactions from STITCH with combined
#'     confidence score larger than this.
#' @param ... Further parameters to \code{\link{import_omnipath_interactions}}.
#'
#' @return A data frame of binary causal interations with effect signs,
#'     resource specific attributes and translated to the desired identifiers.
#'     The ``record_id`` column identifies the original records within each
#'     resource. If one ``record_id`` yields multiple records in the final data
#'     frame, it is the result of one-to-many ID translation or other
#'     processing steps. Before use, it is recommended to select one pair of ID
#'     type columns (by combining the preferred ones) and perform ``distinct``
#'     by the identifier columns and sign.
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M,
#'     Cholley PE, et al. Genome-scale metabolic network reconstruction of model
#'     animals as a platform for translational research. Proceedings of the
#'     National Academy of Sciences. 2021 Jul 27;118(30):e2102344118.
#'
#'     Türei D, Valdeolivas A, Gul L, Palacio‐Escat N, Klein M, Ivanova O, et al.
#'     Integrated intra‐ and intercellular signaling knowledge for multicellular
#'     omics analysis. Molecular Systems Biology. 2021 Mar;17(3):e9923.
#'
#' @examples
#' \dontrun{
#'     human_cosmos <- cosmos_pkn(organism = "human")
#' }
#'
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom rlang exec !!!
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_network}}}
#'     \item{\code{\link{stitch_network}}}
#'     \item{\code{\link{omnipath_for_cosmos}}}
#'     \item{\code{\link{import_omnipath_interactions}}}
#' }
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
#' The prior knowledge network (PKN) used by COSMOS is a network of
#' heterogenous causal interactions: it contains protein-protein,
#' reactant-enzyme and enzyme-product interactions. It is a combination of
#' multiple resources:
#' \itemize{
#'     \item Genome-scale metabolic model (GEM) from Chalmers Sysbio (Wang et
#'     al., 2021.)
#'     \item Network of chemical-protein interactions from STITCH
#'     (\url{http://stitch.embl.de/})
#'     \item Protein-protein interactions from Omnipath (Türei
#'     et al., 2021)
#' }
#' This function downloads, processes and combines the resources above. With
#' all downloads and processing the build might take 30-40 minutes. Data is
#' cached at various levels of processing, shortening processing times. With
#' all data downloaded and HMDB ID translation data preprocessed, the build
#' takes 3-4 minutes; the complete PKN is also saved in the cache, if this is
#' available, loading it takes only a few seconds.
#'
#' @param organism Character or integer: name or NCBI Taxonomy ID of an
#'      organism. Supported organisms are human, mouse and rat.
#' @param protein_ids Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the
#'     "source" and "target" sides of the interaction, respectively. The
#'     default ID type for proteins depends on the resource, hence the "source"
#'     and "target" columns are heterogenous. By default UniProt IDs and Gene
#'     Symbols are included. The Gene Symbols used in the COSMOS PKN are
#'     provided by Ensembl, and do not completely agree with the ones provided
#'     by UniProt and used in OmniPath data by default.
#' @param metabolite_ids Character: translate the metabolite identifiers to
#'     these ID types. Each ID type results two extra columns in the output,
#'     for the "source" and "target" sides of the interaction, respectively.
#'     The default ID type for metabolites depends on the resource, hence the
#'     "source" and "target" columns are heterogenous. By default HMDB IDs and
#'     KEGG IDs are included.
#' @param chalmers_gem_metab_max_degree Numeric: remove metabolites from the
#'     Chalmers GEM network with defgrees larger than this. Useful to remove
#'     cofactors and over-promiscuous metabolites.
#' @param stitch_score Confidence cutoff used for STITCH connections (700 by
#'     default).
#' @param ... Further parameters to \code{\link{import_omnipath_interactions}}.
#'
#' @return A data frame of binary causal interations with effect signs,
#'     resource specific attributes and translated to the desired identifiers.
#'     The ``record_id`` column identifies the original records within each
#'     resource. If one ``record_id`` yields multiple records in the final data
#'     frame, it is the result of one-to-many ID translation or other
#'     processing steps. Before use, it is recommended to select one pair of ID
#'     type columns (by combining the preferred ones) and perform ``distinct``
#'     by the identifier columns and sign.
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M,
#'     Cholley PE, et al. Genome-scale metabolic network reconstruction of model
#'     animals as a platform for translational research. Proceedings of the
#'     National Academy of Sciences. 2021 Jul 27;118(30):e2102344118.
#'
#'     Türei D, Valdeolivas A, Gul L, Palacio‐Escat N, Klein M, Ivanova O, et al.
#'     Integrated intra‐ and intercellular signaling knowledge for multicellular
#'     omics analysis. Molecular Systems Biology. 2021 Mar;17(3):e9923.
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

    stitch <- stitch_network(
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
#' @seealso \itemize{
#'     \item{\code{\link{cosmos_pkn}}}
#'     \item{\code{\link{import_omnipath_interactions}}}
#' }
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
#'     \{code{\link{stitch_network}}.
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
