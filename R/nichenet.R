#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2021
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


#' Builds a NicheNet signaling network
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


#' Builds a NicheNet ligand-receptor network
#'
#' Builds ligand-receptor network prior knowledge for NicheNet using multiple
#' resources.
#'
#' @param omnipath List with paramaters to be passed to
#'     \code{\link{nichenet_lr_network_omnipath}}
#' @param guide2pharma List with paramaters to be passed to
#'     \code{\link{nichenet_lr_network_guide2pharma}}
#' @param ramilowski List with paramaters to be passed to
#'     \code{\link{nichenet_lr_network_ramilowski}}
#' @export
#'
#' @seealso \code{\link{nichenet_lr_network_omnipath},
#'     \link{nichenet_lr_network_guide2pharma},
#'     \link{nichenet_lr_network_ramilowski}}
nichenet_lr_network <- function(
    omnipath = list(),
    guide2pharma = NULL,
    ramilowski = NULL
){



}


#' Builds a NicheNet gene regulatory network
#'
#' Builds gene regulatory network prior knowledge for NicheNet using multiple
#' resources.
#'
#' @param omnipath List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_omnipath}}
#' @export
#'
#' @seealso \code{\link{nichenet_gr_network_omnipath}}
nichenet_gr_network <- function(
    omnipath = list()
){



}


#' Builds signaling network prior knowledge for NicheNet using OmniPath
#'
#' This method never downloads the `ligrecextra` dataset because the
#' ligand-receptor interactions are supposed to come from \code{
#' \link{nichenet_lr_network_omnipath}}.
#'
#' @param min_curation_effort Lower threshold for curation effort
#' @param ... Passed to \code{\link{import_post_translational_interactions}}
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @seealso
nichenet_signaling_network_omnipath <- function(
    min_curation_effort = 0,
    ...
){

    args <- as.list(...)
    args$exclude %<>% union('ligrecextra')
    args$entity_types <- 'protein'

    do.call(import_post_translational_interactions, args) %>%
    omnipath_interactions_postprocess()

}


#' Builds ligand-receptor network prior knowledge for NicheNet using OmniPath
#'
#' This method never downloads the `ligrecextra` dataset because the
#' ligand-receptor interactions are supposed to come from \code{
#' \link{nichenet_lr_network_omnipath}}.
#'
#' @param min_curation_effort Lower threshold for curation effort
#' @param ... Passed to \code{\link{import_intercell_network}}
#' @importFrom magrittr %>%
#' @export
#'
#' @seealso
nichenet_lr_network_omnipath <- function(
    min_curation_effort = 0,
    ...
){

    import_intercell_network(...) %>%
    omnipath_interactions_postprocess()

}


#' Builds gene regulatory network prior knowledge for NicheNet using OmniPath
#'
#' This method never downloads the `ligrecextra` dataset because the
#' ligand-receptor interactions are supposed to come from \code{
#' \link{nichenet_lr_network_omnipath}}.
#'
#' @param min_curation_effort Lower threshold for curation effort
#' @param ... Passed to \code{\link{import_transcriptional_interactions}}
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @seealso
nichenet_signaling_network_omnipath <- function(
    min_curation_effort = 0,
    ...
){

    args <- as.list(...)
    args$exclude %<>% union('ligrecextra')
    args$entity_types <- 'protein'

    do.call(import_transcriptional_interactions, args) %>%
    omnipath_interactions_postprocess()

}


#' Processes OmniPath interactions table into NicheNet format
#'
#' @importFrom dplyr select mutate distinct
#' @importFrom magrittr %>%
omnipath_interactions_postprocess <- function(interactions){

    interactions %>%
    select(from = source_genesymbol, to = target_genesymbol, is_directed) %>%
    mutate(
        source = ifelse(
            is_directed,
            'omnipath_directed',
            'omnipath_undirected'
        ),
        database = 'omnipath'
    ) %>%
    distinct() %>%
    select(-is_directed)

}


#' NicheNet signaling network from PathwayCommons
#'
#' Builds signaling network prior knowledge for NicheNet using PathwayCommons.
#'
#' @param interaction_types Character vector with PathwayCommons interaction
#'     types. Please refer to the default value and the PathwayCommons
#'     webpage.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter
#' @importFrom readr read_tsv cols
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

    harmonizome_nichenet(datasets, dataset_names)

}


#' Combines multiple Harmonizome datasets and converts them to NicheNet format
#'
#' @importFrom dplyr mutate bind_rows
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @seealso \code{\link{harmonizome_download},
#'     \link{harmonizome_nichenet_process}}
harmonizome_nichenet <- function(datasets, dataset_names){

    datasets %>%
    map(harmonizome_nichenet_process) %>%
    bind_rows() %>%
    mutate(
        source = sprintf('harmonizome_%s', dataset_names[source]),
        database = 'harmonizome'
    )

}


#' Processes a table downloaded from Harmonizome to NicheNet format
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang sym
#' @importFrom dplyr select mutate
#' @importFrom stringr str_split_fixed
#' @seealso \code{\link{harmonizome_download}, \link{harmonizome_nichenet}}
harmonizome_nichenet_process <- function(dataset){

    target_desc_col <- c('geotf', 'geokinase', 'geogene')
    to_col <- `if`(
        dataset %in% target_desc_col,
        sym('target_desc'),
        sym('target')
    )
    target_proc <- list(
        geotf = toupper,
        msigdbonc = function(x){
            x %>% str_split_fixed('[._]', 2) %>% `[`(,1)
        }
    )

    dataset %>%
    harmonizome_download() %>%
    select(from = source, to = !!to_col) %>%
    {`if`(
        dataset %in% names(target_proc),
        mutate(., to = target_proc[[dataset]](to)),
        .
    )} %>%
    mutate(source = dataset)

}


#' NicheNet signaling network from Vinayagam
#'
#' Builds signaling network prior knowledge for NicheNet using Vinayagam 2011
#' Supplementary Table S6
#'
#' Find out more at https://doi.org/10.1126/scisignal.2001699
#'
#' @importFrom dplyr %>% select mutate
#' @importFrom readxl read_xls
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
#' @importFrom dplyr select mutate distinct
#' @importFrom magrittr %>%
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
#' @importFrom magrittr %>% `n'est pas`
#' @importFrom dplyt select mutate filter
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


#' NicheNet signaling network from InWeb InBioMap
#'
#' Builds signaling network prior knowledge for NicheNet from the InWeb
#' InBioMap database.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
#'
#' @seealso \code{\link{nichenet_signaling_network}, \link{inbiomap}}
nichenet_signaling_network_inbiomap <- function(...){

    inbiomap(...) %>%
    select(
        from = genesymbol_a,
        to = genesymbol_b,
        source = 'inweb_interaction',
        database = 'inweb_inbiomap'
    )

}


#' Ligand-receptor network from Guide to Pharmacology
#'
#' Downloads ligand-receptor interactions from the Guide to Pharmacology
#' database and converts it to a format suitable for NicheNet.
#'
#' @return Data frame with ligand-receptor interactions in NicheNet format.
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @seealso \code{\link{nichenet_lr_network}}
nichenet_lr_network_guide2pharma <- function(){

    guide2pharma_download() %>%
    filter(
        target_species == 'Human' &
        ligand_species == 'Human'
    ) %>%
    nichenet_common_postprocess(
        source = 'pharmacology',
        database = 'guide2pharmacology',
        from_col = ligand_gene_symbol,
        to_col = target_gene_symbol
    )

}


#' Ligand-receptor network from Ramilowski 2015
#'
#' Downloads ligand-receptor interactions from Supplementary Table 2 of the
#' paper 'A draft network of ligand–receptor-mediated multicellular signalling
#' in human' (Ramilowski et al. 2015,
#' https://www.nature.com/articles/ncomms8866). It converts the downloaded
#' table to a format suitable for NicheNet.
#'
#' @param evidences Character: evidence types, "literature supported",
#' "putative" or both.
#'
#' @return Data frame with ligand-receptor interactions in NicheNet format.
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @seealso \code{\link{nichenet_lr_network}}
nichenet_lr_network_ramilowski <- function(
    evidences = c('literature supported', 'putative')
){

    ramilowski_download() %>%
    filter(Pair.Evidence %in% evidences) %>%
    nichenet_common_postprocess(
        source = 'ramilowski_known',
        database = 'ramilowski',
        from_col = Ligand.ApprovedSymbol,
        to_col = Receptor.ApprovedSymbol
    )

}


#' NicheNet gene regulatory network from Harmonizome
#'
#' Builds gene regulatory network prior knowledge for NicheNet using
#' Harmonizome
#'
#' @param datasets The datasets to use. For possible values please refer to
#'     default value and the Harmonizome webpage.
#' @importFrom magrittr %>%
#' @importFrom dplyr rename
#' @export
nichenet_gr_network_harmonizome <- function(
    datasets = c(
        'cheappi',
        'encodetfppi',
        'jasparpwm',
        'transfac',
        'transfacpwm',
        'motifmap',
        'geotf',
        'geokinase',
        'geogene'
    ),
    ...
){

    dataset_names <- list(
        cheappi = 'CHEA',
        encodetfppi = 'ENCODE',
        jaspar = 'JASPAR',
        transfac = 'TRANSFAC_CUR',
        transfacpwm = 'TRANSFAC',
        motifmap = 'MOTIFMAP',
        geotf = 'GEO_TF',
        geokinase = 'GEO_KINASE',
        geogene = 'GEO_GENE',
        msigdbonc = 'MSIGDB_GENE'
    )

    harmonizome_nichenet(datasets, dataset_names) %>%
    rename(from = to, to = from)

}


#' NicheNet gene regulatory network from RegNetwork
#'
#' Builds a gene regulatory network using data from the RegNetwork database
#' and converts it to a format suitable for NicheNet.
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @seealso \code{\link{regnetwork_download}}
nichenet_gr_network_regnetwork <- function(){

    regnetwork_download() %>%
    filter(
        source_type == 'protein' &
        target_type == 'protein'
    ) %>%
    nichenet_common_postprocess(
        source = 'regnetwork_source',
        database = 'regnetwork'
    )

}


#' NicheNet gene regulatory network from TRRUST
#'
#' Builds a gene regulatory network using data from the TRRUST database
#' and converts it to a format suitable for NicheNet.
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \code{\link{trrust_download}}
nichenet_gr_network_trrust <- function(){

    trrust_download() %>%
    nichenet_common_postprocess(
        source = 'trrust',
        database = 'trrust'
    )

}


#' NicheNet gene regulatory network from HTRIdb
#'
#' Builds a gene regulatory network using data from the HTRIdb database
#' and converts it to a format suitable for NicheNet.
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \code{\link{htridb_download}}
nichenet_gr_network_htridb <- function(){

    htridb_download() %>%
    nichenet_common_postprocess(
        source = 'HTRIDB',
        database = 'HTRIDB',
        from_col = SYMBOL_TF,
        to_col = SYMBOL_TG
    )

}


#' Common postprocessing from building a NicheNet format network table
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !! enquo
#' @importFrom dplyr select distinct mutate
nichenet_common_postprocess <- function(
    data,
    source,
    database,
    from_col = source_genesymbol,
    to_col = target_genesymbol
){

    from_col <- enquo(from_col)
    to_col <- enquo(to_col)

    data %>%
    select(
        from = !!from_col,
        to = !!to_col
    ) %>%
    distinct() %>%
    mutate(
        source = source,
        database = database
    )

}