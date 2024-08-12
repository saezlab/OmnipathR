#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
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
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#

## Functions for importing interactions from OmniPath.
## The interactions database of OminPath consists of several datastes.

PPI_DATASETS <- c('omnipath', 'pathwayextra', 'kinaseextra', 'ligrecextra')
GRN_DATASETS <- c('dorothea', 'tf_target', 'collectri')


#' Interactions from OmniPath
#'
#' Interactions from the \url{https://omnipathdb.org/interactions} endpoint of
#' the OmniPath web service.
#' By default, it downloads only the "omnipath" dataset, which corresponds to
#' the curated causal interactions described in Turei et al. 2016.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -query_type
#'
#' @return A dataframe of molecular interactions.
#'
#' @examples
#' interactions = import_omnipath_interactions(
#'     resources = c('SignaLink3'),
#'     organism = 9606
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{omnipath}}}
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom magrittr %<>%
#' @importFrom rlang exec !!!
#' @export
import_omnipath_interactions <- function(...){

    args <- omnipath_args(list(...), query_type = 'interactions')
    defaults <- list(datasets = 'omnipath')
    args %<>% modifyList(defaults, .)

    exec(omnipath_query, !!!args)

}


#' Literature curated signaling pathways
#'
#' Imports interactions from the `omnipath` dataset of OmniPath, a dataset
#' that inherits most of its design and contents from the original OmniPath
#' core from the 2016 publication. This dataset consists of about 40k
#' interactions.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe of literature curated, post-translational signaling
#'     interactions.
#'
#' @param ... optional additional arguments, passed to
#' \code{\link{import_omnipath_interactions}}.
#'
#' @examples
#' pathways <- omnipath()
#' pathways
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_omnipath_interactions}}}
#'     \item{\code{\link{import_post_translational_interactions}}}
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
omnipath <- function(...){

    args <- omnipath_args(list(...), datasets = 'omnipath')

    exec(import_omnipath_interactions, !!!args)

}


#' Interactions from the `pathway extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=pathwayextra},
#' which contains activity flow interactions without literature reference.
#' The activity flow interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#'
#' @return A dataframe containing activity flow interactions between proteins
#' without literature reference
#'
#' @examples
#' interactions <-
#'     import_pathwayextra_interactions(
#'         resources = c('BioGRID', 'IntAct'),
#'         organism = 9606
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
import_pathwayextra_interactions <- function(...){

    args <- omnipath_args(list(...), datasets = 'pathwayextra')

    exec(import_omnipath_interactions, !!!args)

}


#' Interactions from the `kinase extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=kinaseextra},
#' which contains enzyme-substrate interactions without literature reference.
#' The enzyme-substrate interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe containing enzyme-substrate interactions without
#' literature reference
#'
#' @examples
#' interactions <-
#'    import_kinaseextra_interactions(
#'        resources = c('PhosphoPoint', 'PhosphoSite'),
#'        organism = 9606
#'    )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
import_kinaseextra_interactions <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    ...
){

    args <- omnipath_args(list(...), datasets = 'kinaseextra')

    exec(import_omnipath_interactions, !!!args)

}


#' Interactions from the `ligrec extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=ligrecextra},
#' which contains ligand-receptor interactions without literature reference.
#' The ligand-receptor interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe containing ligand-receptor interactions including
#' the ones without literature references
#'
#' @examples
#' interactions <- import_ligrecextra_interactions(
#'     resources = c('HPRD', 'Guide2Pharma'),
#'     organism = 9606
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
import_ligrecextra_interactions <- function(...){

    args <- omnipath_args(list(...), datasets = 'ligrecextra')

    exec(import_omnipath_interactions, !!!args)

}


#' All post-translational interactions from OmniPath
#'
#' Imports interactions from all post-translational datasets of OmniPath.
#' The datasets are "omnipath", "kinaseextra", "pathwayextra" and
#' "ligrecextra".
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -query_type
#'
#' @return A dataframe containing post-translational interactions
#'
#' @examples
#' interactions <-
#'     import_post_translational_interactions(
#'         resources = c('BioGRID')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#' @importFrom rlang %||% exec !!!
#' @export
import_post_translational_interactions <- function(...){


    args <- omnipath_args(list(...), query_type = 'interactions')
    args$datasets %<>% {. %||% PPI_DATASETS}

    exec(omnipath_query, !!!args)

}


#' TF-target interactions from DoRothEA
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=dorothea}
#' which contains transcription factor (TF)-target interactions from DoRothEA
#' \url{https://github.com/saezlab/DoRothEA}
#' DoRothEA is a comprehensive resource of transcriptional regulation,
#' consisting of 16 original resources, in silico TFBS prediction, gene
#' expression signatures and ChIP-Seq binding site analysis.
#'
#' @param dorothea_levels Vector detailing the confidence levels of the
#' interactions to be downloaded. In dorothea, every TF-target interaction
#' has a confidence score ranging from A to E, being A the most reliable
#' interactions.
#' By default we take A and B level interactions (\code{c(A, B)}).
#' It is to note that E interactions are not available in OmnipathR.
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A data frame of TF-target interactions from DoRothEA.
#'
#' @examples
#' dorothea_grn <- dorothea(
#'     resources = c('DoRothEA', 'ARACNe-GTEx_DoRothEA'),
#'     organism = 9606,
#'     dorothea_levels = c('A', 'B', 'C')
#' )
#' dorothea_grn
#'
#' @seealso \itemize{
#'     \item{\code{\link{collectri}}}
#'     \item{\code{\link{import_transcriptional_interactions}}}
#'     \item{\code{\link{import_tf_target_interactions}}}
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
#' @aliases import_dorothea_interactions
dorothea <- function(dorothea_levels = c('A', 'B'), ...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'dorothea'
    )

    exec(omnipath_query, !!!args)

}


#' Interactions from the TF-target dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_target},
#' which contains transcription factor-target protein coding gene
#' interactions. Note: this is not the only TF-target dataset in OmniPath,
#' `dorothea` is the other one and the `tf_mirna` dataset provides
#' TF-miRNA gene interactions.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe containing TF-target interactions
#'
#' @examples
#' interactions <-
#'     import_tf_target_interactions(
#'         resources = c('DoRothEA', 'SIGNOR')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#' @importFrom rlang exec !!!
#' @export
import_tf_target_interactions <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    ...
){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'tf_target'
    )

    exec(omnipath_query, !!!args)

}


#' All TF-target interactions from OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_target,dorothea},
#' which contains transcription factor-target protein coding gene
#' interactions.
#'
#' @param dorothea_levels Vector detailing the confidence levels of the
#' interactions to be downloaded. In dorothea, every TF-target interaction
#' has a confidence score ranging from A to E, being A the most reliable
#' interactions.
#' By default we take A and B level interactions (\code{c(A, B)}).
#' It is to note that E interactions are not available in OmnipathR.
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -query_type
#'
#' @return A dataframe containing TF-target interactions.
#'
#' @examples
#' grn <-
#'     import_transcriptional_interactions(
#'         resources = c('PAZAR', 'ORegAnno', 'DoRothEA')
#'     )
#' grn
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang %||% exec !!!
#' @export
import_transcriptional_interactions <- function(
    dorothea_levels = c('A', 'B'),
    ...
){

    args <- omnipath_args(list(...), query_type = 'interactions')
    args$datasets %<>% {. %||% GRN_DATASETS} %>% intersect(GRN_DATASETS)

    exec(omnipath_query, !!!args)

}


#' TF-target interactions from CollecTRI
#'
#' CollecTRI is a comprehensive resource of transcriptional regulation,
#' published in 2023, consisting of 14 resources and original literature
#' curation.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe of TF-target interactions.
#'
#' @examples
#' collectri_grn <- collectri()
#' collectri_grn
#'
#' @seealso \itemize{
#'     \item{\code{\link{import_transcriptional_interactions}}}
#'     \item{\code{\link{dorothea}}}
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @export
#' @importFrom rlang exec !!!
collectri <- function(...){

    args <-
        omnipath_args(
            list(...),
            query_type = 'interactions',
            datasets = 'collectri'
        )

    exec(omnipath_query, !!!args)

}


#' Interactions from the miRNA-target dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=mirnatarget},
#' which contains miRNA-mRNA interactions.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe containing miRNA-mRNA interactions
#'
#' @examples
#' interactions <-
#'     import_mirnatarget_interactions(
#'         resources = c('miRTarBase', 'miRecords')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
import_mirnatarget_interactions <- function(...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'mirnatarget'
    )

    exec(omnipath_query, !!!args)

}


#' Interactions from the TF-miRNA dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_mirna},
#' which contains transcription factor-miRNA gene interactions
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe containing TF-miRNA interactions
#'
#' @examples
#' interactions <-
#'     import_tf_mirna_interactions(
#'         resources = c('TransmiR')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @export
#' @importFrom rlang exec !!!
import_tf_mirna_interactions <- function(...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'tf_mirna'
    )

    exec(omnipath_query, !!!args)

}


#' lncRNA-mRNA interactions from OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=lncrna_mrna},
#' which contains lncRNA-mRNA interactions
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe containing lncRNA-mRNA interactions
#'
#' @examples
#' interactions <-
#'     import_lncrna_mrna_interactions(
#'         resources = c('ncRDeathDB')
#'     )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
import_lncrna_mrna_interactions <- function(...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'lncrna_mrna'
    )

    exec(omnipath_query, !!!args)

}


#' Small molecule-protein interactions from OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=small_molecule},
#' which contains small molecule-protein interactions. Small molecules
#' can be metabolites, intrinsic ligands or drug compounds.
#'
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @return A dataframe of small molecule-protein interactions
#'
#' @examples
#' # What are the targets of aspirin?
#' interactions <-
#'     import_small_molecule_protein_interactions(
#'         sources = 'ASPIRIN'
#'     )
#' # The prostaglandin synthases:
#' interactions
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
import_small_molecule_protein_interactions <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    ...
){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'small_molecule'
    )

    exec(omnipath_query, !!!args)

}


#' Imports all interaction datasets available in OmniPath
#'
#' The interaction datasets currently available in OmniPath:
#'
#' \itemize{
#'     \item{omnipath: the OmniPath data as defined in the 2016 paper, an
#'     arbitrary optimum between coverage and quality}
#'     \item{pathwayextra: activity flow interactions without literature
#'     references}
#'     \item{kinaseextra: enzyme-substrate interactions without literature
#'     reference}
#'     \item{ligrecextra: ligand-receptor interactions without
#'     literature reference}
#'     \item{collectri: transcription factor (TF)-target
#'     interactions from CollecTRI}
#'     \item{dorothea: transcription factor (TF)-target
#'     interactions from DoRothEA}
#'     \item{tf_target: transcription factor
#'     (TF)-target interactions from other resources}
#'     \item{mirnatarget: miRNA-mRNA interactions}
#'     \item{tf_mirna: TF-miRNA interactions}
#'     \item{lncrna_mrna: lncRNA-mRNA interactions}
#' }
#'
#' @return A dataframe containing all the datasets in the interactions query
#'
#' @param dorothea_levels The confidence levels of the dorothea
#' interactions (TF-target) which range from A to D. Set to A and B by
#' default.
#' @param types Character: interaction types, such as "transcriptional",
#'     "post_transcriptional", "post_translational", etc.
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @examples
#' interactions <- import_all_interactions(
#'     resources = c('HPRD', 'BioGRID'),
#'     organism = 9606
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @importFrom magrittr %<>% %>% extract2
#' @export
import_all_interactions <- function(
    dorothea_levels = c('A', 'B'),
    types = NULL,
    ...
){

    q_info <- query_info('interactions')
    datasets <- q_info %>% extract2('datasets') %>% setdiff(exclude)
    types %<>% if_null(q_info %>% extract2('types'))
    # it does not make sense without the type field
    fields %<>% c('type', 'dorothea_level') %>% unique

    args <- omnipath_args(
        list(...),
        query_type = 'ineractions',
        types = types,
        datasets = datasets,
        fields = fields
    )

    exec(omnipath_query, !!!args)

}


#' Retrieve a list of interaction resources available in Omnipath
#'
#' Gets the names of the resources from
#' \url{https://omnipath.org/interactions}.
#'
#' @param dataset a dataset within the interactions query type. Currently
#' available datasets are `omnipath`, `kinaseextra`, `pathwayextra`,
#' `ligrecextra`, `dorothea`, `tf_target`, `tf_mirna`, `mirnatarget` and
#' `lncrna_mrna`
#'
#' @return character vector with the names of the interaction databases
#'
#' @examples
#' get_interaction_resources()
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{import_omnipath_interactions}}}
#'     \item{\code{\link{import_pathwayextra_interactions}}}
#'     \item{\code{\link{import_kinaseextra_interactions}}}
#'     \item{\code{\link{import_ligrecextra_interactions}}}
#'     \item{\code{\link{import_mirnatarget_interactions}}}
#'     \item{\code{\link{import_dorothea_interactions}}}
#' }
#' @export
get_interaction_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'interactions', datasets = dataset))

}


#' Create a vector with dataset names from an environment with logical
#' variables.
#'
#' @param envir Environment from the calling function where dataset names
#'     present as logical variables.
#'
#' @importFrom magrittr %>% extract
#' @importFrom purrr keep
#'
#' @noRd
select_interaction_datasets <- function(envir){

    envir %>%
    as.list %>%
    extract(PPI_DATASETS) %>%
    keep(identity) %>%
    names

}


#' Interactions having references
#'
#' @param data An interaction data frame.
#' @param resources Character: consider only these resources. If `NULL`,
#'     records with any reference will be accepted.
#'
#' @return A subset of the input interaction data frame.
#'
#' @examples
#' cc <- import_post_translational_interactions(resources = 'CellChatDB')
#' with_references(cc, 'CellChatDB')
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows distinct across
#' @importFrom purrr map
#' @importFrom tidyselect everything
#' @export
with_references <- function(data, resources = NULL){

    resources %>%
    ensure_list %>%
    map(
        .with_references,
        data = data
    ) %>%
    bind_rows() %>%
    distinct(across(everything()))

}


#' Interactions having references
#'
#' @param data An interaction data frame.
#' @param resource Character: consider only this resource. If `NULL`, records
#'     with any reference will be accepted.
#'
#' @return A subset of the input interaction data frame.
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom dplyr filter
#' @noRd
.with_references <- function(data, resource = NULL){

    predicate <- `if`(
        is.null(resource),
        function(references){
            !is.na(references)
        },
        function(references){
            str_detect(
                references,
                sprintf('%s[^:]*:', resource)
            )
        }
    )

    data %>%
    filter(predicate(references))

}


#' Datasets in the OmniPath Interactions database
#'
#' @return Character: labels of interaction datasets.
#'
#' @examples
#' interaction_datasets()
#'
#' @export
interaction_datasets <- function() {

    query_info('interactions')$datasets

}


#' Interaction types in the OmniPath Interactions database
#'
#' @return Character: labels of interaction types.
#'
#' @examples
#' interaction_types()
#'
#' @export
interaction_types <- function() {

    query_info('interactions')$types

}


#' Create a column with dataset names listed
#'
#' From logical columns for each dataset, here we create a column that is
#' a list of character vectors, containing dataset labels.
#'
#' @param data Interactions data frame with dataset columns (i.e. queried
#'     with the option `fields = "datasets"`).
#' @param remove_logicals Logical: remove the per dataset logical columns.
#'
#' @return The input data frame with the new column "datasets" added.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang syms
#' @importFrom dplyr rowwise mutate c_across select
#' @importFrom tidyselect all_of
#' @export
datasets_one_column <- function(data, remove_logicals = TRUE) {

    dataset_cols <- data %>% colnames %>% intersect(interaction_datasets())

    data %>%
    rowwise() %>%
    mutate(datasets = list(dataset_cols[c_across(all_of(dataset_cols))])) %>%
    `if`(remove_logicals, select(., -dataset_cols), .)

}
