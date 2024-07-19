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


#' Imports interactions from the `omnipath` dataset of Omnipath
#'
#' Imports the database from \url{https://omnipathdb.org/interactions}, which
#' contains only interactions supported by literature references.
#' This part of the interaction database compiled a similar way as it has
#' been presented in the first paper describing OmniPath (Turei et al. 2016).
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param datasets Names of the interaction datasets to download: omnipath
#' (by default). Other possiblites are: pathwayextra, kinaseextra,
#' ligrecextra, dorothea,tf_target, mirnatarget, tf_mirna, lncrna_mrna.
#' The user can select multiple datasets as for example: c('omnipath',
#' 'pathwayextra', 'kinaseextra')
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
#'
#' @return A dataframe of protein-protein interactions
#'
#' @examples
#' interactions = import_omnipath_interactions(
#'     resources = c('SignaLink3'),
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
#' @aliases import_Omnipath_Interactions import_OmniPath_Interactions
#' @importFrom rlang exec !!!
#' @export
import_omnipath_interactions <- function(
    resources = NULL,
    organism = 'human',
    datasets = 'omnipath',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(query_type = 'interactions')

    exec(import_omnipath, !!!args)

}


#' Literature curated signaling pathways
#'
#' Imports interactions from the `omnipath` dataset of Omnipath, a dataset
#' that inherits most of its design and contents from the original OmniPath
#' core from the 2016 publication. This dataset consists of about 40k
#' interactions.
#'
#' @return A dataframe of literature curated, post-translational signaling
#'     interactions.
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
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
omnipath <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(datasets = 'omnipath')

    exec(import_omnipath_interactions, !!!args)

}


# Aliases (old names) to be deprecated
#' @rdname import_omnipath_interactions
#' @param ... Passed to \code{import_omnipath_interactions}.
#' @export
#'
#' @noRd
import_Omnipath_Interactions <- function(...){
    .Deprecated("import_omnipath_interactions")
    import_omnipath_interactions(...)
}


# Aliases (old names) to be deprecated
#' @rdname import_omnipath_interactions
#' @param ... Passed to \code{import_omnipath_interactions}.
#' @export
#'
#' @noRd
import_OmniPath_Interactions <- function(...){
    .Deprecated("import_omnipath_interactions")
    import_omnipath_interactions(...)
}


#' Imports interactions from the `pathway extra` dataset of Omnipath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=pathwayextra},
#' which contains activity flow interactions without literature reference.
#' The activity flow interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
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
#' @aliases import_PathwayExtra_Interactions
#' @importFrom rlang exec !!!
#' @export
import_pathwayextra_interactions <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(datasets = 'pathwayextra')

    exec(import_omnipath_interactions, !!!args)

}


# Aliases (old names) to be deprecated
#'
#' @rdname import_pathwayextra_interactions
#' @param ... Passed to \code{import_pathwayextra_interactions}.
#' @export
#'
#' @noRd
import_PathwayExtra_Interactions <- function(...){
    .Deprecated("import_pathwayextra_interactions")
    import_pathwayextra_interactions(...)
}


#' Imports interactions from the `kinase extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=kinaseextra},
#' which contains enzyme-substrate interactions without literature reference.
#' The enzyme-substrate interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... Optional additional arguments.
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
#' @aliases import_KinaseExtra_Interactions
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
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(datasets = 'kinaseextra')

    exec(import_omnipath_interactions, !!!args)

}


# Aliases (old names) to be deprecated
#' @rdname import_kinaseextra_interactions
#' @param ... Passed to \code{import_kinaseextra_interactions}.
#' @export
#'
#' @noRd
import_KinaseExtra_Interactions <- function(...){
    .Deprecated("import_kinaseextra_interactions")
    import_kinaseextra_interactions(...)
}


#' Imports interactions from the `ligrec extra` dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=ligrecextra},
#' which contains ligand-receptor interactions without literature reference.
#' The ligand-receptor interactions supported by literature references
#' are part of the `omnipath` dataset.
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
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
#' @aliases import_LigrecExtra_Interactions
#' @importFrom rlang exec !!!
#' @export
import_ligrecextra_interactions <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(datasets = 'ligrecextra')

    exec(import_omnipath_interactions, !!!args)

}


# Aliases (old names) to be deprecated
#' @rdname import_ligrecextra_interactions
#' @param ... Passed to \code{import_ligrecextra_interactions}.
#' @export
#'
#' @noRd
import_LigrecExtra_Interactions <- function(...){
    .Deprecated("import_ligrecextra_interactions")
    import_ligrecextra_interactions(...)
}


#' All post-translational interactions from OmniPath
#'
#' Imports interactions from all post-translational datasets of OmniPath.
#' The datasets are "omnipath", "kinaseextra", "pathwayextra" and
#' "ligrecextra".
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param exclude Character: datasets or resources to exclude
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
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
import_post_translational_interactions <- function(
    resources = NULL,
    organism = 'human',
    exclude = NULL,
    references_by_resource = TRUE,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){


    args <- omnipath_args(query_type = 'interactions')
    args$datasets %<>% {. %||% PPI_DATASETS}

    exec(import_omnipath, !!!args)

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
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param dorothea_levels Vector detailing the confidence levels of the
#' interactions to be downloaded. In dorothea, every TF-target interaction
#' has a confidence score ranging from A to E, being A the most reliable
#' interactions.
#' By default we take A and B level interactions (\code{c(A, B)}).
#' It is to note that E interactions are not available in OmnipathR.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources. In case of
#' DoRothEA this is not desirable for most of the applications. For most of
#' the other interaction querying functions it is `FALSE` by default.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
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
#' @aliases import_tfregulons_interactions import_dorothea_interactions
#' @importFrom rlang exec !!!
#' @export
dorothea <- function(
    resources = NULL,
    organism = 'human',
    dorothea_levels = c('A', 'B'),
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = TRUE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(query_type = 'interactions', datasets = 'dorothea')

    exec(import_omnipath, !!!args)

}


# Aliases (old names) to be deprecated
#' @rdname dorothea
#' @param ... Passed to \code{\link{dorothea}}.
#' @export
#'
#' @noRd
import_dorothea_interactions <- function(...){
    .Deprecated("dorothea")
    dorothea(...)
}


#' @rdname dorothea
#' @param ... Passed to \code{\link{dorothea}}.
#' @export
#'
#' @noRd
import_tfregulons_interactions <- function(...){
    .Deprecated("dorothea")
    dorothea(...)
}


#' Imports interactions from the TF-target dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_target},
#' which contains transcription factor-target protein coding gene
#' interactions. Note: this is not the only TF-target dataset in OmniPath,
#' `dorothea` is the other one and the `tf_mirna` dataset provides
#' TF-miRNA gene interactions.
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... Optional additional arguments
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
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(query_type = 'interactions', datasets = 'tf_target')

    exec(import_omnipath, !!!args)

}


#' Imports all TF-target interactions from OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_target,dorothea},
#' which contains transcription factor-target protein coding gene
#' interactions.
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param dorothea_levels Vector detailing the confidence levels of the
#' interactions to be downloaded. In dorothea, every TF-target interaction
#' has a confidence score ranging from A to E, being A the most reliable
#' interactions.
#' By default we take A and B level interactions (\code{c(A, B)}).
#' It is to note that E interactions are not available in OmnipathR.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... Optional additional arguments.
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
    resources = NULL,
    organism = 'human',
    dorothea_levels = c('A', 'B'),
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(query_type = 'interactions')
    args$datasets %<>% {. %||% GRN_DATASETS} %>% intersect(GRN_DATASETS)

    exec(import_omnipath, !!!args)

}


#' TF-target interactions from CollecTRI
#'
#' CollecTRI is a comprehensive resource of transcriptional regulation,
#' published in 2023, consisting of 14 resources and original literature
#' curation.
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources. In case of
#' CollecTRI this is not desirable for most of the applications. For most of
#' the other interaction querying functions it is `FALSE` by default.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... Optional additional arguments, passed to
#'     \code{\link{import_transcriptional_interactions}}.
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
collectri <- function(
    resources = NULL,
    organism = 'human',
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = TRUE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(datasets = 'collectri')

    exec(import_transcriptional_interactions, !!!args)

}


#' Imports interactions from the miRNA-target dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=mirnatarget},
#' which contains miRNA-mRNA interactions.
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
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
#' @aliases import_miRNAtarget_Interactions
#' @importFrom rlang exec !!!
#' @export
import_mirnatarget_interactions <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(
        query_type = 'interactions',
        datasets = 'mirnatarget'
    )

    exec(import_omnipath, !!!args)

}


# Aliases (old names) to be deprecated
#' @rdname import_mirnatarget_interactions
#' @param ... Passed to \code{import_mirnatarget_interactions}.
#' @export
#'
#' @noRd
import_miRNAtarget_Interactions <- function(...){
    .Deprecated("import_mirnatarget_interactions")
    import_mirnatarget_interactions(...)
}


#' Imports interactions from the TF-miRNA dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=tf_mirna},
#' which contains transcription factor-miRNA gene interactions
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
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
import_tf_mirna_interactions <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(query_type = 'interactions', datasets = 'tf_mirna')

    exec(import_omnipath, !!!args)

}


#' Imports interactions from the lncRNA-mRNA dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=lncrna_mrna},
#' which contains lncRNA-mRNA interactions
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
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
import_lncrna_mrna_interactions <- function(
    resources = NULL,
    organism = 'human',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(
        query_type = 'interactions',
        datasets = 'lncrna_mrna'
    )

    exec(import_omnipath, !!!args)

}


#' Interactions from the small molecule-protein dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=small_molecule},
#' which contains small molecule-protein interactions. Small molecules
#' can be metabolites, intrinsic ligands or drug compounds.
#'
#' @param resources interactions not reported in these databases are
#'     removed. See \code{\link{get_interaction_resources}} for more
#'     information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat.
#' @param fields Optional fields to be added.
#' @param default_fields whether to include the default fields (columns) for
#'     the query type. If FALSE, only the fields defined by the user in the
#'     `fields` argument will be added.
#' @param references_by_resource If \code{FALSE}, removes the resource name
#'     prefixes from the references (PubMed IDs); this way the information
#'     which reference comes from which resource will be lost and the PubMed
#'     IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param ... optional additional arguments
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
    genesymbol_resource = NULL,
    ...
){

    args <- omnipath_args(
        query_type = 'interactions',
        datasets = 'small_molecule'
    )

    exec(import_omnipath, !!!args)

}


#' Imports all interaction datasets available in OmniPath
#'
#' The interaction datasets currently available in OmniPath:
#'
#' omnipath: the OmniPath data as defined in the paper, an arbitrary optimum
#' between coverage and quality
#' pathwayextra: activity flow interactions without literature reference
#' kinaseextra: enzyme-substrate interactions without literature reference
#' ligrecextra: ligand-receptor interactions without literature reference
#' dorothea: transcription factor (TF)-target interactions from DoRothEA
#' tf_target: transcription factor (TF)-target interactions from other
#' resources
#' mirnatarget: miRNA-mRNA interactions
#' tf_mirna: TF-miRNA interactions
#' lncrna_mrna: lncRNA-mRNA interactions
#'
#' @return A dataframe containing all the datasets in the interactions query
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Character or integer: Name or NCBI Taxonomy ID of one or
#'     organisms. The web service currently provides interactions for
#'     human, mouse and rat. For other organisms, the data will be translated
#'     by orthologous gene pairs from human. In this case, only one organism
#'     can be provided. If miRNA, lncRNA or small molecule datasets included,
#'     orthology translation is not possible and will remove the interactions
#'     with non-protein partners.
#' @param dorothea_levels The confidence levels of the dorothea
#' interactions (TF-target) which range from A to D. Set to A and B by
#' default.
#' @param exclude Character: datasets or resources to exclude.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param strict_evidences Logical: restrict the evidences to the queried
#' datasets and resources. If set to FALSE, the directions and effect signs
#' and references might be based on other datasets and resources.
#' @param genesymbol_resource Character: either "uniprot" or "ensembl". The
#'     former leaves intact the gene symbols returned by the web service,
#'     originally set from UniProt. The latter updates the gene symbols from
#'     Ensembl, which uses a slightly different gene symbol standard. In this
#'     case a few records will be duplicated, where Ensembl provides ambiguous
#'     translation.
#' @param types Character: interaction types, such as "transcriptional",
#' "post_transcriptional", "post_translational", etc.
#' @param ... optional additional arguments
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
#' @aliases import_AllInteractions
#'
#' @importFrom rlang exec !!!
#' @importFrom magrittr %<>% %>% extract2
#' @export
import_all_interactions <- function(
    resources = NULL,
    organism = 'human',
    dorothea_levels = c('A', 'B'),
    exclude = NULL,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    strict_evidences = FALSE,
    genesymbol_resource = NULL,
    types = NULL,
    ...
){

    q_info <- query_info('interactions')
    datasets <- q_info %>% extract2('datasets') %>% setdiff(exclude)
    types %<>% if_null(q_info %>% extract2('types'))
    # it does not make sense without the type field
    fields %<>% c('type', 'dorothea_level') %>% unique

    args <- omnipath_args(
        query_type = 'ineractions',
        types = types,
        datasets = datasets,
        fields = fields
    )

    exec(import_omnipath, !!!args)

}


# Aliases (old names) to be deprecated
#' @rdname import_all_interactions
#' @param ... Passed to \code{import_all_interactions}.
#' @export
import_AllInteractions <- function(...){
    .Deprecated("import_all_interactions")
    import_all_interactions(...)
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
#' @export
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
#'
#' @aliases get_interaction_databases
get_interaction_resources <- function(dataset = NULL){

    return(get_resources(query_type = 'interactions', datasets = dataset))

}

# Aliases (old names) to be deprecated
#' @rdname get_interaction_resources
#' @param ... Passed to \code{get_interaction_resources}.
#' @export
#'
#' @noRd
get_interaction_databases <- function(...){
    .Deprecated("get_interaction_resources")
    get_interaction_resources(...)
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
