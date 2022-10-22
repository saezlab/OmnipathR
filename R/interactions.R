#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2022
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

## Functions for importing interactions from OmniPath.
## The interactions database of OminPath consists of several datastes.

#' Imports interactions from the `omnipath` dataset of Omnipath
#'
#' Imports the database from \url{https://omnipathdb.org/interactions}, which
#' contains only interactions supported by literature references.
#' This part of the interaction database compiled a similar way as it has
#' been presented in the first paper describing OmniPath (Turei et al. 2016).
#'
#' @return A dataframe of protein-protein interactions
#' @export
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
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
#' @param ... optional additional arguments
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
import_omnipath_interactions <- function(
    resources = NULL,
    organism = 9606,
    datasets = 'omnipath',
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = datasets,
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

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
#' @return A dataframe containing activity flow interactions between proteins
#' without literature reference
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose one of those: 9606 human (default), 10116 rat or 10090 Mouse.
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param ... optional additional arguments
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
import_pathwayextra_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'pathwayextra',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

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
#' @return A dataframe containing enzyme-substrate interactions without
#' literature reference
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param ... Optional additional arguments.
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
import_kinaseextra_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'kinaseextra',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

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
#' @return A dataframe containing ligand-receptor interactions including
#' the ones without literature references
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param ... optional additional arguments
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
import_ligrecextra_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'ligrecextra',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

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
#' @return A dataframe containing post-translational interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param exclude Character: datasets or resources to exclude
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param ... optional additional arguments
#'
#' @importFrom rlang %||% exec !!!
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
import_post_translational_interactions <- function(
    resources = NULL,
    organism = 9606,
    exclude = NULL,
    references_by_resource = TRUE,
    ...
){

    args <- list(...)
    args$datasets <-
        args$datasets %||%
        c('omnipath', 'pathwayextra', 'kinaseextra', 'ligrecextra')

    args %<>% merge_lists(
        list(
            query_type = 'interactions',
            resources = resources,
            organism = organism,
            references_by_resource = references_by_resource,
            exclude = exclude
        )
    )


    exec(import_omnipath, !!!args)

}


#' From the OmniPath webservice imports interactions from the
#' DoRothEA dataset
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=dorothea}
#' which contains transcription factor (TF)-target interactions from DoRothEA
#' \url{https://github.com/saezlab/DoRothEA}
#'
#' @return A dataframe containing TF-target interactions from DoRothEA
#' @export
#'
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
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
#' @param ... optional additional arguments
#'
#' @examples
#' interactions <- import_dorothea_interactions(
#'     resources = c('DoRothEA', 'ARACNe-GTEx_DoRothEA'),
#'     organism = 9606,
#'     dorothea_levels = c('A', 'B', 'C')
#' )
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_TFregulons_Interactions import_tfregulons_interactions
import_dorothea_interactions <- function(
    resources = NULL,
    organism = 9606,
    dorothea_levels = c('A', 'B'),
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        dorothea_levels = dorothea_levels,
        datasets = 'dorothea',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

}


# Aliases (old names) to be deprecated
#' @rdname import_dorothea_interactions
#' @param ... Passed to \code{import_dorothea_interactions}.
#' @export
#'
#' @noRd
import_TFregulons_Interactions <- function(...){
    .Deprecated("import_dorothea_interactions")
    import_dorothea_interactions(...)
}


#' @rdname import_dorothea_interactions
#' @param ... Passed to \code{import_dorothea_interactions}.
#' @export
#'
#' @noRd
import_tfregulons_interactions <- function(...){
    .Deprecated("import_dorothea_interactions")
    import_dorothea_interactions(...)
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
#' @return A dataframe containing TF-target interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param ... Optional additional arguments
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
import_tf_target_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'tf_target',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

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
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
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
#' @param ... Optional additional arguments.
#'
#' @return A dataframe containing TF-target interactions.
#'
#' @examples
#' interactions <-
#'     import_transcriptional_interactions(
#'         resources = c('PAZAR', 'ORegAnno', 'DoRothEA')
#'     )
#'
#'
#' @export
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang %||% exec !!!
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
import_transcriptional_interactions <- function(
    resources = NULL,
    organism = 9606,
    dorothea_levels = c('A', 'B'),
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    is_directed <- NULL

    args <- list(...)
    tr_datasets <- c('dorothea', 'tf_target')
    args$datasets %<>% {. %||% tr_datasets} %>% intersect(tr_datasets)

    result <-
        exec(
            import_omnipath,
            query_type = 'interactions',
            exclude = exclude,
            dorothea_levels = dorothea_levels,
            organism = organism,
            resources = resources,
            references_by_resource = references_by_resource,
            !!!args
        )

    return(result)

}


#' Imports interactions from the miRNA-target dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=mirnatarget},
#' which contains miRNA-mRNA interactions.
#'
#' @return A dataframe containing miRNA-mRNA interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param ... optional additional arguments
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
import_mirnatarget_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'mirnatarget',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

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
#' @return A dataframe containing TF-miRNA interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param ... optional additional arguments
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
import_tf_mirna_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'tf_mirna',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

}


#' Imports interactions from the lncRNA-mRNA dataset of OmniPath
#'
#' Imports the dataset from:
#' \url{https://omnipathdb.org/interactions?datasets=lncrna_mrna},
#' which contains lncRNA-mRNA interactions
#'
#' @return A dataframe containing lncRNA-mRNA interactions
#' @export
#' @param resources interactions not reported in these databases are
#' removed. See \code{\link{get_interaction_resources}} for more information.
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param fields The user can define here the fields to be added. If used, set
#' the next argument, `default_fields`, to FALSE.
#' @param default_fields whether to include the default fields (columns) for
#' the query type. If FALSE, only the fields defined by the user in the
#' `fields` argument will be added.
#' @param references_by_resource if FALSE, removes the resource name prefixes
#' from the references (PubMed IDs); this way the information which reference
#' comes from which resource will be lost and the PubMed IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
#' @param ... optional additional arguments
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
import_lncrna_mrna_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'lncrna_mrna',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

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
#' @param organism Interactions are available for human, mouse and rat.
#'     Choose among: 9606 human (default), 10116 rat and 10090 Mouse.
#' @param fields Optional fields to be added.
#' @param default_fields whether to include the default fields (columns) for
#'     the query type. If FALSE, only the fields defined by the user in the
#'     `fields` argument will be added.
#' @param references_by_resource If \code{FALSE}, removes the resource name
#'     prefixes from the references (PubMed IDs); this way the information
#'     which reference comes from which resource will be lost and the PubMed
#'     IDs will be unique.
#' @param exclude Character: datasets or resources to exclude.
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
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{import_all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
import_small_molecule_protein_interactions <- function(
    resources = NULL,
    organism = 9606,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    exclude = NULL,
    ...
){

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        datasets = 'small_molecule',
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        exclude = exclude,
        ...
    )

    return(result)

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
#' @param organism Interactions are available for human, mouse and rat.
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
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
#' @param ... optional additional arguments
#'
#' @examples
#' interactions <- import_all_interactions(
#'     resources = c('HPRD', 'BioGRID'),
#'     organism = 9606
#' )
#'
#' @importFrom magrittr %<>% %>% extract2
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{get_interaction_resources}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @aliases import_AllInteractions
import_all_interactions <- function(
    resources = NULL,
    organism = 9606,
    dorothea_levels = c('A', 'B'),
    exclude = NULL,
    fields = NULL,
    default_fields = TRUE,
    references_by_resource = TRUE,
    ...
){

    all_datasets <-
        omnipath_url('queries/interactions?format=json') %>%
        safe_json(path = .) %>%
        extract2('datasets') %>%
        setdiff(exclude)

    # it does not make sense without the type field
    fields %<>% c('type', 'dorothea_level') %>% unique

    result <- import_omnipath(
        query_type = 'interactions',
        resources = resources,
        organism = organism,
        dorothea_levels = dorothea_levels,
        exclude = exclude,
        datasets = all_datasets,
        fields = fields,
        default_fields = default_fields,
        references_by_resource = references_by_resource,
        ...
    )
    return(result)

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
#' @importFrom magrittr %>%
#' @importFrom purrr keep
#'
#' @noRd
select_interaction_datasets <- function(envir){

    envir %>%
    as.list %>%
    `[`(
        c(
            'omnipath',
            'pathwayextra',
            'kinaseextra',
            'ligrecextra'
        )
    ) %>%
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
