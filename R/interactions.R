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


#' Molecular interactions from OmniPath
#'
#' The functions listed here all download pairwise, causal molecular
#' interactions from the \url{https://omnipathdb.org/interactions} endpoint of
#' the OmniPath web service. They are different only in the type of
#' interactions and the kind of resources and data they have been compiled
#' from. A complete list of these functions is available below, these cover
#' the interaction datasets and types  currently available in OmniPath:
#'
#' **Post-translational (protein-protein, PPI) interactions**
#'
#' \itemize{
#'     \item{\code{omnipath}: the OmniPath data as defined in the 2016 paper,
#'         an arbitrary optimum between coverage and quality. This dataset
#'         contains almost entirely causal (stimulatory or inhibitory; i.e.
#'         activity flow , according to the SBGN standard), physical
#'         interactions between pairs of proteins, curated by experts
#'         from the literature.}
#'     \item{\code{pathwayextra}: activity flow interactions without literature
#'         references.}
#'     \item{\code{kinaseextra}: enzyme-substrate interactions without
#'         literature references.}
#'     \item{\code{ligrecextra}: ligand-receptor interactions without
#'         literature references.}
#'     \item{\code{post_translational}: all post-translational
#'         (protein-protein, PPI) interactions; this is the combination of the
#'         *omnipath*, *pathwayextra*, *kinaseextra* and *ligrecextra*
#'         datasets.}
#' }
#'
#' **TF-target (gene regulatory, GRN) interactions**
#'
#' \itemize{
#'     \item{\code{collectri}: transcription factor (TF)-target
#'         interactions from CollecTRI.}
#'     \item{\code{dorothea}: transcription factor (TF)-target
#'         interactions from DoRothEA}
#'     \item{\code{tf_target}: transcription factor
#'         (TF)-target interactions from other resources}
#'     \item{\code{transcriptional}: all transcription factor
#'         (TF)-target interactions; this is the combination of the
#'         *collectri*, *dorothea* and *tf_target* datasets.}
#' }
#'
#' **Post-transcriptional (miRNA-target) and other RNA related interactions**
#'
#' In these datasets we intend to collect the literature curated resources,
#' hence we don't include some of the most well known large databases if those
#' are based on predictions or high-throughput assays.
#'
#' \itemize{
#'     \item{\code{mirna_target}: miRNA-mRNA interactions}
#'     \item{\code{tf_mirna}: TF-miRNA interactions}
#'     \item{\code{lncrna_mrna}: lncRNA-mRNA interactions}
#' }
#'
#' **Other interaction access functions**
#'
#' \itemize{
#'     \item{\code{small_molecule}: interactions between small molecules and
#'         proteins. Currently this is a small, experimental dataset that
#'         includes drug-target, ligand-receptor, enzyme-metabolite and other
#'         interactions. In the future this will be largely expanded and
#'         divided into multiple datasets.}
#'     \item{\code{all_interactions}: all the interaction datasets combined.}
#' }
#'
#' @examples
#' op <- omnipath(resources = c("CA1", "SIGNOR", "SignaLink3"))
#' op
#'
#' @seealso \itemize{
#'     \item{\code{\link{interaction_resources}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#'     \item{\code{\link{annotated_network}}}
#' }
#'
#' @name omnipath-interactions
#' @md
NULL


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
#' interactions = omnipath_interactions(
#'     resources = "SignaLink3",
#'     organism = 9606
#' )
#'
#' @importFrom magrittr %<>%
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_omnipath_interactions
omnipath_interactions <- function(...){

    args <- omnipath_args(list(...), query_type = 'interactions')
    defaults <- list(datasets = 'omnipath')
    args %<>% modifyList(defaults, .)

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{omnipath_interactions}.
#' @export
#'
#' @noRd
import_omnipath_interactions <- function(...){
    .Deprecated('import_omnipath_interactions')
    omnipath_interactions(...)
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
#'     \code{\link{omnipath_interactions}}.
#'
#' @examples
#' pathways <- omnipath()
#' pathways
#'
#' @seealso \itemize{
#'     \item{\code{\link{omnipath_interactions}}}
#'     \item{\code{\link{post_translational}}}
#'     \item{\code{\link{interaction_resources}}}
#'     \item{\code{\link{all_interactions}}}
#'     \item{\code{\link{interaction_graph}}}
#'     \item{\code{\link{print_interactions}}}
#' }
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
omnipath <- function(...){

    args <- omnipath_args(list(...), datasets = 'omnipath')

    exec(omnipath_interactions, !!!args)

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
#'     pathwayextra(
#'         resources = c("BioGRID", "IntAct"),
#'         organism = 9606
#'     )
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_pathwayextra_interactions
pathwayextra <- function(...){

    args <- omnipath_args(list(...), datasets = 'pathwayextra')

    exec(omnipath_interactions, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{pathwayextra}.
#' @export
#'
#' @noRd
import_pathwayextra_interactions <- function(...){
    .Deprecated('import_pathwayextra_interactions')
    pathwayextra(...)
}


#' Interactions from the `kinaseextra` dataset of OmniPath
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
#' kinase_substrate <-
#'    kinaseextra(
#'        resources = c('PhosphoPoint', 'PhosphoSite'),
#'        organism = 9606
#'    )
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_kinaseextra_interactions
kinaseextra <- function(...){

    args <- omnipath_args(list(...), datasets = 'kinaseextra')

    exec(omnipath_interactions, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{kinaseextra}.
#' @export
#'
#' @noRd
import_kinaseextra_interactions <- function(...){
    .Deprecated('import_kinaseextra_interactions')
    kinaseextra(...)
}


#' Interactions from the `ligrecextra` dataset of OmniPath
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
#' ligand_receptor <- ligrecextra(
#'     resources = c('HPRD', 'Guide2Pharma'),
#'     organism = 9606
#' )
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_ligrecextra_interactions
ligrecextra <- function(...){

    args <- omnipath_args(list(...), datasets = 'ligrecextra')

    exec(omnipath_interactions, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{ligrecextra}.
#' @export
#'
#' @noRd
import_ligrecextra_interactions <- function(...){
    .Deprecated('import_ligrecextra_interactions')
    ligrecextra(...)
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
#' interactions <- post_translational(resources = "BioGRID")
#'
#' @importFrom rlang %||% exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_post_translational_interactions
post_translational <- function(...){


    args <- omnipath_args(list(...), query_type = 'interactions')
    args$datasets %<>% {. %||% PPI_DATASETS}

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{post_translational}.
#' @export
#'
#' @noRd
import_post_translational_interactions <- function(...){
    .Deprecated('import_post_translational_interactions')
    post_translational(...)
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
#' @importFrom rlang exec !!!
#' @export
#' @aliases import_dorothea_interactions
#' @rdname omnipath-interactions
dorothea <- function(dorothea_levels = c('A', 'B'), ...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'dorothea'
    )

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{dorothea}.
#' @export
#'
#' @noRd
import_dorothea_interactions <- function(...){
    .Deprecated('import_dorothea_interactions')
    dorothea(...)
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
#' interactions <- tf_target(resources = c("DoRothEA", "SIGNOR"))
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_tf_target_interactions
tf_target <- function(...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'tf_target'
    )

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{tf_target}.
#' @export
#'
#' @noRd
import_tf_target_interactions <- function(...){
    .Deprecated('import_tf_target_interactions')
    tf_target(...)
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
#' grn <- transcriptional(resources = c("PAZAR", "ORegAnno", "DoRothEA"))
#' grn
#'
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang %||% exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_transcriptional_interactions
transcriptional <- function(dorothea_levels = c('A', 'B'), ...){

    args <- omnipath_args(list(...), query_type = 'interactions')
    args$datasets %<>% {. %||% GRN_DATASETS} %>% intersect(GRN_DATASETS)

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{transcriptional}.
#' @export
#'
#' @noRd
import_transcriptional_interactions <- function(...){
    .Deprecated('import_transcriptional_interactions')
    transcriptional(...)
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
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
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
#' interactions <- mirna_target( resources = c("miRTarBase", "miRecords"))
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_mirnatarget_interactions
mirna_target <- function(...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'mirnatarget'
    )

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{mirna_target}.
#' @export
#'
#' @noRd
import_mirnatarget_interactions <- function(...){
    .Deprecated('import_mirnatarget_interactions')
    mirna_target(...)
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
#' interactions <- tf_mirna(resources = "TransmiR")
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_tf_mirna_interactions
tf_mirna <- function(...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'tf_mirna'
    )

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{tf_mirna}.
#' @export
#'
#' @noRd
import_tf_mirna_interactions <- function(...){
    .Deprecated('import_tf_mirna_interactions')
    tf_mirna(...)
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
#' interactions <- lncrna_mrna(resources = c("ncRDeathDB"))
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_lncrna_mrna_interactions
lncrna_mrna <- function(...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'lncrna_mrna'
    )

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{lncrna_mrna}.
#' @export
#'
#' @noRd
import_lncrna_mrna_interactions <- function(...){
    .Deprecated('import_lncrna_mrna_interactions')
    lncrna_mrna(...)
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
#' interactions <- small_molecule(sources = "ASPIRIN")
#' # The prostaglandin synthases:
#' interactions
#'
#' @importFrom rlang exec !!!
#' @export
#' @rdname omnipath-interactions
#' @aliases import_small_molecule_protein_interactions
small_molecule <- function(...){

    args <- omnipath_args(
        list(...),
        query_type = 'interactions',
        datasets = 'small_molecule'
    )

    exec(omnipath_query, !!!args)

}


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{small_molecule}.
#' @export
#'
#' @noRd
import_small_molecule_protein_interactions <- function(...){
    .Deprecated('import_small_molecule_protein_interactions')
    small_molecule(...)
}


#' Imports all interaction datasets available in OmniPath
#'
#' @return A dataframe containing all the datasets in the interactions query
#'
#' @param dorothea_levels The confidence levels of the dorothea
#' interactions (TF-target) which range from A to D. Set to A and B by
#' default.
#' @param types Character: interaction types, such as "transcriptional",
#'     "post_transcriptional", "post_translational", etc.
#' @param fields Character: additional fields (columns) to be included in the
#'     result. For a list of available fields, see \code{\link{query_info}}.
#' @param exclude Character: names of datasets or resource to be excluded from
#'     the result. By deafult, the records supported by only these resources or
#'     datasets will be removed from the output. If \code{strict_evidences =
#'     TRUE}, the resource, reference and causality information in the data
#'     frame will be reconstructed to remove all information coming from the
#'     excluded resources.
#' @param ... Arguments passed to \code{\link{omnipath_query}}.
#' @inheritDotParams omnipath_query -datasets -query_type
#'
#' @examples
#' interactions <- all_interactions(
#'     resources = c("HPRD", "BioGRID"),
#'     organism = 9606
#' )
#'
#' @importFrom rlang exec !!!
#' @importFrom magrittr %<>% %>% extract2
#' @export
#' @rdname omnipath-interactions
#' @aliases import_all_interactions
all_interactions <- function(
    dorothea_levels = c('A', 'B'),
    types = NULL,
    fields = NULL,
    exclude = NULL,
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


# Aliases (old names) to be Deprecated
#' @rdname omnipath-interactions
#' @param ... Passed to \code{all_interactions}.
#' @export
#'
#' @noRd
import_all_interactions <- function(...){
    .Deprecated('import_all_interactions')
    all_interactions(...)
}


#' Interaction resources available in Omnipath
#'
#' Names of the resources available in \url{https://omnipathdb.org/interactions}.
#'
#' @param dataset a dataset within the interactions query type. Currently
#' available datasets are `omnipath`, `kinaseextra`, `pathwayextra`,
#' `ligrecextra`, `collectri`, `dorothea`, `tf_target`, `tf_mirna`,
#' `mirnatarget`, `lncrna_mrna` and `small_molecule_protein`.
#'
#' @return Character: names of the interaction resources.
#'
#' @examples
#' interaction_resources()
#'
#' @seealso \itemize{
#'     \item{\code{\link{resources}}}
#'     \item{\code{\link{omnipath}}}
#'     \item{\code{\link{pathwayextra}}}
#'     \item{\code{\link{kinaseextra}}}
#'     \item{\code{\link{ligrecextra}}}
#'     \item{\code{\link{post_translational}}}
#'     \item{\code{\link{dorothea}}}
#'     \item{\code{\link{collectri}}}
#'     \item{\code{\link{tf_target}}}
#'     \item{\code{\link{transcriptional}}}
#'     \item{\code{\link{mirna_target}}}
#'     \item{\code{\link{tf_mirna}}}
#'     \item{\code{\link{small_molecule}}}
#'     \item{\code{\link{all_interactions}}}
#' }
#' @export
#' @aliases get_interaction_resources
interaction_resources <- function(dataset = NULL){

    return(resources(query_type = 'interactions', datasets = dataset))

}


# Aliases (old names) to be Deprecated
#' @rdname complexes
#' @param ... Passed to \code{interaction_resources}.
#' @export
#'
#' @noRd
get_interaction_resources <- function(...){
    .Deprecated('get_interaction_resources')
    interaction_resources(...)
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
