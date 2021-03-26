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


#' Downloads TF-target interactions from ReMap
#'
#' ReMap (\url{http://remap.univ-amu.fr/}) is a database of ChIP-Seq
#' experiments. It provides raw and merged peaks and CRMs (cis regulatory
#' motifs) with their associations to regulators (TFs). TF-target
#' relationships can be derived as it is written in Garcia-Alonso et al.
#' 2019: "For ChIP-seq, we downloaded the binding peaks from ReMap and scored
#' the interactions between each TF and each gene according to the distance
#' between the TFBSs and the genes’ transcription start sites. We evaluated
#' different filtering strategies that consisted of selecting only the
#' top-scoring 100, 200, 500, and 1000 target genes for each TF."
#' (\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/#s1title}).
#' This function returns the top TF-target relationships as used in DoRothEA:
#' \url{
#' https://github.com/saezlab/dorothea/blob/master/inst/scripts/02_chip_seq.R
#' }).
#'
#' @return Data frame with TF-target relationships.
#'
#' @examples
#' remap_interactions <- remap_dorothea_download()
#' remap_interactions
#' # # A tibble: 136,988 x 2
#' #    tf    target
#' #    <chr> <chr>
#' #  1 ADNP  ABCC1
#' #  2 ADNP  ABCC6
#' #  3 ADNP  ABHD5
#' #  4 ADNP  ABT1
#' #  5 ADNP  AC002066.1
#' # # . with 136,978 more rows
#'
#' @export
#' @importFrom utils download.file
#' @seealso \code{\link{remap_tf_target_download}}
remap_dorothea_download <- function(){

    url <- url_parser(url_key = 'omnipath.remap_dorothea_url')

    version <- omnipath_cache_latest_or_new(url = url)

    if(version$status != CACHE_STATUS$READY){

        download_base(
            url = url,
            fun = download.file,
            destfile = version$path,
            quiet = TRUE
        )
        omnipath_cache_download_ready(version)

    }

    readRDS(version$path)

}


#' Downloads TF-target interactions from ReMap
#'
#' ReMap (\url{http://remap.univ-amu.fr/}) is a database of ChIP-Seq
#' experiments. It provides raw and merged peaks and CRMs (cis regulatory
#' motifs) with their associations to regulators (TFs). TF-target
#' relationships can be derived as it is written in Garcia-Alonso et al.
#' 2019: "For ChIP-seq, we downloaded the binding peaks from ReMap and scored
#' the interactions between each TF and each gene according to the distance
#' between the TFBSs and the genes’ transcription start sites. We evaluated
#' different filtering strategies that consisted of selecting only the
#' top-scoring 100, 200, 500, and 1000 target genes for each TF."
#' (\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/#s1title}).
#' This function retrieves the full processed TF-target list from the data
#' deposited in \url{https://zenodo.org/record/3713238}.
#'
#' @return Data frame with TF-target relationships.
#'
#' @examples
#' \donttest{
#' remap_interactions <- remap_tf_target_download()
#' remap_interactions
#' # # A tibble: 9,546,470 x 4
#' #    source_genesymbol target_genesymbol target_ensembl     score
#' #    <chr>             <chr>             <chr>              <dbl>
#' #  1 ADNP              PTPRS             ENSG00000105426.16  1000
#' #  2 AFF4              PRKCH             ENSG00000027075.14  1000
#' #  3 AHR               CTNND2            ENSG00000169862.18  1000
#' #  4 AR                PDE4D             ENSG00000113448.18  1000
#' #  5 ARID1A            PLEC              ENSG00000178209.14  1000
#' # # . with 9,546,460 more rows
#' }
#'
#'
#' @export
#' @importFrom readr read_tsv cols
#' @importFrom magrittr %T>%
#' @seealso \itemize{
#'     \item{\code{\link{remap_dorothea_download}}}
#'     \item{\code{\link{remap_filtered}}}
#' }
remap_tf_target_download <- function(){

    zenodo_download(
        zenodo_record = 3713238,
        zenodo_fname = 'tf_target_sources.zip',
        path = (
            'tf_target_sources/chip_seq/remap/gene_tf_pairs_genesymbol.txt'
        ),
        reader = read_tsv,
        reader_param = list(
            col_names = c(
                'source_genesymbol',
                'target_genesymbol',
                'target_ensembl',
                'score'
            ),
            col_types = cols(),
            progress = FALSE
        ),
        resource = 'ReMap'
    ) %T>%
    load_success()

}


#' Downloads TF-target interactions from ReMap
#'
#' Downloads the ReMap TF-target interactions as processed by Garcia-Alonso et
#' al. (\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/#s1title})
#' and filters them based on a score threshold, the top targets and whether
#' the TF is included in the TF census (Vaquerizas et al. 2009). The code for
#' filtering is adapted from DoRothEA, written by Christian Holland.
#'
#' @param score Numeric: a minimum score between 0 and 1000, records with
#'     lower scores will be excluded. If NULL no filtering performed.
#' @param top_targets Numeric: the number of top scoring targets for each
#'     TF. Essentially the maximum number of targets per TF. If NULL the
#'     number of targets is not restricted.
#' @param only_known_tfs Logical: whether to exclude TFs which are not in
#'     TF census.
#'
#' @return Data frame with TF-target relationships.
#'
#' @examples
#' \donttest{
#' remap_interactions <- remap_filtered()
#' nrow(remap_interactions)
#' # [1] 145680
#'
#' remap_interactions <- remap_filtered(top_targets = 100)
#' remap_interactions
#' # # A tibble: 30,330 x 2
#' #    source_genesymbol target_genesymbol
#' #    <chr>             <chr>
#' #  1 ADNP              ABCC1
#' #  2 ADNP              ABT1
#' #  3 ADNP              AC006076.1
#' #  4 ADNP              AC007792.1
#' #  5 ADNP              AC011288.2
#' # # . with 30,320 more rows
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter group_by slice_max mutate ungroup arrange
#' @seealso \itemize{
#'     \item{\code{\link{remap_tf_target_download}}}
#'     \item{\code{\link{remap_filtered}}}
#'     \item{\code{\link{tfcensus_download}}}
#' }
remap_filtered <- function(
    score = 100,
    top_targets = 500,
    only_known_tfs = TRUE
){

    # NSE vs. R CMD check workaround
    `HGNC symbol` <- min_score <- NULL

    score_threshold <- score

    remap_tf_target_download() %>%
    {`if`(
        is.null(score_threshold),
        .,
        filter(., score >= score_threshold)
    )} %>%
    {`if`(
        only_known_tfs,
        inner_join(
            .,
            tfcensus_download() %>%
            select(source_genesymbol = `HGNC symbol`),
            by = 'source_genesymbol'
        ),
        .
    )} %>%
    {`if`(
        is.null(top_targets),
        .,
        group_by(., source_genesymbol) %>%
        slice_max(score, n = top_targets) %>%
        mutate(min_score = min(score)) %>%
        filter(score > min_score) %>%
        ungroup() %>%
        select(source_genesymbol, target_genesymbol) %>%
        arrange(source_genesymbol, target_genesymbol)
    )}

}