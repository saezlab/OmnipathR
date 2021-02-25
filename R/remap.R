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
#' ReMap (http://remap.univ-amu.fr/) is a database of ChIP-Seq experiments.
#' It provides raw and merged peaks and CRMs (cis regulatory motifs) with
#' their associations to regulators (TFs). TF-target relationships can be
#' derived as it is written in Garcia-Alonso et al. 2019: "For ChIP-seq, we
#' downloaded the binding peaks from ReMap and scored the interactions
#' between each TF and each gene according to the distance between the TFBSs
#' and the genes’ transcription start sites. We evaluated different filtering
#' strategies that consisted of selecting only the top-scoring 100, 200, 500,
#' and 1000 target genes for each TF."
#' (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/#s1title).
#' This function returns the top TF-target relationships as used in DoRothEA:
#' https://github.com/saezlab/dorothea/blob/master/inst/scripts/02_chip_seq.R
#' ).
#'
#' @return Data frame with TF-target relationships.
#'
#' @export
#' @seealso \code{\link{remap_tf_target_download}}
remap_dorothea_download <- function(){

    url <- url_parser(url_key = 'omnipath.remap_url')

    version <- omnipath_cache_latest_or_new(url = url)

    if(version$status != CACHE_STATUS$READY){

        download.file(url = url, destfile = version$path, quiet = TRUE)
        omnipath_cache_download_ready(version)

    }

    readRDS(version$path)

}


#' Downloads TF-target interactions from ReMap
#'
#' ReMap (http://remap.univ-amu.fr/) is a database of ChIP-Seq experiments.
#' It provides raw and merged peaks and CRMs (cis regulatory motifs) with
#' their associations to regulators (TFs). TF-target relationships can be
#' derived as it is written in Garcia-Alonso et al. 2019: "For ChIP-seq, we
#' downloaded the binding peaks from ReMap and scored the interactions
#' between each TF and each gene according to the distance between the TFBSs
#' and the genes’ transcription start sites. We evaluated different filtering
#' strategies that consisted of selecting only the top-scoring 100, 200, 500,
#' and 1000 target genes for each TF."
#' (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/#s1title). This
#' function retrieves the full processed TF-target list from the data
#' deposited in https://zenodo.org/record/3713238.
#'
#' @return Data frame with TF-target relationships.
#'
#' @export
remap_tf_target_download <- function(){



}