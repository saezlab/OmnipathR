#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2020
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


#' Default values for the package options
#'
#' These options describe the default settings for OmnipathR so you do not
#' need to pass these parameters at each function call.
#' Currently the only option useful for the public web service at
#' omnipathdb.org is ``omnipath.license``. If you are a for-profit
#' user set it to ``'commercial'`` to make sure all the data you download
#' from OmniPath is legally allowed for commercial use. Otherwise just leave
#' it as it is: ``'academic'``.
#' If you don't use omnipathdb.org but within your organization you deployed
#' your own pypath server and want to share data whith a limited availability
#' to outside users, you may want to use a password. For this you can use
#' the ``omnipath.password`` option.
#' Also if you want the R package to work from another pypath server instead
#' of omnipathdb.org, you can change the option ``omnipath.url``.
.omnipath_options_defaults <- list(
    omnipath.url = 'https://omnipathdb.org/',
    omnipath.license = 'academic',
    omnipath.password = NULL,
    omnipath.print_urls = FALSE,
    omnipath.pathwaycommons_url = paste0(
        'https://www.pathwaycommons.org/archives/PC2/v12/',
        'PathwayCommons12.All.hgnc.sif.gz'
    ),
    omnipath.harmonizome_url = paste0(
        'https://maayanlab.cloud/static/hdfs/harmonizome/data/',
        '%s/gene_attribute_edges.txt.gz'
    ),
    omnipath.vinayagam_url = paste0(
        'https://stke.sciencemag.org/content/sigtrans/suppl/2011/09/01/',
        '4.189.rs8.DC1/4_rs8_Tables_S1_S2_and_S6.zip'
    ),
    omnipath.cpdb_url = paste0(
        'http://cpdb.molgen.mpg.de/download/',
        'ConsensusPathDB_human_PPI.gz'
    )
)

.onLoad <- function(libname, pkgname){

    do.call(options, .omnipath_options_defaults)

}