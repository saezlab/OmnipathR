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


.onLoad <- function(libname, pkgname){

    omnipath_init_config()
    patch_logger()
    omnipath_init_log(pkgname = pkgname)

    if(Sys.info()['user'] == 'biocbuild'){

        omnipath_set_console_loglevel('trace')

    }

    omnipath_init_cache()
    omnipath_init_db(pkgname)
    .load_magic_bytes(pkgname)
    .load_urls(pkgname)

    logger::log_info('Welcome to OmnipathR!')

}