#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
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
#
#  This file is from https://github.com/klmr/decorator
#  Author: Kondrad Rudolph
#  License: Apache-2.0
#


#' Decorator for trying UniProt subdomains
#'
#' This has any relevance only in rare cases with OS networking issues.
#' 10/2022: deprecated since the uniprot api change, should be removed
#'
#' @importFrom logger log_trace
#' @noRd
uniprot_domains <- decorator %@% function(FUN){

    function(...){

        for(subd in c('legacy')){

            result <- tryCatch(
                FUN(..., .subdomain = subd),
                error = identity
            )

            if(!inherits(result, 'error')){

                break

            }else{

                log_trace(
                    'Failed download attempt: `%s`.',
                    conditionMessage(result)
                )

            }

        }

        if(inherits(result, 'error')){

            stop(conditionMessage(result))

        }

        return(result)

    }

}
