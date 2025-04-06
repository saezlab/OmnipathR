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


#' Acquire a cookie if necessary
#'
#' @param url Character. URL to download to get the cookie.
#' @param http_headers List: HTTP headers.
#' @param default_headers Logical: if TRUE, add the default headers.
#' @param post List: HTTP POST parameters.
#' @param payload Data to send as payload.
#' @param init_post List: HTTP POST parameters for ``init_url``.
#' @param init_payload Data to send as payload with ``init_url``.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv cols
#' @importFrom curl new_handle handle_setopt curl_fetch_memory handle_cookies
#' @importFrom logger log_trace
#' @importFrom purrr pluck
#' @noRd
#'
#' @return A list with cache file path, cookies and response headers.
omnipath_init_cookies <- function(
        url,
        http_headers = list(),
        default_headers = FALSE,
        post = NULL,
        payload = NULL,
        curlopt = list()
    ){

    cookies <- NULL

    while (TRUE) {

        curlopt$followlocation <- FALSE

        req <-
            url %>%
            omnipath_httr2_req(
                post = post,
                payload = payload,
                http_headers = http_headers,
                default_headers = default_headers,
                curlopt = curlopt
            )

        resp <-
            req %>%
            omnipath_httr2_perform()

        cookies %<>% append(
            resp %>%
            omnipath_resp_headers('set-cookie') %T>%
            log_trace('Cookies: %s', .)
        )

        if (resp %>% is_http_redirect) {

            url <- resp %>% redirect_location

        } else {

            break

        }

    }

    cookies

}
