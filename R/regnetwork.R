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


#' Interactions from RegNetwork
#'
#' Downloads transcriptional and post-transcriptional regulatory interactions
#' from the RegNetwork database (\url{http://www.regnetworkweb.org/}). The
#' information about effect signs (stimulation or inhibition), provided by
#' \code{\link{regnetwork_directions}} are included in the result.
#'
#' @param organism Character: either human or mouse.
#'
#' @return Data frame with interactions.
#'
#' @examples
#' regn_interactions <- regnetwork_download()
#' regn_interactions
#' # # A tibble: 372,778 x 7
#' #    source_genesymb. source_entrez target_genesymb. target_entrez
#' #    <chr>            <chr>         <chr>            <chr>
#' #  1 USF1             7391          S100A6           6277
#' #  2 USF1             7391          DUSP1            1843
#' #  3 USF1             7391          C4A              720
#' #  4 USF1             7391          ABCA1            19
#' #  5 TP53             7157          TP73             7161
#' # # . with 372,768 more rows, and 3 more variables: effect <dbl>,
#' # #   source_type <chr>, target_type <chr>
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom dplyr left_join mutate select
#' @importFrom readr read_tsv cols col_character
#' @importFrom tibble tibble
regnetwork_download <- function(organism = 'human'){

    # NSE vs. R CMD check workaround
    source_entrez <- target_entrez <- NULL

    columns <-
        c(
            'source_genesymbol',
            'source_entrez',
            'target_genesymbol',
            'target_entrez'
        )

    path <-
        tryCatch(
            {download_to_cache(
                url_key = 'regnetwork',
                url_param = list(organism)
            )},
            error = regnetwork_error
        )

    if(is.null(path)){

        return(empty_char_tibble(columns))

    }

    con <- unz(path, sprintf('%s.source', organism), open = 'rb')

    on.exit(close(con))

    con %>%
    read_tsv(
        col_names = columns,
        col_types = cols(
            .default = col_character()
        ),
        progress = FALSE
    ) %>%
    left_join(
        regnetwork_directions(organism = organism) %>%
        select(-source_entrez, -target_entrez),
        by = c('source_genesymbol', 'target_genesymbol')
    ) %>%
    mutate(
        source_type = mirna_or_protein(source_entrez),
        target_type = mirna_or_protein(target_entrez)
    ) %>%
    copy_source_attrs(path, resource = 'RegNetwork') %T>%
    load_success()

}


#' Transcription factor effects from RegNetwork
#'
#' @param organism Character: either human or mouse.
#'
#' @return A data frame (tibble) of TF-target interactions with effect signs.
#'
#' @examples
#' regn_dir <- regnetwork_directions()
#' regn_dir
#' # # A tibble: 3,954 x 5
#' #    source_genesymb. source_entrez target_genesymb. target_entrez
#' #    <chr>            <chr>         <chr>            <chr>
#' #  1 AHR              196           CDKN1B           1027
#' #  2 APLNR            187           PIK3C3           5289
#' #  3 APLNR            187           PIK3R4           30849
#' #  4 AR               367           KLK3             354
#' #  5 ARNT             405           ALDOA            226
#' # # . with 3,944 more rows, and 1 more variable: effect <dbl>
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_delim cols col_character
#' @importFrom dplyr mutate recode
#'
#' @export
regnetwork_directions <- function(organism = 'human'){

    # NSE vs. R CMD check workaround
    effect <- NULL

    columns <-
        c(
            'source_genesymbol',
            'source_entrez',
            'target_genesymbol',
            'target_entrez',
            'effect'
        )

    path <-
        tryCatch(
            {download_to_cache(
                url_key = 'regnetwork',
                url_param = list('RegulatoryDirections')
            )},
            error = regnetwork_error
        )

    if(is.null(path)){
        return(empty_char_tibble(columns))
    }

    con <-
        path %>%
        unz(sprintf('new_kegg.%s.reg.direction.txt', organism), open = 'rb')

    on.exit(close(con))

    con %>%
    read_delim(
        delim = ' ',
        skip = 1,
        col_names = columns,
        col_types = cols(
            .default = col_character()
        ),
        progress = FALSE
    ) %>%
    mutate(
        effect = recode(
            effect,
            `--|`    = -1,
            `-->`    =  1,
            .default =  0
        )
    ) %>%
    copy_source_attrs(path, resource = 'RegNetwork') %T>%
    load_success()

}


#' For a vector of miRBase IDs and various other IDs mixed, returns a vector
#' of molecule types "mirna" and "protein"
#'
#' @noRd
mirna_or_protein <- function(ids){

    ifelse(
        grepl('^MI(?:MAT)?[0-9]', ids),
        'mirna',
        'protein'
    )

}


#' Handling the RegNetwork SSL error
#'
#' This error happens because the RegNetwork server is not configured properly
#' to send its intermediate certificate to clients. The missing intermediate
#' certificate is the "InCommon RSA Server CA" one.
#'
#' @importFrom stringr str_detect
#' @importFrom logger log_error
#' @noRd
regnetwork_error <- function(e){

    cond <- conditionMessage(e)

    if(str_detect(cond, 'SSL')){

        msg <-
            paste0(
                'The server at regnetworkweb.org missed to send ',
                'the InCommon RSA Server CA certificate that is ',
                'required for the verification of the server`s ',
                'identity and establishment of a TLS connection. ',
                'You can work around this issue by adding the ',
                'InCommon RSA Server CA certificate to the ',
                'trusted certificate store which is used by the ',
                'R `curl` package. Now we raise only a warning ',
                'and return an empty data frame.'
            ) %>%
            sprintf('%s The original error: %s', ., cond)

        log_error(msg)
        warning(msg)
        return(NULL)

    }else{

        stop(e)

    }

}


#' An empty tibble with character type columns
#'
#' @param header Character vector: column names.
#'
#' @importFrom rlang set_names exec !!!
#' @importFrom tibble tibble
#' @noRd
empty_char_tibble <- function(header){

    character() %>%
    list() %>%
    rep(length(header)) %>%
    set_names(header) %>%
    exec(tibble, !!!.)

}
