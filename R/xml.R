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


#' Read an XML in chunks and parse by event driven parser
#'
#' @param con A connection to an XML file.
#' @param record The name of the record that must be complete within chunks;
#'     incomplete records will be added to the next chunk.
#' @param parser Callable: a function for parsing each chunk. It should accept
#'     the XML chunk as character vector as its first argument, and ... will be
#'     passed as additional arguments.
#' @param post_parse A function to extract and post process the parsed data
#'     structure.
#' @param chunk_size Number of lines to read in each chunk.
#' @param header_lines Number of lines to read in the header. These lines will
#'     be added in front of each chunk.
#' @param closing Character: an XML fragment to be appended to the end of each
#'     chunk: typically it closes the root node. If `NULL` (the default), we
#'     will attempt to guess it from the root tag found in the beginning of the
#'     XML.
#' @param ... All arguments for the `parser`.
#'
#' @return A list of the outputs of `post_parse`, one item for each chunk.
#'
#' @importFrom magrittr %>% %<>% extract2 add
#' @importFrom stringi stri_locate_last_fixed
#' @importFrom progress progress_bar
#' @importFrom purrr reduce map2
#' @noRd
parse_in_chunks <- function(
        con,
        record,
        parser = xml_event_parse,
        post_parse = identity,
        chunk_size = 2e7L,
        header_lines = NULL,
        closing = NULL,
        ...
    ) {

    on.exit(close(con))

    size <- con %>% file_size
    log_trace(
        'Processing XML from `%s` in chunks of %i characters; total size: %i.',
        summary(con)$description,
        chunk_size,
        size
    )
    endtag <- sprintf('</%s>', record)
    incomplete <- ''
    header <- NULL
    result <- list()

    pb <- progress_bar$new(
        total = size,
        format = paste0(
            '  Processing XML: ',
            '[:bar] :percent, :rate, eta: :eta'
        )
    )

    repeat {

        chunk <- con %>% readChar(n = chunk_size)

        if (length(chunk) == 0L) break

        if(is.null(header)) {
            parts <- chunk %>% extract_header(header_lines, record)
            header <- parts$header
            chunk <- parts$chunk
            closing %<>% if_null(sprintf('</%s>', parts$root))
            log_trace('XML closing part: %s', closing)
        }

        chunk %<>% paste(incomplete, ., sep = '')

        i_endtag <-
            chunk %>%
            stri_locate_last_fixed(endtag) %>%
            extract(1L, 2L) %>%
            add(1L)

        if(i_endtag == -1L) {

            incomplete <- chunk

        } else {

            incomplete <- chunk %>% substr(i_endtag, nchar(.))

            result %<>%
                c(
                    chunk %>%
                    substr(1L, i_endtag - 1L) %>%
                    paste(header, ., sep = '') %>%
                    paste(closing, sep = '\n') %>%
                    parser(...) %>%
                    post_parse %>%
                    list
                )
            }

            if(!pb$finished) {
                pb$tick(nchar(chunk))
            }

        }

    result %<>% reduce(~map2(.x, .y, ~c(.x, .y)))

    return(result)

}


#' Wrapper around xmlEventParse
#'
#' @importFrom XML xmlEventParse
#' @noRd
xml_event_parse <- function(file, ...) {

    file %>%
    xmlEventParse(
        asText = TRUE,
        error = xml_error_report,
        ...
    )

}


#' Handle XML parsing errors from the XML package
#'
#' @noRd
xml_error_report <- function(msg, code = -1, domain = 'unknown', ...) {

    log_warn('XML error %i (in %s): %s', code, domain, msg)

}


#' Attempt to separate the header
#'
#' Every XML has a part at its beginning which contains the xml opening tag,
#' the opening of the root node, and potentially other nodes. Here we look for
#' the opening of the first record tag, and consider everything before that to
#' be a header.
#'
#' @return A list with the header, the rest of the chunk, and the name of the
#'     root tag.
#'
#' @importFrom stringi stri_locate_first_fixed stri_locate_all_fixed
#' @importFrom magrittr %>% extract extract2 subtract
#' @noRd
extract_header <- function(xml, header_lines = NULL, record = NULL) {

    log_trace(
        'Extracting XML header; header lines: %s; record: %s.',
        header_lines %>% if_null('NULL'),
        record %>% if_null('NULL')
    )

    `if`(
        is.null(header_lines) && !is.null(record),
        xml %>%
        stri_locate_first_fixed(sprintf('<%s>', record)) %>%
        extract(1L, 1L) %>%
        subtract(1L) %>%
        {list(
            header = substr(xml, 1L, .),
            chunk = substr(xml, . + 1L, nchar(xml))
        )},
        NA
    ) %>%
    {`if`(
        !is.list(.),
        xml %>%
        stri_locate_all_fixed('\n') %>%
        extract2(1L) %>%
        extract(header_lines %>% if_null(2L), 1L) %>%
        substr(xml, 1L, .) %>%
        list(
            header = .,
            chunk = substr(xml, nchar(.) + 1L, nchar(xml))
        ),
        .
    )} %>%
    c(list(root = root_from_header(.$header)))

}


#' Attempt to guess the root tag from the header
#'
#' If we process a huge XML iteratively, we have to add the closing tag of the
#' root node to the end of each chunk. Here we read the header part which is
#' supposed to contain the opening tag of the root node, and extract the name
#' of this tag from the error message: as this tag is not closed, we expect an
#' error about "Premature end of data in tag ...".
#'
#' @return Character: name of the first unclosed tag within the header; NA if
#'     no such tag was found.
#'
#' @importFrom xml2 read_xml
#' @importFrom purrr safely
#' @importFrom magrittr %>% extract2
#' @importFrom stringr str_match
#' @noRd
root_from_header <- function(header) {

    safely(read_xml)(header) %>%
    extract2('error') %>%
    as.character %>%
    str_match('end of data in tag (\\w+)') %>%
    extract2(1L, 2L)

}
