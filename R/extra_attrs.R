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


#' Converts the extra_attrs column from JSON encoded to list
#'
#' @param data A data frame from the OmniPath web service.
#'
#' @return The input data frame with the extra_attrs column converted
#'     to list.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom jsonlite fromJSON
#' @importFrom purrr map
#' @noRd
deserialize_extra_attrs <- function(data){

    data %>%
    {`if`(
        has_extra_attrs(.),
        mutate(., extra_attrs = map(extra_attrs, fromJSON)),
        .
    )}

}


#' Extra attribute names in an interaction data frame
#'
#' Interaction data frames might have an `extra_attrs` column if this field
#' has been requested in the query by passing the `fields = 'extra_attrs'
#' argument. This column contains resource specific attributes for the
#' interactions. The names of the attributes consist of the name of the
#' resource and the name of the attribute, separated by an underscore.
#' This function returns the names of the extra attributes available in
#' the provided data frame.
#'
#' @return Character: the names of the extra attributes in the data frame.
#'
#' @examples
#' i <- import_omnipath_interactions(fields = 'extra_attrs')
#' extra_attrs(i)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @export
extra_attrs <- function(data){

    data %>%
    {`if`(
        has_extra_attrs(data),
        pull(., extra_attrs) %>% map(names) %>% unlist() %>% unique(),
        character(0)
    )}

}


#' New columns from extra attributes
#'
#' @param data An interaction data frame.
#' @param ... The names of the extra attributes; NSE is supported.
#' @param flatten Logical: unnest the list column even if some records have
#'     multiple values for the attributes; these will yield multiple records
#'     in the resulted data frame.
#' @param keep_empty Logical: if `flatten` is `TRUE`, shall we keep the
#'     records which do not have the attribute?
#'
#' @return Data frame with the new column created; the new column is list
#'     type if one interaction might have multiple values of the attribute,
#'     or character type if
#'
#' @examples
#' i <- import_omnipath_interactions(fields = 'extra_attrs')
#' extra_attrs_to_cols(i, Cellinker_type, Macrophage_type)
#' extra_attrs_to_cols(
#'     i,
#'     Cellinker_type,
#'     Macrophage_type,
#'     flatten = TRUE,
#'     keep_empty = FALSE
#' )
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang enquos !!!
#' @importFrom purrr reduce map_chr
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{has_extra_attrs}}}
#'     \item{\code{\link{with_extra_attrs}}}
#' }
extra_attrs_to_cols <- function(
    data,
    ...,
    flatten = FALSE,
    keep_empty = TRUE
){

    if(!keep_empty){

        data %<>% with_extra_attrs(!!!enquos(...))
    }

    map_chr(enquos(...), .nse_ensure_str) %>%
    reduce(
        .extra_attr_to_col,
        .init = data,
        flatten = flatten
    )

}


#' New column from one extra attribute
#'
#' @param data An interaction data frame.
#' @param attr The name of an extra attribute; NSE is supported.
#' @param flatten Logical: unnest the list column even if some records have
#'     multiple values for the attributes; these will yield multiple records
#'     in the resulted data frame.
#' @param keep_empty Logical: if `flatten` is `TRUE`, shall we keep the
#'     records which do not have the attribute?
#'
#' @return Data frame with the new column created; the new column is list
#'     type if one interaction might have multiple values of the attribute,
#'     or character type if
#'
#' @importFrom magrittr %>% is_less_than
#' @importFrom rlang sym !! := enquo
#' @importFrom purrr map map_int pluck
#' @importFrom dplyr first mutate pull
#' @importFrom tidyr unnest
#' @noRd
.extra_attr_to_col <- function(
    data,
    attr,
    flatten = FALSE,
    keep_empty = TRUE
){

    # NSE vs. R CMD check workaround
    extra_attrs <- NULL

    attr_str <- .nse_ensure_str(!!enquo(attr))
    attr <- as.symbol(attr_str)

    data %>%
    {`if`(
        has_extra_attrs(data),
        mutate(
            .,
            !!attr := map(extra_attrs, pluck, attr_str)
        ) %>%
        {`if`(
            flatten,
            unnest(., !!attr, keep_empty = keep_empty),
            {`if`(
                pull(., !!attr) %>%
                    map_int(length) %>%
                    magrittr::is_less_than(2) %>%
                    all,
                mutate(
                    .,
                    !!attr := map(!!attr, first) %>%
                        null_to_na %>%
                        unlist
                ),
                .
            )}
        )},
        .
    )}

}


#' Tells if an interaction data frame has an extra_attrs column
#'
#' @param data An interaction data frame.
#'
#' @return Logical: TRUE if the data frame has the "extra_attrs" column.
#'
#' @examples
#' i <- import_omnipath_interactions(fields = 'extra_attrs')
#' has_extra_attrs(i)
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{extra_attrs_to_cols}}}
#'     \item{\code{\link{with_extra_attrs}}}
#' }
has_extra_attrs <- function(data){

    'extra_attrs' %in% colnames(data)

}


#' Interaction records having certain extra attributes
#'
#' @param data An interaction data frame.
#' @param ... The name(s) of the extra attributes; NSE is supported.
#'
#' @return The data frame filtered to the records having the extra attribute.
#'
#' @examples
#' i <- import_omnipath_interactions(fields = 'extra_attrs')
#' with_extra_attrs(i, Macrophage_type)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull
#' @importFrom purrr map_lgl map_chr
#' @importFrom rlang enquos
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{has_extra_attrs}}}
#'     \item{\code{\link{extra_attrs_to_cols}}}
#' }
with_extra_attrs <- function(data, ...){

    attrs_str <- map_chr(enquos(...), .nse_ensure_str)

    data %>%
    {`if`(
        has_extra_attrs(.),
        filter(
            .,
            pull(., 'extra_attrs') %>%
            map_lgl(
                function(x){
                    attrs_str %>%
                    intersect(names(x)) %>%
                    length %>%
                    as.logical
                }
            )
        ),
        .
    )}

}
