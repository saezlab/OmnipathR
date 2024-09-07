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


#' Converts the extra_attrs column from JSON encoded to list
#'
#' @param data A data frame from the OmniPath web service.
#' @param ... Passed to `jsonlite::fromJSON`.
#'
#' @return The input data frame with the extra_attrs column converted
#'     to list.
#'
#' @importFrom magrittr %>%
#' @noRd
deserialize_extra_attrs <- function(data, ...){

    data %>%
    deserialize_json_col('extra_attrs', ...)

}


#' Converts a JSON encoded column to list
#'
#' @param data A data frame.
#' @param col Character or symbol: name of the JSON encoded column.
#' @param ... Passed to `jsonlite::fromJSON`.
#'
#' @return The input data frame with the JSON column converted to list.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom dplyr mutate
#' @importFrom jsonlite fromJSON
#' @importFrom purrr map
#' @importFrom rlang !! enquo sym := !!! list2
#' @importFrom logger log_trace
#' @noRd
deserialize_json_col <- function(data, col, ...) {

    col <- .nse_ensure_str(!!enquo(col))

    fromjson_args <-
        list2(...) %>%
        add_defaults(fromJSON, list(simplifyVector = FALSE))

    data %>%
    {`if`(
        has_column(., col),
        identity(.) %T>%
        {log_trace('Converting JSON column `%s` to list.', col)} %>%
        mutate(!!sym(col) := map(!!sym(col), fromJSON, !!!fromjson_args)),
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
#' @param data An interaction data frame, as provided by any of the
#'     \code{\link{omnipath-interactions}} functions.
#'
#' @return Character: the names of the extra attributes in the data frame.
#'
#' @examples
#' i <- omnipath(fields = "extra_attrs")
#' extra_attrs(i)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs_to_cols}}}
#'     \item{\code{\link{has_extra_attrs}}}
#'     \item{\code{\link{with_extra_attrs}}}
#'     \item{\code{\link{filter_extra_attrs}}}
#'     \item{\code{\link{extra_attr_values}}}
#' }
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
#' @param ... The names of the extra attributes; NSE is supported. Custom
#'     column names can be provided as argument names.
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
#' i <- omnipath(fields = "extra_attrs")
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
#' @importFrom purrr reduce2 map_chr
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{has_extra_attrs}}}
#'     \item{\code{\link{with_extra_attrs}}}
#'     \item{\code{\link{filter_extra_attrs}}}
#'     \item{\code{\link{extra_attr_values}}}
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
    reduce2(
        if_null_len0(names(.), .),
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
#' @importFrom magrittr %>% %<>% is_less_than extract
#' @importFrom rlang sym !! := enquo
#' @importFrom purrr map map_int pluck
#' @importFrom dplyr  mutate pull
#' @importFrom tidyr unnest
#' @noRd
.extra_attr_to_col <- function(
    data,
    attr,
    col_name = NULL,
    flatten = FALSE,
    keep_empty = TRUE
){

    # NSE vs. R CMD check workaround
    extra_attrs <- NULL

    attr_str <- .nse_ensure_str(!!enquo(attr))
    attr <- as.symbol(attr_str)
    col_name %<>% if_null_len0(attr_str)
    col <- as.symbol(col_name)

    data %>%
    {`if`(
        has_extra_attrs(data),
        mutate(
            .,
            !!col := map(extra_attrs, pluck, attr_str)
        ) %>%
        {`if`(
            flatten,
            unnest(., !!col, keep_empty = keep_empty),
            {`if`(
                pull(., !!col) %>%
                    map_int(length) %>%
                    magrittr::is_less_than(2) %>%
                    all,
                mutate(
                    .,
                    !!col := map(!!col, ~extract(.x, 1L)) %>%
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
#' i <- omnipath(fields = "extra_attrs")
#' has_extra_attrs(i)
#'
#' @importFrom magrittr %>%
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{extra_attrs_to_cols}}}
#'     \item{\code{\link{with_extra_attrs}}}
#'     \item{\code{\link{filter_extra_attrs}}}
#'     \item{\code{\link{extra_attr_values}}}
#' }
has_extra_attrs <- function(data){

    data %>% has_column('extra_attrs')

}


#' Tells if a column exists in the data frame
#'
#' @noRd
has_column <- function(data, col) {

    col %in% colnames(data)

}

#' Interaction records having certain extra attributes
#'
#' @param data An interaction data frame.
#' @param ... The name(s) of the extra attributes; NSE is supported.
#'
#' @return The data frame filtered to the records having the extra attribute.
#'
#' @examples
#' i <- omnipath(fields = "extra_attrs")
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
#'     \item{\code{\link{filter_extra_attrs}}}
#'     \item{\code{\link{extra_attr_values}}}
#' }
with_extra_attrs <- function(data, ...){

    must_have_extra_attrs(data)

    attrs_str <- map_chr(enquos(...), .nse_ensure_str)

    data %>%
    filter(
        pull(., 'extra_attrs') %>%
        map_lgl(
            function(x){
                attrs_str %>%
                intersect(names(x)) %>%
                length %>%
                as.logical
            }
        )
    )

}


#' Filter interactions by extra attribute values
#'
#' @param data An interaction data frame with \emph{extra_attrs} column.
#' @param ... Extra attribute names and values. The contents of the extra
#'     attribute \emph{name} for each record will be checked against the
#'     values provided. The check by default is a set intersection: if any
#'     element is common between the user provided values and the values of
#'     the extra attribute for the record, the record will be matched.
#'     Alternatively, any value can be a custom function which accepts
#'     the value of the extra attribute and returns a single logical
#'     value. Finally, if the extra attribute name starts with a dot,
#'     the result of the check will be negated.
#' @param na_ok Logical: keep the records which do not have the extra
#'     attribute. Typically these are the records which are not from the
#'     resource providing the extra attribute.
#'
#' @return The input data frame with records removed according to the
#'     filtering criteria.
#'
#' @examples
#' cl <- post_translational(
#'     resources = "Cellinker",
#'     fields = "extra_attrs"
#' )
#' # Only cell adhesion interactions from Cellinker
#' filter_extra_attrs(cl, Cellinker_type = "Cell adhesion")
#'
#' op <- omnipath(fields = "extra_attrs")
#' # Any mechanism except phosphorylation
#' filter_extra_attrs(op, .SIGNOR_mechanism = "phosphorylation")
#'
#' @importFrom magrittr %>% %<>% not
#' @importFrom purrr reduce2 map_lgl compose pluck
#' @importFrom stringr str_replace str_starts
#' @importFrom dplyr first
#' @importFrom rlang list2
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs}}}
#'     \item{\code{\link{has_extra_attrs}}}
#'     \item{\code{\link{extra_attrs_to_cols}}}
#'     \item{\code{\link{with_extra_attrs}}}
#'     \item{\code{\link{extra_attr_values}}}
#' }
filter_extra_attrs <- function(data, ..., na_ok = TRUE){

    must_have_extra_attrs(data)

    list2(...) %>%
    reduce2(

        names(.),

        function(d, val, key){

            negate <- `if`(str_starts(key, '\\.'), not, identity)
            key %<>% str_replace('^\\.', '')

            check <-
                `if`(
                    is.function(val),
                    val,
                    compose(
                        ~ .x > 0,
                        length,
                        ~ intersect(.x, val)
                    )
                ) %>%
                compose(negate, .) %>%
                {`if`(na_ok, ~ first(is.na(.x)) || .(.x), .)}

            d %>%
            filter(
                map_lgl(
                    extra_attrs,
                    compose(
                        check,
                        ~ pluck(.x, key)
                    )
                )
            )

        },

        .init = data

    )

}


#' Possible values of an extra attribute
#'
#' Extracts all unique values of an extra attribute occuring in this data
#' frame.
#'
#' @details
#' Note, at the end we unlist the result, which means it works well
#' for attributes which are atomic vectors but gives not so useful result
#' if the attribute values are more complex objects. As the time of
#' writing this, no such complex extra attribute exist in OmniPath.
#'
#' @param data An interaction data frame with \emph{extra_attrs} column.
#' @param key The name of an extra attribute.
#'
#' @return A vector, most likely character, with the unique values of the
#'     extra attribute occuring in the data frame.
#'
#' @examples
#' op <- omnipath(fields = "extra_attrs")
#' extra_attr_values(op, SIGNOR_mechanism)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom purrr map pluck discard compose
#' @importFrom rlang !! enquo
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{extra_attrs_to_cols}}}
#'     \item{\code{\link{has_extra_attrs}}}
#'     \item{\code{\link{with_extra_attrs}}}
#'     \item{\code{\link{filter_extra_attrs}}}
#'     \item{\code{\link{extra_attrs}}}
#' }
extra_attr_values <- function(data, key){

    # NSE vs. R CMD check workaround
    extra_attrs <- NULL

    key <- .nse_ensure_str(!!enquo(key))

    data %>%
    must_have_extra_attrs %>%
    pull(extra_attrs) %>%
    map(pluck, key) %>%
    discard(compose(all, is.na)) %>%
    unlist %>%
    unique

}


#' Raise an error if the data frame has no extra_attrs column
#'
#' @param data A data frame.
#'
#' @return Invisibly returns the input data frame.
#'
#' @noRd
must_have_extra_attrs <- function(data){

    if(!has_extra_attrs(data)){

        stop('Data frame has no "extra_attrs" column.')

    }

    invisible(data)

}
