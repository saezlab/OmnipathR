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

EVIDENCES_KEYS <- c('positive', 'negative', 'directed', 'undirected')


#' Converts the evidences column from JSON encoded to list
#'
#' @param data A data frame from the OmniPath web service.
#' @param ... Passed to `jsonlite::fromJSON`.
#'
#' @return The input data frame with the "evidences" column converted
#'     to list.
#'
#' @importFrom magrittr %>%
#' @noRd
deserialize_evidences <- function(data, ...){

    data %>%
    deserialize_json_col('evidences', ...)

}


#' Tells if a data frame has a column "evidences"
#'
#' @importFrom magrittr %>% or
#' @noRd
has_evidences <- function(data, wide_ok = FALSE){

    data %>%
    has_column('evidences') %>%
    or(wide_ok && has_evidences_wide(data))

}


#' Tells if a data frame has evidence columns by direction and sign
#'
#' @importFrom magrittr %>% is_in
#' @noRd
has_evidences_wide <- function(data) {

    EVIDENCES_KEYS %>% is_in(colnames(data)) %>% all

}


#' Check for "evidences" column, throw an error if it is missing
#'
#' @importFrom magrittr %>% extract
#' @noRd
must_have_evidences <- function(data, wide_ok = FALSE, env = parent.frame()){

    if(!has_evidences(data, wide_ok = wide_ok)) {

        parent_call <- sys.call(-1L)

        m <-
            paste(
                '`%s` can be called only on data',
                'frames with `evidences` column.'
            ) %>%
            sprintf(parent_call %>% as.character %>% extract(1L))

        log_error(m)
        do.call(stop, list(m), envir = env)

    }

}


#' Separate evidences by direction and effect sign
#'
#' @param data An interaction data frame with "evidences" column.
#' @param longer Logical: If TRUE, the "evidences" column is split into rows.
#' @param .keep Logical: keep the "evidences" column. When unnesting to longer
#'     data frame, the "evidences" column will contain the unnested evidences,
#'     while the original column will be retained under the "all_evidences"
#'     name (if `.keep = TRUE`).
#'
#' @return The data frame with new columns or new rows by direction and sign.
#'
#' @examples
#' \dontrun{
#' omnipath <- import_omnipath_interactions(fields = 'evidences')
#' omnipath <- unnest_evidences(omnipath)
#' colnames(omnipath)
#' }
#'
#' @importFrom magrittr %>% extract
#' @importFrom dplyr mutate rename select
#' @importFrom purrr map
#' @importFrom tidyr unnest_longer unnest_wider
#' @importFrom rlang exec !!! !! := sym
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{only_from}}}
#'     \item{\code{\link{filter_evidences}}}
#'     \item{\code{\link{from_evidences}}}
#' }
unnest_evidences <- function(data, longer = FALSE, .keep = FALSE) {

    must_have_evidences(data)

    unnest_method <- `if`(longer, unnest_longer, unnest_wider)
    unnest_args <- `if`(
        longer,
        list(indices_to = 'direction'),
        list(simplify = FALSE)
    )
    tmp_col <- 'evs_tmp' %>% sym
    evs_col <- `if`(longer, 'all_evidences', 'evidences') %>% sym

    data %>%
    mutate(
        !!tmp_col := map(evidences, extract, EVIDENCES_KEYS),
        .after = evidences
    ) %>%
    exec(unnest_method, ., as.character(tmp_col), !!!unnest_args) %>%
    {`if`(
        as.character(tmp_col) %in% colnames(.),
        rename(., all_evidences = evidences, evidences = !!tmp_col),
        .
    )} %>%
    {`if`(.keep, ., select(., -!!evs_col))}

}


#' Filter evidences by dataset, resource and license
#'
#' @param data An interaction data frame with some columns containing evidences
#'      as nested lists.
#' @param ... The evidences columns to filter: tidyselect syntax is supported.
#'      By default the columns "evidences", "positive", "negative", "directed"
#'      and "undirected" are filtered, if present.
#' @param datasets A character vector of dataset names.
#' @param resources A character vector of resource names.
#'
#' @return The input data frame with the evidences in the selected columns
#'    filtered.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang expr
#' @importFrom tidyselect eval_select
#' @importFrom dplyr mutate across
#' @importFrom logger log_trace
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{only_from}}}
#'     \item{\code{\link{unnest_evidences}}}
#'     \item{\code{\link{from_evidences}}}
#' }
filter_evidences <- function(
        data,
        ...,
        datasets = NULL,
        resources = NULL,
        exclude = NULL
    ) {

    # NSE vs. R CMD check workaround
    .x <- NULL

    columns <-
        expr(...) %>%
        eval_select(data) %>%
        names %>%
        if_null_len0(
            EVIDENCES_KEYS %>%
            c('evidences') %>%
            intersect(colnames(data))
        )

    log_trace(
        paste0(
            'Filtering evidence columns: %s; to datasets: %s; ',
            'and resources: %s; excluding resources: %s'
        ),
        paste(columns, collapse = ', '),
        paste(if_null(datasets, 'any'), collapse = ', '),
        paste(if_null(resources, 'any'), collapse = ', '),
        paste(if_null(exclude, 'none'), collapse = ', ')
    )

    data %>%
    mutate(across(
        all_of(columns),
        ~filter_evs_lst(.x, datasets, resources, exclude)
    ))

}


#' Filter evidences within a nested list
#'
#' @param lst A list extracted from the JSON "evidences" column returned by
#'     OmniPath interactions queries.
#'
#' @importFrom magrittr %>% is_in not
#' @importFrom purrr map2 keep map
#' @noRd
filter_evs_lst <- function(
        lst,
        datasets = NULL,
        resources = NULL,
        exclude = NULL
    ) {

    filter_evs <- function(x) {

        x %>%
        {`if`(
            names(.) %>% is.null %>% not,
            map2(
                 names(.),
                 .,
                 function(name, value) {
                     `if`(
                          name %>% is_in(EVIDENCES_KEYS),
                          value %>% filter_evs,
                          value
                     )
                 }
            ),
            keep(., ~match_evidence(.x, datasets, resources, exclude))
        )} %>%
        {`if`(length(.) == 0L, NULL, .)}

    }


    lst %>%
    map(filter_evs)

}


#' Check if a single evidence matches the conditions
#'
#' @param ev An evidence data structure, which is a nested list, as it is
#'     extracted from the JSON in the "evidences" column from OmniPath queries.
#'
#' @importFrom magrittr %>%
#' @noRd
match_evidence <- function(
        ev,
        datasets = NULL,
        resources = NULL,
        exclude = NULL
    ) {

    (is.null(datasets) || ev$dataset %in% datasets) &&
    (
        is.null(resources) ||
        ev$resource %>% is_in(resources) ||
        ev$via %>% if_null_len0('') %>% is_in(resources)
    ) &&
    (
        is.null(exclude) ||
        (
            ev$resource %>% is_in(exclude) %>% not &&
            ev$via %>% if_null_len0('') %>% is_in(., exclude) %>% not
        )
    )

}

#' Recreate interaction data frame based on certain datasets and resources
#'
#' @param data An interaction data frame from the OmniPath web service with
#'     evidences column.
#' @param datasets Character: a vector of dataset labels. Only evidences from
#'     these datasets will be used.
#' @param resources Character: a vector of resource labels. Only evidences
#'     from these resources will be used.
#' @param .keep Logical: keep the "evidences" column.
#'
#' @return A copy of the interaction data frame restricted to the given
#'     datasets and resources.
#'
#' @details
#' The OmniPath interactions database fully integrates all attributes from all
#' resources for each interaction. This comes with the advantage that
#' interaction data frames are ready for use in most of the applications;
#' however, it makes it impossible to know which of the resources and
#' references support the direction or effect sign of the interaction. This
#' information can be recovered from the "evidences" column. The "evidences"
#' column preserves all the details about interaction provenances. In cases
#' when you want to use a faithful copy of a certain resource or dataset, this
#' function will help you do so. Still, in most of the applications the best is
#' to use the interaction data as it is returned by the web service.
#'
#' \strong{Note:} This function is automatically applied if the
#' `strict_evidences` argument is passed to any function querying interactions
#' (e.g. \code{\link{import_omnipath_interactions}}).
#'
#' @examples
#' \dontrun{
#' ci <- collectri(evidences = TRUE)
#' ci <- only_from(ci, datasets = 'collectri')
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang as_function
#' @importFrom purrr keep
#' @importFrom dplyr select
#' @importFrom tidyselect any_of
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{filter_evidences}}}
#'     \item{\code{\link{unnest_evidences}}}
#'     \item{\code{\link{from_evidences}}}
#' }
only_from <- function(
        data,
        datasets = NULL,
        resources = NULL,
        exclude = NULL,
        .keep = FALSE
    ) {

    if(is.null(datasets) && is.null(resources) && is.null(exclude)) {
        return(data)
    }

    must_have_evidences(data, wide_ok = TRUE)
    has_wide <- data %>% has_evidences_wide

    log_trace(
        'Restricting interaction records to datasets: %s; and resources: %s',
        paste(if_null(datasets, 'any'), collapse = ', '),
        paste(if_null(resources, 'any'), collapse = ', ')
    )

    data %>%
    {`if`(has_wide, ., unnest_evidences(., .keep = .keep))} %>%
    filter_evidences(
        datasets = datasets,
        resources = resources,
        exclude = exclude
    ) %>%
    from_evidences(.keep = .keep) %>%
    filter(if_any(EVIDENCES_KEYS, ~not(map_lgl(.x, is.null)))) %>%
    {`if`(has_wide, ., select(., -any_of(EVIDENCES_KEYS)))}

}


#' Recreate interaction records from evidences columns
#'
#' @param data An interaction data frame from the OmniPath web service with
#'     evidences column.
#' @param .keep Logical: keep the original "evidences" column when unnesting to
#'     separate columns by direction.
#'
#' @return A copy of the input data frame with all the standard columns
#'     describing the direction, effect, resources and references of the
#'     interactions recreated based on the contents of the nested list
#'     evidences column(s).
#'
#' @details
#' The OmniPath interaction data frames specify interactions primarily by
#' three columns: "is_directed", "is_stimulation" and "is_inhibition".
#' Besides these, there are the "sources" and "references" columns that are
#' always included in data frames created by OmnipathR and list the resources
#' and literature references for each interaction, respectively. The optional
#' "evidences" column is required to find out which of the resources and
#' references support the direction or effect sign of the interaction. To
#' properly recover information for arbitrary subsets of resources or
#' datasets, the evidences can be filtered first, and then the standard
#' data frame columns can be reconstructed from the selected evidences.
#' This function is able to do the latter. It expects either an "evidences"
#' column or evidences in their wide format 4 columns layout. It overwrites
#' the standard columns of interaction records based on data extracted from
#' the evidences, including the "curation_effort" and "consensus..." columns.
#'
#' \strong{Note:} The "curation_effort" might be calculated slightly
#' differently from the version included in the OmniPath web service. Here
#' we count the resources and the also add the number of references for each
#' resource. E.g. a resource without any literatur reference counts as 1,
#' while a resource with 3 references adds 4 to the value of the curation
#' effort.
#'
#' \strong{Note:} If the "evidences" column has been already unnested to
#' multiple columns ("positive", "negative", etc.) by
#' \code{\link{unnest_evidences}}, then these will be used;
#' otherwise, the column will be unnested within this function.
#'
#' \strong{Note:} This function (or rather its wrapper,
#' \code{\link{only_from}}) is automatically applied if the `strict_evidences`
#' argument is passed to any function querying interactions (e.g.
#' \code{\link{import_omnipath_interactions}}).
#'
#' @examples
#' \dontrun{
#' ci <- collectri(evidences = TRUE)
#' ci <- unnest_evidences(ci)
#' ci <- filter_evidences(datasets = 'collectri')
#' ci <- from_evidences(ci)
#' # the three lines above are equivalent to only_from(ci)
#' # and all the four lines above is equivalent to:
#' # collectri(strict_evidences = TRUE)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom purrr map_lgl
#' @importFrom dplyr select mutate left_join
#' @importFrom tidyr replace_na
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{filter_evidences}}}
#'     \item{\code{\link{unnest_evidences}}}
#'     \item{\code{\link{only_from}}}
#' }
from_evidences <- function(data, .keep = FALSE) {

    # NSE vs. R CMD check workaround
    target <- ce_directed <- ce_positive <- positive <- negative <- .x <-
    directed <- ce_negative <- ce_directed.x <- ce_directed.y <- NULL

    must_have_evidences(data, wide_ok = TRUE)

    prefix <- data$references[1L] %>% str_detect('^\\w+:')

    data %>%
    {`if`(has_evidences_wide(.), ., unnest_evidences(., .keep = .keep))} %>%
    mutate(
        is_directed = as.integer(lgl_from(., positive, negative, directed)),
        is_stimulation = as.integer(map_lgl(positive, ~length(.x) > 0L)),
        is_inhibition = as.integer(map_lgl(negative, ~length(.x) > 0L)),
        sources = resources_from(., EVIDENCES_KEYS),
        references = references_from(., EVIDENCES_KEYS, prefix = prefix),
        ce_positive = curation_effort_from(., positive),
        ce_negative = curation_effort_from(., negative),
        ce_directed = curation_effort_from(., directed),
        curation_effort = curation_effort_from(., EVIDENCES_KEYS),
        consensus_stimulation = as.integer(ce_positive >= ce_negative),
        consensus_inhibition = as.integer(ce_positive <= ce_negative)
    ) %>%
    {`if`(
        has_evidences(.),
        evidences_from(.),
        .
    )} %>%
    left_join(
        select(., source = target, target = source, ce_directed),
        by = c('source', 'target')
    ) %>%
    mutate(
        consensus_direction =
            (ce_directed.x >= ce_directed.y) %>%
            replace_na(TRUE) %>%
            as.integer
    ) %>%
    select(-ce_positive, -ce_negative, -ce_directed.x, -ce_directed.y) %>%
    {`if`('n_references' %in% names(.), count_references(.), .)} %>%
    {`if`('n_resources' %in% names(.), count_resources(.), .)}

}


#' Update the contents of the original "evidences" column
#'
#' Once the evidences have been split by direction and sign and undergone some
#' processing or filtering, we might wish to keep the original, composite
#' "evidences" column up to date, ensuring consistency. Here we copy the
#' processed evidences from the columns "positive", "negative", "directed" and
#' "undirected" back to the original "evidences" column.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom dplyr mutate pick
#' @importFrom tidyselect all_of
#' @importFrom purrr map2 transpose
#' @importFrom logger log_warn
#' @noRd
evidences_from <- function(data) {

    data %>%
    {`if`(
        has_evidences_wide(.),
        mutate(
            .,
            evidences = map2(
                evidences,
                transpose(pick(all_of(EVIDENCES_KEYS))),
                c
            )
        ),
        identity(.) %T>%
        log_warn(
            'Unable to update "evidences" column: split evidence columns ',
            '("positive", "directed", ...) are missing!'
        )
    )}

}


#' Logical vector from a set of list columns
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom purrr pmap_lgl
#' @noRd
lgl_from <- function(data, ...) {

    data %>%
    select(...) %>%
    pmap_lgl(~length(c(...)) > 0L)

}


#' Extract resources from a set of evidence list columns
#'
#' @importFrom magrittr %>%
#' @noRd
resources_from <- function(data, ..., collapse = ';') {

    data %>%
    chr_from(..., fn = resource_name, collapse = collapse)

}


#' Resource name from evidence
#'
#' @param ev Evidence data structure (list).
#'
#' @noRd
resource_name <- function(ev) {

    paste0(
        ev$resource,
        `if`(is.null(ev$via), '', '_'),
        ev$via
    )

}

#' Extract references from a set of evidence list columns
#'
#' @importFrom magrittr %>%
#' @noRd
references_from <- function(data, ..., prefix = TRUE, collapse = ';') {

    extract_references <- function(ev) {
        `if`(
            prefix,
            paste(ev$resource, ev$references, sep = ':'),
            ev$references, collapse = collapse
        ) %>%
        unique %>%
        sort %>%
        paste(collapse = collapse)
    }

    data %>%
    chr_from(..., fn = extract_references, collapse = collapse)

}


#' Character vector from a set of list columns
#'
#' @importFrom magrittr %>%
#' @importFrom purrr pmap_chr map_chr
#' @importFrom dplyr select
#' @noRd
chr_from <- function(data, ..., fn, collapse = ';') {

    map_fn <- function(...) {
        c(...) %>%
        map_chr(fn) %>%
        unique %>%
        sort %>%
        paste(collapse = collapse) %>%
        if_null_len0('')
    }

    data %>%
    select(...) %>%
    pmap_chr(map_fn)

}


#' Curation effort from a set of evidence list columns
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom purrr pmap_int map_int
#' @noRd
curation_effort_from <- function(data, ...) {

    ce <- function(...) {
        c(...) %>%
        map_int(~length(.x$references) + 1L) %>%
        sum
    }

    data %>%
    select(...) %>%
    pmap_int(ce)

}
