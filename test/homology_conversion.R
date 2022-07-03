#' Helper Function to 'decomplexify' ligands and receptors into individual subunits
#'
#' @param resource a ligrec resource
#'
#' @param columns columns to separate and pivot long (e.g. genesymbol or uniprot),
#' `source_genesymbol` and `target_genesymbol` by default
#'
#' @return returns a longer tibble with complex subunits on seperate rows
#'
#' @details takes any number of columns, and assumes `_` as sep.
#'
#' @export
decomplexify <- function(resource,
                         columns = c("source_genesymbol",
                                     "target_genesymbol")){
    columns %>%
        map(function(col){
            sep_cols <- c(str_glue("col{rep(1:5)}"))
            col.complex <- str_glue("{col}_complex")

            resource <<- resource %>%
                mutate({{ col.complex }} :=
                           resource[[str_glue("{col}")]]) %>%
                separate(col,
                         into = sep_cols,
                         sep = "_",
                         extra = "drop",
                         fill = "right") %>%
                pivot_longer(cols = all_of(sep_cols),
                             values_to = col,
                             names_to = NULL) %>%
                tidyr::drop_na(col) %>%
                distinct() %>%
                mutate_at(.vars = c(col),
                          ~str_replace(., "COMPLEX:", ""))
        })
    return(resource)
}

#' Modified `dplyr::recode` function
#'
#' @inheritParams dplyr::recode
#'
#' @param .missing_fun Function to modify any missing homologs/strings
#'  `NULL` by default and any missing values will be discarded.
#'  For example, one could be set it to `str_to_title` to format all symbols, or
#'  any other format in a scenario where a homolog dictionary is not available
#'  for the organism of interest.
#'
#' @details enables to modify unmatched genesymbols
#' @import stringr
#'
#' @keywords internal
recode.character2 <- function(.x,
                              ...,
                              .default = NULL,
                              .missing = NULL,
                              .missing_fun) {
    .x <- as.character(.x)
    values <- rlang::list2(...)
    if (!all(rlang::have_name(values))) {
        bad <- which(!rlang::have_name(values)) + 1
        msg <- glue::glue("{dplyr:::fmt_pos_args(bad)} must be named.")
        rlang::abort(msg)
    }

    n <- length(.x)
    template <- dplyr:::find_template(values, .default, .missing)
    out <- template[rep(NA_integer_, n)]
    replaced <- rep(FALSE, n)

    for (nm in names(values)) {
        out <- dplyr:::replace_with(out,
                                    .x == nm,
                                    values[[nm]],
                                    paste0("`", nm, "`"))
        replaced[.x == nm] <- TRUE
    }

    .default <- dplyr:::validate_recode_default(.default, .x, out, replaced)

    if(!is.null(.missing_fun)){
        out <- dplyr:::replace_with(out,
                                    !replaced & !is.na(.x),
                                    exec(.missing_fun, .default),
                                    "`.default`")
    }

    out <- dplyr:::replace_with(out,
                                is.na(.x),
                                .missing,
                                "`.missing`")
    out
}

#' Helper function to generate a dictionary also for the complexes
#'
#' @param translated_subunits decomplexified op_resource with already recoded
#' subunits
#'
#' @param entity column of interest (target or source)
#'
#' @noRd
.generate_complex_dict <- function(translated_subunits, entity = "target"){
    genesymbol_complex <- str_glue("{entity}_genesymbol_complex")
    genesymbol <- str_glue("{entity}_genesymbol")

    filter(translated_subunits, str_detect(.data[[genesymbol_complex]], "_")) %>%
        dplyr::select(.data[[genesymbol_complex]], .data[[genesymbol]]) %>%
        distinct() %>%
        group_by(.data[[genesymbol_complex]]) %>%
        group_nest(keep = FALSE, .key = 'subunits') %>%
        mutate(translated_complex = map(subunits, function(sub){
            glue::glue_collapse(sub[[genesymbol]], sep = "_") %>%
                as.character()
        })) %>%
        dplyr::select(.data[[genesymbol_complex]], translated_complex) %>%
        unnest(translated_complex)

}

#' Helper function to get homologene dictionary
#'
#' @param entities genes to be converted
#'
#' @param target_organism target organism (obtain tax id from `show_homologene`)
#'
#' @keywords internal
#'
#' @importFrom OmnipathR homologene_download
get_homologene_dict <- function(entities,
                                target_organism){

    # Load homology geneset
    hg_gs <- homologene_download(target = !!target_organism,
                                 source = 9606L, # always human
                                 id_type = "genesymbol") %>%
        select(-hgroup) %>%
        # Limit to the universe of the resource
        filter(genesymbol_source %in% entities)

    # Convert to dictionary
    return(hg_gs %>% deframe())
}

#' Function to translate Human complexes to organism X homologs
#'
#' @param op_resource resource in the format of OmniPath/LIANA
#'
#' @param symbols_dict dictionary (named list) with human genesymbols and ortholog
#'  in a second species
#'
#' @param columns columns relevant for homology conversion
#'
#' @details Complexes cannot be joined directly, thus this function will
#' translate each subunit one by one.
#'
#' @noRd
.handle_complexes <- function(op_resource,
                              symbols_dict,
                              columns=columns,
                              .missing_fun = NULL){
    # decomplexify
    op_resource_decomplex <- op_resource %>%
        decomplexify(columns=columns)

    # translate subunits
    translated_subunits <- op_resource_decomplex %>%
        mutate(across(columns,
                      ~recode.character2(.x, !!!symbols_dict,
                                         .missing_fun = .missing_fun)))

    # Generate Dictionaries for complexes
    target_complex_dict <- .generate_complex_dict(translated_subunits,
                                                  entity="target") %>%
        deframe()
    source_complex_dict <- .generate_complex_dict(translated_subunits,
                                                  entity="source") %>%
        deframe()

    # Bind all dictionaries
    dict <- pmap( # append multiple lists
        list(
            list(
                symbols_dict,
                target_complex_dict,
                source_complex_dict
            )
        ), c) %>%
        flatten %>%
        flatten()

    # get orthologous resource
    op_ortholog <- op_resource %>%
        mutate(across(ends_with("genesymbol"),
                      ~recode.character2(.x, !!!dict,
                                         .missing_fun = .missing_fun)))

    return(op_ortholog)
}




require(tidyverse)

# Test with SignaLink Pathways
op_resource <-
    OmnipathR::import_omnipath_annotations(
        resources = 'SignaLink_pathway',
        wide = TRUE
    ) %>%
    select(genesymbol, pathway)
columns <- c("genesymbol")

op_resource <- liana::select_resource("Consensus")[[1]]

# Params
columns <- c("target_genesymbol", "source_genesymbol")
target_organism <- 10090
max_homologs <- 5
.missing_fun <- NULL

### FUN
# Minimum column set resource
minres <- op_resource %>% select(!!columns)

# Get decomplexified resource
decomp <- decomplexify(minres, columns = columns)

# Get union of symbols
entities <- purrr::reduce(map(columns, function(col) decomp[[col]]), union)

# generate homology geneset
symbols_dict <- get_homologene_dict(entities = entities,
                                    target_organism = target_organism)


# Remove any missing antities
if(is.null(.missing_fun)){

    # All missing entities
    missing_entities <- setdiff(entities,
                                names(symbols_dict))

    liana_message(
        str_glue("Entries without homologs:
                     {paste(missing_entities, collapse = '; ')}"),
        verbose = verbose
    )

    # Keep only interactions for which all proteins (incl. subunits) have a matching homologue
    missing <- decomp %>%
        # check if neither is in missing entities
        mutate(lr_present = !if_any(columns, function(x) x %in% missing_entities)) %>%
        group_by(across(ends_with("complex"))) %>%
        summarise(all_present = mean(lr_present), .groups = "keep") %>%
        # only keep those that are present
        filter(all_present < 1) %>%
        # remove _complex for join
        rename_with(~gsub("_complex", "", .x), ends_with("complex"))

    # Remove interactions without matches
    minres <- anti_join(minres,
                        missing,
                        by = columns)
}


# Obtain genes with multiple matches
entity_2many <- symbols_dict %>%
    enframe(name = "genesymbol_source",
            value = "genesymbol_target") %>%
    group_by(genesymbol_source) %>%
    count(name = "n_match") %>%
    filter(n_match > 1 & n_match <= max_homologs) %>%
    pull(genesymbol_source)

### IF NOT many2many .handle_complexes for all
op_notmany <- minres %>%
    decomplexify(columns = columns) %>%
    filter(!if_any(!ends_with("complex"), function(x) x %in% entity_2many))

# homologous omnipath resource (with 1to1 alone)
or_notmany <- op_notmany %>%
    # recomplexify
    select(-columns) %>%
    rename_with(~gsub("_complex", "", .x), ends_with("complex")) %>%
    distinct() %>%
    .handle_complexes(symbols_dict = symbols_dict,
                      .missing_fun = .missing_fun,
                      columns=columns)


