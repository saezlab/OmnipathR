####
# Test with SignaLink Pathways
signalink <-
  OmnipathR::import_omnipath_annotations(
    resources = 'SignaLink_pathway',
    wide = TRUE
  ) %>%
  select(genesymbol, pathway)

or_signalink <- generate_homologs(op_resource = signalink,
                                  target_organism = 10090,
                                  columns = c("genesymbol"))


# test with LIANA
op_resource <- liana::select_resource("Consensus")[[1]]
or_resource <- generate_homologs(op_resource = op_resource,
                                  target_organism = 10090,
                                  columns = c("target_genesymbol", "source_genesymbol"))


#  FUNNN step by step -----
require(tidyverse)

# Params
target_organism <- 10090
max_homologs <- 5
.missing_fun <- NULL
verbose <- TRUE

### FUN
op_resource %<>% mutate(across(all_of(columns),
                               ~str_replace(., "COMPLEX:", "")))

# Minimum column set resource
minres <- op_resource %>% select(!!columns)

# Get decomplexified resource
decomp <- decomplexify(minres, columns = columns)

# Get union of symbols
entities <- purrr::reduce(map(columns, function(col) decomp[[col]]), union)

# generate homology geneset
symbols_dict <- get_homologene_dict(entities = entities,
                                    target_organism = target_organism,
                                    id_type = "genesymbol")


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
        group_by(across(all_of(ends_with("complex")))) %>%
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
    filter(!if_any(!ends_with("complex"), function(x) x %in% entity_2many)) %>%
    # recomplexify
    select(-columns) %>%
    rename_with(~gsub("_complex", "", .x), ends_with("complex")) %>%
    distinct()

# homologous omnipath resource (with 1to1 alone)
or_notmany <- op_notmany %>%
    # join back missing cols
    left_join(op_resource, by = columns) %>%
    .handle_complexes(symbols_dict = symbols_dict,
                      .missing_fun = .missing_fun,
                      columns=columns)

### Get interactions with many-many
op_1many <- anti_join(minres,
                      op_notmany,
                      by = columns)

# Recursively translate the 1many resource
op_one2_many <- op_1many %>%
    liana::decomplexify(columns = columns) %>%
    group_by(across(ends_with("complex"))) %>%
    group_split()

# On all one2many !!!
suppressWarnings(pb <- dplyr::progress_estimated(length(op_one2_many)))
or_many <- op_one2_many %>%
    map(function(op_row.decomp){
        if(verbose) pb$tick()$print()

        # create dictonary by row
        dicts_row <- .create_row_dict(symbols_dict,
                                      op_row.decomp,
                                      entity_2many,
                                      columns)

        # The tibble to be translated to all homologs
        op_row <- inner_join(op_row.decomp,
                             op_1many,
                             by = columns) %>%
            # recomplexify
            select(-columns) %>%
            rename_with(~gsub("_complex", "", .x),
                        ends_with("complex")) %>%
            distinct() %>%
            # join back remainder of cols
            left_join(op_resource,
                      by = columns)

        # Return all matches for all homologs
        or_row <- map(dicts_row, function(d){
            .handle_complexes(op_row,
                              symbols_dict = d,
                              .missing_fun = .missing_fun,
                              columns = columns)
        }) %>%
            bind_rows() %>%
            distinct()
    }) %>%
    bind_rows()

# Bind 1to1 and 1tomany
or_resource <- bind_rows(or_notmany, or_many) %>%
    select(-ends_with("complex")) %>%
    select(!!columns, everything())





### create row dict revamp ----
op_row.decomp <- op_one2_many[[7]]

dicts_row <- .create_row_dict(symbols_dict,
                              op_row.decomp,
                              entity_2many)






logical_row <- map_lgl(columns, function(col){
    any(op_row.decomp[[col]] %in% entity_2many)
})

dicts_row <- map(columns, function(col){
    symbols_dict[names(symbols_dict) %in% op_row.decomp[[col]]] %>%
        enframe(name = "genesymbol_source", value = "genesymbol_target")
    })

# if only 1 col -> return dict
if(length(columns)==1){
    symbols_dict[names(symbols_dict) %in% op_row.decomp[[columns]]] %>%
        enframe(name = "genesymbol_source", value = "genesymbol_target")

} else if(all(logical_row)){ # if all columns have 1-to-many
    c(.bind_dicts(dicts_row[[1]], dicts_row[[2]]),
      .bind_dicts(dicts_row[[2]], dicts_row[[1]]))
} else{ # main entity is the one with 1-to-many
    .bind_dicts(main_entity = keep(dicts_row, logical_row) %>% pluck(1),
                secondary_entity = discard(dicts_row, logical_row) %>% pluck(1))
}








if(length(dicts_row)==1){

}


# do the ligand or receptor posses homologs
is.l.2many <- any(op_row.decomp$source_genesymbol %in% entity_2many)
is.r.2many <- any(op_row.decomp$target_genesymbol %in% entity_2many)

# Ligands
dicts_row.l <-
    symbols_dict[names(symbols_dict) %in% op_row.decomp$source_genesymbol] %>%
    enframe(name = "genesymbol_source", value = "genesymbol_target")
# Receptors
dicts_row.r <-
    symbols_dict[names(symbols_dict) %in% op_row.decomp$target_genesymbol] %>%
    enframe(name = "genesymbol_source", value = "genesymbol_target")

if(all(is.l.2many, is.r.2many)){
    c( # both
        .bind_dicts(dicts_row.l, dicts_row.r),
        .bind_dicts(dicts_row.r, dicts_row.l)
    )

} else if(is.l.2many){
    .bind_dicts(dicts_row.l, dicts_row.r)
} else if(is.r.2many){
    .bind_dicts(dicts_row.r, dicts_row.l)
}

