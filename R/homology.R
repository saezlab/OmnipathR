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



#' Function to generate a homologous OmniPath resource
#'
#' @param op_resource a resource in the format of OmniPath/LIANA
#'
#' @param target_organism `ncbi_taxid` or `name` of the target organism.
#' See `show_homologene` for available organisms via OmnipathR's `HomoloGene`
#'
#' @param max_homologs Determines the max number of homologs to be translated.
#' Certain genes will have multiple homolog matches, with some having also
#' certain isoforms considered. To exclude cases in which the number of
#' matched homologs is too high, one can adjust the homologs parameter.
#' Setting this to `1` would mean that one-to-many homolog matches are discarded
#'
#' @param .missing_fun approach to handle missing interactions. By default
#' set to `NULL` which would mean that any interactions without a homology
#' match will be filtered. This can be set to e.g. `str_to_title` when working
#' with murine symbols. Then if a gene has no matched homolog, instead of
#' discarding it, the `.missing_fun` will be used to format the name from human.
#' Hence, increasing the number of matches, but likely introducing some
#' mismatches.
#'
#' @param symbols_dict `NULL` by default, then `get_homologene_dict` is called
#' to generate a dictionary from OmniPathR's homologene resource. Alternatively,
#' one can pass their own symbols_dictionary.
#'
#' @param source name of the source (ligand) column
#'
#' @param target name of the target (receptor) column
#'
#' @param verbose logical for verbosity
#'
#' @return a converted ligand-receptor resource
#'
#' @export
generate_homologs <- function(op_resource,
                              target_organism,
                              max_homologs = 5,
                              .missing_fun = NULL,
                              symbols_dict = NULL,
                              columns = c("source_genesymbol",
                                          "target_genesymbol"),
                              verbose = TRUE){

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
                                        target_organism = target_organism)


    # Remove any missing antities
    if(is.null(.missing_fun)){

        # All missing entities
        missing_entities <- setdiff(entities,
                                    names(symbols_dict))

        if(verbose) message(
            stringr::str_glue("Entries without homologs:
                     {paste(missing_entities, collapse = '; ')}")
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
        minres <- dplyr::anti_join(minres,
                                   missing,
                                   by = columns)
    }


    # Obtain genes with multiple matches
    entity_2many <- symbols_dict %>%
        tibble::enframe(name = "genesymbol_source",
                       value = "genesymbol_target") %>%
        group_by(genesymbol_source) %>%
        dplyr::count(name = "n_match") %>%
        filter(n_match > 1 & n_match <= max_homologs) %>%
        pull(genesymbol_source)

    if(verbose) message(
        stringr::str_glue("One-to-many homolog matches: {paste(entity_2many, collapse = '; ')}"),
        verbose = verbose
    )

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
        dplyr::left_join(op_resource, by = columns) %>%
        .handle_complexes(symbols_dict = symbols_dict,
                          .missing_fun = .missing_fun,
                          columns=columns)

    ### Get interactions with many-many
    op_1many <- dplyr::anti_join(minres,
                                 op_notmany,
                                 by = columns)

    # Recursively translate the 1many resource
    op_one2_many <- op_1many %>%
        decomplexify(columns = columns) %>%
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
            op_row <- dplyr::inner_join(op_row.decomp,
                                        op_1many,
                                        by = columns) %>%
                # recomplexify
                select(-columns) %>%
                rename_with(~gsub("_complex", "", .x),
                            ends_with("complex")) %>%
                distinct() %>%
                # join back remainder of cols
                dplyr::left_join(op_resource,
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

    return(or_resource)

}


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
            sep_cols <- c(stringr::str_glue("col{rep(1:5)}"))
            col.complex <- stringr::str_glue("{col}_complex")

            resource <<- resource %>%
                mutate({{ col.complex }} :=
                           resource[[stringr::str_glue("{col}")]]) %>%
                separate(col,
                         into = sep_cols,
                         sep = "_",
                         extra = "drop",
                         fill = "right") %>%
                tidyr::pivot_longer(cols = all_of(sep_cols),
                                    values_to = col,
                                    names_to = NULL) %>%
                tidyr::drop_na(all_of(col)) %>%
                distinct() %>%
                mutate(across(all_of(c(col, col.complex)),
                              ~str_replace(., "COMPLEX:", "")))
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
                                target_organism,
                                id_type = "genesymbol"){

    # Load homology geneset
    hg_gs <- homologene_download(target = !!target_organism,
                                 source = 9606L, # always human
                                 id_type = !!id_type) %>%
        select(-hgroup) %>%
        # Limit to the universe of the resource
        filter(.data[[stringr::str_glue("{id_type}_source")]] %in% entities)

    # Convert to dictionary
    return(hg_gs %>% tibble::deframe())
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
    complex_dict <- map(columns,
                        ~.generate_complex_dict(translated_subunits, col=.x))


    # Bind all dictionaries
    dict <- pmap( # append multiple lists
        list(
            list(
                symbols_dict,
                purrr::flatten(complex_dict)
            )
        ), c) %>%
        purrr::flatten() %>%
        purrr::flatten()

    # get orthologous resource
    op_ortholog <- op_resource %>%
        mutate(across(columns,
                      ~recode.character2(.x,
                                         !!!dict,
                                         .missing_fun = .missing_fun)))

    return(op_ortholog)
}

#' Helper function to generate a dictionary also for the complexes
#'
#' @param translated_subunits decomplexified op_resource with already recoded
#' subunits
#'
#' @param entity column of interest (target or source)
#'
#' @noRd
.generate_complex_dict <- function(translated_subunits, col){
    col_complex <- stringr::str_glue("{col}_complex")

    filter(translated_subunits, str_detect(.data[[col_complex]], "_")) %>%
        dplyr::select(.data[[col_complex]], .data[[col]]) %>%
        distinct() %>%
        group_by(.data[[col_complex]]) %>%
        dplyr::group_nest(keep = FALSE, .key = 'subunits') %>%
        mutate(translated_complex = map(subunits, function(sub){
            glue::glue_collapse(sub[[col]], sep = "_") %>%
                as.character()
        })) %>%
        dplyr::select(.data[[col_complex]], translated_complex) %>%
        tidyr::unnest(translated_complex) %>%
        tibble::deframe()

}


#' Function to create a dictionary for one2many maps
#'
#' @param symbols_dict the human to mouse dict (should be a vect)
#' @param op_row.decomp decomp omni (we need all genes)
#' @param entity_2many entities that match to many homologs
#' @param columns columns
#'
#' @return a dictionary for the provided interaction
#'
#' @noRd
.create_row_dict <- function(symbols_dict,
                             op_row.decomp,
                             entity_2many,
                             columns){

    logical_row <- map_lgl(columns, function(col){
        any(op_row.decomp[[col]] %in% entity_2many)
    })

    dicts_row <- map(columns, function(col){
        symbols_dict[names(symbols_dict) %in% op_row.decomp[[col]]] %>%
            tibble::enframe(name = "genesymbol_source", value = "genesymbol_target")
    })

    # if only 1 col -> return dict
    if(length(columns)==1){
        symbols_dict[names(symbols_dict) %in% op_row.decomp[[columns]]] %>%
            tibble::enframe(name = "genesymbol_source", value = "genesymbol_target") %>%
            dplyr::group_by_all() %>%
            group_split() %>%
            map(~tibble::deframe(.x))

    } else if(all(logical_row)){ # if all columns have 1-to-many
        c(.bind_dicts(dicts_row[[1]], dicts_row[[2]]),
          .bind_dicts(dicts_row[[2]], dicts_row[[1]]))
    } else{ # main entity is the one with 1-to-many
        .bind_dicts(main_entity = keep(dicts_row, logical_row) %>% pluck(1),
                    secondary_entity = discard(dicts_row, logical_row) %>% pluck(1))
    }
}



#' Helper function to bind dictionaries
#'
#' @param main_entity entity with one-to-many mapping (ligand or receptor)
#' @param secondary_entity other entity (ligand or receptor)
#'
#' @details split by interaction, deframe and bind the secondary entity.
#' For example, if a ligand (`main_entity`) has one-to-many mapping to
#' multiple homologs, then main_entity will contain all homologs that match to
#' the ligand's subunits (or protein if not a complex). To this, we also attach
#' the genes of the receptor (`secondary_entity`) so that `generate_orthologs`
#' can translate both the ligand and receptor, and vice versa.
#'
#' @noRd
.bind_dicts <- function(main_entity, secondary_entity){
    main_entity %>%
        dplyr::group_by_all() %>%
        group_split() %>%
        map(~tibble::deframe(bind_rows(.x, secondary_entity)))
}


#' Homology translation
#'
#' Translates identifiers between organisms using orthology data from Ensembl.
#'
#' @param d Data frame or character vector.
#' @param ... Column specification: from zero to up to three arguments, with
#'     or without names. NSE is supported. Arguments beyond the third one
#'     will be ignored.
#'     \itemize{
#'         \item{
#'             The name of the arguments should be column names, the
#'             values identifier types, either as character or as symbols.
#'         }
#'         \item{
#'             Arguments without names assumed to be both column names and
#'             identifier types, e.g. a column called "uniprot" containing
#'             UniProt IDs.
#'         }
#'         \item{
#'             The first column spefication describes the source column,
#'             with identifiers of the source organism. This column must
#'             exist in the data and this will be the input of the homology
#'             translation. This column will be removed from the returned
#'             data frame.
#'         }
#'         \item{
#'             In case of "uniprot", the source column name can be anything,
#'             if it contains only UniProt IDs it will be handled accordingly.
#'         }
#'         \item{
#'             In case of "genesymbol", is enough if the source column name
#'             contains the word "genesymbol", e.g. "ligand_genesymbol".
#'         }
#'         \item{
#'             The second column spefication describes the target column,
#'             with its name and identifier type. If not provided, both
#'             the column name and type will be the same as the source
#'         }
#'         \item{
#'             Optionally a third column can be specified with another
#'             identifier type. This is convenient if you want, for example
#'             also Gene Symbols along with UniProt IDs.
#'         }
#'         \item{
#'             If no specification provided, the input assumed to have
#'             a column named either "uniprot" or "genesymbol", or be
#'             a character vector of UniProt IDs or Gene Symbols.
#'         }
#'     }
#' @param target Character or integer: name or identifier of the target
#'     organism (the one we translate to). The default target organism is
#'     mouse.
#' @param source Character or integer: name of identifier of the source
#'     organism (the one the IDs in the input data belong to). The default
#'     source organism is human.
#' @param ensembl_orthology_types Character vector: use only this orthology
#'     relationship types. Possible values are "one2one", "one2many" and
#'     "many2many".
#' @param ensembl_min_orthology_confidence Integer: use only orthology
#'     relations with at least this level of confidence. In Ensembl the
#'     confidence can be either 0 or 1, so only these values make sense.
#'     If 0, all the orthology records will be used, if 1, only the ones
#'     with higher confidence.
#'
#' @return Data frame with the translated columns or character vector with
#'     translated identifiers.
#'
#' @examples
#' \dontrun{
#' # these proteins are ULK1, IFNG, EGFR, TGFB1, IL1R1
#' human_uniprots <- c("O75385", "P01579", "P00533", "P01137", "P14778")
#' homology_translate(human_uniprots)
#' }
#'
#' @importFrom magrittr %<>% %>% extract
#' @importFrom purrr map_chr
#' @importFrom rlang enquos !! enquo sym %||% :=
#' @importFrom dplyr rename_with inner_join rename filter select mutate pull
#' @importFrom tidyselect ends_with
#' @export
homology_translate <- function(
    d,
    ...,
    target = 10090,
    source = 9606,
    ensembl_orthology_types = c('one2one', 'one2many'),
    ensembl_min_orthology_confidence = 1L
){

    UNIPROT_DEFAULTS <- c('uniprot', 'uniprot', 'genesymbol')
    GENESYMBOL_DEFAULTS <- c('genesymbol')

    ensembl_orthology_types %<>% map_chr(~sprintf('ortholog_%s', .x))

    ids <-
        enquos(...) %>%
        map(.nse_ensure_str) %>%
        if_null_len0(
            `if`(
                'uniprot' %in% colnames(d) || is_uniprot(d),
                UNIPROT_DEFAULTS,
                GENESYMBOL_DEFAULTS
            )
        ) %>%
        set_names(names(.) %||% unlist(.)) %>%
        set_names(ifelse(nchar(names(.)), names(.), unlist(.)))

    id_cols <- names(ids)
    id_types <- unlist(ids)

    if(length(id_cols) == 1L || id_cols[2] == '.replace'){

        id_cols[2] <- id_cols[1]
        id_types[2] <- id_types[1]

    }

    target_col <- id_cols[2]
    target_id_type <- id_types[2]
    target_ens_org <- ensembl_name(target)
    target_up_col <- sym(sprintf('%s_uniprot', target_ens_org))
    target_col <- sym(sprintf('%s_%s', target_ens_org, target_id_type))

    target2_col <- id_cols[3]
    target2_id_type <- sym(id_types[3])
    target2_col_tmp <- sym(sprintf('%s__tmp', target2_col))


    source_col_str <- id_cols[1]
    source_col <- sym(source_col_str)
    source_col_pos <- d %>% colnames %>% {which(. == source_col_str)}

    # to handle single vectors as inputs, we convert them to data frames
    vector_input <- !is.data.frame(d)

    if(vector_input){

        d <- tibble(!!source_col := d)

    }

    # have to keep these to move the new column
    # to the position of the original column
    before_col <-
        d %>% colnames %>%
        extract(source_col_pos - 1L) %>%
        {`if`(length(.) == 0L, NULL, .)}
    after_col <-
        d %>% colnames %>%
        extract(source_col_pos + 1L) %>%
        if_null_len0(NA) %>% # this happens if the input is vector
        {`if`(is.na(.) || !is.null(before_col), NULL, .)}

    # trying to be smart: if the name is not a registered ID type, but
    # "genesymbol" is in the name then we assume it's genesymbol, otherwise
    # if the column contains UniProt IDs, we know it's uniprot, and finally
    # we fall back to the provided value, it still might be a non registered
    # ID type
    source_id_type <-
        id_types[1] %>%
        {`if`(
            is_id_type(.),
            .,
            `if`(
                str_detect(., 'genesymbol'),
                'genesymbol',
                `if`(
                    is_uniprot(d$.),
                    'uniprot',
                    .
                )
            )
        )}

    source_ncbi <- ncbi_taxid(source)
    target_ncbi <- ncbi_taxid(target)

    # if we start from anything else than UniProt,
    # first we have to translate it to UniProt:
    if(source_id_type != 'uniprot'){

        source_uniprot <- sprintf('%s_uniprot', source_col_str)

        d %<>%
            translate_ids(
                !!source_col := source_id_type,
                !!sym(source_uniprot) := 'uniprot',
                organism = source_ncbi
            ) %>%
            select(-!!source_col)

        source_col_str <- source_uniprot
        source_col <- sym(source_col_str)

    }

    # homology table from Ensembl, filtered:
    ensembl_orthology(
        organism_a = source,
        organism_b = target
    ) %>%
    filter(
        b_orthology_type %in% ensembl_orthology_types &
        b_orthology_confidence >= ensembl_min_orthology_confidence
    ) %>%
    # translating source Ensembl peptide IDs to UniProts
    translate_ids(
        ensembl_peptide_id = ensp,
        uniprot,
        ensembl = TRUE
    ) %>%
    # removing the untranslated ones and
    # the ones without orthologous peptide
    filter(
        !is.na(b_ensembl_peptide) &
        !is.na(uniprot)
    ) %>%
    # adding an ugly label so we avoid any duplicate column names
    rename_with(~sprintf('%s__homology', .x)) %>%
    # joining the input data frame by source side uniprots
    inner_join(
        d,
        by = c(uniprot__homology = source_col_str)
    ) %>%
    # we want to keep the name from the original data frame
    # but that was on the left side, so we rename
    rename(!!source_col := uniprot__homology) %>%
    # translating the target Ensembl peptide IDs to UniProts
    translate_ids(
        b_ensembl_peptide__homology := ensp,
        !!target_up_col := uniprot,
        organism = target_ncbi,
        ensembl = TRUE
    ) %>%
    # if the final ID type is not UniProt, translating the target
    # UniProts to the desired ID type
    {`if`(
        target_id_type != 'uniprot',
        translate_ids(
            .,
            !!target_up_col := uniprot,
            !!target_col := !!sym(target_id_type),
            organism = target_ncbi
        ),
        .
    )} %>%
    # now we have the translation done
    # let's get rid of the redundant homology columns
    select(-ends_with('__homology')) %>%
    select(., -!!source_col) %>%
    # if a secondary identifier type is requested for the target
    # organism, let's translate it from UniProt:
    {`if`(
        is.na(target2_col),
        .,
        translate_ids(
            .,
            !!target_up_col := uniprot,
            !!target2_col_tmp := !!target2_id_type,
            organism = target_ncbi
        ) %>%
        mutate(
            !!sym(target2_col) := !!target2_col_tmp
        ) %>%
        select(-!!target2_col_tmp)
    )} %>%
    # if the target ID type is not UniProt,
    # we can remove the target UniProt column
    {`if`(
        target_id_type != 'uniprot',
        select(., -!!target_up_col),
        .
    )} %>%
    # remove the untranslated records
    filter(!is.na(!!target_col)) %>%
    # replace the original column with the translated one
    rename(!!source_col := !!target_col) %>%
    {`if`(
        vector_input,
        # output is a vector
        pull(., !!source_col) %>% unique,
        # move it to its original position
        {`if`(
            is.null(after_col),
            {`if`(
                is.null(before_col),
                # this happens to single column data frames only
                .,
                relocate(., !!source_col, .after = before_col)
            )},
            relocate(., !!source_col, .before = after_col)
        )}
    )}


}


#' Helper function to show available organisms via OmnipathR's homologene resource
#'
#' @importFrom OmnipathR homologene_raw
#'
#' @export
show_homologene <- function(){
    homologene_raw() %>%
        pull(ncbi_taxid) %>%
        unique %>%
        tibble(ncbi_taxid = .,
               name = common_name(.),
               latin = latin_name(.)) %>%
        na.omit() %>%
        arrange(name) %>%
        print(n=nrow(.))
}

