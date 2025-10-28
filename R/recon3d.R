
#' Recon3D model from BiGG
#'
#' Returns the content extracted from Matlab file.
#'
#' @importFrom magrittr %>%
#' @importFrom R.utils gunzip
#' @export
recon3d_raw_matlab <- function() {

    .slow_doctest()

    'recon3d_bigg' %>%
    download_to_cache() %>%
    gunzip(temporary = TRUE, remove = FALSE) %>%
    R.matlab::readMat()

}


#' Recon-3D from JSON
#'
#' @examples
#' recon3d_raw()
#'
#' @importFrom magrittr %>%
#' @export
recon3d_raw <- function() {

    .slow_dowtest()

    download_to_cache('recon3d_bigg_json') %>%
    safe_json()

}


#' Move Recon-3D from Matlab object to tibble
#'
#' @param name Character: top level name in the Matlab object.
#'
#' @noRd
recon3d_table <- function(name){

    .slow_dowtest()

    recon3d_raw() %>%
    extract2(name) %>%
    tibble()

}


#' Metabolites from Recon-3D
#'
#' @rdname recon3d
#'
#' @param extra_hmdb Logical: add extra HMDB IDs from Virtual Metabolic Human.
#'
#' @examples
#' recon3d_metabolites()
#'
#' @return Data frame: tibble of metabolites.
#' @importFrom tidyr unnest
#' @export
recon3d_metabolites <- function(extra_hmdb = TRUE){

    .slow_dowtest()

    recon3d_table("metabolites") %>%
    unnest(notes) %>%
    unnest(original_bigg_ids) %>%
    unnest(annotation) %>%
    {`if`(
        extra_hmdb,
        left_join(
            .,
            recon3d_raw_vmh() %>%
            gem_matlab_tibble("mets", "metHMDBID") %>%
            unnest(cols = c('mets', 'metHMDBID')) %>%
            unnest(cols = c('mets')),
            by = c("original_bigg_ids" = "mets")
        ) %>%
        rowwise() %>%
        mutate(hmdb = list(if_null(unique(c(hmdb, metHMDBID)), NA_character_))) %>%
        select(-metHMDBID)
        ,
        .
    )}

}

#' Reactions from Recon-3D
#'
#' @rdname recon3d
#'
#' @examples
#' recon3d_reactions()
#'
#' @return Data frame: tibble of reactions.
#' @export
recon3d_reactions <- function(){

    .slow_dowtest()

    recon3d_table("reactions") %>%
    unnest(notes) %>%
    unnest(metabolites) %>%
    unnest(original_bigg_ids) %>%
    relocate(original_bigg_ids, .after = 2L)

}


#' Genes from Recon-3D
#'
#' @rdname recon3d
#'
#' @examples
#' recon3d_genes()
#'
#' @return Data frame: tibble of genes.
#' @export
recon3d_genes <- function(){

    .slow_dowtest()

    recon3d_table("genes") %>%
    unnest(notes) %>%
    unnest(annotation) %>%
    unnest(original_bigg_ids) %>%
    slice(-1)

}


#' Compartments from Recon-3D
#'
#' @rdname recon3d
#'
#' @examples
#' recon3d_compartments()
#'
#' @return Data frame: tibble of compartments.
#' @export
recon3d_compartments <- function(){

    .slow_dowtest()

    recon3d_raw() %>%
    extract2('compartments') %>%
    {tibble(code = names(.), name = unlist(.))}

}


#' Recon-3D model from Virtual Metabolic Human
#'
#' @examples
#' recon3d_raw_vmh()
#'
#' @importFrom magrittr %>%
#' @export
recon3d_raw_vmh <- function() {

    .slow_dowtest()

    'recon3d_model' %>%
    generic_downloader(
        reader = R.matlab::readMat,
        url_key_param = list(),
        reader_param = list(),
        resource = 'Recon-3D',
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}
