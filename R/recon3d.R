
#' @export
recon3d_raw_matlab <- function() {
    'recon3d_bigg' %>%
    generic_downloader(
      reader = R.matlab::readMat,
      url_key_param = list(),
      reader_param = list(),
      resource = 'recon3d_bigg',
      post = NULL,
      use_httr = FALSE
    ) %T>%
    load_success()
}


#' @export
recon3d_raw <- function() {
    download_to_cache('recon3d_bigg_json') %>%
    safe_json()
}


#' @noRd
recon3d_table <- function(name){
   recon3d_raw() %>% extract2(name) %>% tibble()
}


#' Metabolites from Recon-3D
#'
#' @rdname recon3d
#'
#' @examples
#' recon3d_metabolites()
#'
#' @return Data frame: tibble of metabolites.
#' @importFrom tidyr unnest
#' @export
recon3d_metabolites <- function(){

    recon3d_table("metabolites") %>%
    unnest(notes) %>%
    unnest(original_bigg_ids) %>%
    unnest(annotation)

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
    recon3d_table("compartments") %>% unnest() 
}
