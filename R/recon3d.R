
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

recon3d_metabolites <- function(){
  recon3d_table("metabolites")
}
recon3d_reactions <- function(){
  recon3d_table("reactions")
}
recon3d_genes <- function(){
  recon3d_table("genes")
}
recon3d_compartments <- function(){
  recon3d_table("compartments")
}