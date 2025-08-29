
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

    path <- download_to_cache('recon3d_bigg_json')

    safe_json(path)

}
