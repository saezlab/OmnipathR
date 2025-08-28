recon3d_raw <- function(){
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


