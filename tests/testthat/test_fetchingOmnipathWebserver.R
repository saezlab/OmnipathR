library(OmnipathR)
library(dplyr)

################################################################################
## Test of the functions fetching Omnipath Webserver
################################################################################

########## ########## ########## ##########
########## Resource Queries      ##########   
########## ########## ########## ##########
## This function is convenient for appropriate resource retrieval. Following:
## http://bioconductor.org/developers/how-to/web-query/
## It tries to retrieve the resource one or several times before failing.
omnipath_download <- function(URL, FUN, ..., N.TRIES=1L) {
  N.TRIES <- as.integer(N.TRIES)
  stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))
  
  while (N.TRIES > 0L) {
    result <- tryCatch(FUN(URL, ...), error=identity)
    if (!inherits(result, "error"))
      break
    N.TRIES <- N.TRIES - 1L
  }
  
  if (N.TRIES == 0L) {
    stop("'omnipath_download()' failed:",
         "\n  URL: ", URL,
         "\n  error: ", conditionMessage(result))
  }
  
  return(result) 
}



################################################################################
## Test of the functions getting the list of Omnipath available databases
################################################################################

## get_ptms_databases
op <- Sys.info()['sysname']

if (op != "Windows"){
    url_ptms <- 
        paste0("http://omnipathdb.org/ptms/?fields=sources&fields=references&", 
        "fields=curation_effort&genesymbols=1")

    ptms <- omnipath_download(url_ptms, read.table, sep = '\t', header = TRUE, 
       stringsAsFactors = FALSE)
    ptms_resources <- 
       unique(unlist(strsplit(x = as.character(ptms$sources),split = ";")))

    ## get_interaction_databases
    url_interactions <- paste0('http://omnipathdb.org/interactions?',
        'datasets=omnipath,pathwayextra,kinaseextra,ligrecextra,dorothea,',
        'tfregulons,mirnatarget,tf_target,tf_mirna,lncrna_mrna',
        '&fields=sources,references&genesymbols=1&dorothea_levels=A,B,C,D')
    interactions <- omnipath_download(url_interactions, read.table, sep = '\t', 
         header = TRUE, stringsAsFactors = FALSE)
    interaction_resources  <- 
        unique(unlist(strsplit(x = as.character(interactions$sources), 
            split = ";")))

    ## get_complexes_databases
    url_complexes <- 'http://omnipathdb.org/complexes?&fields=sources'
    complexes <- omnipath_download(url_complexes, read.csv, sep = '\t', 
        header = TRUE, stringsAsFactors = FALSE)
    complexes_resources <-
        unique(unlist(strsplit(x = as.character(complexes$sources), 
            split = ";")))

    ## get_annotation_databases
    url_annotations <- 'http://omnipathdb.org/annotations_summary'
    annotations <- omnipath_download(url_annotations, read.table, sep = '\t', 
        header = TRUE, stringsAsFactors = FALSE)
    annotations_db <- unique(annotations$source)

    ## get_intercell_categories
    url_intercell <- 'http://omnipathdb.org/intercell_summary'
    intercell <- omnipath_download(url_intercell, read.csv, sep = '\t', 
        header = TRUE, stringsAsFactors = FALSE)
    intercell_categories <- unique(intercell$category)
    intercell_parents <- unique(intercell$parent)

    ## Check the results between simulations and original functions
    test_that("Check the databases/categories available in Omnipath", {
        expect_equal(get_enzsub_resources(), sort(ptms_resources))
        # expect_equal(get_interaction_resources(), sort(interaction_resources))
        expect_equal(get_complex_resources(), sort(complexes_resources))
        expect_equal(get_annotation_resources(), sort(annotations_db))
        expect_equal(get_intercell_categories(), intercell_categories)
        expect_equal(get_intercell_generic_categories(),intercell_parents)
    })
}
