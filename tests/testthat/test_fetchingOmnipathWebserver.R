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

    ptms <- getURL(url_ptms, read.table, sep = '\t', header = TRUE, 
       stringsAsFactors = FALSE)
    ptms_databases <- 
       unique(unlist(strsplit(x = as.character(ptms$sources),split = ";")))

    ## get_interaction_databases
    url_interactions <- paste0('http://omnipathdb.org/interactions?',
        'datasets=omnipath,pathwayextra,kinaseextra,ligrecextra',
        ',tfregulons,mirnatarget&fields=sources,references&genesymbols=1')
    interactions <- getURL(url_interactions, read.table, sep = '\t', 
         header = TRUE, stringsAsFactors = FALSE)
    interaction_databases  <- 
        unique(unlist(strsplit(x = as.character(interactions$sources), 
            split = ";")))

    ## get_complexes_databases
    url_complexes <- 'http://omnipathdb.org/complexes?&fields=sources'
    complexes <- getURL(url_complexes, read.csv, sep = '\t', header = TRUE,
        stringsAsFactors = FALSE)
    complexes_databases <-
        unique(unlist(strsplit(x = as.character(complexes$sources), 
            split = ";")))

    ## get_annotation_databases
    url_annotations <- 'http://omnipathdb.org/annotations_summary'
    annotations <- getURL(url_annotations, read.table, sep = '\t', 
        header = TRUE, stringsAsFactors = FALSE)
    annotations_db <- unique(annotations$source)

    ## get_intercell_categories
    url_intercell <- 'http://omnipathdb.org/intercell_summary'
    intercell <- getURL(url_intercell, read.csv, sep = '\t', header = TRUE,
        stringsAsFactors = FALSE)
    intercell_categories <- unique(intercell$category)
    intercell_parents <- unique(intercell$parent)

    ## Check the results between simulations and original functions
    test_that("Check the databases/categories available in Omnipath", {
        expect_equal(get_enzsub_resources(), ptms_resources)
        expect_equal(get_interaction_resources(), interaction_resources)
        expect_equal(get_complex_resources(), complexes_resources)
        expect_equal(get_annotation_resources(), annotations_db)
        expect_equal(get_intercell_categories(), intercell_categories)
        expect_equal(get_intercell_generic_categories(), intercell_parents)
    })

################################################################################
## Test the import of the databases with given filters
################################################################################

    ## import_Omnipath_PTMS
    ptms_filter_func <- import_omnipath_enzsub(resources=c("PhosphoSite"), 
        organism = 9606)
    ptms_filter_test <-OmnipathR:::filter_by_resource(ptms,"PhosphoSite")
    ## We replicate the functionality of the function
    ptms_filter_test$residue_offset <- 
        as.character(as.numeric(ptms_filter_test$residue_offset))
    ptms_filter_test$sources <- as.character(ptms_filter_test$sources)
    ptms_filter_test$references <- as.character(ptms_filter_test$references)
    ptms_filter_test$references <- 
        unlist(lapply(strsplit(ptms_filter_test$references,split = ";"),
            function(x)paste(unique(x),collapse=";")))
    ptms_filter_test$nsources <-
        unlist(lapply(strsplit(ptms_filter_test$sources,split = ";"),length))
    ptms_filter_test$nrefs <-
        unlist(lapply(strsplit(ptms_filter_test$references,split = ";"),length))

    ## import_AllInteractions
    interactions_filter_func <- 
        import_all_interactions(resources=c("SignaLink3"), 
             organism = 9606)
    interactions_filter_test <-
        OmnipathR:::filter_by_resource(interactions,"SignaLink3")
    ## We replicate the functionality of the function
    interactions_filter_test$sources <- 
        as.character(interactions_filter_test$sources)
    interactions_filter_test$nsources <-
        unlist(lapply(strsplit(interactions_filter_test$sources,split = ";"), 
            length))
    interactions_filter_test$references <- 
        as.character(interactions_filter_test$references)
    interactions_filter_test$references <- 
        unlist(lapply(strsplit(interactions_filter_test$references,split = ";"),
            function(x)paste(unique(x),collapse=";")))
    interactions_filter_test$nrefs <-
        unlist(lapply(strsplit(interactions_filter_test$references,split = ";"),
            length))

    ## import_Omnipath_complexes
    complexes_filter_func <- 
        import_omnipath_complexes(resources=c("CORUM", "hu.MAP"))
    complexes_filter_test <- 
        OmnipathR:::filter_by_resource(complexes,c("CORUM","hu.MAP"))
    ## We replicate the functionality of the function
    complexes_filter_test$sources <- as.character(complexes_filter_test$sources)
    complexes_filter_test$references <- 
        as.character(complexes_filter_test$references)
    # we remove references mentioned multiple times:
    complexes_filter_test$references <-
        unlist(lapply(strsplit(complexes_filter_test$references,split = ";"),
            function(x)paste(unique(x),collapse=";")))
    complexes_filter_test$nsources <-
        unlist(lapply(strsplit(complexes_filter_test$sources,split = ";"), 
            length))
    complexes_filter_test$nrefs <-
        unlist(lapply(strsplit(complexes_filter_test$references,split = ";"),
            length))

    ## import_Omnipath_annotations 
    annotations_filter_func <- 
        import_omnipath_annotations(proteins=c("TP53","LMNA"),
            resources=c("HPA_subcellular"))
    url_annotations <- 'http://omnipathdb.org/annotations?&proteins=TP53,LMNA'
    annotations_genes <- getURL(url_annotations, read.csv, sep = '\t', 
        header = TRUE, stringsAsFactors = FALSE)
    annotations_filter_test <- 
        dplyr::filter(annotations_genes, source=="HPA_subcellular")

    ## Check the results between simulations and original functions
    test_that("Check the fecthing of Omnipath webserver with filters", {
        expect_equal(ptms_filter_func, ptms_filter_test)
        expect_equal(interactions_filter_func, interactions_filter_test)
        expect_equal(complexes_filter_func, complexes_filter_test)
        expect_equal(annotations_filter_func, annotations_filter_test)
   })
}
