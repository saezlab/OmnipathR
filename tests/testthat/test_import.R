library(OmnipathR)
library(dplyr)
library(tidyr)
library(jsonlite)

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


## Retrieves a list of resources for a certain query type and dataset.
## We test each query with downloading only a few resources.
get_resources_test <- function(
        query_type,
        dataset = NULL,
        license = NULL,
        top = NULL
    ){

    license_levels <- list(
        academic = 1,
        commercial = 2
    )

    resources_url <- 'https://omnipathdb.org/resources'
    resources <- jsonlite::fromJSON(txt = resources_url)

    resources <- Filter(
        function(res_name){

            res <- resources[[res_name]]

            (
                query_type %in% names(res$queries) && (
                    is.null(dataset) ||
                    dataset %in% res$queries[[query_type]]$datasets
                ) && (
                    is.null(license) || (
                        license_levels[[res$license$purpose]] >=
                        license_levels[[license]]
                    )
                )
            )

        },
        names(resources)
    )

    resources <- sort(unique(resources))

    `if`(is.null(top), resources, head(resources, top))

}


## Here we build a list with parameters that we later iterate through
## in order to test all the various download methods.
test_items <- list(
    omnipath_complexes = NULL,
    omnipath_enzsub = NULL,
    omnipath_intercell = NULL,
    omnipath_annotations = NULL,
    omnipath_interactions = list(datasets = 'omnipath'),
    kinaseextra_interactions = list(datasets = 'kinaseextra'),
    pathwayextra_interactions = list(datasets = 'pathwayextra'),
    ligrecextra_interactions = list(datasets = 'ligrecextra'),
    dorothea_interactions = list(datasets = 'dorothea'),
    tf_target_interactions = list(datasets = 'tf_target'),
    tf_mirna_interactions = list(datasets = 'tf_mirna'),
    mirnatarget_interactions = list(datasets = 'mirnatarget')
)

for(item in names(test_items)){

    resource_col <- list(
        complexes = 'sources',
        interactions = 'sources',
        enzsub = 'sources',
        annotations = 'source',
        intercell = 'database'
    )

    method <- get(sprintf('import_%s', item))
    dataset <- `if`(
        'datasets' %in% names(test_items[[item]]),
        test_items[[item]]$datasets,
        NULL
    )
    query_type <- tail(strsplit(item, '_')[[1]], 1)

    resources <- get_resources_test(
        query_type = query_type,
        dataset = dataset,
        license = 'academic',
        top = 5
    )

    response <- method(resources = resources)

    if(query_type %in% c('enzsub', 'interactions', 'complexes')){
        response <- response %>% separate_rows(sources, sep = ';')
    }

    resources_in_response <- response %>%
        pull(resource_col[[query_type]]) %>%
        unique()

    test_that(
        expect_true(
            all(resources %in% resources_in_response)
        ),
        expect_gt(nrow(response), 10)
    )

}