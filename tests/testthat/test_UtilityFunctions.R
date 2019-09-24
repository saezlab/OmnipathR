library(OmnipathR)
library(dplyr)
library(igraph)

################################################################################
## Test utlity functions from the OmnipathR package
################################################################################

## get_complex_genes. We simulate a result for this function.
complexes <- import_Omnipath_complexes()
query_genes <- c("LMNA","BANF1")
complexes_geneset_test <- 
    complexes[which(unlist(lapply(strsplit(complexes$components_genesymbols,
    '-'), function(x){any(x %in% query_genes)}))),]

complexes_geneset_func <- get_complex_genes(select_genes=query_genes)

## get_signed_ptms. 
ptms <- import_Omnipath_PTMS()
interactions <- import_Omnipath_Interactions()
query_gene <- "Q9Y6N7"
ptms_gene <- dplyr::filter(ptms, enzyme ==query_gene)
interactions_gene <- dplyr::filter(interactions, 
    target==query_gene | source==query_gene)

get_signed_ptms_test <- merge(ptms_gene,interactions_gene[,
    c("source","target","is_stimulation","is_inhibition")],
    by.x = c("enzyme","substrate"),by.y = c("source","target"),all.x = TRUE)

get_signed_ptms_test <-  
    get_signed_ptms_test[with(get_signed_ptms_test, 
    order(as.numeric(residue_offset))), ]
rownames(get_signed_ptms_test)  <- seq(nrow(get_signed_ptms_test)) 

get_signed_ptms_func <- get_signed_ptms() %>% 
    dplyr::filter(enzyme == query_gene)   %>% 
    dplyr::arrange(as.numeric(residue_offset))

## Check the results between simulations and original functions. 
test_that("Check that the original functions return the expected result", {
	expect_equal(complexes_geneset_func, complexes_geneset_test)
	expect_equal(get_signed_ptms_func, get_signed_ptms_test)
	expect_equal(is_igraph(ptms_graph(ptms_gene)), TRUE)
	expect_equal(is_igraph(interaction_graph(interactions)), TRUE)
})
