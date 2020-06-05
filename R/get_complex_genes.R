#' Get all the molecular complexes for a given gene(s)
#'
#' This function returns all the molecular complexes where an input set
#' of genes participate. User can choose to retrieve every complex where
#' any of the input genes participate or just retrieve these complexes where
#' all the genes in input set participate together.
#'
#' @param complexes complexes data frame (obtained using
#'    \code{\link{import_omnipath_complexes}})
#' @param select_genes vector containing the genes for whom complexes will be
#' retrieved (hgnc format).
#' @param total_match [default=FALSE] logical indicating if the user wants to
#' get all the complexes where any of the input genes participate (FALSE) or
#' to get only the complexes where all the input genes participate together
#' (TRUE)
#' @export
#' @return data.frame of complexes
#' @examples
#' complexes = import_omnipath_complexes(filter_databases=c("CORUM", "hu.MAP"))
#' query_genes = c("LMNA","BANF1")
#' complexes_query_genes = get_complex_genes(complexes,query_genes)
#' @seealso \code{\link{import_omnipath_complexes}})
get_complex_genes = function(complexes = import_Omnipath_complexes(),
    select_genes, total_match = FALSE){

    if(is.null(select_genes)){
        stop("A vector of genes should be provided")
    }

    if (!is.logical(total_match)){
        stop("total_match parameter should be logical")
    }

    if (total_match){
        complexes_geneset <-
        complexes[which(unlist(lapply(strsplit(complexes$components_genesymbols,
            '_'), function(x){
            sum(x %in% select_genes) == length(select_genes)}))),]
    } else {
        complexes_geneset <-
        complexes[which(unlist(lapply(strsplit(complexes$components_genesymbols,
            '_'), function(x){any(x %in% select_genes)}))),]
    }
    return(complexes_geneset)
}
