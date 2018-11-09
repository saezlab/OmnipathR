#' get signs for ptms interactions
#'
#' ptms data does not contain sign (activation/inhibition), we generate this information
#' based on the interaction network
#'
#' @param ptms ptms data frame
#' @param interactions interaction data frame
#' @return data.frame of ptms with is_inhibition and is_stimulation columns
#' @examples
#' ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))
#' interactions = import_Omnipath_Interactions()
#' ptms = get_signed_ptms(ptms,interactions)
get_signed_ptms <- function(ptms = import_Omnipath_PTMS(), interactions = import_Omnipath_Interactions()){

	signed.ptms = merge(ptms,interactions[,c("source","target","is_stimulation","is_inhibition")],by.x = c("enzyme","substrate"),by.y = c("source","target"),all.x = T)

}
