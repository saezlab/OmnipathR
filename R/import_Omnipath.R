

#' Import Omnipath PTMS
#'
#' imports the ptms database from http://omnipathdb.org/ptms
#'
#' @return data.frame
#' @export
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are removed.
#' @examples
#' ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))
import_Omnipath_PTMS = function (from_cache_file=NULL,
								 filter_databases = .get_ptms_databases()){
	url_ptms <- 'http://omnipathdb.org/ptms/?fields=sources&fields=references&genesymbols=1'

	if(is.null(from_cache_file)){
		ptms = read.table(url_ptms, sep = '\t', header = TRUE,stringsAsFactors = F)
		print(paste0("Downloaded ", nrow(ptms), " interactions"))
	}else{
		load(from_cache_file)
	}

	if(!is.null(filter_databases)){
		filteredPTMS = .filter_sources(ptms,databases = filter_databases)
	}else{
		filteredPTMS = ptms
	}
	filteredPTMS$residue_offset = as.character(as.numeric(filteredPTMS$residue_offset))


	filteredPTMS$sources = as.character(filteredPTMS$sources)
	filteredPTMS$references = as.character(filteredPTMS$references)
	filteredPTMS$nsources =	unlist(lapply(strsplit(filteredPTMS$sources,split = ";"),length))
	filteredPTMS$nrefs =	unlist(lapply(strsplit(filteredPTMS$references,split = ";"),length))


	return(filteredPTMS)
}

#' Get ptms databases
#'
#' get the names of the databases from omnipath.org/ptms
#' @return character vector with the databases
.get_ptms_databases = function(){
	url_ptms <- 'http://omnipathdb.org/ptms/?fields=sources'
	ptms = read.table(url_ptms, sep = '\t', header = TRUE,stringsAsFactors = F)
	return(unique(unlist(strsplit(x = as.character(ptms$sources),split = ";"))))
}

#' Import Omnipath interaction database
#'
#' imports the database from http://omnipathdb.org/interactions, which contains
#' only interactions with references
#'
#' @return data.frame
#' @export
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are removed.
#' @examples
#' interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite", "Signor"))
import_Omnipath_Interactions = function (from_cache_file=NULL,
										 filter_databases = .get_interaction_databases()){
	url <- 'http://omnipathdb.org/interactions?fields=sources&fields=references&genesymbols=1'
	if(is.null(from_cache_file)){
		interactions = read.table(url, sep = '\t', header = TRUE, stringsAsFactors = F)
		print(paste0("Downloaded ", nrow(interactions), " interactions"))
	}else{
		load(from_cache_file)
	}

	if(!is.null(filter_databases)){
		filteredInteractions = .filter_sources(interactions,databases = filter_databases)
	}else{
		filteredInteractions = interactions
	}

	filteredInteractions$sources = as.character(filteredInteractions$sources)
	filteredInteractions$references = as.character(filteredInteractions$references)
	filteredInteractions$nsources =	unlist(lapply(strsplit(filteredInteractions$sources,split = ";"),length))
	filteredInteractions$nrefs =	unlist(lapply(strsplit(filteredInteractions$references,split = ";"),length))


	return(filteredInteractions)
}

#' Get interaction databases
#'
#' get the names of the databases from omnipath.org/interactions
#' @return character vector with the databases
#'
.get_interaction_databases = function(){
	url_interactions <- 'http://omnipathdb.org/interactions/?fields=sources'
	interactions = read.table(url_interactions, sep = '\t', header = TRUE,stringsAsFactors = F)
	return(unique(unlist(strsplit(x = as.character(interactions$sources),split = ";"))))
}




# filtering the PTMS and interactions which are reported in Signor and PhosphoSite
#' .filter_sources
.filter_sources = function(interactions, databases){
	# takes either the interactions or the PTMS and removes interactions which are
	# not reported by the asked databases.

	nInter = nrow(interactions)
	subsetInteractions = plyr::adply(interactions,1,function(interaction){
		if (any(strsplit(interaction[["sources"]],split = ";")[[1]] %in% databases)){
			return(interaction)
		}else
			return(NULL)
	}
	)
	nInterPost = nrow(subsetInteractions)

	print(paste0("removed ",nInter-nInterPost, " interactions during database filtering."))
	return(subsetInteractions)
}
