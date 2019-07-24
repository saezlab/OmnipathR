########## ########## ########## ##########
########## PTMS                  ##########   
########## ########## ########## ##########

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
	# we remove references mentioned multiple times:
	filteredPTMS$references = unlist(lapply(strsplit(filteredPTMS$references,split = ";"),function(x)paste(unique(x),collapse=";")))
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


########## ########## ########## ##########
########## INTERACTIONS          ##########   
########## ########## ########## ##########


## Import Omnipath interaction Database. The new version of Ominpath contains 
## several different datastes. 

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
	# we remove references mentioned multiple times:
	filteredInteractions$references = unlist(lapply(strsplit(filteredInteractions$references,split = ";"),function(x)paste(unique(x),collapse=";")))
	filteredInteractions$nsources =	unlist(lapply(strsplit(filteredInteractions$sources,split = ";"),length))
	filteredInteractions$nrefs =	unlist(lapply(strsplit(filteredInteractions$references,split = ";"),length))
  
	return(filteredInteractions)
}



#' Imports from Omnipath webservice the interactions from Pathwayextra dataset
#'
#' Imports the dataset from: 
#' 'http://omnipathdb.org/interactions?datasets=pathwayextra', 
#' which contains activity flow interactions without literature reference
#' 
#' @return data.frame
#' @export
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are removed.
#' @examples
#' interactions <- import_PathwayExtra_Interactions(filter_databases=c("BioGRID","IntAct"))
import_PathwayExtra_Interactions = function (from_cache_file=NULL,
        filter_databases = .get_interaction_databases()){
  url <- 
    'http://omnipathdb.org/interactions?datasets=pathwayextra&fields=sources&genesymbols=1'
  if(is.null(from_cache_file)){
    interactions = read.table(url, sep = '\t', header = TRUE, stringsAsFactors = F)
    print(paste0("Downloaded ", nrow(interactions), " interactions"))
  }else{
    load(from_cache_file)
  }
  
  if(!is.null(filter_databases)){
    filteredInteractions <- .filter_sources(interactions,databases = filter_databases)
  }else{
    filteredInteractions <- interactions
  }
  
  filteredInteractions$sources <- as.character(filteredInteractions$sources)
  filteredInteractions$nsources <-	
    unlist(lapply(strsplit(filteredInteractions$sources,split = ";"),length))
  
  
  return(filteredInteractions)
}


#' Imports from Omnipath webservice the interactions from kinaseextra dataset
#'
#' Imports the dataset from: 
#' 'http://omnipathdb.org/interactions?datasets=kinaseextra', 
#' which contains enzyme-substrate interactions without literature reference
#' 
#' @return data.frame
#' @export
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are removed.
#' @examples
#' interactions <- import_KinaseExtra_Interactions(filter_databases=c("PhosphoPoint","PhosphoSite"))
import_KinaseExtra_Interactions = function (from_cache_file=NULL,
      filter_databases = .get_interaction_databases()){
  url <- 
    'http://omnipathdb.org/interactions?datasets=kinaseextra&fields=sources&genesymbols=1'
  if(is.null(from_cache_file)){
    interactions = read.table(url, sep = '\t', header = TRUE, stringsAsFactors = F)
    print(paste0("Downloaded ", nrow(interactions), " interactions"))
  }else{
    load(from_cache_file)
  }
  
  if(!is.null(filter_databases)){
    filteredInteractions <- .filter_sources(interactions,databases = filter_databases)
  }else{
    filteredInteractions <- interactions
  }
  
  filteredInteractions$sources <- as.character(filteredInteractions$sources)
  filteredInteractions$nsources <-	
    unlist(lapply(strsplit(filteredInteractions$sources,split = ";"),length))
  
  
  return(filteredInteractions)
}


#' Imports from Omnipath webservice the interactions from ligrecextra dataset
#'
#' Imports the dataset from: 
#' 'http://omnipathdb.org/interactions?datasets=ligrecextra', 
#' which contains ligand-receptor interactions without literature reference
#' 
#' @return data.frame
#' @export
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are removed.
#' @examples
#' interactions <- import_LigrecExtra_Interactions(filter_databases=c("HPRD","Guide2Pharma"))
import_LigrecExtra_Interactions = function (from_cache_file=NULL,
              filter_databases = .get_interaction_databases()){
  url <- 
    'http://omnipathdb.org/interactions?datasets=ligrecextra&fields=sources&genesymbols=1'
  if(is.null(from_cache_file)){
    interactions = read.table(url, sep = '\t', header = TRUE, stringsAsFactors = F)
    print(paste0("Downloaded ", nrow(interactions), " interactions"))
  }else{
    load(from_cache_file)
  }
  
  if(!is.null(filter_databases)){
    filteredInteractions <- .filter_sources(interactions,databases = filter_databases)
  }else{
    filteredInteractions <- interactions
  }
  
  filteredInteractions$sources <- as.character(filteredInteractions$sources)
  filteredInteractions$nsources <-	
    unlist(lapply(strsplit(filteredInteractions$sources,split = ";"),length))
  
  
  return(filteredInteractions)
}



#' Imports from Omnipath webservice the interactions from Dorothea dataset
#'
#' Imports the dataset from: 
#' 'http://omnipathdb.org/interactions?datasets=tfregulons', 
#' which contains transcription factor (TF)-target interactions from DoRothEA
#' 'https://github.com/saezlab/DoRothEA'
#' 
#' @return data.frame
#' @export
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are removed.
#' @examples
#' interactions <- import_TFregulons_Interactions(filter_databases=c("tfact","ARACNe-GTEx"))
import_TFregulons_Interactions = function (from_cache_file=NULL,
                filter_databases = .get_interaction_databases()){
  url <- 
    'http://omnipathdb.org/interactions?datasets=tfregulons&fields=sources,tfregulons_level&genesymbols=1'
  if(is.null(from_cache_file)){
    interactions = read.table(url, sep = '\t', header = TRUE, stringsAsFactors = F)
    print(paste0("Downloaded ", nrow(interactions), " interactions"))
  }else{
    load(from_cache_file)
  }
  
  if(!is.null(filter_databases)){
    filteredInteractions <- .filter_sources(interactions,databases = filter_databases)
  }else{
    filteredInteractions <- interactions
  }
  
  filteredInteractions$sources <- as.character(filteredInteractions$sources)
  filteredInteractions$nsources <-	
    unlist(lapply(strsplit(filteredInteractions$sources,split = ";"),length))
  
  
  return(filteredInteractions)
}


#' Imports from Omnipath webservice the interactions from miRNAtarget dataset
#'
#' Imports the dataset from: 
#' 'http://omnipathdb.org/interactions?datasets=mirnatarget', 
#' which contains miRNA-mRNA and TF-miRNA interactions
#' 
#' @return data.frame
#' @export
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are removed.
#' @examples
#' interactions <- import_miRNAtarget_Interactions(filter_databases=c("miRTarBase","miRecords"))
import_miRNAtarget_Interactions = function (from_cache_file=NULL,
              filter_databases = .get_interaction_databases()){
  url <- 
    'http://omnipathdb.org/interactions?datasets=mirnatarget&fields=sources,references&genesymbols=1'
  if(is.null(from_cache_file)){
    interactions = read.table(url, sep = '\t', header = TRUE, stringsAsFactors = F)
    print(paste0("Downloaded ", nrow(interactions), " interactions"))
  }else{
    load(from_cache_file)
  }
  
  if(!is.null(filter_databases)){
    filteredInteractions <- .filter_sources(interactions,databases = filter_databases)
  }else{
    filteredInteractions <- interactions
  }
  
  filteredInteractions$sources <- as.character(filteredInteractions$sources)
  filteredInteractions$nsources <-	
    unlist(lapply(strsplit(filteredInteractions$sources,split = ";"),length))
  
  
  return(filteredInteractions)
}

#' Get the different interaction databases
#'
#' get the names of the databases from omnipath.org/interactions
#' @return character vector with the databases
#'
.get_interaction_databases = function(){
  url_interactions <- 'http://omnipathdb.org/interactions?datasets=omnipath,pathwayextra,kinaseextra,ligrecextra,tfregulons,mirnatarget&fields=sources'
  interactions = read.table(url_interactions, sep = '\t', header = TRUE,
                            stringsAsFactors = F)
  return(unique(unlist(strsplit(x = as.character(interactions$sources),
                                split = ";"))))
}

# filtering the PTMS and interactions which are reported in the import functions
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
  
  print(paste0("removed ",nInter-nInterPost, 
               " interactions during database filtering."))
  return(subsetInteractions)
}