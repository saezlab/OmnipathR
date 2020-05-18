########## ########## ########## ##########
########## PTMS                  ##########   
########## ########## ########## ##########

#' Import Omnipath post-translational modifications (PTMs)
#'
#' imports the PTMs database from \url{http://omnipathdb.org/ptms}
#'
#' @return A data frame containing the information about ptms
#' @export
#' @importFrom utils read.table
#' @param from_cache_file path to an earlier data file
#' @param filter_databases PTMs not reported in these databases are 
#' removed. See \code{\link{get_ptms_databases}} for more information
#' @param select_organism PTMs are available for human, mouse and rat. 
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @examples
#' ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"),
#'        select_organism=9606)
#'        
#' @seealso \code{\link{get_ptms_databases}, 
#'   \link{import_Omnipath_Interactions}}         
import_Omnipath_PTMS = function (from_cache_file=NULL,
    filter_databases = get_ptms_databases(),select_organism = 9606){

    url_ptms_common <- paste0(
        'http://omnipathdb.org/ptms/?',
        'fields=sources&fields=references&fields=curation_effort'
    )

    url_ptms <- organism_url(url_ptms_common, select_organism)
    
    if(is.null(from_cache_file)){
        ptms <- getURL(url_ptms, read.table, sep = '\t', header = TRUE,
            stringsAsFactors = FALSE)
        message("Downloaded ", nrow(ptms), " PTMs")
    } else {
        load(from_cache_file)
    }

    filteredPTMS <- filter_format_inter(ptms,filter_databases)
    return(filteredPTMS)
}

#' Get Post-translational modification (PTMs) databases
#'
#' get the names of the different databases available for ptms databases 
#' \url{http://omnipath.org/ptms}
#' 
#' @return character vector with the names of the PTMs databases
#' @export
#' @importFrom utils read.table
#' @examples
#' get_ptms_databases()
#' @seealso  \code{\link{import_Omnipath_PTMS}} 
get_ptms_databases = function(){
    url_ptms <- 'http://omnipathdb.org/ptms/?fields=sources'
    ptms <- getURL(url_ptms, read.table, sep = '\t', header = TRUE,
        stringsAsFactors = FALSE)
    return(unique(unlist(strsplit(x = as.character(ptms$sources),split = ";"))))
}

########## ########## ########## ##########
########## INTERACTIONS          ##########   
########## ########## ########## ##########

## Import Omnipath interaction Database. The new version of Ominpath contains 
## several different datastes. 

#' Import Omnipath interaction database
#'
#' imports the database from \url{http://omnipathdb.org/interactions}, which 
#' contains only interactions with references. These interactions are the 
#' original ones from the first Omnipath version. 
#' 
#' @return A dataframe containing information about protein-protein interactions
#' @export
#' @importFrom utils read.table
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are 
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param select_organism Interactions are available for human, mouse and rat. 
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' 
#' @examples
#' interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3"),
#'                select_organism = 9606)
#'                
#' @seealso \code{\link{get_interaction_databases}, 
#'   \link{import_AllInteractions}}                 
import_Omnipath_Interactions = function (from_cache_file=NULL,
    filter_databases = get_interaction_databases(),select_organism = 9606){

    url_interactions_common <- paste0(
        'http://omnipathdb.org/interactions?',
        'fields=sources&fields=references&fields=curation_effort'
    )
    
    url_interactions <- organism_url(url_interactions_common, select_organism)
    
    if(is.null(from_cache_file)){
        interactions <- getURL(url_interactions,read.table,sep = '\t', 
            header = TRUE, stringsAsFactors = FALSE)
        message("Downloaded ", nrow(interactions), " interactions")
    } else {
        load(from_cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions,filter_databases)
    return(filteredInteractions)
}    
    

#' Imports from Omnipath webservice the interactions from 
#' Pathwayextra dataset
#'
#' Imports the dataset from: 
#' \url{http://omnipathdb.org/interactions?datasets=pathwayextra}, 
#' which contains activity flow interactions without literature reference
#' 
#' @return A dataframe containing activity flow interactions between proteins
#' without literature reference
#' @export
#' @importFrom utils read.table
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are 
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param select_organism Interactions are available for human, mouse and rat. 
#' Choose one of those: 9606 human (default), 10116 rat or 10090 Mouse
#' @examples
#' interactions <- 
#'     import_PathwayExtra_Interactions(filter_databases=c("BioGRID","IntAct"),
#'      select_organism = 9606)
#' @seealso \code{\link{get_interaction_databases}, 
#'   \link{import_AllInteractions}}
import_PathwayExtra_Interactions = function (from_cache_file=NULL,
    filter_databases = get_interaction_databases(),select_organism = 9606){

    url_pathwayextra_common <- 
        paste0('http://omnipathdb.org/interactions?datasets=',
        'pathwayextra&fields=sources') 

    url_pathwayextra <- organism_url(url_pathwayextra_common, select_organism)

    if(is.null(from_cache_file)){
        interactions <- getURL(url_pathwayextra,read.table,sep = '\t', 
            header = TRUE, stringsAsFactors = FALSE)
        message("Downloaded ", nrow(interactions), " interactions")
    } else {
        load(from_cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions,filter_databases)
    return(filteredInteractions)
}

#' Imports from Omnipath webservice the interactions from 
#' kinaseextra dataset
#'
#' Imports the dataset from: 
#' \url{http://omnipathdb.org/interactions?datasets=kinaseextra}, 
#' which contains enzyme-substrate interactions without literature reference
#' 
#' @return A dataframe containing enzyme-substrate interactions without 
#' literature reference
#' @export
#' @importFrom utils read.table
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are 
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param select_organism Interactions are available for human, mouse and rat. 
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' 
#' @examples
#' interactions <- 
#'    import_KinaseExtra_Interactions(filter_databases=c("PhosphoPoint",
#'    "PhosphoSite"), select_organism = 9606)
#' @seealso \code{\link{get_interaction_databases}, 
#'   \link{import_AllInteractions}}          
import_KinaseExtra_Interactions = function (from_cache_file=NULL,
    filter_databases = get_interaction_databases(),select_organism = 9606){

    url_kinaseextra_common <- 
        paste0('http://omnipathdb.org/interactions?datasets=',
            'kinaseextra&fields=sources') 
    
    url_kinaseextra <- organism_url(url_kinaseextra_common, select_organism)
    
    if(is.null(from_cache_file)){
        interactions <- getURL(url_kinaseextra,read.table,sep = '\t', 
            header = TRUE, stringsAsFactors = FALSE)
        message("Downloaded ", nrow(interactions), " interactions")
    } else {
        load(from_cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions,filter_databases)
    return(filteredInteractions)
}


#' Imports from Omnipath webservice the interactions from 
#' ligrecextra dataset
#'
#' Imports the dataset from: 
#' \url{http://omnipathdb.org/interactions?datasets=ligrecextra}, 
#' which contains ligand-receptor interactions without literature reference
#' 
#' @return A dataframe containing ligand-receptor interaction without literature
#' reference
#' @export
#' @importFrom utils read.table
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are 
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param select_organism Interactions are available for human, mouse and rat. 
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @examples
#' interactions <- import_LigrecExtra_Interactions(filter_databases=c("HPRD",
#'       "Guide2Pharma"), select_organism=9606)
#' @seealso \code{\link{get_interaction_databases}, 
#'   \link{import_AllInteractions}}
import_LigrecExtra_Interactions = function (from_cache_file=NULL,
    filter_databases = get_interaction_databases(), select_organism=9606){

    url_ligrecextra_common <- 
        paste0('http://omnipathdb.org/interactions?datasets=',
        'ligrecextra&fields=sources') 

    url_ligrecextra <- organism_url(url_ligrecextra_common, select_organism)  

    if(is.null(from_cache_file)){
        interactions <- getURL(url_ligrecextra, read.table, sep = '\t', 
            header = TRUE, stringsAsFactors = FALSE)
        message("Downloaded ", nrow(interactions), " interactions")
    } else {
        load(from_cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions,filter_databases)
    return(filteredInteractions)
}


#' Imports from Omnipath webservice the interactions from 
#' Dorothea dataset
#'
#' Imports the dataset from: 
#' \url{http://omnipathdb.org/interactions?datasets=tfregulons} 
#' which contains transcription factor (TF)-target interactions from DoRothEA
#' \url{https://github.com/saezlab/DoRothEA}
#' 
#' @return A dataframe containing TF-target interactions from DoRothEA
#' @export
#' @importFrom utils read.table
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are 
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param select_organism Interactions are available for human, mouse and rat. 
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @param confidence_level Vector detailing the confidence levels of the 
#' interactions to be downloaded. In dorothea, every TF-target interaction has a 
#' confidence score ranging from A to E, being A the most reliable interactions.
#' By default we take A and B level interactions (\code{c(A,B)}). It is to note 
#' that E interactions are not available in OmnipathR.   
#' @examples
#' interactions <- import_TFregulons_Interactions(filter_databases=c("DoRothEA_A",
#'     "ARACNe-GTEx_DoRothEA"), select_organism=9606)
#' @seealso \code{\link{get_interaction_databases}, 
#'   \link{import_AllInteractions}}
import_TFregulons_Interactions = function (from_cache_file=NULL,
    filter_databases = get_interaction_databases(),select_organism=9606, 
    confidence_level = c('A','B')){

    url_tfregulons_common <- 
        paste0('http://omnipathdb.org/interactions?datasets=tfregulons&',
        'fields=sources,tfregulons_level')

    url_tfregulons <- organism_url(url_tfregulons_common, select_organism)  

    confidence_level <- as.vector(confidence_level)
    
    if (length(confidence_level) > 4 | length(confidence_level) < 1){
        stop("The confidence levels vector is not correct")
    } else {
        if (all(confidence_level %in% c("A","B","C","D"))){
            url_tfregulons <- paste0(url_tfregulons, "&tfregulons_levels=",
                 paste0(confidence_level,collapse = ","))    
        } else {
            stop("Your confident levels are not correct, they should range from 
                A to D.")
        }
    }
        
    if(is.null(from_cache_file)){
        interactions <- getURL(url_tfregulons, read.table, sep = '\t', 
            header = TRUE, stringsAsFactors = FALSE)
        message("Downloaded ", nrow(interactions), " interactions")
    } else {
        load(from_cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions,filter_databases)
    return(filteredInteractions)
}


#' Imports from Omnipath webservice the interactions from 
#' miRNAtarget dataset
#'
#' Imports the dataset from: 
#' \url{http://omnipathdb.org/interactions?datasets=mirnatarget}, 
#' which contains miRNA-mRNA and TF-miRNA interactions
#' 
#' @return A dataframe containing miRNA-mRNA and TF-miRNA interactions
#' @export
#' @importFrom utils read.table
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are 
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @examples
#' interactions <- 
#'   import_miRNAtarget_Interactions(filter_databases=c("miRTarBase",
#'   "miRecords"))
#' @seealso \code{\link{get_interaction_databases}, 
#'   \link{import_AllInteractions}}
import_miRNAtarget_Interactions = function (from_cache_file=NULL,
    filter_databases = get_interaction_databases()){

    url <- paste0('http://omnipathdb.org/interactions?datasets=mirnatarget',
        '&fields=sources,references&genesymbols=1') 

    if(is.null(from_cache_file)){
        interactions <- getURL(url, read.table, sep = '\t', header = TRUE, 
            stringsAsFactors = FALSE)
        message("Downloaded ", nrow(interactions), " interactions")
    } else {
        load(from_cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions,filter_databases)
    return(filteredInteractions)
}


#' Imports from Omnipath webservice all the available interactions 
#' from the different datasets
#'
#' Imports the dataset from: 
#' \url{http://omnipathdb.org/interactions?datasets=omnipath,pathwayextra,kinaseextra,ligrecextra,tfregulons,mirnatarget&fields=sources,references&genesymbols=1}, 
#' which contains all the different interactions available in the webserver:
#' 
#' omnipath: the OmniPath data as defined in the paper, an arbitrary optimum 
#' between coverage and quality
#' pathwayextra: activity flow interactions without literature reference
#' kinaseextra: enzyme-substrate interactions without literature reference
#' ligrecextra: ligand-receptor interactions without literature reference
#' tfregulons: transcription factor (TF)-target interactions from DoRothEA
#' mirnatarget: miRNA-mRNA and TF-miRNA interactions
#' 
#' @return A dataframe containing all the datasets in the interactions query
#' @export
#' @importFrom utils read.table
#' @param from_cache_file path to an earlier data file
#' @param filter_databases interactions not reported in these databases are 
#' removed. See \code{\link{get_interaction_databases}} for more information.
#' @param select_organism Interactions are available for human, mouse and rat. 
#' Choose among: 9606 human (default), 10116 rat and 10090 Mouse
#' @examples
#' interactions <- import_AllInteractions(filter_databases=c("HPRD","BioGRID"),
#'     select_organism = 9606)
#' @seealso \code{\link{get_interaction_databases}}
import_AllInteractions = function (from_cache_file=NULL,
    filter_databases = get_interaction_databases(),select_organism = 9606){

    url_allinteractions_common <- 
        paste0('http://omnipathdb.org/interactions?datasets=omnipath',
            ',pathwayextra,kinaseextra,ligrecextra,tfregulons,mirnatarget', 
            '&fields=sources,references')

    url_allinteractions <- organism_url(url_allinteractions_common, 
        select_organism)    

    if(is.null(from_cache_file)){
        interactions <- getURL(url_allinteractions, read.table, sep = '\t', 
            header = TRUE, stringsAsFactors = FALSE)
        message("Downloaded ", nrow(interactions), " interactions")
    } else {
        load(from_cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions,filter_databases)
    return(filteredInteractions)
}

#' Get the different interaction databases
#'
#' get the names of the databases from \url{http://omnipath.org/interactions}
#' 
#' @return character vector with the names of the interaction databases
#' @export
#' @importFrom utils read.table
#' @examples
#' get_interaction_databases()
#' @seealso \code{\link{import_AllInteractions}, 
#' \link{import_Omnipath_Interactions}, \link{import_PathwayExtra_Interactions},
#' \link{import_KinaseExtra_Interactions}, 
#' \link{import_LigrecExtra_Interactions},
#' \link{import_miRNAtarget_Interactions}, 
#' \link{import_TFregulons_Interactions}}
get_interaction_databases = function(){
    url_interactions <- paste0('http://omnipathdb.org/interactions?',
        'datasets=omnipath,pathwayextra,kinaseextra,ligrecextra',
        ',tfregulons,mirnatarget&fields=sources')
    interactions <- getURL(url_interactions, read.table, sep = '\t', 
        header = TRUE,stringsAsFactors = FALSE)
    return(unique(unlist(strsplit(x = as.character(interactions$sources),
        split = ";"))))
}

########## ########## ########## ##########
########## Complexes             ##########   
########## ########## ########## ##########

#' Import Omnipath Complexes
#'
#' imports the complexes stored in Omnipath database from 
#' \url{http://omnipathdb.org/complexes}
#'
#' @return A dataframe containing information about complexes
#' @export
#' @importFrom utils read.csv
#' @param from_cache_file path to an earlier data file
#' @param filter_databases complexes not reported in these databases are 
#' removed. See \code{\link{get_complexes_databases}} for more information.
#' @examples
#' complexes = import_Omnipath_complexes(filter_databases=c("CORUM", "hu.MAP"))
#' @seealso \code{\link{get_complexes_databases}}
import_Omnipath_complexes = function (from_cache_file=NULL,
    filter_databases = get_complexes_databases()){

    url_complexes <- 'http://omnipathdb.org/complexes?&fields=sources'

    if(is.null(from_cache_file)){
        complexes <- getURL(url_complexes, read.csv, sep = '\t', header = TRUE,
            stringsAsFactors = FALSE)
        message("Downloaded ", nrow(complexes), " complexes")
    } else {
        load(from_cache_file)
    }

    filteredcomplexes <- filter_format_inter(complexes,filter_databases)
    return(filteredcomplexes)
}


#' Get the different complexes databases integrated in Omnipath
#'
#' get the names of the databases from \url{http://omnipath.org/complexes}
#' @return character vector with the names of the databases
#' @export
#' @importFrom utils read.csv
#' @examples
#' get_complexes_databases()
#' @seealso \code{\link{import_Omnipath_complexes}}
get_complexes_databases = function(){
    url_complexes <- 'http://omnipathdb.org/complexes?&fields=sources'
    complexes <- getURL(url_complexes, read.csv, sep = '\t', header = TRUE,
        stringsAsFactors = FALSE)
    return(unique(unlist(strsplit(x = as.character(complexes$sources),
        split = ";"))))
}

########## ########## ########## ##########
########## Annotations           ##########
########## ########## ########## ##########

#' Import Omnipath Annotations
#'
#' imports the annotations stored in Omnipath database from
#' \url{http://omnipathdb.org/annotations}
#'
#' @return A data.frame containing different gene/complex annotations
#' @export
#' @importFrom utils read.csv
#' @param from_cache_file Path to an earlier data file
#' @param select_genes Vector containing the genes or proteins for whom
#' annotations will be retrieved (UniProt IDs or HGNC Gene Symbols or
#' miRBase IDs). It is also possible to donwload annotations for protein
#' complexes. To do so, write "COMPLEX:" right before the genesymbols of
#' the genes integrating the complex. Check the vignette for examples.
#' @param filter_databases Load the annotations only from these databases.
#' See \code{\link{get_annotation_databases}} for possible values.
#' @param force_full_download Force the download of the entire annotations
#' dataset. This is disabled by default because the size of this data is
#' around 1GB. We recommend to retrieve the annotations for a set of proteins
#' or only from a few databases, depending on your interest.
#' @param ... Additional arguments. 
#' @examples
#' annotations = import_Omnipath_annotations(select_genes=c("TP53","LMNA"),
#'      filter_databases=c("HPA_subcellular"))
#' @seealso \code{\link{get_annotation_databases}}
import_Omnipath_annotations = function (from_cache_file=NULL,
    select_genes = NULL, filter_databases = NULL,
    force_full_download = FALSE, ...){

    if(
        !force_full_download &&
        is.null(from_cache_file) &&
        is.null(select_genes) &&
        is.null(filter_databases)
    ){
        
        stop(
            paste(
                'Downloading the entire annotations database is not allowed',
                'by default because of its huge size (>1GB). If you really',
                'want to do this use the `force_full_download` parameter.',
                'However we recommend to query a set of proteins or a few',
                'databases, depending on your interest.'
            )
        )
        
    }

    url_annotations <- 'http://omnipathdb.org/annotations?'
    
    annotations <- NULL
    
    if(!is.null(from_cache_file)){
        
        annotations <-
            filter_sources_annotations(
                load(from_cache_file),
                databases = filter_databases
            )
        
        if(is.null(annotations) && !is.null(filter_databases)){
            warning(
                paste(
                    'Empty result from cache file.',
                    'Might be the cache file does not contain data from',
                    'the requested databases? Trying without cache.'
                )
            )
        }
        
    }
    
    if(is.null(annotations)){
        
        databases_part <- ''
        
        if(!is.null(filter_databases)){
            
            all_databases <- get_annotation_databases()
            unknown_databases <- setdiff(filter_databases, all_databases)
            
            if(length(unknown_databases) != 0){
                
                warning(
                    sprintf(
                        paste(
                            'The following databases are not available: %s.',
                            'Check the database names for spelling mistakes.'
                        ),
                        paste0(unknown_databases, collapse = ', ')
                    )
                )
                
            }
            
            databases_part <- sprintf(
                'databases=%s',
                paste0(filter_databases, collapse = ',')
            )
            
        }
        
        if(length(select_genes) > 600){
            
            annotations <- list()
            
            proteins_chunks <-
                split(
                    select_genes,
                    cut(
                        seq_along(select_genes),
                        breaks = length(select_genes) / 500,
                        labels = FALSE
                    )
                )
            
            cat('Downloading ')
            
            for(proteins_chunk in proteins_chunks){
                
                cat('.')
                
                annotations <- append(
                    annotations,
                    list(
                        import_Omnipath_annotations(
                            select_genes = proteins_chunk,
                            filter_databases = filter_databases,
                            recursive_call = TRUE
                        )
                    )
                )
                
            }
            
            cat(' ready.\n')
            
            annotations <- do.call(rbind, annotations)
            
        }else{
            
            proteins_part <- `if`(
                is.null(select_genes),
                '',
                sprintf(
                    '&proteins=%s',
                    paste0(select_genes, collapse = ',')
                )
            )
            
            url_annotations <- paste0(
                url_annotations,
                databases_part,
                proteins_part
            )
            
            annotations <- getURL(
                url_annotations,
                read.csv,
                sep = '\t',
                header = TRUE,
                stringsAsFactors = FALSE
            )
            
        }
        
        args <- list(...)
        
        if(!'recursive_call' %in% names(args) || !args$recursive_call){
            
            message(sprintf('Downloaded %d annotations.', nrow(annotations)))
            
        }
        
    }
    
    return(annotations)
    
}


#' Get the different annotation databases integrated in Omnipath
#'
#' get the names of the databases from \url{http://omnipath.org/annotation}
#' 
#' @return character vector with the names of the annotation databases
#' @export
#' @examples
#' get_annotation_databases()
#' @seealso \code{\link{import_Omnipath_annotations}}
get_annotation_databases = function(){
    url_annotations <- 'http://omnipathdb.org/annotations_summary'
    annotations <- getURL(url_annotations, read.table, sep = '\t', 
        header = TRUE,stringsAsFactors = FALSE)

    annotations_db <- unique(annotations$source)
    return(annotations_db)
}


########## ########## ########## ##########
########## Intercell             ##########   
########## ########## ########## ##########

#' Import Omnipath Intercell Data
#'
#' imports the intercell data stored in Omnipath database from 
#' \url{http://omnipathdb.org/intercell}. Intercell provides 
#' information on the roles in inter-cellular signaling. E.g. if a protein is 
#' a ligand, a receptor, an extracellular matrix (ECM) component, etc.
#'
#' @return A dataframe cotaining information about roles in inter-cellular
#' signaling. 
#' @export
#' @importFrom utils read.csv
#' @param from_cache_file path to an earlier data file
#' @param select_categories vector containing the categories to be retrieved.
#' All the genes belonging to those categories will be returned. For furter 
#' information about the categories see \code{\link{get_intercell_categories}} 
#' @param select_classes vector containing the main classes to be retrieved.
#' All the genes belonging to those classes will be returned. For furter 
#' information about the main classes see \code{\link{get_intercell_classes}} 
#' @examples
#' intercell = import_Omnipath_intercell(select_categories=c("ecm"))
#' @seealso \code{\link{get_intercell_categories}} 
import_Omnipath_intercell = function (from_cache_file=NULL,
    select_categories = get_intercell_categories(), 
    select_classes = get_intercell_classes()){

    url_intercell <- 'http://omnipathdb.org/intercell'

    if(is.null(from_cache_file)){
        intercell <- getURL(url_intercell, read.csv, sep = '\t', header = TRUE,
            stringsAsFactors = FALSE)
        message("Downloaded ", nrow(intercell), " intercell records")
    } else {
        load(from_cache_file)
    }

    if (!(all(select_categories == get_intercell_categories())) | 
        !(all(select_classes == get_intercell_classes()))){
            filteredintercell <- 
               filter_intercell(intercell,select_categories,select_classes)
            return(filteredintercell) 
    } else {
        return(intercell)
    }
}

#' Imports an intercellular network combining annotations and interactions
#'
#' Imports an intercellular network by mapping intercellular annotations 
#' and protein interactions. It first imports the PPI interactions from the
#' different datasets here described. Then, it takes proteins with the 
#' intercellular roles defined by the user. Some proteins should be defined 
#' as transmiters (eg. ligand) and other as receivers (receptor). We find the 
#' interactions which source is a transmiter and its target a receiver. 
#'
#' @return A dataframe containing information about protein-protein 
#' interactions and the inter-cellular roles of the protiens involved in those
#' interactions. 
#' @export
#' @importFrom utils read.csv
#' @param from_cache_file path to an earlier data file
#' @param filter_databases vector containing interactions databases. 
#' Interactions not reported in these databases are removed. 
#' See \code{\link{get_interaction_databases}} for more information.
#' @param classes_source A list containing two vectors. The first one with 
#' the main classes to be considered as transmiters and the second with the 
#' main classes to be considered as receivers. For furter information
#' about the main classes see \code{\link{get_intercell_classes}} 
#' @examples
#' intercellNetwork <- import_intercell_network(
#' classes_source = list(transmiters=c('ligand'),receivers=c('receptor')))
#' @seealso \code{\link{get_intercell_categories}} 
import_intercell_network = function (from_cache_file=NULL,
    filter_databases = get_interaction_databases(), 
    classes_source = list(transmiters=c('ligand'),receivers=c('receptor'))) {

    mainclass <- genesymbol <- NULL
    AllClasses <- unlist(classes_source)
    
    if (!all(AllClasses %in% get_intercell_classes())){
        stop("Some all the classes are not correct. 
            Check get_intercell_classes()")
    }
    
    url_allinteractions_common <- 
        paste0('http://omnipathdb.org/interactions?datasets=omnipath',
            ',pathwayextra,kinaseextra,ligrecextra', 
            '&fields=sources,references')

    url_allinteractions <- organism_url(url_allinteractions_common, 9606)    
    
    if(is.null(from_cache_file)){
        interactions <- getURL(url_allinteractions, read.table, sep = '\t', 
            header = TRUE, stringsAsFactors = FALSE)
        message("Downloaded ", nrow(interactions), " interactions")
    } else {
        load(from_cache_file)
    }

    filteredInteractions <- filter_format_inter(interactions,filter_databases)
    
    intercellAnnotations <- 
        import_Omnipath_intercell(select_classes = AllClasses)

    genesTransmiters <- intercellAnnotations %>%
        dplyr::filter(mainclass %in% classes_source$transmiters) %>%
        dplyr::distinct(genesymbol,mainclass)
    genesReceivers <- intercellAnnotations %>%
        dplyr::filter(mainclass %in% classes_source$receivers) %>%
        dplyr::distinct(genesymbol,mainclass)
    
    intercelNetwork <- 
        dplyr::inner_join(filteredInteractions, genesTransmiters, 
            by=c("source_genesymbol"="genesymbol")) %>% 
        dplyr::rename(class_source = mainclass) %>%
        dplyr::inner_join(genesReceivers, 
            by=c("target_genesymbol"="genesymbol")) %>% 
        dplyr::rename(class_target = mainclass)
        
    return(intercelNetwork)
}

#' Get the different intercell categories described in Omnipath
#'
#' get the names of the categories from \url{http://omnipath.org/intercell}
#' @return character vector with the different intercell categories
#' @export
#' @importFrom utils read.csv
#' @examples
#' get_intercell_categories()
#' @seealso \code{\link{import_Omnipath_intercell}, 
#' \link{get_intercell_classes}}
get_intercell_categories = function(){

    url_intercell <- 'http://omnipathdb.org/intercell_summary'
    intercell <- getURL(url_intercell, read.csv, sep = '\t', header = TRUE,
        stringsAsFactors = FALSE)

    return(unique(intercell$category))
}

#' Get the different intercell main classes described in Omnipath
#'
#' get the names of the main classes from \url{http://omnipath.org/intercell}
#' @return character vector with the different intercell main classes
#' @export
#' @importFrom utils read.csv
#' @examples
#' get_intercell_classes()
#' @seealso \code{\link{import_Omnipath_intercell},
#' \link{get_intercell_categories}}
get_intercell_classes = function(){

    url_intercell <- 'http://omnipathdb.org/intercell_summary'
    intercell <- getURL(url_intercell, read.csv, sep = '\t', header = TRUE,
        stringsAsFactors = FALSE)

    return(unique(intercell$mainclass))
}

########## ########## ########## ##########
########## SOURCE FILTERING      ##########   
########## ########## ########## ##########
## Non exported functions (package internal functions) to filter PTMs,
## interactions, complexes and annotations according to the databases passed
## to the main functions

## Filtering Interactions, PTMs and complexes
filter_sources = function(interactions, databases){

    nInter = nrow(interactions)

    subsetInteractions <- 
        interactions[which(unlist(lapply(strsplit(interactions$sources,";"),
        function(x){any(x %in% databases)}))),]

    nInterPost = nrow(subsetInteractions)

    message("removed ",nInter-nInterPost,
        " interactions during database filtering.")
    return(subsetInteractions)
}


## Filtering Annotations
filter_sources_annotations = function(annotations, databases){
## takes annotations and removes those which are
## not reported by the given databases.

    if(is.null(databases)){
        return(annotations)
    }

    nAnnot = nrow(annotations)
    subsetAnnotations <- dplyr::filter(annotations, source %in% databases)
    nAnnotPost = nrow(subsetAnnotations)

    message("removed ",nAnnot-nAnnotPost,
        " annotations during database filtering.")

    if (nAnnotPost > 0){
        return(subsetAnnotations)
    } else {
        return(NULL)
    }
}

## Filtering intercell records according to the categories and/or classes 
## selected
filter_intercell = function(intercell, categories, classes){
## takes intercell removes and removes those not reported by the given 
## databases
    nIntercell = nrow(intercell)
    subsetIntercell <- dplyr::filter(intercell, .data$category %in% categories)    
    subsetIntercell <- 
        dplyr::filter(subsetIntercell, .data$mainclass %in% classes)    
    
    nIntercellPost = nrow(subsetIntercell)

    message("removed ",nIntercell-nIntercellPost, 
        " intercell records during category/class filtering.")

    if (nIntercellPost > 0){
        return(subsetIntercell)
    } else {
        return(NULL)
    }
}

########## ########## ########## ##########
########## Resource Queries      ##########   
########## ########## ########## ##########
## This function is convenient for appropriate resource retrieval. Following:
## http://bioconductor.org/developers/how-to/web-query/
## It tries to retrieve the resource one or several times before failing.
getURL <- function(URL, FUN, ..., N.TRIES=1L) {
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

    while (N.TRIES > 0L) {
        result <- tryCatch(FUN(URL, ...), error=identity)
        if (!inherits(result, "error"))
            break
            N.TRIES <- N.TRIES - 1L
        }

    if (N.TRIES == 0L) {
        stop("'getURL()' failed:",
            "\n  URL: ", URL,
            "\n  error: ", conditionMessage(result))
    }

    return(result) 
}

########## ########## ########## ##########
########## Queries Format        ##########
########## ########## ########## ##########
## This function format de url for the queries to the Omnipath webserver 
## according to the selected organism
organism_url <- function(url,organism){

    if (organism %in% c(9606, 10116, 10090)){
        if (organism == 9606){
            url_final <- paste0(url,'&genesymbols=1') 
        } else {
            if (organism == 10116){
                url_final <- paste0(url,'&genesymbols=1&organisms=10116') 
            }
            if (organism == 10090){
                url_final <- paste0(url,'&genesymbols=1&organisms=10090') 
            }
        }     
    } else {
        stop("The selected organism is not correct")
    }
    return(url_final)
}

########## ########## ########## ########## ##########
########## Format and filter of interactions #########
########## ########## ########## ########## ##########
## This function calls to the filtering functions and gives format
## to the data frames containing the interactions (For instance, it generates
## a new field with the number of sources/references reporting a given
## interaction)

filter_format_inter <- function(interaction_df,databases){

    if(!is.null(databases)){
        interaction_df <- 
            filter_sources(interaction_df,databases = databases)
    } else {
        interaction_df <- interaction_df
    }

    if (nrow(interaction_df)==0){
        stop("Try another database filtering: No records found")
    }

    if ("residue_offset" %in% colnames(interaction_df)){
        interaction_df$residue_offset <- 
            as.character(as.numeric(interaction_df$residue_offset))      
    }

    interaction_df$sources <- as.character(interaction_df$sources)
    interaction_df$nsources <-
        unlist(lapply(strsplit(interaction_df$sources,split = ";"),length))

    if ("references" %in% colnames(interaction_df)){
        interaction_df$references <- as.character(interaction_df$references)
        ## we remove references mentioned multiple times:
        interaction_df$references <- 
            unlist(lapply(strsplit(interaction_df$references,split = ";"),
            function(x)paste(unique(x),collapse=";")))
        interaction_df$nrefs <- 
            unlist(lapply(strsplit(interaction_df$references,split = ";"),
            length))
    }
    return(interaction_df)        
}
