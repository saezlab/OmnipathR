#' print interactions
#'
#' prints the interactions in a nice format
#'
#' @param interDF data.frame with the interacgtions
#' @param writeRefs [FALSE] writes also the PubMed IDs
#' @export
#' @return NULL
#' @examples
#' ptms = import_Omnipath_PTMS()
#' print_interactions(head(ptms))
#' print_interactions(tail(ptms),writeRefs=T)
#' print_interactions(dplyr::filter(ptms,enzyme_genesymbol=="MAP2K1",substrate_genesymbol=="MAPK3"))
print_interactions = function(interDF,writeRefs=F){

	if(nrow(interDF)==0) {
		print("no interactions")
		return(invisible(NULL))
	}

	if("enzyme" %in% colnames(interDF)){  #PTMS
		interDF = interDF[order(interDF$nrefs,interDF$nsources,decreasing = T),]
		interDF$enzyme = paste0(interDF$enzyme_genesymbol, " (", interDF$enzyme ,")")
		interDF$substrate = paste0(interDF$substrate_genesymbol,"_",interDF$residue_type,interDF$residue_offset, " (", interDF$substrate ,")")


		signs = ifelse(interDF$is_stimulation==1,
					   ifelse(interDF$is_inhibition==1,"(+/-)","( + )"),   # stimulation: T, is also inhibition?
					   ifelse(interDF$is_inhibition==1,"( - )","( ? )"))     # stimulation: F, is inhibition?
		interDF$interaction = paste0("==", signs,"==>")
		if(writeRefs){
			print(interDF[,c('enzyme',"interaction","substrate","modification","nsources","nrefs","references")])
		}else{
			print(interDF[,c('enzyme',"interaction","substrate","modification","nsources","nrefs")])
		}


	}else{ # print interactons

		interDF = interDF[order(interDF$nrefs,interDF$nsources,decreasing = T),]
		interDF$source = paste0(interDF$source_genesymbol, " (", interDF$source ,")")
		interDF$target = paste0(interDF$target_genesymbol, " (", interDF$target ,")")


		signs = ifelse(interDF$is_stimulation==1,
					   ifelse(interDF$is_inhibition==1,"(+/-)","( + )"),   # stimulation: T, is also inhibition?
					   ifelse(interDF$is_inhibition==1,"( - )","( ? )"))     # stimulation: F, is inhibition?

		direction = ifelse(interDF$is_directed==1, ">","")
		interDF$interaction = paste0("==", signs,"==",direction)

		if(writeRefs){
			print(interDF[,c('source',"interaction","target","nsources","nrefs","references")])
		}else{
			print(interDF[,c('source',"interaction","target","nsources","nrefs")])
		}



	}
}

#' print paths given by edge sequence
#'
#' prints the interactions in the path in a nice format
#'
#' @param edgeSeq edge sequence
#' @param G igraph object (from ptms or interactions)
#' @export
#' @return NULL
#' @examples
#' printPath_es(shortest_paths(OPI_g,from = "TYRO3",to = "STAT3", output = 'epath')$epath[[1]],OPI_g)
printPath_es = function(edgeSeq,G){

	# outDF = data.frame()
	if(length(edgeSeq)==0) {
		cat("\nempty path\n")
		return(NULL)
	}

	signs = ifelse(edgeSeq$is_stimulation==1,
				   ifelse(edgeSeq$is_inhibition==1,"(+/-)","( + )"),   # stimulation: T, is also inhibition?
				   ifelse(edgeSeq$is_inhibition==1,"( - )","( ? )"))     # stimulation: F, is inhibition?
	interaction = paste0("==", signs,"==>")

	if(! is.null(edgeSeq$residue_type)){
		edgeSeq$residue_type

		df = data.frame(source = paste(tail_of(G, edgeSeq)$name," (", tail_of(G, edgeSeq)$up_ids,")",sep = ""),
						interaction = interaction,
						target = paste(
							paste0(head_of(G, edgeSeq)$name, "_", edgeSeq$residue_type , edgeSeq$residue_offset),
									   " (", head_of(G, edgeSeq)$up_ids,")",sep = ""),
						nsources = edgeSeq$nsources,
						nrefs = edgeSeq$nrefs)
	}else{

		df = data.frame(source = paste(tail_of(G, edgeSeq)$name," (", tail_of(G, edgeSeq)$up_ids,")",sep = ""),
						interaction = interaction,
						target = paste(head_of(G, edgeSeq)$name," (", head_of(G, edgeSeq)$up_ids,")",sep = ""),
						nsources = edgeSeq$nsources,
						nrefs = edgeSeq$nrefs)
	}

	print(df)
}



# convert vertex sequence to named sequence to find unique
unique_nodeSeq = function(nodeSeq_list){
	# takes a list of nodeSequences converts them to names and takes the unique paths
	name_path = list()
	for(i in 1 : length(nodeSeq_list)){

		path1 = nodeSeq_list[[i]]

		name_seq = c()
		for( j in 1:length(path1)){
			name_seq = c(name_seq, path1[j]$name)
		}
		name_path[[i]] = name_seq
	}
	unique(name_path)
}

#unique_nodeSeq(nodeSeq)


####


#' print paths given by node sequence
#'
#' prints the interactions in the path in a nice format
#'
#' @param nodeSeq node sequence
#' @param G igraph object (from ptms or interactions)
#' @export
#' @return NULL
#' @examples
#' interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite", "Signor"))
#' OPI_g = interaction_graph(interactions = interactions )
#' printPath_vs(all_shortest_paths(OPI_g,from = "DYRK2",to = "MAPKAPK2")$res,OPI_g)
#' printPath_vs(shortest_paths(OPI_g,from = "TYRO3",to = "STAT3")$vpath,OPI_g)
#'
#' ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))
#' ptms_g = ptms_graph(ptms)
#' printPath_vs(all_shortest_paths(ptms_g,from = "SRC",to = "STAT1")$res,ptms_g)
printPath_vs = function(nodeSeq,G){

	if(length(nodeSeq)==0){
		print("empty path")
		return(invisible(NULL))
	}
	nodeSeq_names = unique_nodeSeq(nodeSeq)
	#browser()
	for(i in 1:length(nodeSeq_names)){
		print(paste0("pathway ", i, ": ", paste(nodeSeq_names[[i]],collapse = " -> ")))
		edgeSet = c()
		for(j in 2:length(nodeSeq_names[[i]])){
			edgeSet = c(edgeSet, E(G)[nodeSeq_names[[i]][[j-1]]  %->% nodeSeq_names[[i]][[j]]])
		}
		printPath_es(E(G)[edgeSet],G)
		#
		# for(j in 2:length(nodeSeq_names[[i]])){
		#
		# 	printPath_es(E(G)[nodeSeq_names[[i]][[j-1]]  %->% nodeSeq_names[[i]][[j]]],G)
		#
		# }
	}
}

#printPath_vs(nodeSeq,G)
