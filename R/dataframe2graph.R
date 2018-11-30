
#' PTMS graph
#'
#' transforms the ptmds interactions data.frame to igraph object
#'
#' @return igraph object
#' @export
#' @param interactions data.frame created by \code{\link{import_Omnipath_PTMS}}
#' @examples
#' ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))
#' ptms_g = ptms_graph(ptms = ptms )
ptms_graph <- function(ptms){
	# This is a gene_name based conversion to igraph, i.e. the vertices are identified
	# by genenames, and not by uniprot IDs.
	# This might cause issue when a gene name encodes multiple uniprot IDs.


	library(dplyr)

	# This does not preserve the p-Site info
	# edges = raw %>% select(c(enzyme_genesymbol, substrate_genesymbol)) %>% distinct()

	# keep only edge attributes
	edges = ptms %>% select(- c(enzyme, substrate))

	# lets try to do it at once
	# edgesB = raw %>% select(- c(enzyme, substrate)) %>%  # remove UniprotID
	#	group_by(enzyme_genesymbol, substrate_genesymbol) %>%
	#	summarise("residues" = paste(paste(residue_type, residue_offset, ifelse(is_stimulation=="1","+","?"),ifelse(is_inhibition=="1","-","?"),sep="_"),collapse=",")) %>% ungroup()

	# build vertices: gene_names and gene_uniprotIDs
	nodesA = select(ptms, c(enzyme_genesymbol, enzyme))
	nodesB = select(ptms, c(substrate_genesymbol, substrate))
	colnames(nodesA) = colnames(nodesB) = c("genesymbol", "up_id")
	nodes = rbind(nodesA,nodesB)
	nodes = unique(nodes)
	nodes = nodes %>% group_by(genesymbol) %>% summarise("up_ids" = paste0(up_id,collapse=",")) %>% ungroup()


	op_dfs = list(edges = edges,
				  nodes = nodes)

	directed = TRUE

	op_g   <- igraph::graph_from_data_frame(d = op_dfs$edges,
											directed = directed,
											vertices = op_dfs$nodes)

	igraph::E(op_g)$sources    <- strsplit(igraph::E(op_g)$sources,    ';')
	igraph::E(op_g)$references <- strsplit(igraph::E(op_g)$references, ';')

	return(op_g)

}

#' Build Omnipath interaction graph
#'
#' transforms the interactions data.frame to igraph object
#'
#' @return igraph object
#' @export
#' @param interactions data.frame created by \code{\link{import_Omnipath_Interactions}}
#' @examples
#' interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite", "Signor"))
#' OPI_g = interaction_graph(interactions = interactions )
interaction_graph <- function(interactions){
	# This is a gene_name based conversion to igraph, i.e. the vertices are identified
	# by genenames, and not by uniprot IDs.
	# This might cause issue when a gene name encodes multiple uniprot IDs.

	require(plyr)
	require(dplyr)

	# keep only edge attributes
	edges = interactions %>% select(- c(source, target))


	# build vertices: gene_names and gene_uniprotIDs
	nodesA = select(interactions, c(source_genesymbol, source))
	nodesB = select(interactions, c(target_genesymbol, target))
	colnames(nodesA) = colnames(nodesB) = c("genesymbol", "up_id")
	nodes = rbind(nodesA,nodesB)
	nodes = unique(nodes)
	nodes = nodes %>% group_by(genesymbol) %>% summarise("up_ids" = paste0(up_id,collapse=",")) %>% ungroup()
	op_dfs = list(edges = edges,
				  nodes = nodes)

	directed = TRUE

	op_g <- igraph::graph_from_data_frame(d = op_dfs$edges,
										  directed = directed,
										  vertices = op_dfs$nodes)

	igraph::E(op_g)$sources    <- strsplit(igraph::E(op_g)$sources,    ';')
	igraph::E(op_g)$references <- strsplit(igraph::E(op_g)$references, ';')

	return(op_g)

}

