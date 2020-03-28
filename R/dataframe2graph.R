
#' Post-translational modifications (PTMs) graph
#'
#' transforms the ptms interactions data.frame to igraph object
#'
#' @return An igraph object
#' @export
#' @import igraph
#' @importFrom dplyr %>% select group_by summarise ungroup
#' @importFrom rlang .data
#' @param ptms data.frame created by \code{\link{import_Omnipath_PTMS}}
#' @examples
#' ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))
#' ptms_g = ptms_graph(ptms = ptms )
#' @seealso  \code{\link{import_Omnipath_PTMS}}
ptms_graph = function(ptms){
    # This is a gene_name based conversion to igraph, i.e. the vertices are 
    # identified by genenames, and not by uniprot IDs.
    # This might cause issue when a gene name encodes multiple uniprot IDs.

    # We check that the input dataframe contain the required columns.
    if (!all(c("enzyme", "substrate", "enzyme_genesymbol",
        "substrate_genesymbol","sources") %in% colnames(ptms))) {
    stop("The input data frame does not contain the required columns")
    }

    # We give the proper format to the edges by calling to the function below
    output_graph <- format_graph_edges(ptms,flag = "ptms_dataset")
    return(output_graph)
}

#' Build Omnipath interaction graph
#'
#' transforms the interactions data.frame to an igraph object
#'
#' @return An igraph object
#' @export
#' @import igraph
#' @importFrom dplyr %>% select group_by summarise ungroup
#' @param interactions data.frame created by 
#' \code{\link{import_Omnipath_Interactions}},
#' \code{\link{import_PathwayExtra_Interactions}}, 
#' \code{\link{import_KinaseExtra_Interactions}},
#' \code{\link{import_LigrecExtra_Interactions}}, 
#' \code{\link{import_TFregulons_Interactions}},
#' \code{\link{import_miRNAtarget_Interactions}} or 
#' \code{\link{import_AllInteractions}} 
#' @examples
#' interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3"))
#' OPI_g = interaction_graph(interactions)
#' @seealso \code{\link{import_Omnipath_Interactions}},
#' \code{\link{import_PathwayExtra_Interactions}}, 
#' \code{\link{import_KinaseExtra_Interactions}},
#' \code{\link{import_LigrecExtra_Interactions}}, 
#' \code{\link{import_TFregulons_Interactions}},
#' \code{\link{import_miRNAtarget_Interactions}} or 
#' \code{\link{import_AllInteractions}} 
interaction_graph <- function(interactions = interactions){
    # This is a gene_name based conversion to igraph, i.e. the vertices are 
    # identified by genenames, and not by uniprot IDs.
    # This might cause issue when a gene name encodes multiple uniprot IDs.

    # We check that the input dataframe contain the required columns.
    if (!all(c("source", "target", "source_genesymbol","target_genesymbol", 
        "sources") %in% colnames(interactions))) {
        stop("The input data frame does not contain the required columns")
    }
    
    # We give the proper format to the edges by calling to the function below
    output_graph <- 
        format_graph_edges(interactions,flag = "interactions_dataset")
    return(output_graph)  
}

########## ########## ########## ##########
########## ########## ########## ##########
########## ########## ########## ##########
## Non exported function to give the proper format to the edges in order to
## transform the interactions o PTMs data frame into a graph. The input 
## parameters are the data frame containing interactions or PTMs and a flag 
## indicating if we are dealing with the interactions dataset or the PTMs 
## dataset.

format_graph_edges <- function(df_interact, flag){

    if (flag=="ptms_dataset"){
        # keep only edge attributes
        edges <- df_interact %>% 
            dplyr::select(- c(.data$enzyme, .data$substrate))

        # build vertices: gene_names and gene_uniprotIDs
        nodesA <- 
            dplyr::select(df_interact, c(.data$enzyme_genesymbol, .data$enzyme))
        nodesB <- 
            dplyr::select(df_interact, c(.data$substrate_genesymbol, 
            .data$substrate))
    } else {
        if (flag=="interactions_dataset"){
            # keep only edge attributes
            edges <- df_interact %>% 
                dplyr::select(- c(.data$source, .data$target))      
            # build vertices: gene_names and gene_uniprotIDs
            nodesA <- 
                dplyr::select(df_interact,
                c(.data$source_genesymbol,.data$source))
            nodesB <- 
                dplyr::select(df_interact, 
                c(.data$target_genesymbol,.data$target))
        } else {
            stop("Incorrect Input Flag")
        }
    }  

    colnames(nodesA) = colnames(nodesB) = c("genesymbol", "up_id")
    nodes <- rbind(nodesA,nodesB)
    nodes <- unique(nodes)
    nodes <- nodes %>% dplyr::group_by(.data$genesymbol) %>% 
        dplyr::summarise("up_ids" = paste0(.data$up_id,collapse=",")) %>% 
        dplyr::ungroup()

    op_dfs <- list(edges = edges, nodes = nodes)
    directed <- TRUE
    op_g <- igraph::graph_from_data_frame(d = op_dfs$edges,directed = directed,
        vertices = op_dfs$nodes)

    igraph::E(op_g)$sources    <- strsplit(igraph::E(op_g)$sources,    ';')
    if ("references" %in% colnames(df_interact)){
        igraph::E(op_g)$references <- strsplit(igraph::E(op_g)$references, ';')
    }

    return(op_g)
}
