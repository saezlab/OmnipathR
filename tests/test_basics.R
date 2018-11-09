# test basic function in import_Omnipath.R

.get_ptms_databases()
.get_interaction_databases()
ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))
interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite", "Signor"))

ptms_g = omnipath_graph(ptms = ptms )
OPI_g = omnipath_graph(interactions = interactions )
