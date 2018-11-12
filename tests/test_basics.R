# test basic function in import_Omnipath.R
library(OmnipathR)

.get_ptms_databases()
.get_interaction_databases()
ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))
interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite", "Signor"))

ptms_g = ptms_graph(ptms = ptms )
OPI_g = interaction_graph(interactions = interactions )
