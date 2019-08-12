# test basic function in import_Omnipath.R
library(OmnipathR)

.get_ptms_databases()
.get_interaction_databases()
.get_complexes_databases()
.get_annotation_databases()
.get_intercell_categories()


ptms <- import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))
interactions <- import_Omnipath_Interactions(filter_databases=c("SignaLink3",
  "PhosphoSite", "Signor"))
interactions_PE <- import_PathwayExtra_Interactions(filter_databases=c("BioGRID",
  "IntAct"), select_organism = 9606)
interactions_KE <- import_KinaseExtra_Interactions(filter_databases=c("PhosphoPoint",
  "PhosphoSite"),select_organism = 9606)
interactions_LE <- import_LigrecExtra_Interactions(filter_databases=c("HPRD",
  "Guide2Pharma"),select_organism=9606)
interactions_TF <- import_TFregulons_Interactions(filter_databases=c("tfact",
  "ARACNe-GTEx"),select_organism=9606)
interactions_mi <- import_miRNAtarget_Interactions(filter_databases=c("miRTarBase",
  "miRecords"))
interactions_ALL <- import_AllInteractions(filter_databases=c("HPRD","BioGRID"),
                                       select_organism = 9606)


complexes = import_Omnipath_complexes(filter_databases=c("CORUM", "hu.MAP"))

annotations = import_Omnipath_annotations(select_genes=c("TP53","LMNA"),
                                          filter_databases=c("HPA"))

intercell = import_Omnipath_intercell(select_categories=c("ecm"))


ptms = get_signed_ptms(ptms,interactions)

ptms_g = ptms_graph(ptms = ptms )
OPI_g = interaction_graph(interactions = interactions )


query_genes = c("LMNA","BANF1")
complexes_query_genes = get_complex_genes(complexes,query_genes)
