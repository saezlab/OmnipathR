# OmnipathR
utility functions to work with Omnipath in R


### Description

basic utility functions to download and interact with data from Omnipath webservice (www.omnipathdb.org).

### Install
you can use the `devtools` package to install from the GitHub in one line: 
```{r}
if(!require(devtools)) install.packages("devtools")
devtools::install_github("saezlab/omnipathR")  
```
Or download, unzip and install the usual way:
`install.packages('./OmnipathR',repo=NULL)`

### Examples

Download post-translational modifications:  
`ptms = import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "Signor"))`

Download protein-protein interactions:  
`interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite", "Signor"))`

Convert to igraph objects:  
`ptms_g = ptms_graph(ptms = ptms )`
`OPI_g = interaction_graph(interactions = interactions )`

Print some interactions:  
```{r}
print_interactions(head(ptms))

          enzyme interaction            substrate    modification nsources nrefs
1 PRKCA (P17252)  ==( + )==> NPHS1_T1120 (O60500) phosphorylation        3     0
2 PRKCA (P17252)  ==( + )==> NPHS1_T1125 (O60500) phosphorylation        3     0
3 PRKCA (P17252)  ==( + )==>  PDE3A_S465 (Q14432) phosphorylation        3     0
4 PRKCA (P17252)  ==( + )==>  PDE3A_S428 (Q14432) phosphorylation        3     0
5 PRKCA (P17252)  ==( + )==>  PDE3A_S438 (Q14432) phosphorylation        3     0
6 PRKCA (P17252)  ==( + )==>  PDE3A_S312 (Q14432) phosphorylation        3     0
```

Interactions with references:  
`print_interactions(tail(ptms),writeRefs=T)`

Find interactions between kinase and substrate:  
```{r}
print_interactions(dplyr::filter(ptms,enzyme_genesymbol=="MAP2K1",substrate_genesymbol=="MAPK3"))

           enzyme interaction           substrate    modification nsources nrefs
2 MAP2K1 (Q02750)  ==( + )==> MAPK3_Y204 (P27361) phosphorylation        6    13
1 MAP2K1 (Q02750)  ==( + )==> MAPK3_T202 (P27361) phosphorylation        6     8
3 MAP2K1 (Q02750)  ==( + )==>  MAPK3_T80 (P27361) phosphorylation        1     0
4 MAP2K1 (Q02750)  ==( + )==> MAPK3_Y222 (P27361) phosphorylation        1     0
5 MAP2K1 (Q02750)  ==( + )==> MAPK3_Y210 (P27361) phosphorylation        1     0
6 MAP2K1 (Q02750)  ==( + )==> MAPK3_T207 (P27361) phosphorylation        1     0
```

Find shortest paths on the directed network between proteins:  
```{r}
> printPath_es(shortest_paths(OPI_g,from = "TYRO3",to = "STAT3", output = 'epath')$epath[[1]],OPI_g)

      source     interaction       target nsources nrefs
1  TYRO3 (Q06418)  ==( + )==> PIK3R1 (P27986)        4     3
2 PIK3R1 (P27986)  ==( ? )==>     AR (P10275)        3     3
3     AR (P10275)  ==( ? )==>  STAT3 (P40763)        3     4`
```

Find all shortest paths between proteins:  
```{r}
printPath_vs(all_shortest_paths(OPI_g,from = "DYRK2",to = "MAPKAPK2")$res,OPI_g)
[1] "pathway 1: DYRK2 -> TP53 -> MAPK3 -> MAPKAPK2"
          source interaction            target nsources nrefs
1 DYRK2 (Q92630)  ==( + )==>     TP53 (P04637)        6   102
2  TP53 (P04637)  ==( ? )==>    MAPK3 (P27361)        4     6
3 MAPK3 (P27361)  ==( + )==> MAPKAPK2 (P49137)        3     4
[1] "pathway 2: DYRK2 -> TP53 -> MAPK14 -> MAPKAPK2"
           source interaction            target nsources nrefs
1  DYRK2 (Q92630)  ==( + )==>     TP53 (P04637)        6   102
2   TP53 (P04637)  ==( ? )==>   MAPK14 (Q16539)        3     8
3 MAPK14 (Q16539)  ==( + )==> MAPKAPK2 (P49137)       19    40
[1] "pathway 3: DYRK2 -> TP53 -> MAPK1 -> MAPKAPK2"
          source interaction            target nsources nrefs
1 DYRK2 (Q92630)  ==( + )==>     TP53 (P04637)        6   102
2  TP53 (P04637)  ==( ? )==>    MAPK1 (P28482)        5    11
3 MAPK1 (P28482)  ==( + )==> MAPKAPK2 (P49137)        9    11
```
