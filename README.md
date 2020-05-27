<img src='man/figures/logo_omnipath.png' align="right" height="140">

# OmnipathR

Utility functions to work with **Omnipath** in `R`. 

## Description

*OmnipathR* is an `R` package built to provide easy access to the data stored 
in the Omnipath webservice: 
    
  <http://omnipathdb.org/>
    
The webservice implements a very simple REST style API. This package make 
requests by the HTTP protocol to retreive the data. Hence, fast Internet 
access is required for a proper use of *OmnipathR*. 

The package also provides some utility functions to filter, analyse and 
visualize the data.
    
## Query types

We provide here a brief summary about the data available through *OmnipathR*.
*OmnipathR* provides access to 5 types of queries:  

1. **Interactions**: protein-protein interactions from different datasets.
2. **Post-translational modifications (PTMs)**: enzyme-substrate interactions. 
3. **Complexes**: comprehensive database of more than 22000 protein complexes.
4. **Annotations**: large variety of data about proteins and complexes features.
5. **Intercell**: information on the roles in inter-cellular signaling.

For a more detailed information, we recommend you to visit the following sites:

  <http://omnipathdb.org/>
  
  <http://omnipathdb.org/info>
  
  <https://github.com/saezlab/pypath/blob/master/webservice.rst>
  
  <https://bioconductor.org/packages/devel/bioc/vignettes/OmnipathR/inst/doc/OmnipathR.pdf>
  

## Installation
First of all, you need a current version of `R` (<http://www.r-project.org>).
*OmnipathR* is a freely available package deposited on *Bioconductor* and 
*Github*: 
(<http://bioconductor.org/>, <https://github.com/saezlab/OmnipathR>).

You can install it by running the following commands on a `R` console:
 
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## Last release in Bioconductor
BiocManager::install("OmnipathR")
## Development version with the lastest updates
BiocManager::install(version='devel')
```

## Getting started and some usage examples
To get started, we strongly recommend to read our vignette in order to deal with 
the different types of queries and handle the data they return:

  <https://bioconductor.org/packages/devel/bioc/vignettes/OmnipathR/inst/doc/OmnipathR.pdf>
  
You can also check the manual:

  <https://bioconductor.org/packages/devel/bioc/manuals/OmnipathR/man/OmnipathR.pdf>
   
In addition, we provide here some examples for a quick start: 

```{r}
library(OmnipathR)
```

Download human protein-protein interactions for some source databases:  

```{r}
interactions <- 
  import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite", 
  "SIGNOR"))
```


Download human post-translational modifications for some source databases:  

```{r}
ptms <- import_Omnipath_PTMS(filter_databases=c("PhosphoSite", "SIGNOR"))
```

Convert both data frames into networks (`igraph` objects)
```{r}
ptms_g = ptms_graph(ptms = ptms )
OPI_g = interaction_graph(interactions = interactions )
```

Print some interactions in a nice format:  
```{r}
print_interactions(head(interactions))

            source interaction           target nsources nrefs
1     UBC (P0CG48)  ==( ? )==>   IKBKG (Q9Y6K9)        5    41
2   IKBKG (Q9Y6K9)  ==( ? )==>     UBC (P0CG48)        5    41
4   PINK1 (Q9BXM7)  ==( + )==>     UBC (P0CG48)        5    19
3     UBC (P0CG48)  ==( + )==>    PRKN (O60260)        3    16
5     UBC (P0CG48)  ==( ? )==> TNFAIP3 (P21580)        5    14
6 TNFAIP3 (P21580)  ==( ? )==>     UBC (P0CG48)        5    14
```

Find interactions between a specific kinase and a specific substrate:  
```{r}
print_interactions(dplyr::filter(ptms,enzyme_genesymbol=="MAP2K1",
  substrate_genesymbol=="MAPK3"))

           enzyme interaction           substrate    modification nsources
4 MAP2K1 (Q02750)       ====> MAPK3_Y204 (P27361) phosphorylation        6
3 MAP2K1 (Q02750)       ====> MAPK3_T202 (P27361) phosphorylation        6
1 MAP2K1 (Q02750)       ====> MAPK3_Y210 (P27361) phosphorylation        1
2 MAP2K1 (Q02750)       ====> MAPK3_T207 (P27361) phosphorylation        1
           
```

Find shortest paths on the directed network between proteins:  
```{r}
printPath_es(shortest_paths(OPI_g,from = "TYRO3",to = "STAT3", 
    output = 'epath')$epath[[1]],OPI_g)

          source interaction         target nsources nrefs
1 TYRO3 (Q06418)  ==( + )==>  GRB2 (P62993)        1     1
2  GRB2 (P62993)  ==( + )==>  EGFR (P00533)       11    63
3  EGFR (P00533)  ==( + )==> STAT3 (P40763)       10    21
```

Find all shortest paths between proteins:  
```{r}
printPath_vs(all_shortest_paths(OPI_g,from = "DYRK2",to = "MAPKAPK2")$res,OPI_g)
[1] "pathway 1: DYRK2 -> TP53 -> MAPK3 -> MAPKAPK2"
          source interaction            target nsources nrefs
1 DYRK2 (Q92630)  ==( + )==>     TP53 (P04637)        7   100
2  TP53 (P04637)  ==( - )==>    MAPK3 (P27361)        5     3
3 MAPK3 (P27361)  ==( + )==> MAPKAPK2 (P49137)        4     4
[1] "pathway 2: DYRK2 -> TP53 -> MAPK14 -> MAPKAPK2"
           source interaction            target nsources nrefs
1  DYRK2 (Q92630)  ==( + )==>     TP53 (P04637)        7   100
2   TP53 (P04637)  ==( ? )==>   MAPK14 (Q16539)        5     7
3 MAPK14 (Q16539)  ==( + )==> MAPKAPK2 (P49137)       21    27
[1] "pathway 3: DYRK2 -> TP53 -> MAPK1 -> MAPKAPK2"
          source interaction            target nsources nrefs
1 DYRK2 (Q92630)  ==( + )==>     TP53 (P04637)        7   100
2  TP53 (P04637)  ==( ? )==>    MAPK1 (P28482)        7     8
3 MAPK1 (P28482)  ==( + )==> MAPKAPK2 (P49137)       10     9
```

## Feedbacks, bug reports, features
Feedbacks and bugreports are always very welcomed!  

Please use the Github issue page to report bugs or for questions: 

  <https://github.com/saezlab/OmnipathR/issues>

Many thanks for using *OmnipathR*!
