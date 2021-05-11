## OmnipathR v2.99.19

+ NicheNet pipeline works
+ Fixed an error which resulted value 1 in the `n_references` columns
even for records without references

## OmnipathR v2.99.17

+ Improved quality filtering of intercell networks

## OmnipathR v2.99.16

+ Improved API for translate_ids
+ Quality filtering of intercell networks

## OmnipathR v2.99.13

+ Functions to access KEGG Pathways
+ More robust access to UniProt (in case of network failures)

## OmnipathR v2.99.11

+ Improved downloader backends

## OmnipathR v2.99.8

+ OBO parser
+ Gene Ontology access, functions to work with ontology trees
+ Database manager

## OmnipathR v2.99.7

+ New interactions query parameter: loops
+ Fixed many caching bugs

## OmnipathR v2.99.6

+ All downloaders attempt 3 times
+ New resources: Human Phenotype Ontology and Gene Ontology annotations

## OmnipathR v2.99.0 (2021-03-08)

+ Many new resources besides OmniPath:
BioPlex, ConsensusPathDB, EVEX, Guide to Pharmacology (IUPHAR/BPS),
Harmonizome, HTRIdb, InWeb InBioMap, Pathway Commons,
Ramilowski et al. 2015, RegNetwork, ReMap, TF census,
TRRUST and Vinayagam et al. 2011
+ Methods for converting network elements to Bio Model
Analyzer (BMA) motifs
+ NicheNet workflow
+ Some igraph related methods: finding all paths up to
certain length; extracting the giant component of a graph
+ Caching
+ Logging
+ Configuration handling

## OmnipathR v1.99.0 (2020-10-03)

+ New tests
+ Many little bugfixes
+ Updated package metadata
+ Preparation for Bioconductor 3.12

## OmnipathR v1.3.4 (2020-08-04)

+ License and password handling
+ Can be directed to different server by options

## OmnipathR v1.3.1 (2020-06-05)

+ Addition of the package website.
+ Modification of the pdf vignette for an html one.
+ Major refactoring of functions.
+ Addition of the intercellular network function.
+ Modification for the intercellular categories.

## OmnipathR v1.1.4 (2020-03-28)
+ Refactored the `import_annotations` method

## OmnipathR v0.99.12 (2019-10-10)

+ Modification in the separation between genes within a complex (from dash
to underscore)

## OmnipathR v0.99.0 (2019-10-10)

+ Submitted to Bioconductor

