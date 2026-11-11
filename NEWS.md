# OmnipathR v4.1.6

## Bug fixes

+ `consensuspathdb_download`: the UniProt ID mapping query used the ID type
  labels of the retired UniProt API ("ACC", "GENENAME"), and UniProt
  rejected it with HTTP 400
+ `uniprot_id_mapping_table`: if the UniProt ID mapping job does not finish
  within `omnipathr.uniprot_idmapping_timeout`, it stops with an error that
  names the job and its status. Before, it went on, and failed with "object
  'result' not found"

# OmnipathR v4.1.5

## Bug fixes

+ Cache lock: the cache database is locked by an operating system lock on
  the file `cache.flock`, instead of the existence of the file `cache.lock`.
  The lock is released even if the R process crashes or is killed, so a
  stale lock can not block other processes any more. Before, a killed
  process left the lock behind, and every later process failed with "The
  cache directory is locked for longer than 5 seconds". Locking is also
  atomic now, and the lock is released if an error occurs while it is held
+ The cache database is written to a temporary file and then renamed, so
  it is never read half written
+ The default of the `omnipathr.cache_timeout` option is 60 seconds instead
  of 5: the wait applies only while another live process holds the lock
+ `omnipath_unlock_cache_db` releases the lock of the current session, and
  removes the legacy `cache.lock` file of older versions
+ On build servers, the cache is not wiped any more when the package is
  loaded. Each R process of the build or check of OmnipathR wiped the cache
  and removed its lock: all data was downloaded again several times, which
  exceeded the time limit of the checks, and packages built at the same
  time by the same user failed with a locked or wiped cache

## New features

+ Errors raised by OmnipathR print debugging information to the console:
  the OmnipathR version and the library it is loaded from, R and the
  platform, the versions of the main dependencies, libcurl, and the cache
  directory. The full information is printed at the first error of the
  session, later errors print only the version line. The errors have the
  class `omnipathr_error`

## Technical

+ New dependency: filelock
+ `database_summary` is not run in the examples on build servers

# OmnipathR v4.1.4

## Bug fixes

+ Organism names: the Ensembl species page is not available any more, and
  the functions that translate organism names (e.g. `ncbi_taxid`,
  `latin_name`) failed (#130). The Ensembl organisms now come from the
  `species.json` file at https://ftp.ebi.ac.uk/pub/ensemblorganisms/. It
  has more than 11,000 organisms, before there were about 340
+ `ensembl_organisms_raw` returns a different table: one row for each
  organism in `species.json`, with the NCBI Taxonomy IDs of the organism and
  of its species, the scientific and the common name
+ The Ensembl organism IDs (e.g. "hsapiens") are the names of the BioMart
  datasets, from a table in the package. Before, they were derived from the
  scientific names, and for some organisms (e.g. gorilla) they were wrong
+ If one of the resources of organism names (Ensembl, OMA or UniProt) is
  not available, the names from the other resources are still available.
  The incomplete table is not saved to the cache
+ The database cache keys contain the package version: after an upgrade,
  the databases are built again from their sources

## Technical

+ rvest is not a dependency any more
+ On build servers, the verbose console log and the cache wipe at package
  load happen only when R builds or checks OmnipathR itself. Before, the
  log output gave a NOTE in the checks of the packages that depend on
  OmnipathR

# OmnipathR v4.1.3

## Technical

+ Downloads from rest.uniprot.org use HTTP/1.1. Over HTTP/2, UniProt
  sometimes resets slow transfers, and the download fails with "HTTP/2
  stream 1 was not closed cleanly". The new option `omnipathr.http1_hosts`
  lists the hosts that use HTTP/1.1
+ UniProt queries (`all_uniprots` and the functions built on it): on build
  servers, if UniProt fails, they give a warning and an empty table, and
  skip UniProt for the rest of the session
+ Faster examples on build servers: `kinasephos` and the deprecated
  NicheNet network functions do not download in the checks of OmnipathR

# OmnipathR v4.1.2

## Deprecated

+ The NicheNet functions (`nichenet_*`) are deprecated, and the NicheNet
  vignette is removed. A new version will retrieve the NicheNet prior
  knowledge built by the OmniPath server. The resource clients written for
  NicheNet (e.g. EVEX, Harmonizome, ReMap) stay
+ `static_table` and `static_tables` are deprecated. The regular query
  functions try the plain HTTP server if the HTTPS server fails
+ The packages mlrMBO, parallelMap, ParamHelpers and smoof are not in
  Suggests any more. `nichenet_optimization` asks for them if they are not
  installed

## Bug fixes

+ `static_table`: works again for interaction tables. Before, it stopped
  with "argument is of length zero" in `unnest_evidences`
+ Slow functions return their full result when another package calls them
  on a build server. Before, they returned `NULL` in the checks of
  dependent packages, such as wppi

## Technical

+ The downloader accepts fallback URLs (mirrors). Each round of attempts
  tries all URLs once, then the next round starts after a short wait
+ HTTP errors: the log shows the status, the final URL and the response
  headers. After HTTP 429 or 503, the wait follows the `Retry-After`
  header, but it is never longer than the new option
  `omnipathr.retry_downloads_max_wait` (30 s)
+ `biomart_query`: tries the US East and Asia mirrors of Ensembl if the
  main site fails. On build servers, if all mirrors fail, it gives a
  warning and an empty table, and skips Ensembl for the rest of the session

# OmnipathR v4.1.1

## Bug fixes

+ KEGG: the organism list comes from the `list/genome` endpoint, because
  `list/organism` returns 400
+ `oma_organisms`: reads the new OMA header, and downloads the table again
  when the cached copy is malformed
+ `subnetwork`: selects the vertices with `igraph::V()`
+ `omnipath_cache_remove`: with `max_age` or `min_age`, it removes only the
  matching versions. Before, it removed all versions when `status` was not
  given
+ `resources`: the category filter works again

## Technical

+ `resources`: downloads with retries, and keeps the list in the cache for
  one hour at most
+ The environment variable `OMNIPATHR_BUILDSERVER` marks a build server
+ Vignettes build when a third-party server does not answer; the chunk
  shows the error
+ Tests of the live web service do not run on Bioconductor build machines
+ Runnable examples for 32 more manual pages
+ CITATION file in ASCII

# OmnipathR v4.0.0

## Features

+ ID translation ambiguity analysis
+ KEGG REST API client
+ New metabolomics resources: HMDB, RaMP, Chalmers Sysbio GEMs, STITCH,
+ Metabolite-protein interactions from MetalinksDB
+ Metabolite ID translation using RaMP and HMDB data

## Technical

+ Rewritten and improved parameter processing for OmniPath queries
+ Rewritten downloaders based on `curl` and `httr2`
+ Fine control over `curl` handlers
+ Detailed log messages about HTTP requests
+ Diagnostic log messages: session info, libraries, versions, platform info,
  curl options, HTTP timings, headers
+ Option to include `curl` debug log in the log file
+ Robust TCP keep-alive parameters that hopefully fix rest.uniprot.org dowloads

# OmnipathR v3.2.0

+ New resource: PrePPI
+ Configuration option to completely disable logging to file
+ New vignettes: Path reconstruction, and Bioconductor 2021 Workshop
+ Many minor bugfixes and improvements

# OmnipathR v2.99.19 (v3.0.0)

+ NicheNet pipeline works
+ Fixed an error which resulted value 1 in the `n_references` columns
even for records without references

# OmnipathR v2.99.17

+ Improved quality filtering of intercell networks

# OmnipathR v2.99.16

+ Improved API for translate_ids
+ Quality filtering of intercell networks

# OmnipathR v2.99.13

+ Functions to access KEGG Pathways
+ More robust access to UniProt (in case of network failures)

# OmnipathR v2.99.11

+ Improved downloader backends

# OmnipathR v2.99.8

+ OBO parser
+ Gene Ontology access, functions to work with ontology trees
+ Database manager

# OmnipathR v2.99.7

+ New interactions query parameter: loops
+ Fixed many caching bugs

# OmnipathR v2.99.6

+ All downloaders attempt 3 times
+ New resources: Human Phenotype Ontology and Gene Ontology annotations

# OmnipathR v2.99.0 (2021-03-08)

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

# OmnipathR v1.99.0 (2020-10-03)

+ New tests
+ Many little bugfixes
+ Updated package metadata
+ Preparation for Bioconductor 3.12

# OmnipathR v1.3.4 (2020-08-04)

+ License and password handling
+ Can be directed to different server by options

# OmnipathR v1.3.1 (2020-06-05)

+ Addition of the package website.
+ Modification of the pdf vignette for an html one.
+ Major refactoring of functions.
+ Addition of the intercellular network function.
+ Modification for the intercellular categories.

# OmnipathR v1.1.4 (2020-03-28)
+ Refactored the `import_annotations` method

# OmnipathR v0.99.12 (2019-10-10)

+ Modification in the separation between genes within a complex (from dash
to underscore)

# OmnipathR v0.99.0 (2019-10-10)

+ Submitted to Bioconductor

