# CLAUDE.md: OmnipathR - Your Comprehensive Guide

## Package Overview

**OmnipathR** is a sophisticated R client for molecular biology databases and the flagship client for the OmniPath web service. This Bioconductor package provides seamless access to curated molecular interaction data, pathway information, and biological annotations from more than 30 public resources. With its robust download machinery, intelligent caching system, and extensive utility functions, OmnipathR serves as an essential tool for systems biology, network analysis, and bioinformatics research.

## Core Features

### 🔗 OmniPath Database Integration
- **Primary Focus**: Full-featured client for [OmniPath](https://omnipathdb.org/) - a comprehensive database combining molecular signaling data from 100+ resources
- **Data Types**:
  - Protein-protein interactions (PPIs) with directionality and effect signs
  - Enzyme-PTM (post-translational modification) relationships
  - Gene regulatory interactions
  - Protein complexes
  - Intercellular communication roles (ligands, receptors, mediators, etc.)
  - Functional annotations (GO, localization, structure, disease associations)

### 🌐 Multi-Resource Database Access
Direct clients for 20-30 major molecular biology databases:
- **Interaction Networks**: BioPlex, InWeb InBioMap, ConsensusPathDB, Pathway Commons
- **Literature Mining**: EVEX, Ramilowski et al. 2015, Vinayagam et al. 2011
- **Regulatory Data**: TRRUST, RegNetwork, TF census, HTRIdb
- **Pathways**: KEGG, Guide to Pharmacology (IUPHAR/BPS)
- **Annotations**: Gene Ontology, Human Phenotype Ontology, Harmonizome

### 🔄 Advanced ID Mapping & Translation
The **`translate_ids`** function exemplifies the package's sophisticated approach to biological data integration:
- **Multi-source translation**: UniProt, Ensembl, HMDB, RAMP, Chalmers GEMs
- **Complex handling**: Supports protein complexes with configurable expansion rules
- **Ambiguity management**: Advanced handling of one-to-many mappings with quantification and qualification options
- **Tidyverse integration**: Heavy use of `purrr`, `tidyr`, and `dplyr` for functional programming workflows

### ⚙️ Robust Infrastructure

#### Download Machinery (`R/internals.R`)
- **Resilient HTTP/HTTPS handling** with retry logic and error management
- **Configurable chunking** for large downloads
- **Compression support** (gzip, bzip2, xz)
- **Parallel processing** capabilities
- **Progress tracking** and detailed logging

#### Caching System (`R/cache.R`)
- **Intelligent local caching** with SQLite-based metadata management
- **Cache versioning** and automatic updates
- **Configurable cache directories** and cleanup routines
- **Status tracking** (READY, STARTED, FAILED, DELETED)
- **Metadata preservation** for download provenance

#### Configuration & Logging
- **Comprehensive options system** for customization
- **Structured logging** with multiple levels (TRACE, DEBUG, INFO, WARN, ERROR)
- **Session management** and configuration persistence

## Package Architecture

### Directory Structure
```
R/
├── internals.R          # Core download and URL handling machinery
├── cache.R              # Local caching system implementation
├── id_mapping.R         # Advanced ID translation functions
├── omnipath.R           # Main OmniPath API client
├── interactions.R       # Interaction network functions
├── annotations.R        # Biological annotation utilities
├── complexes.R          # Protein complex handling
├── intercell.R          # Intercellular communication data
├── orthology.R          # Cross-species translation
├── taxonomy.R           # Species and taxonomy utilities
├── nichenet.R           # NicheNet integration
├── cosmos.R             # Causal reasoning integration
└── [20+ other R files]  # Database-specific clients and utilities
```

### Key Design Patterns
- **Tidyverse-first approach**: Consistent use of pipe operators (`%>%`), functional programming with `purrr`, and data manipulation with `dplyr`/`tidyr`
- **Error handling**: Comprehensive validation with `checkmate` and custom error messages
- **Modular design**: Each database resource has dedicated modules with standardized interfaces
- **Caching by default**: All downloads are cached by default with intelligent cache invalidation
- **Configuration-driven**: Extensive options system for customization without code modification

## Integration Highlights

### NicheNet Integration
- **Seamless integration** with the NicheNet ligand activity prediction pipeline
- **Direct data access** for NicheNet's weighted networks and target genes
- **Pre-built workflows** for common ligand-receptor analysis tasks

### Bioconductor Standards Compliance
- **R CMD check**: Passes all Bioconductor quality standards
- **Documentation**: Comprehensive roxygen2 documentation with examples
- **Vignettes**: Extensive tutorials including Bioconductor workshop materials
- **Testing**: Robust test suite with `testthat`
- **Licensing**: MIT license with proper attribution

### Causal Reasoning Support
- **COSMOS integration**: Automatic generation of prior knowledge networks for causal reasoning
- **CARNIVAL compatibility**: Data formatted for causal network reconstruction
- **DoRothEA & PROGENy**: Integrated workflows for transcriptomics analysis

## Usage Patterns

### Basic Workflow
```r
# Initialize OmnipathR
library(OmnipathR)

# Download protein-protein interactions
ppi <- import_omnipath_interactions()

# Translate gene IDs
genes_hgnc <- translate_ids(
    ppi$source,
    from = "uniprot",
    to = "genesymbol"
)

# Access intercellular communication data
intercell <- import_intercell_network()
```

### Advanced Features
```r
# Custom download with caching
interactions <- import_post_translational_interactions(
    resources = c("PhosphoSite", "SIGNOR"),
    organism = 9606,
    reviewed = TRUE
)

# Multi-resource ID translation with ambiguity handling
id_mapping <- translate_ids(
    data_frame,
    uniprot,
    genesymbol,
    ensembl,
    uploadlists = TRUE,
    quantify_ambiguity = TRUE,
    keep_untranslated = TRUE
)
```

## Developer Notes

### Code Style
- **Tidyverse conventions**: Consistent with RStudio's tidyverse style guide
- **Functional programming**: Heavy use of `purrr::map()` and anonymous functions
- **Pipe-heavy workflows**: Chaining operations with `magrittr` pipes
- **Type safety**: Extensive input validation with `checkmate`

### Performance Considerations
- **Lazy loading**: Functions load data only when needed
- **Memory efficiency**: Streaming downloads for large datasets
- **Parallel processing**: Support for parallel execution where beneficial
- **Chunked operations**: Large datasets processed in configurable chunks

### Extensibility
- **Plugin architecture**: Easy addition of new database resources
- **Standardized interfaces**: Consistent patterns across all resource clients
- **Configuration system**: Runtime customization without code changes
- **Hook system**: Extension points for custom processing

## Getting Started

For new developers working with OmnipathR:

1. **Package Development Setup**: Since we're actively developing the package, don't use `library(OmnipathR)`. Instead use:
   ```r
   library(devtools)
   load_all()
   ```

2. **Read the Bioconductor workshop vignette** (`vignettes/bioc_workshop.Rmd`) for a comprehensive overview
3. **Examine `translate_ids`** in `R/id_mapping.R` to understand the package's sophisticated approach to data integration
4. **Review `R/internals.R`** to understand the robust download machinery (see summary below)
5. **Check the package tests** for usage patterns and expected behaviors
6. **Explore the extensive vignettes** covering specific features and use cases

## Development Troubleshooting

### Cache Issues
- **Common problem**: Corrupted cache files can cause issues during development
- **Warning about incomplete final line**: This comes from jsonlite when JSON lacks trailing newline - it doesn't necessarily mean corruption
- **Solution**: Remove corrupted files and rebuild cache inventory:
   ```bash
   rm /path/to/cache/filename
   ```
   Then in R:
   ```r
   omnipath_cache_clean_db()
   ```

### Error Handling Guidelines
- **Proper logging**: Always use this pattern for errors/warnings:
   1. Create message with `sprintf()`
   2. Log with `log_error()` or `log_warn()`
   3. Raise condition with `stop()` or `warning()`

### Package Architecture

#### Download Machinery (`R/internals.R`)
The package provides several sophisticated download functions that integrate caching, logging, and reliability:

- **`download_base()`**: Core low-level downloader with retry logic and HTTP header handling
- **`download_to_cache()`**: Downloads files with automatic caching - use for most downloads
- **`generic_downloader()`**: Downloads and reads tabular data with reader functions
- **`archive_downloader()`**: Downloads and extracts archive files (zip, tar, etc.)
- **`archive_extractor()`**: Downloads archives and extracts specific files with optional reader callbacks

**Key patterns**:
- Always use package download functions instead of direct httr2/httr calls
- Use `get_url()` or `url_parser()` to construct URLs from `inst/internal/urls.json`
- URLs are registered in `inst/internal/urls.json` and loaded via `.load_urls()`
- Archive extraction supports callbacks for immediate reading (e.g., `reader = SBMLR::readSBML`)
- All downloads integrate with the caching system automatically

#### URL Management
- **URL registry**: All resource URLs stored in `inst/internal/urls.json`
- **Parameter substitution**: URLs support `%s` placeholders via `url_param` arguments
- **Base URL construction**: Use `url_parser()` for consistent URL building
- **Resource-specific patterns**: Each database follows standardized URL keys

## Quality Assurance

- **Continuous integration**: Automated testing on multiple R versions
- **Bioconductor compliance**: Regular checks against Bioconductor standards
- **Code coverage**: Comprehensive test coverage for core functionality
- **Documentation**: Complete roxygen2 documentation with examples
- **Performance monitoring**: Regular benchmarking of download and processing speeds

This package represents years of development in computational systems biology, combining rigorous database integration with user-friendly R programming patterns to create a powerful tool for molecular biology research.

## Saezverse Reference

Use this Saezverse entry as the cross-repo reference for OmnipathR. When
working on OmnipathR here, check the Saezverse OmnipathR spec for the
high-level package description, architecture summary, dependencies, and
known issues to ensure changes stay aligned.

- Saezverse repository: https://github.com/saezlab/saezverse
- OmnipathR spec: human/packages/omnipathr.md
- Last reviewed: 2026-03-03
