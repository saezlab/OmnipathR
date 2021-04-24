#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2021
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Alberto Valdeolivas
#                  Dénes Türei (turei.denes@gmail.com)
#                  Attila Gábor
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Website: https://saezlab.github.io/omnipathr
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Gene annotations from Gene Ontology
#'
#' Gene Ontology is an ontology of gene subcellular localizations, molecular
#' functions and involvement in biological processes. Gene products across
#' many organisms are annotated with the ontology terms. This function
#' downloads the gene-ontology term associations for certain model organisms
#' or all organisms. For a description of the columns see
#' \url{http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/}.
#'
#' @param organism Character: either "chicken", "cow", "dog", "human", "pig"
#'     or "uniprot_all".
#' @param aspects Character vector with some of the following elements:
#'     "C" (cellular component), "F" (molecular function) and "P" (biological
#'     process). Gene Ontology is three separate ontologies called as three
#'     aspects. By this parameter you can control which aspects to include
#'     in the output.
#' @param slim Character: if not \code{NULL}, the name of a GOsubset (slim).
#'     instead of the full GO annotation, the slim annotation will be
#'     returned. See details at \code{\link{go_annot_slim}}.
#'
#' @return A tibble (data frame) of annotations as it is provided by the
#'     database
#'
#' @examples
#' goa_data <- go_annot_download()
#' goa_data
#' # # A tibble: 606,840 x 17
#' #    db       db_object_id db_object_symbol qualifier go_id   db_ref
#' #    <fct>    <chr>        <chr>            <fct>     <chr>   <chr>
#' #  1 UniProt. A0A024RBG1   NUDT4B           NA        GO:000. GO_REF:00.
#' #  2 UniProt. A0A024RBG1   NUDT4B           NA        GO:000. GO_REF:00.
#' #  3 UniProt. A0A024RBG1   NUDT4B           NA        GO:004. GO_REF:00.
#' #  4 UniProt. A0A024RBG1   NUDT4B           NA        GO:005. GO_REF:00.
#' #  5 UniProt. A0A024RBG1   NUDT4B           NA        GO:005. GO_REF:00.
#' # # . with 606,830 more rows, and 11 more variables:
#' # #   evidence_code <fct>, with_or_from <chr>, aspect <fct>,
#' # #   db_object_name <chr>, db_object_synonym <chr>,
#' # #   db_object_type <fct>, taxon <fct>, date <date>,
#' # #   assigned_by <fct>, annotation_extension <chr>,
#' # #   gene_product_from_id <chr>
#'
#' @export
#' @importFrom magrittr %>% %T>%
#' @importFrom readr cols col_factor col_character col_date
#' @importFrom dplyr filter
go_annot_download <- function(
    organism = 'human',
    aspects = c('C', 'F', 'P'),
    slim = NULL
){

    # NSE vs. R CMD check workaround
    aspect <- NULL

    if(!is.null(slim)){

        exec(go_annot_slim, !!!as.list(environment()))

    }

    'go_annot' %>%
    generic_downloader(
        url_param = list(organism),
        reader_param = list(
            skip = 40,
            col_names = c(
                'db',
                'db_object_id',
                'db_object_symbol',
                'qualifier',
                'go_id',
                'db_ref',
                'evidence_code',
                'with_or_from',
                'aspect',
                'db_object_name',
                'db_object_synonym',
                'db_object_type',
                'taxon',
                'date',
                'assigned_by',
                'annotation_extension',
                'gene_product_from_id'
            ),
            col_types = cols(
                db = col_factor(),
                qualifier = col_factor(),
                evidence_code = col_factor(),
                aspect = col_factor(),
                db_object_type = col_factor(),
                taxon = col_factor(),
                date = col_date('%Y%m%d'),
                assigned_by = col_factor(),
                annotation_extension = col_character(),
                gene_product_from_id = col_character()
            )
        ),
        resource = 'Gene Ontology'
    ) %>%
    filter(aspect %in% aspects) %T>%
    load_success()

}


#' The Gene Ontology tree
#'
#' @param basic Logical: use the basic or the full version of GO. As written
#'     on the GO home page: "the basic version of the GO is filtered such
#'     that the graph is guaranteed to be acyclic and annotations can be
#'     propagated up the graph. The relations included are is a, part of,
#'     regulates, negatively regulates and positively regulates. This
#'     version excludes relationships that cross the 3 GO hierarchies.
#'     This version should be used with most GO-based annotation tools."
#' @param tables In the result return data frames or nested lists. These
#'     later can be converted to each other if necessary. However converting
#'     from table to list is faster.
#' @param subset Character: the GO subset (GO slim) name. GO slims are
#'     subsets of the full GO which "give a broad overview of the ontology
#'     content without the detail of the specific fine grained terms". This
#'     option, if not \code{NULL}, overrides the \code{basic} parameter.
#'     Available GO slims are: "agr" (Alliance for Genomics Resources),
#'     "generic", "aspergillus", "candida", "drosophila", "chembl",
#'     "metagenomic", "mouse", "plant", "pir" (Protein Information Resource),
#'     "pombe" and "yeast".
#' @param relations Character vector: the relations to include in the
#'     processed data.
#'
#' @return A list with the following elements: 1) "names" a list with
#'     terms as names and names as values; 2) "namespaces" a list with
#'     terms as names and namespaces as values; 3) "relations" a list with
#'     relations between terms: terms are keys, values are lists with
#'     relations as names and character vectors of related terms as
#'     values; 4) "subsets" a list with terms as keys and character
#'     vectors of subset names as values (or \code{NULL} if the term
#'     does not belong to any subset); 5) "obsolete" character vector
#'     with all the terms labeled as obsolete. If the \code{tables}
#'     parameter is \code{TRUE}, "names", "namespaces", "relations"
#'     and "subsets" will be data frames (tibbles).
#'
#' @examples
#' # retrieve the generic GO slim, a small subset of the full ontology
#' go <- go_ontology_download(subset = 'generic')
#'
#' @importFrom magrittr %>%
#' @export
go_ontology_download <- function(
    basic = TRUE,
    tables = TRUE,
    subset = NULL,
    relations = c(
        'is_a', 'part_of', 'occurs_in', 'regulates',
        'positively_regulates', 'negatively_regulates'
    )
){

    url_key_param <- `if`(is.null(subset), 'full', 'slim')
    url_param <- `if`(is.null(subset), `if`(basic, 'go-basic', 'go'), subset)

    path <-
        'go_%s' %>%
        download_to_cache(
            url_param = list(url_param),
            url_key_param = list(url_key_param)
        )

    path %>%
    obo_parser(relations = relations, tables = tables) %>%
    copy_source_attrs(path, resource = 'Gene Ontology')

}


#' GO slim gene annotations
#'
#' GO slims are subsets of the full GO which "give a broad overview of the
#' ontology content without the detail of the specific fine grained terms".
#' In order to annotate genes with GO slim terms, we take the annotations
#' and search all ancestors of the terms up to the root of the ontology
#' tree. From the ancestors we select the terms which are part of the slim
#' subset.
#'
#' @param organism Character: either "chicken", "cow", "dog", "human", "pig"
#'     or "uniprot_all".
#' @param slim Character: the GO subset (GO slim) name. Available GO
#'     slims are: "agr" (Alliance for Genomics Resources), "generic",
#'     "aspergillus", "candida", "drosophila", "chembl", "metagenomic",
#'     "mouse", "plant", "pir" (Protein Information Resource), "pombe"
#'     and "yeast".
#' @param aspects Character vector with some of the following elements:
#'     "C" (cellular component), "F" (molecular function) and "P" (biological
#'     process). Gene Ontology is three separate ontologies called as three
#'     aspects. By this parameter you can control which aspects to include
#'     in the output.
#' @param cache Logical: Load the result from cache if available.
#'
#' @return A tibble (data frame) of genes annotated with ontology terms in
#'     in the GO slim (subset).
#'
#' @details
#' Building the GO slim is resource intensive in its current implementation.
#' For human annotation and generic GO slim it might take around 20 minutes.
#' The result is saved into the cache so next time loading the data from
#' there is really quick. If the \code{cache} option is \code{FALSE} the
#' data will be built fresh (the annotation and ontology files still might
#' come from cache), and the newly build GO slim will overwrite the cache
#' instance.
#'
#' @examples
#' \donttest{
#' goslim <- go_annot_slim(organism = 'human', slim = 'generic')
#' goslim
#' # # A tibble: 276,371 x 8
#' #    db     db_object_id db_object_symbol go_id aspect db_object_name
#' #    <fct>  <chr>        <chr>            <chr> <fct>  <chr>
#' #  1 UniPr. A0A024RBG1   NUDT4B           GO:0. F      Diphosphoinosito.
#' #  2 UniPr. A0A024RBG1   NUDT4B           GO:0. F      Diphosphoinosito.
#' #  3 UniPr. A0A024RBG1   NUDT4B           GO:0. C      Diphosphoinosito.
#' #  4 UniPr. A0A024RBG1   NUDT4B           GO:0. C      Diphosphoinosito.
#' #  5 UniPr. A0A024RBG1   NUDT4B           GO:0. C      Diphosphoinosito.
#' # # . with 276,366 more rows, and 2 more variables:
#' # #   db_object_synonym <chr>, db_object_type <fct>
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang exec !!!
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{go_annot_download}}}
#'     \item{\code{\link{go_ontology_download}}}
#'     \item{\code{\link{get_db}}}
#' }
go_annot_slim <- function(
    organism = 'human',
    slim = 'generic',
    aspects = c('C', 'F', 'P'),
    cache = TRUE
){

    cache_pseudo_url <- 'go_annot_slim'
    cache_pseudo_post <- list(
        organism = organism,
        slim = slim,
        aspects = aspects
    )

    in_cache <-
        omnipath_cache_get(
            url = cache_pseudo_url,
            post = cache_pseudo_post,
            create = FALSE
        ) %>%
        omnipath_cache_latest_version

    if(is.null(in_cache) || !cache){

        annot <-
            exec(.go_annot_slim, !!!cache_pseudo_post) %>%
            omnipath_cache_save(
                url = cache_pseudo_url,
                post = cache_pseudo_post
            )

    }else{

        annot <- omnipath_cache_load(
            url = cache_pseudo_url,
            post = cache_pseudo_post
        )

    }

    return(annot)

}


#' See \code{\link{go_annot_slim}}
#'
#' @importFrom magrittr %>%
#' @importFrom progress progress_bar
#' @importFrom dplyr left_join select distinct mutate n_distinct
#' @importFrom purrr map
#' @importFrom tidyr unnest
#'
#' @noRd
.go_annot_slim <- function(
    organism = 'human',
    slim = 'generic',
    aspects = c('C', 'F', 'P')
){

    # NSE vs. R CMD check workaround
    go_id <- qualifier <- db_ref <- evidence_code <- with_or_from <-
    assigned_by <- annotation_extension <- taxon <- gene_product_from_id <-
    NULL

    annot <- go_annot_download(organism = organism, aspects = aspects)
    slim_terms <- go_ontology_download(subset = slim)$names$term

    pb <- progress_bar$new(
        total = annot %>% pull(go_id) %>% n_distinct,
        format = '  Looking up ancestors [:bar] :percent eta: :eta'
    )

    annot %>%
    select(
        -qualifier, -db_ref, -evidence_code, -with_or_from,
        -date, -assigned_by, -annotation_extension, -taxon,
        -gene_product_from_id
    ) %>%
    left_join(
        select(., go_id) %>%
        distinct %>%
        mutate(
            ancestors = map(go_id, function(g){
                pb$tick()
                union(ancestors(g), g)
            })
        ),
        by = 'go_id'
    ) %>%
    mutate(
        go_id = map(
            ancestors,
            intersect,
            slim_terms
        )
    ) %>%
    select(-ancestors) %>%
    unnest(go_id) %>%
    distinct

}