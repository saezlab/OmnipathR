#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Diego Mananes
#                  Alberto Valdeolivas
#                  Dénes Türei (turei.denes@gmail.com)
#                  Attila Gábor
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#

#' Genome scale metabolic model by Wang et al. 2021
#'
#' Process the GEMs from Wang et al., 2021
#' (\url{https://github.com/SysBioChalmers}) into convenient tables.
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicus), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#' @param orphans Logical: include orphan reactions (reactions without known
#'     enzyme).
#'
#' @return List containing the following elements: \itemize{
#'     \item{reactions: tibble of reaction data;}
#'     \item{metabolites: tibble of metabolite data;}
#'     \item{reaction_ids: translation table of reaction identifiers;}
#'     \item{metabolite_ids: translation table of metabolite identifiers;}
#'     \item{S: Stoichiometric matrix (sparse).}
#' }
#'
#' @examples
#' gem <- chalmers_gem()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'     PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'     Genome-scale metabolic network reconstruction of model animals as a
#'     platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'     27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>% extract extract2 equals
#' @importFrom dplyr pull across na_if first mutate slice rename n
#' @importFrom vctrs vec_cast
#' @importFrom tibble as_tibble
#' @importFrom purrr map map_chr map_lgl map2 pmap_int
#' @importFrom logger log_info
#' @importFrom tidyselect everything
#' @importFrom rlang are_na
#' @importFrom stringr str_split str_replace_all str_detect
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_network}}}
#'     \item{\code{\link{chalmers_gem_metabolites}}}
#'     \item{\code{\link{chalmers_gem_reactions}}}
#'     \item{\code{\link{chalmers_gem_raw}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#'     \item{\code{\link{cosmos_pkn}}}
#' }
chalmers_gem <- function(organism = 'Human', orphans = TRUE) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    mets <- metNames <- metFormulas <- inchis <- rxns <- lb <- ub <-
    grRules <- reversible <- direction <- reactants <- products <-
    orphan <- NULL

    log_info('Processing Chalmers SysBio GEM.')

    raw <- chalmers_gem_raw(organism = organism)
    metabolite_ids <- chalmers_gem_metabolites(organism = organism)
    reaction_ids <- chalmers_gem_reactions(organism = organism)

    S <- raw %>% extract2(1L) %>% extract(,,1L) %>% extract2('S')

    metabolites <-
        raw %>%
        chalmers_gem_matlab_tibble(mets, metNames, metFormulas, inchis) %>%
        mutate(
            across(
                everything(),
                ~na_if(map_chr(unlist(.x, recursive = FALSE), first), '')
            )
        ) %T>%
        {log_info('Chalmers GEM: %i metabolites.', nrow(.))}

    reactions <-
        raw %>%
        chalmers_gem_matlab_tibble(rxns, lb, ub, rev, grRules)  %>%
        rename(reversible = rev) %>%
        mutate(
            rxns = unlist(rxns),
            direction = ifelse(ub + lb >= 0L, 1L, -1L),
            reversible = vec_cast(reversible, logical()),
            grRules = (
                unlist(grRules, recursive = FALSE) %>%
                map_chr(extract, 1L) %>%
                str_replace_all('[\\(\\)]|_AT\\d+', '') %>%
                str_split(' or ') %>%
                map(str_split, ' and ') %>%
                map(discard, is_empty_2)
            ),
            orphan = are_na(map_lgl(
                grRules,
                ~any(map_lgl(.x, ~any(str_detect(.x, '^\\d+$'))))
            )),
            reactants = map2(
                1L:n(),
                direction,
                ~.x %>% extract(S,,.) %>% equals(-.y) %>% which %>%
                     slice(metabolites, .) %>% pull(mets)
            ),
            products = map2(
                1L:n(),
                direction,
                ~.x %>% extract(S,,.) %>% equals(.y) %>% which %>%
                    slice(metabolites, .) %>% pull(mets)
            ),
            numof_interactions = pmap_int(
                list(reactants, products),
                ~length(c(..1, ..2))
            )
        ) %T>%
        {log_info(
            'Chalmers GEM: %i reactions, %i orphan.',
            nrow(.),
            pull(., orphan) %>% sum(na.rm = TRUE)
        )} %>%
        {`if`(orphans, ., filter(., !orphan))}

    return(
        list(
            organism = organism,
            reactions = reactions,
            metabolites = metabolites,
            reaction_ids = reaction_ids,
            metabolite_ids = metabolite_ids,
            S = S
        )
    )

}


#' Chalmers SysBio GEM in the form of gene-metabolite interactions
#'
#' Processing GEMs from Wang et al., 2021
#' (\url{https://github.com/SysBioChalmers}) to generate PKN for COSMOS
#'
#' @param organism_or_gem Character or integer or list or data frame: either
#'     an organism (taxon) identifier or a list containing the ``reactions``
#'     data frame as it is provided by \code{\link{chalmers_gem}}, or the
#'     reactions data frame itself.
#'     Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'     10116 (Rattus norvegicus), 7955 (Danio rerio), 7227 (Drosophila
#'     melanogaster) and 6239 (Caenorhabditis elegans).
#' @param metab_max_degree Degree cutoff used to prune metabolites with high
#'     degree assuming they are cofactors (400 by default).
#' @param protein_ids Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the "a"
#'     and "b" sides of the interaction, respectively. The default ID type for
#'     proteins is Esembl Gene ID, and by default UniProt IDs and Gene
#'     Symbols are included.
#' @param metabolite_ids Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the "a"
#'     and "b" sides of the interaction, respectively. The default ID type for
#'     metabolites is Metabolic Atlas ID, and HMDB IDs and KEGG IDs are included.
#'
#' @return Data frame (tibble) of gene-metabolite interactions.
#'
#' @examples
#' gem <- chalmers_gem_network()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'     PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'     Genome-scale metabolic network reconstruction of model animals as a
#'     platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'     27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %<>% %T>% extract2
#' @importFrom rlang is_list enquos
#' @importFrom purrr map map_chr
#' @importFrom dplyr mutate across select bind_rows filter
#' @importFrom dplyr left_join row_number arrange
#' @importFrom tidyr unnest_longer
#' @importFrom logger log_info
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem}}}
#'     \item{\code{\link{chalmers_gem_metabolites}}}
#'     \item{\code{\link{chalmers_gem_reactions}}}
#'     \item{\code{\link{chalmers_gem_raw}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#'     \item{\code{\link{cosmos_pkn}}}
#' }
chalmers_gem_network <- function(
        organism_or_gem = 'Human',
        metab_max_degree = 400L,
        protein_ids = c('uniprot', 'genesymbol'),
        metabolite_ids = c('hmdb', 'kegg')
    ) {

    .slow_doctest()

    # NSE vs R CMD check workaround
    orphan <- ri <- grRules <- reactants <- products <- reversible <-
    met_to_gene <- target <- reverse <- ensg <- metabolicatlas <- NULL

    if (metab_max_degree < 1L) {
        '`metab_max_degree` cannot be less than 1.' %T>%
        log_error_with_info %>%
        stop
    }

    organism_or_gem %<>% `if`(is_list(.), ., chalmers_gem(.))
    organism <- organism_or_gem$organism

    log_info(
        paste0(
           'Building Chalmers Sysbio GEM: organism=%s; metab_max_degree=%i; ',
           'protein_ids=%s; metabolite_ids=%s'
        ),
        organism,
        metab_max_degree,
        paste(protein_ids, collapse = ','),
        paste(metabolite_ids, collapse = ',')
    )

    organism_or_gem %>%
    extract2('reactions') %T>%
    {log_info('Chalmers Sysbio GEM: compiling gene-metabolite interactions.')} %>%
    filter(!orphan) %>%
    # once I see better what's the purpose of this labeling,
    # I might move it a bit further down the pipeline - denes
    mutate(ri = row_number()) %>%
    select(ri, grRules, reactants, products, reversible) %>%
    unnest_longer(grRules, simplify = FALSE) %>%
    # note: this represents the AND relationship between genes,
    # i.e. both genes together are required for the reaction
    mutate(ci = row_number()) %>%
    unnest_longer(grRules, indices_include = FALSE) %>%
    {bind_rows(
        binary_from_reaction(., grRules, reactants) %T>%
        {log_trace('Chalmers GEM: %i reactant-enzyme interactions.', nrow(.))},
        binary_from_reaction(., grRules, products, met_to_gene = FALSE) %T>%
        {log_trace('Chalmers GEM: %i enzyme-product interactions.', nrow(.))}
    )} %T>%
    {log_info(
        'Chalmers GEM: %i records in gene-metabolite interactions table.',
        nrow(.)
    )} %>%
    list(
        lo_degree =
            c(
                pull(., source),
                pull(., target)
            ) %>%
            table %>%
            keep(~.x < metab_max_degree) %>%
            names
    ) %>%
    {filter(
        extract2(., 1L),
        met_to_gene & source %in% .$lo_degree |
        !met_to_gene & target %in% .$lo_degree
    )} %>%
    arrange(ri) %>%
    left_join(
        select(., ri, source = target, target = source, reverse) %>%
            mutate(transporter = TRUE),
        by = c('ri', 'source', 'target', 'reverse')
    ) %T>%
    {log_info(
        paste0(
            'Chalmers GEM: %i gene-metabolite interactions after removing ',
            'metabolites with more than %i interactions.'
        ),
        nrow(.),
        metab_max_degree
    )} %>%
    translate_ids_multi(
        source = ensg,
        target,
        !!!syms(protein_ids),
        ensembl = TRUE,
        organism = organism
    ) %>%
    translate_ids_multi(
        source = metabolicatlas,
        target,
        !!!syms(metabolite_ids),
        chalmers = TRUE,
        organism = organism
    ) %T>%
    {log_info(
        'Chalmers GEM: built gene-metabolite network of %i interactions.',
        nrow(.)
    )}

}


#' Binary interactions from reactions
#'
#' These interactions point from reactants to enzymes or from enzymes to
#' products, according to the parameters.
#'
#' @importFrom rlang enquo !! sym set_names
#' @importFrom magrittr %>% inset2
#' @importFrom dplyr filter select bind_rows mutate rename
#' @importFrom tidyr unnest_longer separate_wider_regex
#' @importFrom stringr str_replace_all str_c
#' @noRd
binary_from_reaction <- function(
        reactions,
        gene_col,
        met_col,
        met_to_gene = TRUE
    ) {


    # NSE vs R CMD check workaround
    source_comp <- target_comp <- target <- ri <- ci <-
    source_id <- target_id <- reversible <- NULL

    list(
        s = .nse_ensure_str(!!enquo(gene_col)),
        t = .nse_ensure_str(!!enquo(met_col))
    ) %>%
    {`if`(met_to_gene, rev(.), .)} %>%
    c(list(unnest_longer(reactions, !!sym(.$t), indices_include = FALSE))) %>%
    set_names(c('s', 't', 'r')) %>%
    {bind_rows(
        select(.$r, ri, ci, source = !!sym(.$s), target = !!sym(.$t)) %>%
        mutate(reverse = FALSE, met_to_gene = met_to_gene),
        inset2(., 'r', filter(.$r, reversible)) %>%
        {select(.$r, ri, ci, source = !!sym(.$t), target = !!sym(.$s))} %>%
        mutate(reverse = TRUE, met_to_gene = !met_to_gene)
    )} %>%
    mutate(across(c(source, target), ~str_replace_all(.x, '[\\(\\)]', ''))) %>%
    filter(!is.na(source), !is.na(target)) %>%
    separate_wider_regex(
        c(source, target),
        patterns = c(
            id = '(?:MAM\\d+|[-\\w\\.:\'\\[\\]]+)',
            comp = '[cmxrelgni]??'
        ),
        names_sep = '_'
    ) %>%
    mutate(comp = as.factor(str_c(source_comp, target_comp))) %>%
    select(-source_comp, -target_comp) %>%
    rename(source = source_id, target = target_id)

}


#' Reactions from the Chalmers SysBio GEM (Wang et al., 2021)
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Data frame of reaction identifiers.
#'
#' @examples
#' chalmers_gem_reactions()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'   PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'   Genome-scale metabolic network reconstruction of model animals as a
#'   platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'   27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_network}}}
#'     \item{\code{\link{chalmers_gem_metabolites}}}
#'     \item{\code{\link{chalmers_gem}}}
#'     \item{\code{\link{chalmers_gem_raw}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#'     \item{\code{\link{cosmos_pkn}}}
#' }
chalmers_gem_reactions <- function(organism = 'Human') {

    .slow_doctest()

    organism %<>% organism_for('chalmers-gem')

    'chalmers_gem' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = organism %>% list('reactions.tsv'),
        reader_param = list(),
        resource = 'Chalmers GEM',
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' Metabolites from the Chalmers SysBio GEM (Wang et al., 2021)
#'
#' @param organism Character or integer: an organism (taxon) identifier.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Data frame of metabolite identifiers.
#'
#' @examples
#' chalmers_gem_metabolites()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'   PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'   Genome-scale metabolic network reconstruction of model animals as a
#'   platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'   27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom readr read_tsv
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_network}}}
#'     \item{\code{\link{chalmers_gem_reactions}}}
#'     \item{\code{\link{chalmers_gem}}}
#'     \item{\code{\link{chalmers_gem_raw}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#'     \item{\code{\link{cosmos_pkn}}}
#' }
chalmers_gem_metabolites <- function(organism = 'Human') {

    .slow_doctest()

    organism %<>% organism_for('chalmers-gem')

    generic_downloader(
        'chalmers_gem',
        reader = read_tsv,
        url_key_param = list(),
        url_param = organism %>% list('metabolites.tsv'),
        reader_param = list(),
        resource = 'Chalmers GEM',
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' GEM matlab file from Chalmers Sysbio (Wang et al., 2021)
#'
#' Downloads and imports the matlab file containing the genome scale metabolic
#' models created by Chalmers SysBio.
#'
#' @param organism Character or integer: name or identifier of the organism.
#'   Supported taxons are 9606 (Homo sapiens), 10090 (Mus musculus),
#'   10116 (Rattus norvegicu), 7955 (Danio rerio), 7227 (Drosophila
#'   melanogaster) and 6239 (Caenorhabditis elegans).
#'
#' @return Matlab object containing the GEM.
#'
#' @details The Matlab object is parsed into a nested list containing a number
#'   of vectors and two sparse matrices. The top level contains a single
#'   element under the name "ihuman" for human; under this key there is an
#'   array of 31 elements. These elements are labeled by the row names of the
#'   array.
#'
#' @examples
#' chalmers_gem_raw()
#'
#' @references Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley
#'   PE, Huang S, Gobom J, Svensson T, Uhlen M, Zetterberg H, Nielsen J.
#'   Genome-scale metabolic network reconstruction of model animals as a
#'   platform for translational research. Proc Natl Acad Sci U S A. 2021 Jul
#'   27;118(30):e2102344118. doi: \doi{10.1073/pnas.2102344118}.
#'
#' @importFrom magrittr %>% %T>% %<>%
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{chalmers_gem_network}}}
#'     \item{\code{\link{chalmers_gem_reactions}}}
#'     \item{\code{\link{chalmers_gem}}}
#'     \item{\code{\link{chalmers_gem_reactions}}}
#'     \item{\code{\link{chalmers_gem_id_mapping_table}}}
#'     \item{\code{\link{cosmos_pkn}}}
#' }
chalmers_gem_raw <- function(organism = 'Human') {

    .slow_doctest()

    organism %<>%
        organism_for('chalmers-gem') %T>%
        log_info('Downloading GEM from Chalmers Sysbio for organism `%s`.', .)

    'chalmers_gem' %>%
    generic_downloader(
        reader = R.matlab::readMat,
        url_key_param = list(),
        url_param = organism %>% list(sprintf('%s-GEM.mat', .)),
        reader_param = list(),
        resource = 'Chalmers GEM',
        post = NULL,
        use_httr = FALSE
    ) %T>%
    load_success()

}


#' Tibble from Chalmers GEM Matlab object
#'
#' @param matlab Chalmers GEM in an R object loaded from the Matlab
#'     dump.
#' @param ... Variable names: should contain either only reaction or
#'     metabolite variables, otherwise num of rows won't be uniform.
#'
#' @return Tibble with the requested variables.
#'
#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom dplyr pull
#' @importFrom magrittr %>% extract extract2 equals
#' @importFrom rlang enquos
#'
#' @noRd
chalmers_gem_matlab_tibble <- function(matlab, ...) {

    cols <-
        enquos(...) %>%
        map_chr(.nse_ensure_str)

    matlab %>%
    as_tibble %>%
    pull(1L) %>%
    extract(,,1L) %>%
    map(
        function(x) {
            if (x %>% dim %>% equals(1L) %>% all) x %>% extract(1L,1L) else x
        }
    ) %>%
    extract(cols) %>%
    map(extract,,1L) %>%
    as_tibble

}
