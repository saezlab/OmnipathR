#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
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


#' Retrieve the STITCH links dataset
#'
#' @param organism Character or integer: name or NCBI Taxonomy ID of an
#'     organism. STITCH supports many organisms, please refer to their web site
#'     at \url{https://stitch.embl.de/}.
#'
#' @return Data frame: organism specific STITCH links dataset.
#'
#' @examples
#' stl <- stitch_links()
#'
#' @importFrom magrittr %<>% %>% %T>%
#' @importFrom readr read_tsv
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{stitch_actions}}}
#'     \item{\code{\link{stitch_links}}}
#'     \item{\code{\link{stitch_network}}}
#' }
stitch_links <- function(organism = 'human', prefixes = FALSE) {

    .slow_doctest()

    organism %<>% organism_for('stitch')
    log_trace('Loading STITCH protein-small molecule links.')

    'stitch_links' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = list(organism),
        reader_param = list(trim_ws = TRUE),
        resource = NULL,
        post = NULL,
        use_httr = FALSE
    ) %>%
    stitch_remove_prefixes(chemical, protein, remove = !prefixes) %T>%
    load_success()

}


#' Retrieve the STITCH actions dataset
#'
#' @param organism Character or integer: name or NCBI Taxonomy ID of an
#'     organism. STITCH supports many organisms, please refer to their web site
#'     at \url{https://stitch.embl.de/}.
#'
#' @return Data frame of STITCH actions.
#'
#' @examples
#' sta <- stitch_actions(organism = 'mouse')
#'
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom readr read_tsv
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{stitch_actions}}}
#'     \item{\code{\link{stitch_links}}}
#'     \item{\code{\link{stitch_network}}}
#' }
stitch_actions <- function(organism = 'human', prefixes = FALSE) {

    .slow_doctest()

    organism %<>% organism_for('stitch')
    log_trace('Loading STITCH actions.')

    'stitch_actions' %>%
    generic_downloader(
        reader = read_tsv,
        url_key_param = list(),
        url_param = list(organism),
        reader_param = list(trim_ws = TRUE),
        resource = NULL,
        post = NULL,
        use_httr = FALSE
    ) %>%
    mutate(action = as.factor(action)) %>%
    stitch_remove_prefixes(item_id_a, item_id_b, remove = !prefixes) %T>%
    load_success()

}


#' Remove the prefixes from STITCH identifiers
#'
#' STITCH adds the NCBI Taxonomy ID as a prefix to Ensembl protein identifiers,
#' e.g. "9606.ENSP00000170630", and "CID" followed by "s" or "m"
#' (stereospecific or merged, respectively) in front of PubChem Compound
#' Identifiers. It also pads the CID with zeros. This function removes these
#' prefixes, leaving only the identifiers.
#'
#' @param d Data frame, typically the output of \code{\link{stitch_links}} or
#'     \code{\link{stitch_actions}}.
#' @param ... Names of columns to remove prefixes from. NSE is supported.
#' @param remove Logical: remove the prefixes? If FALSE, this function does
#'     nothing.
#'
#' @examples
#' stitch_remove_prefixes(
#'     tibble(a = c('9606.ENSP00000170630', 'CIDs00012345')),
#'     a
#' )
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang enquos
#' @importFrom purrr map_chr
#' @importFrom dplyr mutate across
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{stitch_actions}}}
#'     \item{\code{\link{stitch_links}}}
#'     \item{\code{\link{stitch_network}}}
#' }
stitch_remove_prefixes <- function(d, ..., remove = TRUE) {

    if(remove) {

        cols <- enquos(...) %>% map_chr(.nse_ensure_str)

        d %<>%
            mutate(
                across(
                    cols,
                    ~str_replace(.x, '^(\\d+\\.|CID[ms]0*)', '')
                )
            )

    }

    d

}


#' Chemical-protein interactions from STITCH
#'
#' @param organism Character or integer: name or NCBI Taxonomy ID of an
#'     organism. STITCH supports many organisms, please refer to their web site
#'     at \url{https://stitch.embl.de/}.
#' @param min_score Confidence cutoff used for STITCH connections
#'     (700 by default).
#' @param protein_ids Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the "a"
#'     and "b" sides of the interaction, respectively. The default ID type for
#'     proteins is Esembl Protein ID, and by default UniProt IDs and Gene
#'     Symbols are included.
#' @param metabolite_ids Character: translate the protein identifiers to these ID
#'     types. Each ID type results two extra columns in the output, for the "a"
#'     and "b" sides of the interaction, respectively. The default ID type for
#'     metabolites is PubChem CID, and HMDB IDs and KEGG IDs are included.
#' @param cosmos Logical: use COSMOS format?
#'
#' @return A data frame of STITCH chemical-protein and protein-chemical
#' interactions with their effect signs, and optionally with identifiers
#' translated.
#'
#' @examples
#' stn <- stitch_network(protein_ids = 'genesymbol', metabolite_ids = 'hmdb')
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr bind_rows select filter mutate rename inner_join
#' @importFrom dplyr row_number
#' @importFrom tidyr unite
#' @importFrom purrr reduce
#' @importFrom rlang syms !!!
#' @importFrom logger log_info
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{stitch_actions}}}
#'     \item{\code{\link{stitch_links}}}
#'     \item{\code{\link{stitch_remove_prefixes}}}
#' }
stitch_network <- function(
        organism = 'human',
        min_score = 700L,
        protein_ids = c('uniprot', 'genesymbol'),
        metabolite_ids = c('hmdb', 'kegg'),
        cosmos = FALSE
    ) {

    .slow_doctest()

    # NSE vs. R CMD check workaround
    chemical <- protein <- item_id_a <- item_id_b <- CID <- HMDB <-
    a_is_acting <- ensp <- genesymbol <- pubchem <- hmdb <-
    hmdb_source <- hmdb_target <- item_id <- action <- NULL

    log_info(
        paste0(
           'Building STITCH GEM: organism=%s, min_score=%i, ',
           'protein_ids=%s, metabolite_ids=%s, cosmos=%s'
        ),
        organism,
        min_score,
        paste0(protein_ids, collapse = ','),
        paste0(metabolite_ids, collapse = ','),
        cosmos
    )

    organism %<>% organism_for('stitch')

    links <-
        stitch_links(organism) %>%
        filter(
            combined_score >= min_score,
            experimental >= min_score | database >= min_score
        ) %>%
        select(
            item_id_a = chemical,
            item_id_b = protein
        ) %>%
        distinct %>%
        bind_rows(
            .,
            rename(
                .,
                item_id_a = item_id_b,
                item_id_b = item_id_a
            )
        )

    stitch_actions(organism) %>%
    filter(
        mode == 'activation' |
        mode == 'inhibition',
        a_is_acting
    ) %>%
    select(-mode, -a_is_acting) %>%
    inner_join(links, by = c('item_id_a', 'item_id_b')) %>%
    mutate(record_id = row_number()) %T>%
    {log_info(
        'STITCH GEM: %i interactions before ID translation.',
        nrow(.)
    )} %>%
    translate_ids_multi(
        item_id = ensp,
        !!!syms(protein_ids),
        suffixes = c('a', 'b'),
        ensembl = TRUE,
        organism = organism
    ) %>%
    translate_ids_multi(
        item_id = pubchem,
        !!!syms(metabolite_ids),
        suffixes = c('a', 'b'),
        entity_type = 'metabolite',
        organism = organism
    ) %>%
    rename(sign = action) %>%
    {`if`(
        cosmos,
        mutate(
            .,
            sign = ifelse(sign == 'inhibition', -1L, 1L)
        ),
        .,
    )} %T>%
    {log_info('STITCH GEM ready: %i interactions.', nrow(.))}

}
