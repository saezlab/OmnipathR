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


#' Downloads a BioPlex interaction dataset
#'
#' @param version Character: either 1.0, 2.0, 3.0 or HCT116_1.0
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv cols
#'
#' @seealso \code{\link{bioplex1}, \link{bioplex2}, \link{bioplex3},
#'     \link{bioplex_hct116_1}}
#'
#' @noRd
bioplex_download <- function(version){

    result <-
        'omnipath.bioplex_%s_url' %>%
        generic_downloader(
            url_key_param = list(version),
            resource = sprintf('BioPlex %s', version)
        )

    names(result) <- c(
        'GeneA',
        'GeneB',
        'UniprotA',
        'UniprotB',
        'SymbolA',
        'SymbolB',
        'p_wrong',
        'p_no_interaction',
        'p_interaction'
    )

    result %T>% load_success

}


#' Downloads the BioPlex version 1.0 interaction dataset
#'
#' This dataset contains ~24,000 interactions detected in HEK293T cells
#' using 2,594 baits.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @export
#'
#' @examples
#' \donttest{
#' bioplex_interactions <- bioplex1()
#' bioplex_interactions %>% nrow
#' # [1] 23744
#' bioplex_interactions %>% colnames
#' # [1] "GeneA"         "GeneB"        "UniprotA"   "UniprotB"
#' # [5] "SymbolA"       "SymbolB"      "p_wrong"    "p_no_interaction"
#' # [9] "p_interaction"
#' }
#'
#' @seealso \itemize{
#'     \item{\code{\link{bioplex2}}}
#'     \item{\code{\link{bioplex3}}}
#'     \item{\code{\link{bioplex_hct116_1}}}
#'     \item{\code{\link{bioplex_all}}}
#' }
bioplex1 <- function(){

    bioplex_download(version = '1.0')

}


#' Downloads the BioPlex version 2.0 interaction dataset
#'
#' This dataset contains ~56,000 interactions detected in HEK293T cells
#' using 5,891 baits.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @export
#'
#' @examples
#' \donttest{
#' bioplex_interactions <- bioplex2()
#' bioplex_interactions %>% nrow
#' # [1] 56553
#' bioplex_interactions %>% colnames
#' # [1] "GeneA"         "GeneB"        "UniprotA"   "UniprotB"
#' # [5] "SymbolA"       "SymbolB"      "p_wrong"    "p_no_interaction"
#' # [9] "p_interaction"
#' }
#'
#' @seealso \itemize{
#'     \item{\code{\link{bioplex1}}}
#'     \item{\code{\link{bioplex3}}}
#'     \item{\code{\link{bioplex_hct116_1}}}
#'     \item{\code{\link{bioplex_all}}}
#' }
bioplex2 <- function(){

    bioplex_download(version = '2.0')

}


#' Downloads the BioPlex version 3.0 interaction dataset
#'
#' This dataset contains ~120,000 interactions detected in HEK293T cells
#' using 10,128 baits.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @export
#'
#' @examples
#' \donttest{
#' bioplex_interactions <- bioplex3()
#' bioplex_interactions %>% nrow
#' # [1] 118162
#' bioplex_interactions %>% colnames
#' # [1] "GeneA"         "GeneB"        "UniprotA"   "UniprotB"
#' # [5] "SymbolA"       "SymbolB"      "p_wrong"    "p_no_interaction"
#' # [9] "p_interaction"
#' }
#'
#' @seealso \itemize{
#'     \item{\code{\link{bioplex1}}}
#'     \item{\code{\link{bioplex2}}}
#'     \item{\code{\link{bioplex_hct116_1}}}
#'     \item{\code{\link{bioplex_all}}}
#' }
bioplex3 <- function(){

    bioplex_download(version = '3.0')

}


#' Downloads the BioPlex HCT116 version 1.0 interaction dataset
#'
#' This dataset contains ~71,000 interactions detected in HCT116 cells
#' using 5,522 baits.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @export
#'
#' @examples
#' \donttest{
#' bioplex_interactions <- bioplex_hct116_1()
#' bioplex_interactions %>% nrow
#' # [1] 70966
#' bioplex_interactions %>% colnames
#' # [1] "GeneA"         "GeneB"        "UniprotA"   "UniprotB"
#' # [5] "SymbolA"       "SymbolB"      "p_wrong"    "p_no_interaction"
#' # [9] "p_interaction"
#' }
#'
#' @seealso \itemize{
#'     \item{\code{\link{bioplex1}}}
#'     \item{\code{\link{bioplex2}}}
#'     \item{\code{\link{bioplex3}}}
#'     \item{\code{\link{bioplex_all}}}
#' }
bioplex_hct116_1 <- function(){

    bioplex_download(version = 'HCT116_1.0')

}


#' Downloads all BioPlex interaction datasets
#'
#' BioPlex provides four interaction datasets: version 1.0, 2.0, 3.0 and
#' HCT116 version 1.0. This function downloads all of them, merges them to
#' one data frame, removes the duplicates (based on unique pairs of UniProt
#' IDs) and separates the isoform numbers from the UniProt IDs.
#' More details at https://bioplex.hms.harvard.edu/interactions.php
#'
#' @param unique Logical. Collapse the duplicate interactions into single
#'     rows or keep them as they are. In case of merging duplicate records
#'     the maximum p value will be choosen for each record.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate first summarize_all group_by bind_rows
#' @importFrom tidyr separate
#' @export
#'
#' @examples
#' \donttest{
#' bioplex_interactions <- bioplex_all()
#' bioplex_interactions
#' # # A tibble: 195,538 x 11
#' #    UniprotA IsoformA UniprotB IsoformB GeneA GeneB SymbolA SymbolB
#' #    <chr>       <int> <chr>       <int> <dbl> <dbl> <chr>   <chr>
#' #  1 A0AV02          2 Q5K4L6         NA 84561 11000 SLC12A8 SLC27A3
#' #  2 A0AV02          2 Q8N5V2         NA 84561 25791 SLC12A8 NGEF
#' #  3 A0AV02          2 Q9H6S3         NA 84561 64787 SLC12A8 EPS8L2
#' #  4 A0AV96          2 O00425          2 54502 10643 RBM47   IGF2BP3
#' #  5 A0AV96          2 O00443         NA 54502  5286 RBM47   PIK3C2A
#' #  6 A0AV96          2 O43426         NA 54502  8867 RBM47   SYNJ1
#' #  7 A0AV96          2 O75127         NA 54502 26024 RBM47   PTCD1
#' #  8 A0AV96          2 O95208          2 54502 22905 RBM47   EPN2
#' #  9 A0AV96          2 O95900         NA 54502 26995 RBM47   TRUB2
#' # 10 A0AV96          2 P07910          2 54502  3183 RBM47   HNRNPC
#' # # . with 195,528 more rows, and 3 more variables: p_wrong <dbl>,
#' # #   p_no_interaction <dbl>, p_interaction <dbl>
#' }
#'
#' @seealso \itemize{
#'     \item{\code{\link{bioplex1}}}
#'     \item{\code{\link{bioplex2}}}
#'     \item{\code{\link{bioplex3}}}
#'     \item{\code{\link{bioplex_hct116_1}}}
#' }
bioplex_all <- function(unique = TRUE){

    # NSE vs. R CMD check workaround
    UniprotA <- UniprotB <- p_wrong <- p_no_interaction <- p_interaction <-
        NULL

    bind_rows(
        bioplex1(),
        bioplex2(),
        bioplex3(),
        bioplex_hct116_1()
    ) %>%
    {`if`(
        unique,
        group_by(., UniprotA, UniprotB) %>%
        mutate(
            p_wrong = max(p_wrong),
            p_no_interaction = max(p_no_interaction),
            p_interaction = max(p_interaction)
        ) %>%
        summarize_all(first),
        .
    )} %>%
    separate(
        UniprotA,
        into = c('UniprotA', 'IsoformA'),
        sep = '-',
        fill = 'right',
        convert = TRUE
    ) %>%
    separate(
        UniprotB,
        into = c('UniprotB', 'IsoformB'),
        sep = '-',
        fill = 'right',
        convert = TRUE
    )

}