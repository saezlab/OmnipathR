#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2022
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


#' Interactions from PrePPI
#'
#' Retrieves predicted protein-protein interactions from the PrePPI
#' database (\url{http://honig.c2b2.columbia.edu/preppi}). The interactions
#' in this table are supposed to be correct with a > 0.5 probability.
#'
#' @param ... Minimum values for the scores. The available scores are:
#'     str, protpep, str_max, red, ort, phy, coexp, go, total, exp and final.
#'     Furthermore, an operator can be passed, either \code{.op = '&'} or
#'     \code{.op = '|'}, which is then used for combined filtering by
#'     multiple scores.
#'
#' @return A data frame (tibble) of interactions with scores, databases
#'     and literature references.
#'
#' @details
#' PrePPI is a combination of many prediction methods, each resulting a
#' score. For an explanation of the scores see
#' \url{https://honiglab.c2b2.columbia.edu/hfpd/help/Manual.html}.
#' The minimum, median and maximum values of the scores:
#'
#'     | Score   | Minimum | Median   | Maximum            |
#'     | ------- | ------- | -------- | ------------------ |
#'     | str     |       0 |     5.5  |           6,495    |
#'     | protpep |       0 |     3.53 |          38,138    |
#'     | str_max |       0 |    17.9  |          38,138    |
#'     | red     |       0 |     1.25 |              24.4  |
#'     | ort     |       0 |     0    |           5,000    |
#'     | phy     |       0 |     2.42 |               2.42 |
#'     | coexp   |       0 |     2.77 |              45.3  |
#'     | go      |       0 |     5.86 |             181    |
#'     | total   |       0 | 1,292    | 106,197,000,000    |
#'     | exp     |       1 |   958    |           4,626    |
#'     | final   |     600 | 1,778    |            4.91e14 |
#'
#' @examples
#' preppi <- preppi_download()
#' preppi
#' # # A tibble: 1,545,710 x 15
#' #    prot1 prot2 str_score protpep_score str_max_score red_score ort_score
#' #    <chr> <chr>     <dbl>         <dbl>         <dbl>     <dbl>     <dbl>
#' #  1 Q131. P146.     18.6           6.45         18.6      4.25      0.615
#' #  2 P064. Q96N.      1.83         14.3          14.3      4.25      0
#' #  3 Q7Z6. Q8NC.      4.57          0             4.57     0         0
#' #  4 P370. P154.    485.            0           485.       1.77      0.615
#' #  5 O004. Q9NR.     34.0           0            34.0      0.512     0
#' # # . with 1,545,700 more rows, and 8 more variables: phy_score <dbl>,
#' # #   coexp_score <dbl>, go_score <dbl>, total_score <dbl>, dbs <chr>,
#' # #   pubs <chr>, exp_score <dbl>, final_score <dbl>
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom readr read_tsv cols
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace str_replace_all
#' @export
#' @md
#' @seealso \code{\link{preppi_filter}}
preppi_download <- function(...){

    # NSE vs. R CMD check workaround
    dbs <- pubs <- NULL

    'preppi' %>%
    archive_extractor(
        path = 'preppi_final600.txt',
        reader = read_tsv,
        reader_param = list(
            col_types = cols(),
            progress = FALSE
        ),
        resource = 'PrePPI'
    ) %T>%
    load_success() %>%
    preppi_filter(...) %>%
    mutate(
        dbs = str_replace(dbs, ',$', ''),
        pubs =
            str_replace_all(pubs, 'pubmed:', '') %>%
            str_replace_all('\\|', ',')
    )


}


#' Filter PrePPI interactions by scores
#'
#' @param data A data frame of PrePPI interactions as provided by
#'     \code{\link{preppi_download}}.
#' @param ... Minimum values for the scores. The available scores are:
#'     str, protpep, str_max, red, ort, phy, coexp, go, total, exp and final.
#'     See more about the scores at \code{\link{preppi_download}}.
#' @param .op The operator to combine the scores with: either \code{'&'} or
#'     \code{'|'}. With the former, only records where all scores are above
#'     the threshold will be kept; with the latter, records where at least
#'     one score is above its threshold will be kept.
#'
#' @return The input data frame (tibble) filtered by the score thresholds.
#'
#' @examples
#' preppi <- preppi_download()
#' preppi_filtered <- preppi_filter(preppi, red = 10, str = 4.5, ort = 1)
#' nrow(preppi_filtered)
#' # [1] 8443
#'
#' @importFrom magrittr %>% %T>% extract
#' @importFrom purrr keep
#' @importFrom rlang !! parse_expr
#' @importFrom dplyr filter
#' @export
#' @seealso \code{\link{preppi_download}}
preppi_filter <- function(data, ..., .op = '&'){

    list(...) %>%
    extract(
        names(.) %>%
        keep(
            function(x){
                sprintf('%s_score', x) %in% colnames(data)
            }
        )
    ) %>%
    {`if`(
        length(.),
        sprintf('%s_score >= %g', names(.), .) %>%
        paste(collapse = .op) %>%
        {filter(data, !!parse_expr(.))},
        data
    )}

}