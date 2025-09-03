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
gem_matlab_tibble <- function(matlab, ...) {

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

