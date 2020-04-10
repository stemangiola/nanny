#' Creates an histogram from a data frame
#'
#' \lifecycle{experimental}
#'
#' @description create_histogram() creates a `ggplot` object from a data frame
#'
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#'
#' @name create_histogram
#'
#' @param .data A numerical vector
#' @param .abundance A column name for the variable to be plotted
#'
#' @details This function creates a ggpot histogram from a data frame
#'
#' @return A `ggplot` object
#'
#' @examples
#'
#' create_histogram(
#'    tibble::tibble(counts = 1:100),
#'    counts
#' )
#'
#' @export
create_histogram = function(.data, .abundance) {

  .abundance = enquo(.abundance)

  .data %>%
    ggplot(aes(!!.abundance)) +
    geom_histogram()

}
