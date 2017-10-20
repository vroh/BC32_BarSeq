#' Initialization function (loading of the required packages).
#'
#' @details This function loads the required packages
#' @keywords BC32 initialize
#' @export
#' @examples
#' initBC()

initBC <- function() {
  library(ShortRead)
  library(Biostrings)
  library(dplyr)
  library(stringdist)
  library(ggplot2)
  library(gridExtra)
  library(msa)
  library(shiny)
  library(reshape)
}
