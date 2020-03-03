#' A function to export QC plots after the pooling step.
#'
#' @param dir The directory containing the statistics from the pooling step.
#' @details This function shows the frequency of the pooled barcodes according to the distance to the reference barcode.
#' @keywords BC32 processing pooling stats
#' @export
#' @examples
#' poolingStatsPlot(dir = getwd())

poolingStatsPlot <- function(dir = getwd()) {
  files <- paste0(dir, "/", dir())
  statsFiles <- grep("pooling", files)
  toPlot <- files[statsFiles]
  for (file in toPlot) {
    current <- read.delim(file, stringsAsFactors = F)
    pdf(paste0(file, '_plot.pdf'), 10)
    plot(current$distance, current$alt_freq)
    dev.off()
  }
}
