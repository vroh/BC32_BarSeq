#' A function to merge multiple count summaries into one table.
#'
#' @param sampname A data frame with 3 columns containing the name of each sample to
#' process, its file name and the expected multiplex index sequence.
#' @param dir The directory containing the summaries.
#' @param min.count Counts below this number are replaced by zero (considered as noise). Defaults to 10.
#' @param cleancol If true, columns that contain only zeros will be removed. Defaults to F.
#' @details This function merge all count summaries in the sampname table together in one data frame. Counts below
#' min.count are replaced by zeros.
#' @keywords BC32 processing update merge
#' @export
#' @examples
#' counts <- mergeSummaries(sampname,
#'                          getwd(),
#'                          10,
#'                          F)

mergeSummaries <- function(sampname,
                           dir,
                           min.count = 10,
                           cleancol = F) {

  # Load data summaries
  summaries <- paste0(dir, '/', sampname$sample, '_barcode_summary_pooled.txt')
  data.list <- list()
  for (i in 1:length(summaries)) {
    data.list[[i]] <- read.table(summaries[i], sep = '\t', header = T)
  }

  # Merge data
  merged <- NULL
  current <- data.list[1]
  for (i in 2:length(summaries)){
    merged <- merge(current, data.list[i], by = 'seq', all = T)
    current <- merged
  }
  colnames(merged) <- c('seq', sampname$sample)

  # Replace NA's with 0
  merged[is.na(merged)] <- 0

  # Replace counts below threshold by 0 (remove noise)
  for (i in 2:length(merged)) {
    merged[merged[,i] < min.count, i] <- 0
  }

  # Remove columns that only contain zeros
  if(cleancol) {
    merged <- merged[, c(T, as.logical(apply(merged[2:length(merged)], 2, sum) != 0))]
  }

  # Return merged table for counts
  merged
}
